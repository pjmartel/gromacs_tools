#!/usr/bin/env python3
"""
gmx_monitor - live web panel for a running GROMACS MD simulation.

Follows the mdrun log files of one run (the segments of gmx_continue_grompp.sh, or the log
and .partNNNN logs of gmx_continue_extend.sh / a single mdrun) and serves a page that
refreshes itself: state, progress of the segment and of the whole run, speed (ns/day), time
left and finish time, and graphs of temperature, pressure, energies, conserved-energy drift,
constraint rmsd and speed over the simulation. With --edr it also plots the box size, volume
and density from the .edr files (through gmx energy); with --mindist the distance of the protein
to its periodic image, against the cutoff, and its size against the box (gmx mindist -pi); with
--rms the RMSD from the start of the run (gmx rms).

Uses only the Python standard library: runs with any python3 >= 3.8, no packages needed.

Usage:
  gmx_monitor.py <stem> [--dir DIR] [--end NS] [--edr] [--mindist] [--rms] [--refresh S] [--port PORT]
                [--once]
"""
from __future__ import annotations

import argparse
import errno
import json
import math
import os
import re
import shlex
import shutil
import socket
import subprocess
import sys
import tempfile
import threading
import time
from datetime import datetime
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

# ---------------------------------------------------------------------------
# mdrun log parsing
# ---------------------------------------------------------------------------

_MDRUN_HEADER = re.compile(r"^\s*:-\) GROMACS - gmx mdrun, (\S+) \(-:\s*$")
_PARAM = re.compile(r"^   (nsteps|dt|tinit|init-step|nstlog|tcoupl|pcoupl|rvdw|rcoulomb)\s+= (\S+)\s*$")
_PME_OPTIMAL = re.compile(r"^\s+optimal pme grid .*, coulomb cutoff (\d+(?:\.\d+)?)\s*$")
_REF_T = re.compile(r"^   ref-t:\s+(.+?)\s*$")
_REF_P_ROW = re.compile(r"^      ref-p\[\s*(\d)\]=\{\s*(\S+),\s*(\S+),\s*(\S+)\}\s*$")
_NSTEPS_OVERRIDE = re.compile(r"^Overriding nsteps with value passed on the command line: (-?\d+) steps")
_STARTED = re.compile(r"^Started mdrun on rank 0 (.+?)\s*$")
_FINISHED = re.compile(r"^Finished mdrun on rank 0 (.+?)\s*$")
_CHECKPOINT = re.compile(r"^Writing checkpoint, step (\d+) at (.+?)\s*$")
_HOST = re.compile(r"^Hardware detected on host (\S+?)(?: \(.*\))?:\s*$")
_NATOMS = re.compile(r"^There are:\s+(\d+) Atoms\s*$")
_DRIFT = re.compile(r"^\s+Conserved energy drift: (\S+) kJ/mol/ps per atom\s*$")
_STEP_HEADER = re.compile(r"^\s+Step\s+Time\s*$")
_STEP_VALUES = re.compile(r"^\s+(-?\d+)\s+(-?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)\s*$")
_PERF_HEADER = re.compile(r"^\s+\(ns/day\)")
_PERF = re.compile(r"^Performance:\s+(.+?)\s*$")
_ISSUES = [
    ("LINCS warnings", re.compile(r"LINCS WARNING")),
    ("SETTLE warnings", re.compile(r"can not be settled")),
    ("atoms moved too far", re.compile(r"moved more than")),
]
_MONTHS = {m: i for i, m in enumerate(
    ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"], 1)}

# Energy terms kept from the log (the rest of each energy block is not needed)
LOG_TERMS = ["Temperature", "Pressure (bar)", "Potential", "Kinetic En.", "Total Energy",
             "Conserved En.", "Constr. rmsd"]

_FIELD = 15   # width of each column of the energy blocks in the log


def parse_log_date(text: str) -> datetime | None:
    """'Thu Jan  8 16:01:11 2026' -> datetime (locale independent); None if not that format."""
    parts = text.split()
    if len(parts) != 5 or parts[1] not in _MONTHS:
        return None
    try:
        hh, mm, ss = (int(v) for v in parts[3].split(":"))
        return datetime(int(parts[4]), _MONTHS[parts[1]], int(parts[2]), hh, mm, ss)
    except ValueError:
        return None


def is_mdrun_log(path: Path) -> bool | None:
    """True if *path* starts like an mdrun log, False if it starts like something else (e.g.
    the combined output log of gmx_continue_grompp.sh, which starts with gmx grompp), None if
    it is still empty."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(2048).decode(errors="replace").splitlines()
    except OSError:
        return None
    lines = [line for line in head if line.strip()]
    if not lines:
        return None
    return bool(_MDRUN_HEADER.match(lines[0]))


class MdLog:
    """One mdrun log file, read incrementally: each update() parses only the new lines."""

    def __init__(self, path: Path):
        self.path = path
        self._reset()

    def _reset(self):
        self.offset = 0
        self._partial = b""
        self.version = None
        self.host = None
        self.natoms = None
        self.params: dict[str, str] = {}
        self.ref_t: list[float] = []
        self._ref_p: dict[int, float] = {}
        self.starts: list[datetime] = []
        self.start_times: list[float] = []      # simulation time (ps) of the first frame after each start
        self._awaiting_start_frame = False
        self.checkpoints: list[tuple[datetime, int]] = []
        self.finished: datetime | None = None
        self.ns_per_day: float | None = None
        self.coulomb_tuned: float | None = None  # Coulomb cutoff chosen by the PME tuning
        self.drift: float | None = None          # kJ/mol/ps per atom, reported at the end
        self.issues: dict[str, int] = {}
        self.fatal: str | None = None
        self._fatal_next = False
        self.steps: list[int] = []
        self.times: list[float] = []
        self.cols: dict[str, list[float]] = {term: [] for term in LOG_TERMS}
        self.bad_frames = 0
        self._state = "idle"
        self._pending = None
        self._names: list[str] = []
        self._values: dict[str, float] = {}
        self._perf_columns: list[str] | None = None
        self.mtime = 0.0

    # -- reading --
    def update(self):
        try:
            st = self.path.stat()
        except OSError:
            return
        if st.st_size < self.offset:          # truncated or replaced: read it again
            self._reset()
        self.mtime = st.st_mtime
        if st.st_size == self.offset:
            return
        with open(self.path, "rb") as fh:
            fh.seek(self.offset)
            data = fh.read()
        self.offset += len(data)
        data = self._partial + data
        lines = data.split(b"\n")
        self._partial = lines.pop()           # an incomplete last line waits for the next update
        for raw in lines:
            self._line(raw.decode(errors="replace").rstrip("\r"))

    # -- derived values --
    def param_float(self, name: str) -> float | None:
        try:
            return float(self.params[name])
        except (KeyError, ValueError):
            return None

    @property
    def dt(self) -> float | None:
        return self.param_float("dt")

    @property
    def seg_start(self) -> float | None:
        """Simulation time (ps) at which this run starts."""
        tinit, dt, init = self.param_float("tinit"), self.dt, self.param_float("init-step")
        if tinit is None or dt is None:
            return None
        return tinit + (init or 0) * dt

    @property
    def seg_end(self) -> float | None:
        """Simulation time (ps) at which this run ends (None for nsteps = -1)."""
        start, nsteps, dt = self.seg_start, self.param_float("nsteps"), self.dt
        if start is None or nsteps is None or dt is None or nsteps < 0:
            return None
        return start + nsteps * dt

    @property
    def cutoff(self) -> float | None:
        """The largest non-bonded cutoff (nm): rvdw, rcoulomb, or the Coulomb cutoff the PME
        tuning chose."""
        values = [v for v in (self.param_float("rvdw"), self.param_float("rcoulomb"),
                              self.coulomb_tuned) if v is not None]
        return max(values) if values else None

    @property
    def ref_p(self) -> float | None:
        if len(self._ref_p) == 3:
            return sum(self._ref_p.values()) / 3
        return None

    # -- one line --
    def _line(self, line: str):
        if self.version is None:
            m = _MDRUN_HEADER.match(line)
            if m:
                self.version = m.group(1)
                return
        if self._energy_line(line):
            return
        m = _PARAM.match(line)
        if m:
            # the last value counts: an extended run (convert-tpr) prints a new nsteps
            self.params[m.group(1)] = m.group(2)
            return
        m = _PME_OPTIMAL.match(line)
        if m:
            self.coulomb_tuned = float(m.group(1))
            return
        m = _NSTEPS_OVERRIDE.match(line)
        if m:
            self.params["nsteps"] = m.group(1)       # mdrun -nsteps
            return
        m = _REF_T.match(line)
        if m:
            try:
                self.ref_t = [float(v) for v in m.group(1).split()]
            except ValueError:
                pass
            return
        m = _REF_P_ROW.match(line)
        if m:
            row = int(m.group(1))
            try:
                self._ref_p[row] = float(m.group(2 + row))
            except (ValueError, IndexError):
                pass
            return
        m = _STARTED.match(line)
        if m:
            when = parse_log_date(m.group(1))
            if when:
                self.starts.append(when)
                self._awaiting_start_frame = True
                self.finished = None
            return
        m = _CHECKPOINT.match(line)
        if m:
            when = parse_log_date(m.group(2))
            if when:
                self.checkpoints.append((when, int(m.group(1))))
            return
        m = _FINISHED.match(line)
        if m:
            self.finished = parse_log_date(m.group(1))
            return
        m = _HOST.match(line)
        if m:
            self.host = m.group(1)
            return
        m = _NATOMS.match(line)
        if m:
            self.natoms = int(m.group(1))
            return
        m = _DRIFT.match(line)
        if m:
            try:
                self.drift = float(m.group(1))
            except ValueError:
                pass
            return
        if _PERF_HEADER.match(line):
            self._perf_columns = re.findall(r"\(([^)]+)\)", line)
            return
        m = _PERF.match(line)
        if m and self._perf_columns:
            values = m.group(1).split()
            if "ns/day" in self._perf_columns and len(values) == len(self._perf_columns):
                try:
                    self.ns_per_day = float(values[self._perf_columns.index("ns/day")])
                except ValueError:
                    pass
            self._perf_columns = None
            return
        if self._fatal_next and line.strip():
            self.fatal = line.strip()
            self._fatal_next = False
            return
        if "Fatal error:" in line:
            self.fatal = line.split("Fatal error:", 1)[1].strip() or None
            self._fatal_next = self.fatal is None
            return
        for name, pattern in _ISSUES:
            if pattern.search(line):
                self.issues[name] = self.issues.get(name, 0) + 1
                return

    def _energy_line(self, line: str) -> bool:
        """Energy frames: '   Step   Time', its values, a blank line, '   Energies (kJ/mol)',
        then pairs of name/value lines up to a blank line. Anything else ends the frame; a frame
        whose names and values don't match is dropped, never guessed."""
        state = self._state
        if _STEP_HEADER.match(line):
            self._state = "step"
            return True
        if state == "step":
            m = _STEP_VALUES.match(line)
            if m:
                self._pending = (int(m.group(1)), float(m.group(2)))
                self._state = "after_step"
                return True
            self._state = "idle"
            return False
        if state == "after_step":
            if not line.strip():
                return True
            if line.strip() == "Energies (kJ/mol)":
                self._state, self._names, self._values = "names", [], {}
                return True
            if _CHECKPOINT.match(line):
                return False      # the last frame: mdrun writes the checkpoint line before the energies
            self._state = "idle"
            return False
        if state == "names":
            if not line.strip():
                self._commit()
                self._state = "idle"
                return True
            text = line.rstrip()
            if len(text) % _FIELD:
                self._state = "idle"      # the frame was cut off (e.g. by a crash): not an error
                return False
            self._names = [text[i:i + _FIELD].strip() for i in range(0, len(text), _FIELD)]
            self._state = "values"
            return True
        if state == "values":
            values = line.split()
            if not values:
                self._state = "idle"      # cut off after a names line
                return True
            if len(values) != len(self._names):
                self._drop()
                return False
            try:
                for name, value in zip(self._names, values):
                    self._values[name] = float(value)
            except ValueError:
                self._drop()
                return False
            self._state = "names"
            return True
        return False

    def _drop(self):
        self.bad_frames += 1
        self._state = "idle"

    def _commit(self):
        if not self._pending or "Temperature" not in self._values:
            return
        step, t = self._pending
        self._pending = None
        # After a restart from a checkpoint (mdrun -append), the frames after the checkpoint are
        # in the log twice: the newer ones replace the older
        if self.times and t <= self.times[-1]:
            keep = 0
            while keep < len(self.times) and self.times[keep] < t:
                keep += 1
            del self.steps[keep:], self.times[keep:]
            for col in self.cols.values():
                del col[keep:]
        self.steps.append(step)
        self.times.append(t)
        for term in LOG_TERMS:
            value = self._values.get(term, math.nan)
            self.cols[term].append(value)
            if term == "Temperature" and not math.isfinite(value):
                self.issues["non-finite energies (NaN)"] = self.issues.get("non-finite energies (NaN)", 0) + 1
        if self._awaiting_start_frame:
            self.start_times.append(t)
            self._awaiting_start_frame = False


# ---------------------------------------------------------------------------
# Finding the run's files
# ---------------------------------------------------------------------------

_VALUE_OPTIONS = {"--template", "--initial", "--timestep", "--title", "--topology", "--plumed",
                  "--posres-ref", "--commands-file", "--tpr", "--cpt"}


def detect_total_end(directory: Path, stem: str) -> tuple[float, str] | None:
    """Total length of the run (ps) from the invocation recorded by gmx_continue_grompp.sh /
    gmx_continue_extend.sh (their *_commands*.sh files): (end_ps, file name), the newest file
    whose <basename> <replica> is *stem*; None if there is none."""
    found = []
    for path in directory.glob("gmx_continue_*_commands*.sh"):
        try:
            with open(path, errors="replace") as fh:
                head = [next(fh, "") for _ in range(10)]
        except OSError:
            continue
        line = next((l for l in head if l.startswith("# Invocation: ")), None)
        if not line:
            continue
        try:
            tokens = shlex.split(line[len("# Invocation: "):])[1:]
        except ValueError:
            continue
        positional, i = [], 0
        while i < len(tokens):
            if tokens[i] in _VALUE_OPTIONS:
                i += 2
            elif tokens[i].startswith("--"):
                i += 1
            else:
                positional.append(tokens[i])
                i += 1
        if len(positional) < 5 or f"{positional[0]}_{positional[1]}" != stem:
            continue
        try:
            end = float(positional[3])
        except ValueError:
            continue
        end_ps = end if "--ps" in tokens else end * 1000
        found.append((path.stat().st_mtime, end_ps, path.name))
    if not found:
        return None
    _, end_ps, name = max(found)
    return end_ps, name


def find_logs(directory: Path, stem: str) -> tuple[str, list[Path]]:
    """The run's mdrun logs, in order, and the layout: 'segments' (<stem>_<start>_<end>.log,
    gmx_continue_grompp.sh) or 'parts' (<stem>.log and <stem>.partNNNN.log, mdrun -noappend /
    gmx_continue_extend.sh, or one log)."""
    seg_re = re.compile(rf"^{re.escape(stem)}_(\d+(?:\.\d+)?)_(\d+(?:\.\d+)?)\.log$")
    part_re = re.compile(rf"^{re.escape(stem)}(?:\.part(\d+))?\.log$")
    segments, parts = [], []
    for path in directory.iterdir():
        m = seg_re.match(path.name)
        if m:
            segments.append(((float(m.group(1)), float(m.group(2))), path))
            continue
        m = part_re.match(path.name)
        if m:
            parts.append((int(m.group(1) or 1), path))
    if segments:
        return "segments", [p for _, p in sorted(segments)]
    return "parts", [p for _, p in sorted(parts)]


# ---------------------------------------------------------------------------
# .edr: box, volume, density (gmx energy)
# ---------------------------------------------------------------------------

EDR_TERMS = ["Box-X", "Box-Y", "Box-Z", "Volume", "Density"]


class EdrReader:
    """Reads the box/volume/density terms of the run's .edr files with gmx energy, in a
    background thread; a file is read again only when it has changed."""

    def __init__(self, gmx: str, interval: float):
        self.gmx = gmx
        self.interval = interval
        self.status = "waiting for the first read"
        self._cache: dict[Path, tuple[tuple[float, int], dict[str, list[float]]]] = {}
        self._lock = threading.Lock()
        self._running = False
        self._last = 0.0

    def maybe_refresh(self, edrs: list[Path]):
        """Start a background read if the interval has passed and none is running."""
        with self._lock:
            if self._running or time.time() - self._last < self.interval:
                return
            self._running = True
        threading.Thread(target=self._refresh, args=(edrs,), daemon=True).start()

    def _refresh(self, edrs: list[Path]):
        try:
            for path in edrs:
                try:
                    st = path.stat()
                except OSError:
                    continue
                key = (st.st_mtime, st.st_size)
                cached = self._cache.get(path)
                if cached and cached[0] == key:
                    continue
                data = self._read(path)
                if data is not None:
                    self._cache[path] = (key, data)
            n = len([p for p in edrs if p in self._cache])
            self.status = f"read {n} of {len(edrs)} .edr file(s) at {datetime.now():%H:%M:%S}"
        except Exception as exc:     # never let a gmx problem stop the monitor
            self.status = f"error: {exc}"
        finally:
            with self._lock:
                self._running = False
                self._last = time.time()

    def _read(self, path: Path) -> dict[str, list[float]] | None:
        path = path.resolve()        # gmx runs in a temporary folder
        with tempfile.TemporaryDirectory(prefix="gmx_monitor_") as tmp:
            # 1. the list of terms in this file (gmx energy prints it, then exits on empty input)
            try:
                listing = subprocess.run([self.gmx, "energy", "-f", str(path), "-o", "x.xvg"],
                                         input="\n", capture_output=True, text=True,
                                         cwd=tmp, timeout=120)
            except (OSError, subprocess.SubprocessError) as exc:
                self.status = f"cannot run {self.gmx} energy: {exc}"
                return None
            numbers = self._term_numbers(listing.stdout + listing.stderr)
            if numbers is None:
                self.status = f"could not read the term list of {path.name} (gmx energy)"
                return None
            wanted = [(t, numbers[t]) for t in EDR_TERMS if t in numbers]
            if not wanted:
                return {}          # constant volume: no box terms
            try:
                subprocess.run([self.gmx, "energy", "-f", str(path), "-o", "edr.xvg"],
                               input="".join(f"{n}\n" for _, n in wanted) + "\n",
                               capture_output=True, text=True, cwd=tmp, timeout=300)
            except (OSError, subprocess.SubprocessError) as exc:
                self.status = f"gmx energy failed on {path.name}: {exc}"
                return None
            return self._parse_xvg(Path(tmp) / "edr.xvg", [t for t, _ in wanted])

    @staticmethod
    def _term_numbers(text: str) -> dict[str, int] | None:
        """{name: number} from the table between the dashed line after 'Select the terms' and
        the next blank line."""
        lines = text.splitlines()
        try:
            start = next(i for i, l in enumerate(lines) if l.startswith("Select the terms"))
            dash = next(i for i in range(start, len(lines)) if re.match(r"^-{20,}\s*$", lines[i]))
        except StopIteration:
            return None
        numbers = {}
        for line in lines[dash + 1:]:
            if not line.strip():
                break
            for number, name in re.findall(r"(\d+)\s+(\S+)", line):
                numbers[name] = int(number)
        return numbers or None

    @staticmethod
    def _parse_xvg(path: Path, terms: list[str]) -> dict[str, list[float]] | None:
        """Columns of an xvg from gmx energy, matched by their legends to *terms*."""
        legends: dict[int, str] = {}
        rows = []
        try:
            with open(path) as fh:
                for line in fh:
                    m = re.match(r'^@ s(\d+) legend "(.*)"\s*$', line)
                    if m:
                        legends[int(m.group(1))] = m.group(2)
                    elif line.strip() and line[0] not in "#@":
                        rows.append(line.split())
        except OSError:
            return None
        if sorted(legends.values()) != sorted(terms):
            return None
        data: dict[str, list[float]] = {"time": []}
        data.update({legends[i]: [] for i in legends})
        for row in rows:
            if len(row) != len(legends) + 1:
                continue
            try:
                values = [float(v) for v in row]
            except ValueError:
                continue
            data["time"].append(values[0])
            for i in legends:
                data[legends[i]].append(values[i + 1])
        return data

    def series(self, edrs: list[Path]) -> dict[str, list[float]]:
        """The cached data of *edrs*, joined in order (a later file replaces overlapping times);
        a term missing from a file is NaN there."""
        parts = [self._cache[p][1] for p in edrs if p in self._cache and self._cache[p][1]]
        terms = [t for t in EDR_TERMS if any(t in d for d in parts)]
        out: dict[str, list[float]] = {"time": [], **{t: [] for t in terms}}
        for data in parts:
            if out["time"] and data["time"]:
                keep = 0
                while keep < len(out["time"]) and out["time"][keep] < data["time"][0]:
                    keep += 1
                for values in out.values():
                    del values[keep:]
            out["time"].extend(data["time"])
            for t in terms:
                out[t].extend(data.get(t, [math.nan] * len(data["time"])))
        return out


# ---------------------------------------------------------------------------
# Trajectory analyses in the background: gmx mindist -pi, gmx rms
# ---------------------------------------------------------------------------

def trajectory_due(reader: "TrajectoryReader", segment) -> bool:
    """True if *segment* still has frames for the reader to compute (used by --once)."""
    xtc, tpr, finished, t_now = segment
    if xtc in reader._complete or t_now is None or not xtc.exists() or not tpr.exists():
        return False
    rows = reader._rows.get(xtc, [])
    return not rows or finished or t_now >= rows[-1][0] + reader.dt_ps * 1.05


def gmx_error_message(stderr: str) -> str:
    """The message of a GROMACS error: the lines after its "Fatal error:" / "Error in user
    input:" heading, up to the dashed line (else the last line of the output)."""
    lines = stderr.splitlines()
    for i, line in enumerate(lines):
        if line.startswith(("Fatal error", "Error in user input", "Inconsistency in user input")):
            message = []
            for follow in lines[i + 1:]:
                if follow.startswith("-----") or follow.startswith("For more information"):
                    break
                if follow.strip():
                    message.append(follow.strip())
            return " ".join(message)[:300] or line.strip()
    rest = [l.strip() for l in lines if l.strip()]
    return rest[-1][:300] if rest else "no error message"


class TrajectoryReader:
    """Runs a gmx trajectory analysis on the run's trajectories in a background thread, one frame
    per `dt_ps`: a first pass over each segment's trajectory, then only the frames added since
    (gmx seeks to -b quickly in an .xtc, so this stays cheap as the trajectory grows). Subclasses
    give the gmx command (_command), what it reads on standard input (_input) and the parser of
    its .xvg (_parse); rows start with the time in ps."""

    tool = ""          # gmx tool name, for messages
    option = ""        # command-line option of the group, for messages

    def __init__(self, gmx: str, group: str, dt_ps: float, index: Path | None = None):
        self.gmx = gmx
        self.group = group
        self.index = index
        self.dt_ps = dt_ps
        self.status = "waiting for the first calculation"
        self.error: str | None = None
        self._rows: dict[Path, list[tuple[float, ...]]] = {}     # trajectory -> rows
        self._complete: set[Path] = set()                         # finished segments, all read
        self._lock = threading.Lock()
        self._running = False

    @property
    def running(self) -> bool:
        return self._running

    def _command(self, xtc: Path, tpr: Path, out: str) -> list[str]:
        raise NotImplementedError

    def _input(self) -> str:
        return self.group + "\n"

    def _parse(self, path: Path):
        raise NotImplementedError

    def maybe_refresh(self, segments):
        """segments: [(xtc, tpr, finished, sim time reached by its log in ps)]. Starts a background
        calculation if one of them has at least dt_ps of new trajectory."""
        todo = []
        for xtc, tpr, finished, t_now in segments:
            if xtc in self._complete or t_now is None:
                continue
            rows = self._rows.get(xtc, [])
            begin = rows[-1][0] + self.dt_ps if rows else None
            # the frame at `begin` is written once the log is past it (one more frame of margin)
            if begin is None or t_now >= begin + 0.05 * self.dt_ps or finished:
                todo.append((xtc, tpr, finished, begin))
        if not todo:
            return
        with self._lock:
            if self._running:
                return
            self._running = True
        threading.Thread(target=self._refresh, args=(todo,), daemon=True).start()

    def _refresh(self, todo):
        try:
            for xtc, tpr, finished, begin in todo:
                if not xtc.exists() or not tpr.exists():
                    self.status = f"waiting for {xtc.name if not xtc.exists() else tpr.name}"
                    continue
                self.status = f"computing {xtc.name}" + (f" from {begin / 1000:g} ns" if begin else "")
                rows = self._run(xtc, tpr, begin)
                if rows is None:
                    return            # error: reported in self.error, tried again next time
                old = self._rows.setdefault(xtc, [])
                old.extend(r for r in rows if not old or r[0] > old[-1][0])
                if finished:
                    self._complete.add(xtc)
            n = sum(len(r) for r in self._rows.values())
            self.error = None
            self.status = f"{n} frame(s), one per {self.dt_ps / 1000:g} ns, updated {datetime.now():%H:%M:%S}"
        except Exception as exc:     # never let a gmx problem stop the monitor
            self.error = f"gmx {self.tool}: {exc}"
        finally:
            with self._lock:
                self._running = False

    def _run(self, xtc: Path, tpr: Path, begin: float | None):
        command = [self.gmx] + self._command(xtc, tpr, "out.xvg") + ["-dt", f"{self.dt_ps:g}"]
        if begin is not None:
            command += ["-b", f"{begin:g}"]
        if self.index is not None:
            command += ["-n", str(self.index.resolve())]
        if shutil.which("nice"):
            command = ["nice", "-n", "19"] + command    # never slow down the simulation
        with tempfile.TemporaryDirectory(prefix="gmx_monitor_") as tmp:
            try:
                run = subprocess.run(command, input=self._input(), capture_output=True,
                                     text=True, cwd=tmp, timeout=7200)
            except (OSError, subprocess.SubprocessError) as exc:
                self.error = f"cannot run {self.gmx} {self.tool}: {exc}"
                return None
            out = Path(tmp) / "out.xvg"
            if run.returncode != 0 or not out.exists():
                if f"No such group '{self.group}'" in run.stdout:
                    groups = re.findall(r"^Group\s+\d+ \(\s*(.+?)\) has", run.stderr, re.M)
                    source = self.index.name if self.index is not None else "the default groups"
                    self.error = (f"gmx {self.tool}: no group '{self.group}' in {source} ({self.option});"
                                  f" groups: " + ", ".join(dict.fromkeys(groups)))
                    return None
                if begin is not None and "doesn't exist" in gmx_error_message(run.stderr):
                    return []      # no frame at -b yet (still in mdrun's output buffer), or all read
                self.error = f"gmx {self.tool} failed on {xtc.name}: {gmx_error_message(run.stderr)}"
                return None
            return self._parse(out)

    def series(self, xtcs: list[Path]) -> list[tuple[float, ...]]:
        """The rows of *xtcs*, joined in order (a later trajectory replaces overlapping times)."""
        out: list[tuple[float, ...]] = []
        for xtc in xtcs:
            rows = self._rows.get(xtc, [])
            if out and rows:
                while out and out[-1][0] >= rows[0][0]:
                    out.pop()
            out.extend(rows)
        return out


def read_xvg(path: Path) -> tuple[dict[str, str], dict[int, str], list[list[str]]]:
    """(header {title, xaxis, yaxis}, legends {set: text}, data rows as fields) of an .xvg."""
    header: dict[str, str] = {}
    legends: dict[int, str] = {}
    rows = []
    with open(path) as fh:
        for line in fh:
            m = re.match(r'^@\s+(title|xaxis\s+label|yaxis\s+label)\s+"(.*)"\s*$', line)
            if m:
                header[m.group(1).split()[0]] = m.group(2)
                continue
            m = re.match(r'^@ s(\d+) legend "(.*)"\s*$', line)
            if m:
                legends[int(m.group(1))] = m.group(2)
            elif line.strip() and line[0] not in "#@":
                rows.append(line.split())
    return header, legends, rows


def float_rows(rows: list[list[str]], width: int) -> list[tuple[float, ...]]:
    """The rows with exactly *width* numeric fields, as floats."""
    out = []
    for fields in rows:
        if len(fields) != width:
            continue
        try:
            out.append(tuple(float(v) for v in fields))
        except ValueError:
            continue
    return out


class MindistReader(TrajectoryReader):
    """gmx mindist -pi: distance of the group to its periodic image, its largest internal
    distance and the box lengths."""

    tool, option = "mindist", "--mindist-group"
    COLUMNS = ["min per.", "max int.", "box1", "box2", "box3"]

    def _command(self, xtc, tpr, out):
        return ["mindist", "-s", str(tpr.resolve()), "-f", str(xtc.resolve()), "-pi", "-od", out]

    def _parse(self, path: Path):
        """Rows (time ps, min per., max int., box1, box2, box3), checked against the legends."""
        _, legends, rows = read_xvg(path)
        if [legends.get(i) for i in range(len(self.COLUMNS))] != self.COLUMNS:
            self.error = f"unexpected gmx mindist output (legends {sorted(legends.values())})"
            return None
        return float_rows(rows, len(self.COLUMNS) + 1)


class RmsReader(TrajectoryReader):
    """gmx rms: RMSD of the group after fitting it (the same group) to a reference structure:
    --rms-ref, else the .tpr of the run's first segment (the start of the run; each segment's own
    .tpr starts where the previous one ended)."""

    tool, option = "rms", "--rms-group"

    def __init__(self, gmx: str, group: str, dt_ps: float, index: Path | None = None,
                 reference: Path | None = None):
        super().__init__(gmx, group, dt_ps, index)
        self.reference = reference
        self.first_tpr: Path | None = None      # set by the run: the default reference

    @property
    def reference_used(self) -> Path | None:
        return self.reference or self.first_tpr

    def _command(self, xtc, tpr, out):
        reference = self.reference_used or tpr
        return ["rms", "-s", str(reference.resolve()), "-f", str(xtc.resolve()), "-o", out, "-tu", "ps"]

    def _input(self):
        return f"{self.group}\n{self.group}\n"     # fit group, RMSD group

    def _parse(self, path: Path):
        """Rows (time ps, RMSD nm), checked against the title and axis labels."""
        header, _, rows = read_xvg(path)
        if header.get("title") != "RMSD" or header.get("xaxis") != "Time (ps)" \
                or header.get("yaxis") != "RMSD (nm)":
            self.error = f"unexpected gmx rms output ({header})"
            return None
        return float_rows(rows, 2)


# ---------------------------------------------------------------------------
# The run: all logs together, status, speed, time left
# ---------------------------------------------------------------------------

def fmt_duration(seconds: float | None) -> str:
    if seconds is None or not math.isfinite(seconds) or seconds < 0:
        return "-"
    minutes = int(round(seconds / 60))
    if minutes < 60:
        return f"{minutes} min"
    hours, minutes = divmod(minutes, 60)
    if hours < 48:
        return f"{hours} h {minutes:02d} min"
    days, hours = divmod(hours, 24)
    return f"{days} d {hours} h"


def decimate(n: int, max_points: int) -> list[int]:
    """Indices of at most *max_points* evenly spread samples of n, always with the last one."""
    if n <= max_points:
        return list(range(n))
    stride = math.ceil(n / max_points)
    idx = list(range(0, n, stride))
    if idx[-1] != n - 1:
        idx.append(n - 1)
    return idx


def running_mean(values: list[float], window: int) -> list[float]:
    out, acc, count, queue = [], 0.0, 0, []
    for v in values:
        finite = math.isfinite(v)
        queue.append(v if finite else None)
        if finite:
            acc += v
            count += 1
        if len(queue) > window:
            old = queue.pop(0)
            if old is not None:
                acc -= old
                count -= 1
        out.append(acc / count if count else math.nan)
    return out


def linear_slope(xs: list[float], ys: list[float]) -> float | None:
    pts = [(x, y) for x, y in zip(xs, ys) if math.isfinite(y)]
    if len(pts) < 3:
        return None
    mx = sum(x for x, _ in pts) / len(pts)
    my = sum(y for _, y in pts) / len(pts)
    sxx = sum((x - mx) ** 2 for x, _ in pts)
    if sxx == 0:
        return None
    return sum((x - mx) * (y - my) for x, y in pts) / sxx


class Run:
    STALL_MIN_S = 300          # never call a run stalled before 5 min without log output,
    STALL_MAX_S = 3600         # always after 1 h
    BETWEEN_SEGMENTS_S = 600   # after a segment finished: "waiting for the next one" this long
    RATE_WINDOW_S = 900        # current speed: over the last 15 min
    PRESSURE_AVG_FRAMES = 50
    MINDIST_MARGIN = 1.2       # warn when closer to the periodic image than 1.2 x the cutoff

    def __init__(self, directory: Path, stem: str, end_ns: float | None, edr: EdrReader | None,
                 max_points: int, mindist: MindistReader | None = None, rms: RmsReader | None = None):
        self.directory = directory
        self.stem = stem
        self.end_override = end_ns * 1000 if end_ns is not None else None
        self.edr = edr
        self.mindist = mindist
        self.rms = rms
        self.max_points = max_points
        self.logs: dict[Path, MdLog] = {}
        self.order: list[Path] = []
        self.layout = "parts"
        self.samples: list[tuple[float, float]] = []    # (wall time, sim time ps) seen while running
        self.lock = threading.Lock()

    def update(self):
        layout, paths = find_logs(self.directory, self.stem)
        order = []
        for path in paths:
            if path not in self.logs:
                kind = is_mdrun_log(path)
                if kind is None:
                    continue           # still empty
                if not kind:
                    continue           # not an mdrun log (e.g. the combined grompp/mdrun output)
                self.logs[path] = MdLog(path)
            order.append(path)
        self.layout, self.order = layout, order
        for path in order:
            self.logs[path].update()
        current = self.current
        if current and current.times and current.finished is None:
            t = current.times[-1]
            if not self.samples or t > self.samples[-1][1]:
                self.samples.append((current.mtime, t))
        if self.edr:
            self.edr.maybe_refresh(self.edr_paths())
        if self.mindist:
            self.mindist.maybe_refresh(self.trajectory_segments())
        if self.rms:
            segments = self.trajectory_segments()
            if segments:
                self.rms.first_tpr = segments[0][1]
            self.rms.maybe_refresh(segments)

    @property
    def current(self) -> MdLog | None:
        return self.logs[self.order[-1]] if self.order else None

    def edr_paths(self) -> list[Path]:
        return [p.with_suffix(".edr") for p in self.order if p.with_suffix(".edr").exists()]

    def trajectory_segments(self):
        """[(xtc, tpr, finished, sim time reached in ps)] per log, for the trajectory analyses
        (--mindist, --rms): <log name>.xtc and .tpr
        (.partNNNN logs have no .tpr of their own: <stem>.tpr)."""
        out = []
        for path in self.order:
            log = self.logs[path]
            tpr = path.with_suffix(".tpr")
            if not tpr.exists():
                tpr = self.directory / f"{self.stem}.tpr"
            out.append((path.with_suffix(".xtc"), tpr, log.finished is not None,
                        log.times[-1] if log.times else None))
        return out

    # -- whole run series --
    def merged(self) -> tuple[list[float], dict[str, list[float]], list[float]]:
        times: list[float] = []
        cols: dict[str, list[float]] = {term: [] for term in LOG_TERMS}
        boundaries = []
        for path in self.order:
            log = self.logs[path]
            if not log.times:
                continue
            if times:
                boundaries.append(log.times[0])
                keep = 0
                while keep < len(times) and times[keep] < log.times[0]:
                    keep += 1
                del times[keep:]
                for col in cols.values():
                    del col[keep:]
            times.extend(log.times)
            for term in LOG_TERMS:
                cols[term].extend(log.cols[term])
        return times, cols, boundaries

    def speed_points(self) -> list[tuple[float, float]]:
        """(sim time ps, ns/day) from consecutive (wall, sim time) pairs at least 60 s apart:
        checkpoint lines and starts of every log, and what the monitor saw itself."""
        pairs = []
        for path in self.order:
            log = self.logs[path]
            dt = log.dt
            if dt is None:
                continue
            for when, t in zip(log.starts, log.start_times):
                pairs.append((when.timestamp(), t, "start"))
            for when, step in log.checkpoints:
                pairs.append((when.timestamp(), step * dt + (log.param_float("tinit") or 0), "cpt"))
        pairs.extend((wall, t, "seen") for wall, t in self.samples)
        pairs.sort(key=lambda p: p[1])
        points, last = [], None
        for wall, t, kind in pairs:
            if kind == "start":
                last = (wall, t)          # a restart: don't measure across the gap
                continue
            if last and wall - last[0] >= 60 and t > last[1]:
                points.append(((t + last[1]) / 2, (t - last[1]) / 1000 / ((wall - last[0]) / 86400)))
                last = (wall, t)
            elif last is None:
                last = (wall, t)
        return points

    def rates(self) -> tuple[float | None, float | None]:
        """(current, average since the last start) speed in ns/day of the current log."""
        log = self.current
        if not log or not log.times or not log.starts or not log.start_times:
            return None, None
        elapsed = log.mtime - log.starts[-1].timestamp()
        done = log.times[-1] - log.start_times[-1]
        avg = done / 1000 / (elapsed / 86400) if elapsed > 30 and done > 0 else None
        if log.finished and log.ns_per_day:
            return log.ns_per_day, log.ns_per_day
        recent = [(w, t) for w, t in self.samples if w >= log.mtime - self.RATE_WINDOW_S
                  and t >= log.start_times[-1]]
        cpts = [(w.timestamp(), s * log.dt + (log.param_float("tinit") or 0))
                for w, s in log.checkpoints if w.timestamp() >= log.mtime - self.RATE_WINDOW_S] \
            if log.dt else []
        pts = sorted(set(recent + cpts))
        if len(pts) >= 2 and pts[-1][0] - pts[0][0] >= 60 and pts[-1][1] > pts[0][1]:
            now = (pts[-1][1] - pts[0][1]) / 1000 / ((pts[-1][0] - pts[0][0]) / 86400)
        else:
            now = avg
        return now, avg

    # -- everything for the page --
    def snapshot(self) -> dict:
        log = self.current
        now = time.time()
        out: dict = {"stem": self.stem, "dir": str(self.directory.resolve()), "layout": self.layout,
                     "host": socket.gethostname().split(".")[0],
                     "now": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                     "files": [], "state": "no log", "detail": ""}
        for i, path in enumerate(self.order):
            lg = self.logs[path]
            out["files"].append({"name": path.name, "start_ns": _ns(lg.seg_start), "end_ns": _ns(lg.seg_end),
                                 "finished": lg.finished is not None, "frames": len(lg.times)})
        if not log:
            out["detail"] = f"no mdrun log for '{self.stem}' in {self.directory}"
            return out
        times, cols, boundaries = self.merged()
        t_now = times[-1] if times else None
        dt = log.dt
        nstlog = log.param_float("nstlog")
        seg_start, seg_end = log.seg_start, log.seg_end
        # After an extension (convert-tpr) and a restart with -append, mdrun doesn't print the
        # parameters again: the log's nsteps is then out of date (the run is already past it)
        stale = seg_end is not None and t_now is not None and dt is not None and t_now > seg_end + dt
        total = (self.end_override, "--end") if self.end_override else None
        if total is None:
            found = detect_total_end(self.directory, self.stem)
            if found:
                total = (found[0], found[1])
        if total is None and self.layout == "parts" and seg_end is not None and not stale:
            total = (seg_end, "the log (nsteps)")
        if stale:
            seg_end = total[0] if total and total[0] >= t_now - dt else None
        rate_now, rate_avg = self.rates()
        rate = rate_now or rate_avg

        def left(target):
            if target is None or t_now is None or not rate:
                return None
            return max(0.0, target - t_now) / 1000 / rate * 86400

        seg_left, total_left = left(seg_end), left(total[0] if total else None)
        if total_left == 0.0:
            total_left = None
        if log.finished:
            seg_left = None
        # state
        age = now - log.mtime
        total_done = total is not None and t_now is not None and dt is not None and t_now >= total[0] - dt
        if log.fatal:
            state, detail = "error", log.fatal
        elif not log.times:
            state, detail = "starting", "mdrun has started; no energy frame yet"
        elif log.finished:
            since = now - log.finished.timestamp()
            ended_early = seg_end is not None and dt is not None and t_now < seg_end - dt
            if ended_early:
                state, detail = "stopped", (f"mdrun finished at {_ns(t_now):g} ns, before the end of "
                                            f"this log's run ({_ns(seg_end):g} ns): -maxh or -nsteps?")
            elif total_done or (total is None and self.layout == "parts"):
                state, detail = "finished", f"finished {log.finished:%Y-%m-%d %H:%M}"
            elif since < self.BETWEEN_SEGMENTS_S:
                state, detail = "between segments", "segment finished; waiting for the next one"
            else:
                state, detail = "stopped", (f"segment finished {fmt_duration(since)} ago and no next "
                                            f"segment has started")
        else:
            expected = (nstlog * dt) / 1000 / rate * 86400 if rate and nstlog and dt else None
            limit = min(self.STALL_MAX_S, max(self.STALL_MIN_S, 10 * expected)) if expected \
                else self.STALL_MIN_S * 2
            if age > limit:
                state, detail = "stalled?", f"no log output for {fmt_duration(age)}"
            else:
                state, detail = "running", f"last log output {int(age)} s ago"
        issues = {}
        for path in self.order:
            for name, count in self.logs[path].issues.items():
                issues[name] = issues.get(name, 0) + count
        bad = sum(self.logs[p].bad_frames for p in self.order)
        if bad:
            issues["energy frames not understood (skipped)"] = bad
        # distance to the periodic image
        mindist = None
        if self.mindist:
            rows = self.mindist.series([x for x, _, _, _ in self.trajectory_segments()])
            cutoff = log.cutoff
            limit = cutoff * self.MINDIST_MARGIN if cutoff else None
            if cutoff and rows:
                under = [r for r in rows if r[1] < cutoff]
                close = [r for r in rows if cutoff <= r[1] < limit]
                if under:
                    issues[f"closer to its periodic image than the cutoff ({cutoff:g} nm), first at "
                           f"{under[0][0] / 1000:g} ns; frames"] = len(under)
                if close:
                    issues[f"closer to its periodic image than {self.MINDIST_MARGIN:g} x the cutoff "
                           f"({limit:.2f} nm), first at {close[0][0] / 1000:g} ns; frames"] = len(close)
            mi = decimate(len(rows), self.max_points)
            lowest = min(rows, key=lambda r: r[1]) if rows else None
            mindist = {
                "status": self.mindist.status, "error": self.mindist.error, "group": self.mindist.group,
                "cutoff": cutoff, "limit": limit,
                "lowest": lowest[1] if lowest else None, "lowest_ns": lowest[0] / 1000 if lowest else None,
                "t": [round(rows[i][0] / 1000, 6) for i in mi],
                "min": [_num(rows[i][1]) for i in mi], "maxint": [_num(rows[i][2]) for i in mi],
                "box1": [_num(rows[i][3]) for i in mi], "box2": [_num(rows[i][4]) for i in mi],
                "box3": [_num(rows[i][5]) for i in mi],
            }
            if self.mindist.error:
                issues[self.mindist.error] = None
        # RMSD
        rms = None
        if self.rms:
            rows = self.rms.series([x for x, _, _, _ in self.trajectory_segments()])
            ri = decimate(len(rows), self.max_points)
            highest = max(rows, key=lambda r: r[1]) if rows else None
            reference = self.rms.reference_used
            rms = {"status": self.rms.status, "error": self.rms.error, "group": self.rms.group,
                   "reference": reference.name if reference else None,
                   "highest": highest[1] if highest else None,
                   "highest_ns": highest[0] / 1000 if highest else None,
                   "t": [round(rows[i][0] / 1000, 6) for i in ri], "rmsd": [_num(rows[i][1]) for i in ri]}
            if self.rms.error:
                issues[self.rms.error] = None
        # conserved-energy drift of the current log
        drift_fit = linear_slope(log.times, log.cols["Conserved En."])
        out.update({
            "state": state, "detail": detail,
            "version": log.version, "host": log.host or out["host"], "natoms": log.natoms,
            "dt_fs": dt * 1000 if dt else None, "nstlog_ps": nstlog * dt if nstlog and dt else None,
            "ref_t": log.ref_t[0] if log.ref_t else None, "ref_p": log.ref_p,
            "pcoupl": log.params.get("pcoupl"), "tcoupl": log.params.get("tcoupl"),
            "restarts": max(0, len(log.starts) - 1),
            "time_ns": _ns(t_now), "step": log.steps[-1] if log.steps else None,
            "seg": {"index": len(self.order), "count": len(self.order), "name": log.path.name,
                    "start_ns": _ns(seg_start), "end_ns": _ns(seg_end),
                    "progress": _progress(t_now, seg_start, seg_end),
                    "left_s": seg_left, "left": fmt_duration(seg_left),
                    "finish": _clock(now, seg_left)},
            "total": {"end_ns": _ns(total[0]) if total else None, "source": total[1] if total else None,
                      "progress": _progress(t_now, 0.0 if self.layout == "segments" else
                                            (self.logs[self.order[0]].seg_start or 0.0),
                                            total[0] if total else None),
                      "left_s": total_left, "left": fmt_duration(total_left),
                      "finish": _clock(now, total_left)},
            "rate": {"now": rate_now, "avg": rate_avg},
            "log_age_s": age,
            "drift": {"fit_per_ns": drift_fit * 1000 if drift_fit is not None else None,
                      "fit_per_atom_ns": drift_fit * 1000 / log.natoms
                      if drift_fit is not None and log.natoms else None,
                      "gromacs_per_atom_ns": log.drift * 1000 if log.drift is not None else None},
            "issues": issues,
        })
        # series (whole run, thinned out)
        idx = decimate(len(times), self.max_points)
        p_avg = running_mean(cols["Pressure (bar)"], self.PRESSURE_AVG_FRAMES)
        pick = lambda values: [_num(values[i]) for i in idx]       # noqa: E731
        out["series"] = {
            "t": [round(times[i] / 1000, 6) for i in idx],
            "T": pick(cols["Temperature"]), "P": pick(cols["Pressure (bar)"]), "P_avg": pick(p_avg),
            "Epot": pick(cols["Potential"]), "Etot": pick(cols["Total Energy"]),
            "Econs": pick(cols["Conserved En."]), "constr": pick(cols["Constr. rmsd"]),
        }
        out["p_avg_ps"] = self.PRESSURE_AVG_FRAMES * nstlog * dt if nstlog and dt else None
        out["boundaries"] = [round(b / 1000, 6) for b in boundaries]
        out["current_start_ns"] = _ns(log.times[0]) if log.times else None
        speed = self.speed_points()
        out["speed"] = {"t": [round(t / 1000, 6) for t, _ in speed], "v": [_num(v) for _, v in speed]}
        if self.edr:
            data = self.edr.series(self.edr_paths())
            ei = decimate(len(data["time"]), self.max_points)
            out["edr"] = {"status": self.edr.status,
                          "t": [round(data["time"][i] / 1000, 6) for i in ei]}
            for term in EDR_TERMS:
                if term in data:
                    out["edr"][term] = [_num(data[term][i]) for i in ei]
        if mindist is not None:
            out["mindist"] = mindist
        if rms is not None:
            out["rms"] = rms
        return out


def _ns(ps):
    return None if ps is None else round(ps / 1000, 6)


def _num(v):
    return v if isinstance(v, (int, float)) and math.isfinite(v) else None


def _progress(t, start, end):
    if t is None or start is None or end is None or end <= start:
        return None
    return min(1.0, max(0.0, (t - start) / (end - start)))


def _clock(now, seconds):
    if seconds is None:
        return None
    when = datetime.fromtimestamp(now + seconds)
    fmt = "%H:%M" if when.date() == datetime.now().date() else "%a %d %b %H:%M"
    return when.strftime(fmt)


def text_status(snap: dict) -> str:
    """The snapshot as a few lines of text (--once)."""
    if not snap.get("seg"):
        return f"{snap['stem']} on {snap['host']}: {snap['state']} - {snap['detail']}"
    seg, total, rate = snap["seg"], snap["total"], snap["rate"]
    pct = lambda p: "-" if p is None else f"{100 * p:.1f}%"           # noqa: E731
    num = lambda v, f: "-" if v is None else format(v, f)             # noqa: E731
    lines = [
        f"{snap['stem']} on {snap['host']}, in {snap['dir']}  ({snap['layout']}, {len(snap['files'])} log file(s))",
        f"  state     : {snap['state']} - {snap['detail']}",
        f"  time      : {num(snap['time_ns'], '.3f')} ns   step {snap['step']}",
        f"  segment   : {seg['name']}  {num(seg['start_ns'], 'g')}-{num(seg['end_ns'], 'g')} ns  "
        f"{pct(seg['progress'])}  left {seg['left']}  (finish {seg['finish'] or '-'})",
        f"  whole run : end {num(total['end_ns'], 'g')} ns (from {total['source'] or '-'})  "
        f"{pct(total['progress'])}  left {total['left']}  (finish {total['finish'] or '-'})",
        f"  speed     : {num(rate['now'], '.1f')} ns/day now, {num(rate['avg'], '.1f')} average",
        f"  GROMACS {snap['version']} on {snap['host']}, {snap['natoms']} atoms, dt {num(snap['dt_fs'], 'g')} fs",
    ]
    rm = snap.get("rms")
    if rm:
        if rm["t"]:
            lines.append(f"  RMSD      : {rm['group']} {rm['rmsd'][-1]:.3f} nm now, highest {rm['highest']:.3f} nm "
                         f"at {rm['highest_ns']:g} ns (fit to {rm['reference']})")
        else:
            lines.append(f"  RMSD      : {rm['error'] or rm['status']}")
    md = snap.get("mindist")
    if md:
        if md["t"]:
            lines.append(f"  image     : {md['group']} {md['min'][-1]:.3f} nm from its periodic image now, lowest "
                         f"{md['lowest']:.3f} nm at {md['lowest_ns']:g} ns (cutoff {num(md['cutoff'], 'g')} nm, "
                         f"warning below {num(md['limit'], '.2f')} nm)")
        else:
            lines.append(f"  image     : {md['error'] or md['status']}")
    if snap["issues"]:
        lines.append("  issues    : " + ", ".join(k if v is None else f"{k}: {v}"
                                                   for k, v in snap["issues"].items()))
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Web page
# ---------------------------------------------------------------------------

PAGE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>gmx_monitor</title>
<style>
:root {
  --bg: #f6f7f9; --card: #ffffff; --text: #1d2330; --muted: #667085; --line: #e3e6eb;
  --grid: #eef0f3; --accent: #2f6fdf; --accent2: #d9480f; --accent3: #2b8a3e; --ref: #9aa1ad;
  --ok: #2b8a3e; --warn: #c77700; --bad: #c92a2a; --bar: #e8ecf2;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg: #12151b; --card: #1b2029; --text: #e6e9ef; --muted: #98a2b3; --line: #2a313d;
    --grid: #232a35; --accent: #6ea0ff; --accent2: #ff8a50; --accent3: #63c174; --ref: #6b7380;
    --ok: #63c174; --warn: #f0a33a; --bad: #ff6b6b; --bar: #2a313d;
  }
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--bg); color: var(--text);
       font: 14px/1.4 system-ui, -apple-system, "Segoe UI", Roboto, sans-serif; }
main { max-width: 1400px; margin: 0 auto; padding: 16px; }
h1 { font-size: 20px; margin: 0; }
h1 .host { font-weight: 400; color: var(--muted); }
.top { display: flex; align-items: center; gap: 12px; flex-wrap: wrap; margin-bottom: 12px; }
.pill { padding: 3px 10px; border-radius: 999px; font-weight: 600; font-size: 13px; color: #fff; }
.muted { color: var(--muted); }
.card { background: var(--card); border: 1px solid var(--line); border-radius: 10px; padding: 14px; }
.stats { display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 10px 18px; }
.stat .k { color: var(--muted); font-size: 12px; }
.stat .v { font-size: 18px; font-weight: 600; font-variant-numeric: tabular-nums; }
.stat .s { color: var(--muted); font-size: 12px; }
.bars { display: grid; gap: 10px; margin-top: 14px; }
.barrow { display: grid; grid-template-columns: 110px 1fr 70px; gap: 10px; align-items: center; }
.bar { height: 12px; background: var(--bar); border-radius: 6px; overflow: hidden; }
.bar > div { height: 100%; background: var(--accent); border-radius: 6px; transition: width .5s; }
.bar.total > div { background: var(--accent3); }
.issues { margin-top: 12px; color: var(--warn); }
.grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: 12px; margin-top: 12px; }
.chart h2 { font-size: 14px; margin: 0 0 6px; display: flex; justify-content: space-between; gap: 8px; }
.chart h2 span { font-weight: 400; color: var(--muted); font-size: 12px; }
canvas { width: 100%; height: 220px; display: block; }
footer { margin-top: 14px; color: var(--muted); font-size: 12px; }
@media (max-width: 520px) { .grid { grid-template-columns: 1fr; } .barrow { grid-template-columns: 80px 1fr 56px; } }
</style>
</head>
<body>
<main>
  <div class="top">
    <h1 id="title">gmx_monitor</h1>
    <span class="pill" id="state">…</span>
    <span class="muted" id="detail"></span>
  </div>
  <div class="card">
    <div class="stats" id="stats"></div>
    <div class="bars">
      <div class="barrow"><span class="muted" id="seglabel">Segment</span><div class="bar"><div id="segbar" style="width:0"></div></div><span id="segpct">-</span></div>
      <div class="barrow"><span class="muted">Whole run</span><div class="bar total"><div id="totbar" style="width:0"></div></div><span id="totpct">-</span></div>
    </div>
    <div class="issues" id="issues"></div>
  </div>
  <div class="grid" id="charts"></div>
  <footer id="footer"></footer>
</main>
<script>
const REFRESH = __REFRESH__;
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const fmt = (v, d) => (v === null || v === undefined || !isFinite(v)) ? "-" : Number(v).toFixed(d);
const pct = p => (p === null || p === undefined) ? "-" : (100 * p).toFixed(1) + "%";

function niceTicks(lo, hi, n) {
  if (!(isFinite(lo) && isFinite(hi))) return [];
  if (hi === lo) { hi = lo + 1; lo = lo - 1; }
  const raw = (hi - lo) / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => s >= raw) || 10 * mag;
  const ticks = [];
  for (let v = Math.ceil(lo / step) * step; v <= hi + step * 1e-9; v += step) ticks.push(+v.toPrecision(12));
  return ticks;
}
function tickLabel(v, span) {
  const a = Math.abs(v);
  if (a !== 0 && (a >= 1e6 || a < 1e-3)) return v.toExponential(2);
  const d = span >= 20 ? 0 : span >= 2 ? 1 : span >= 0.2 ? 2 : 3;
  return v.toFixed(d);
}

const charts = {};
function chart(id, title, unit) {
  if (charts[id]) return charts[id];
  const card = document.createElement("div");
  card.className = "card chart";
  card.innerHTML = `<h2>${title}<span id="${id}-info"></span></h2><canvas id="${id}"></canvas>`;
  document.getElementById("charts").appendChild(card);
  const c = { canvas: card.querySelector("canvas"), info: card.querySelector("span"), unit, data: null, hover: null };
  c.canvas.addEventListener("mousemove", e => { const r = c.canvas.getBoundingClientRect(); c.hover = e.clientX - r.left; draw(c); });
  c.canvas.addEventListener("mouseleave", () => { c.hover = null; draw(c); });
  charts[id] = c;
  return c;
}

function draw(c) {
  const d = c.data; if (!d) return;
  const cv = c.canvas, dpr = window.devicePixelRatio || 1;
  const W = cv.clientWidth, H = cv.clientHeight;
  if (cv.width !== Math.round(W * dpr)) { cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr); }
  const g = cv.getContext("2d");
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, W, H);
  const L = 62, R = 10, T = 10, B = 28;
  const x = d.x;
  if (!x.length) { g.fillStyle = css("--muted"); g.fillText("no data yet", L, T + 20); return; }
  let ylo = Infinity, yhi = -Infinity;
  for (const s of d.series) for (const v of s.y) if (v !== null) { if (v < ylo) ylo = v; if (v > yhi) yhi = v; }
  for (const r of (d.refs || [])) if (r.y !== null) { ylo = Math.min(ylo, r.y); yhi = Math.max(yhi, r.y); }
  if (!isFinite(ylo)) { g.fillStyle = css("--muted"); g.fillText("no data yet", L, T + 20); return; }
  const pad = (yhi - ylo) * 0.06 || Math.abs(yhi) * 0.01 || 1; ylo -= pad; yhi += pad;
  const xlo = x[0], xhi = x[x.length - 1] > x[0] ? x[x.length - 1] : x[0] + 1e-3;
  const X = v => L + (v - xlo) / (xhi - xlo) * (W - L - R);
  const Y = v => T + (yhi - v) / (yhi - ylo) * (H - T - B);
  g.font = "11px system-ui, sans-serif"; g.lineWidth = 1;
  // current segment shading, segment boundaries
  if (d.shadeFrom !== null && d.shadeFrom !== undefined && d.shadeFrom > xlo) {
    g.fillStyle = css("--grid"); g.fillRect(X(d.shadeFrom), T, W - R - X(d.shadeFrom), H - T - B);
  }
  // grid + ticks
  g.strokeStyle = css("--grid"); g.fillStyle = css("--muted");
  const yt = niceTicks(ylo, yhi, 4), xt = niceTicks(xlo, xhi, 6);
  g.textAlign = "right"; g.textBaseline = "middle";
  for (const v of yt) { g.beginPath(); g.moveTo(L, Y(v)); g.lineTo(W - R, Y(v)); g.stroke(); g.fillText(tickLabel(v, yhi - ylo), L - 6, Y(v)); }
  g.textAlign = "center"; g.textBaseline = "top";
  for (const v of xt) { g.beginPath(); g.moveTo(X(v), T); g.lineTo(X(v), H - B); g.stroke(); g.fillText(tickLabel(v, xhi - xlo), X(v), H - B + 6); }
  g.textAlign = "right"; g.fillText("ns", W - R, H - B + 6);
  g.setLineDash([3, 3]); g.strokeStyle = css("--ref");
  for (const b of (d.boundaries || [])) if (b > xlo && b < xhi) { g.beginPath(); g.moveTo(X(b), T); g.lineTo(X(b), H - B); g.stroke(); }
  for (const r of (d.refs || [])) if (r.y !== null && r.y !== undefined) {
    g.setLineDash(r.dash || [3, 3]); g.strokeStyle = css(r.color || "--ref"); g.lineWidth = r.width || 1;
    g.beginPath(); g.moveTo(L, Y(r.y)); g.lineTo(W - R, Y(r.y)); g.stroke();
  }
  g.setLineDash([]); g.lineWidth = 1;
  // series
  g.save(); g.beginPath(); g.rect(L, T, W - L - R, H - T - B); g.clip();
  for (const s of d.series) {
    g.strokeStyle = css(s.color); g.lineWidth = s.width || 1.2; g.globalAlpha = s.alpha || 1;
    g.beginPath(); let pen = false;
    for (let i = 0; i < x.length; i++) {
      const v = s.y[i];
      if (v === null) { pen = false; continue; }
      if (pen) g.lineTo(X(x[i]), Y(v)); else { g.moveTo(X(x[i]), Y(v)); pen = true; }
    }
    g.stroke();
    if (s.dots) { g.fillStyle = css(s.color); for (let i = 0; i < x.length; i++) if (s.y[i] !== null) { g.beginPath(); g.arc(X(x[i]), Y(s.y[i]), 2.5, 0, 7); g.fill(); } }
  }
  g.restore(); g.globalAlpha = 1;
  // legend
  g.textAlign = "left"; g.textBaseline = "middle"; let lx = L + 8;
  for (const s of d.series.filter(s => s.name)) {
    g.fillStyle = css(s.color); g.fillRect(lx, T + 8, 12, 3); g.fillStyle = css("--text");
    g.fillText(s.name, lx + 16, T + 10); lx += 26 + g.measureText(s.name).width;
  }
  for (const r of (d.refs || []).filter(r => r.name && r.y !== null && r.y !== undefined)) {
    g.strokeStyle = css(r.color || "--ref"); g.setLineDash(r.dash || [3, 3]); g.lineWidth = r.width || 1;
    g.beginPath(); g.moveTo(lx, T + 10); g.lineTo(lx + 14, T + 10); g.stroke(); g.setLineDash([]); g.lineWidth = 1;
    g.fillStyle = css("--text"); g.fillText(r.name, lx + 18, T + 10); lx += 28 + g.measureText(r.name).width;
  }
  // hover readout
  if (c.hover !== null && c.hover >= L && c.hover <= W - R) {
    const xv = xlo + (c.hover - L) / (W - L - R) * (xhi - xlo);
    let lo = 0, hi = x.length - 1;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (x[m] < xv) lo = m; else hi = m; }
    const i = Math.abs(x[lo] - xv) < Math.abs(x[hi] - xv) ? lo : hi;
    g.strokeStyle = css("--muted"); g.beginPath(); g.moveTo(X(x[i]), T); g.lineTo(X(x[i]), H - B); g.stroke();
    const parts = [`${x[i].toFixed(3)} ns`].concat(d.series.filter(s => s.name || d.series.length === 1)
      .map(s => `${s.name ? s.name + " " : ""}${s.y[i] === null ? "-" : (+s.y[i]).toPrecision(6)}`));
    c.info.textContent = parts.join("  ·  ") + (c.unit ? " " + c.unit : "");
  } else {
    c.info.textContent = d.info || "";
  }
}

function stat(k, v, s) { return `<div class="stat"><div class="k">${k}</div><div class="v">${v}</div><div class="s">${s || ""}</div></div>`; }

function render(s) {
  const colors = { "running": "--ok", "finished": "--accent", "between segments": "--accent", "starting": "--accent",
                   "stalled?": "--warn", "stopped": "--warn", "error": "--bad", "no log": "--bad" };
  const pill = document.getElementById("state");
  pill.textContent = s.state; pill.style.background = css(colors[s.state] || "--muted");
  const who = s.host ? `${s.stem} @ ${s.host}` : s.stem;
  const h1 = document.getElementById("title");
  h1.textContent = s.stem;
  if (s.host) { const sp = document.createElement("span"); sp.className = "host"; sp.textContent = " on " + s.host; h1.appendChild(sp); }
  document.getElementById("detail").textContent = s.detail || "";
  if (!s.seg) { document.title = `${who} · ${s.state}`; return; }
  const seg = s.seg, tot = s.total, r = s.rate;
  const totp = tot.progress !== null ? tot.progress : seg.progress;
  document.title = s.state === "finished" ? `✔ ${who} finished`
    : `${s.state === "running" ? "" : "⚠ "}${who} ${pct(totp)}${tot.left !== "-" ? " · " + tot.left + " left" : ""}`;
  document.getElementById("stats").innerHTML = [
    stat("Simulation time", fmt(s.time_ns, 3) + " ns", "step " + (s.step === null ? "-" : s.step)),
    stat("Speed", fmt(r.now, 1) + " ns/day", "average " + fmt(r.avg, 1) + " ns/day"),
    stat("Segment left", seg.left, seg.finish ? "finishes " + seg.finish : ""),
    stat("Whole run left", tot.left, tot.finish ? "finishes " + tot.finish : (tot.end_ns === null ? "total unknown: use --end" : "")),
    stat("Run end", tot.end_ns === null ? "-" : tot.end_ns + " ns", tot.source ? "from " + tot.source : ""),
    stat("Temperature", fmt(s.series.T[s.series.T.length - 1], 1) + " K", s.ref_t !== null ? "target " + s.ref_t + " K (" + (s.tcoupl || "") + ")" : ""),
    stat("Conserved-energy drift", fmt(s.drift.fit_per_atom_ns, 5), "kJ/mol/ns per atom" + (s.drift.gromacs_per_atom_ns !== null ? ` (GROMACS: ${fmt(s.drift.gromacs_per_atom_ns, 5)})` : "")),
    stat("System", (s.natoms === null ? "-" : s.natoms) + " atoms", `dt ${fmt(s.dt_fs, 1)} fs · GROMACS ${s.version || "-"}` + (s.restarts ? ` · ${s.restarts} restart(s)` : "")),
  ].join("");
  document.getElementById("seglabel").textContent = s.layout === "segments" ? `Segment ${seg.index}` : "This log";
  document.getElementById("segbar").style.width = (100 * (seg.progress || 0)) + "%";
  document.getElementById("segpct").textContent = pct(seg.progress);
  document.getElementById("totbar").style.width = (100 * (totp || 0)) + "%";
  document.getElementById("totpct").textContent = pct(totp);
  const iss = Object.entries(s.issues || {});
  document.getElementById("issues").textContent = iss.length ? "⚠ " + iss.map(([k, v]) => v === null ? k : `${k}: ${v}`).join(" · ") : "";
  const S = s.series, common = { x: S.t, boundaries: s.boundaries, shadeFrom: s.boundaries.length ? s.current_start_ns : null };
  const last = a => a.length ? a[a.length - 1] : null;
  const set = (id, title, unit, data) => { const c = chart(id, title, unit); c.data = Object.assign({}, common, data); draw(c); };
  set("T", "Temperature", "K", { series: [{ y: S.T, color: "--accent" }], refs: [{ y: s.ref_t }], info: `now ${fmt(last(S.T), 1)} K` });
  set("P", "Pressure", "bar", { series: [{ y: S.P, color: "--accent", alpha: 0.35, width: 1 }, { name: `running mean (${fmt(s.p_avg_ps, 0)} ps)`, y: S.P_avg, color: "--accent2", width: 2 }],
      refs: [{ y: s.ref_p }], info: `mean now ${fmt(last(S.P_avg), 1)} bar` });
  set("Epot", "Potential energy", "kJ/mol", { series: [{ y: S.Epot, color: "--accent" }], info: `now ${fmt(last(S.Epot), 0)} kJ/mol` });
  set("Econs", "Total and conserved energy", "kJ/mol", { series: [{ name: "total", y: S.Etot, color: "--accent" }, { name: "conserved", y: S.Econs, color: "--accent2", width: 2 }],
      info: s.drift.fit_per_ns !== null ? `drift ${fmt(s.drift.fit_per_ns, 1)} kJ/mol/ns (this segment)` : "" });
  set("constr", "Constraint rmsd", "", { series: [{ y: S.constr, color: "--accent3" }], info: `now ${last(S.constr) === null ? "-" : (+last(S.constr)).toExponential(2)}` });
  const sp = chart("speed", "Speed", "ns/day");
  sp.data = { x: s.speed.t, series: [{ y: s.speed.v, color: "--accent3", dots: true }], boundaries: s.boundaries, info: `now ${fmt(r.now, 1)} ns/day` }; draw(sp);
  if (s.rms) {
    const Q = s.rms, rc = chart("rms", `RMSD (${Q.group})`, "nm");
    rc.data = { x: Q.t, boundaries: s.boundaries, series: [{ y: Q.rmsd, color: "--accent", dots: Q.t.length <= 60 }],
      info: Q.error ? "⚠ " + Q.error
        : (Q.t.length ? `now ${fmt(last(Q.rmsd), 3)} nm · highest ${fmt(Q.highest, 3)} nm at ${fmt(Q.highest_ns, 1)} ns · fit to ${Q.reference}` : Q.status) };
    draw(rc);
  }
  if (s.mindist) {
    const M = s.mindist, mc = { x: M.t, boundaries: s.boundaries };
    const md = chart("mindist", `Distance to periodic image (${M.group})`, "nm");
    md.data = Object.assign({}, mc, {
      series: [{ y: M.min, color: "--accent", dots: M.t.length <= 60 }],
      refs: [{ y: M.cutoff, color: "--bad", dash: [6, 4], width: 1.5, name: `cutoff ${fmt(M.cutoff, 3)} nm` },
             { y: M.limit, color: "--warn", dash: [3, 3], width: 1.5, name: `1.2 × cutoff` }],
      info: M.error ? "⚠ " + M.error
        : (M.t.length ? `now ${fmt(last(M.min), 2)} nm · lowest ${fmt(M.lowest, 2)} nm at ${fmt(M.lowest_ns, 1)} ns` : M.status) });
    draw(md);
    const same = M.box1.every((v, i) => v === M.box2[i] && v === M.box3[i]);
    const boxes = same ? [{ name: "box", y: M.box1, color: "--accent" }]
      : [{ name: "box 1", y: M.box1, color: "--accent" }, { name: "box 2", y: M.box2, color: "--accent3" }, { name: "box 3", y: M.box3, color: "--ref" }];
    const sz = chart("size", `${M.group} size vs box`, "nm");
    sz.data = Object.assign({}, mc, { series: [{ name: "largest internal distance", y: M.maxint, color: "--accent2" }].concat(boxes),
      info: M.t.length ? `now ${fmt(last(M.maxint), 2)} of ${fmt(last(M.box1), 2)} nm` : "" });
    draw(sz);
  }
  if (s.edr) {
    const E = s.edr, ec = { x: E.t, boundaries: s.boundaries };
    const eset = (id, title, unit, series, info) => { const c = chart(id, title, unit); c.data = Object.assign({}, ec, { series, info }); draw(c); };
    if (E.Density) eset("density", "Density (.edr)", "kg/m³", [{ y: E.Density, color: "--accent" }], `now ${fmt(last(E.Density), 1)} kg/m³`);
    if (E.Volume) eset("volume", "Volume (.edr)", "nm³", [{ y: E.Volume, color: "--accent" }], `now ${fmt(last(E.Volume), 2)} nm³`);
    if (E["Box-X"]) eset("box", "Box size (.edr)", "nm", [{ name: "x", y: E["Box-X"], color: "--accent" }, { name: "y", y: E["Box-Y"], color: "--accent2" }, { name: "z", y: E["Box-Z"], color: "--accent3" }], "");
    if (!E.Density && !E.Volume && !E["Box-X"]) { const c = chart("edrnote", ".edr", ""); c.data = { x: [], series: [] }; c.info.textContent = E.status; draw(c); }
  }
  document.getElementById("footer").textContent =
    `${s.dir} · ${s.files.length} log file(s): ${s.files.map(f => f.name).join(", ")} · updated ${s.now} · refresh every ${REFRESH} s` +
    (s.edr ? ` · .edr: ${s.edr.status}` : "");
}

async function refresh() {
  try {
    const res = await fetch("data", { cache: "no-store" });
    render(await res.json());
  } catch (e) {
    const pill = document.getElementById("state");
    pill.textContent = "monitor not reachable"; pill.style.background = css("--bad");
  }
}
window.addEventListener("resize", () => Object.values(charts).forEach(draw));
refresh();
setInterval(refresh, REFRESH * 1000);
</script>
</body>
</html>
"""


def make_handler(run: Run, refresh: float):
    page = PAGE.replace("__REFRESH__", str(refresh)).encode()

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            path = self.path.split("?", 1)[0]
            if path in ("/", "/index.html"):
                self._send(200, "text/html; charset=utf-8", page)
            elif path == "/data":
                with run.lock:
                    run.update()
                    body = json.dumps(run.snapshot(), allow_nan=False).encode()
                self._send(200, "application/json", body)
            else:
                self._send(404, "text/plain", b"not found")

        def _send(self, code, ctype, body):
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Cache-Control", "no-store")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def log_message(self, *args):     # no line per request in the terminal
            pass

    return Handler


# ---------------------------------------------------------------------------
# Server start (as in gmx_panel.py)
# ---------------------------------------------------------------------------

_DEFAULT_PORT = 8090
_PORT_TRIES = 20   # without --port: try 8090, 8091, ... up to this many ports


def describe_port_user(port: int) -> str:
    """Return " by process <pid> (<command>)" for the program listening on *port*, or "" if
    it can't be found (ss shows the process only for your own programs)."""
    try:
        out = subprocess.run(["ss", "-H", "-ltnp", f"sport = :{port}"],
                             capture_output=True, text=True, timeout=5).stdout
    except (OSError, subprocess.SubprocessError):
        return ""
    match = re.search(r"pid=(\d+)", out)
    if not match:
        return ""
    pid = match.group(1)
    try:
        argv = Path(f"/proc/{pid}/cmdline").read_bytes().decode(errors="replace").split("\0")
        cmd = " ".join(os.path.basename(a) if a.startswith("/") else a for a in argv if a)
    except OSError:
        cmd = ""
    if len(cmd) > 100:
        cmd = cmd[:97] + "..."
    return f" by process {pid}" + (f" ({cmd})" if cmd else "")


def start_server(host: str, port: int | None, handler) -> ThreadingHTTPServer:
    explicit = port is not None
    first = port if explicit else _DEFAULT_PORT
    for candidate in ([first] if explicit else range(first, first + _PORT_TRIES)):
        try:
            return ThreadingHTTPServer((host, candidate), handler)
        except OSError as exc:
            if exc.errno != errno.EADDRINUSE:
                raise
            holder = describe_port_user(candidate)
            if explicit:
                sys.exit(f"Error: port {candidate} is already in use{holder}.\n"
                         f"Stop that program, or choose another port with --port.")
            print(f"Port {candidate} is already in use{holder}; trying {candidate + 1}.")
    sys.exit(f"Error: ports {first}-{first + _PORT_TRIES - 1} are all in use; choose one with --port.")


def main():
    parser = argparse.ArgumentParser(
        prog="gmx_monitor",
        description="Live web panel for a running GROMACS MD simulation, from its mdrun log files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
The stem is <basename>_<replica>, as given to gmx_continue_grompp.sh / gmx_continue_extend.sh:
  segments  <stem>_<start>_<end>.log            (gmx_continue_grompp.sh)
  parts     <stem>.log, <stem>.partNNNN.log     (gmx_continue_extend.sh, or one mdrun -deffnm <stem>)
The whole-run length comes from --end, else from the invocation recorded in
gmx_continue_*_commands*.sh in the same folder, else (parts) from nsteps in the log.

Examples:
  %(prog)s md_0                          # files in the current folder
  %(prog)s MnMT4_apo_0 --dir ~/projects/MnMT4/apo/run0
  %(prog)s md_0 --end 1000               # the whole run is 1000 ns
  %(prog)s md_0 --edr                    # also box, volume and density (gmx energy)
  %(prog)s md_0 --mindist                # also distance to the periodic image (gmx mindist -pi)
  %(prog)s md_0 --rms                    # also backbone RMSD from the start of the run (gmx rms)
  %(prog)s md_0 --once                   # print the status in the terminal and exit
""")
    parser.add_argument("stem", help="<basename>_<replica> of the run, e.g. md_0")
    parser.add_argument("--dir", "-d", default=".", help="folder with the run's files (default: .)")
    parser.add_argument("--end", type=float, default=None,
                        help="length of the whole run in ns (default: detected, see below)")
    parser.add_argument("--refresh", type=float, default=10,
                        help="page refresh interval in seconds (default: 10)")
    parser.add_argument("--edr", action="store_true",
                        help="also plot box size, volume and density from the .edr files (gmx energy)")
    parser.add_argument("--edr-interval", type=float, default=60,
                        help="seconds between .edr reads (default: 60)")
    parser.add_argument("--mindist", action="store_true",
                        help="also plot the distance to the periodic image and the size vs the box "
                             "(gmx mindist -pi on the trajectories)")
    parser.add_argument("--mindist-dt", type=float, default=1.0,
                        help="--mindist: one frame per this many ns (default: 1)")
    parser.add_argument("--mindist-group", default="Protein",
                        help="--mindist: index group (default: Protein)")
    parser.add_argument("--mindist-index", default=None,
                        help="--mindist: index file with the group (default: the .tpr's default groups)")
    parser.add_argument("--rms", action="store_true",
                        help="also plot the RMSD (gmx rms on the trajectories)")
    parser.add_argument("--rms-dt", type=float, default=1.0,
                        help="--rms: one frame per this many ns (default: 1)")
    parser.add_argument("--rms-group", default="Backbone",
                        help="--rms: group fitted and measured (default: Backbone)")
    parser.add_argument("--rms-ref", default=None,
                        help="--rms: reference structure (.tpr, .gro or .pdb with the trajectory's "
                             "atoms; default: the .tpr of the run's first segment)")
    parser.add_argument("--rms-index", default=None,
                        help="--rms: index file with the group (default: the default groups)")
    parser.add_argument("--gmx", default="gmx",
                        help="GROMACS executable for --edr, --mindist and --rms (default: gmx)")
    parser.add_argument("--max-points", type=int, default=2000,
                        help="points per graph; longer runs are thinned out (default: 2000)")
    parser.add_argument("--port", "-p", type=int, default=None,
                        help=f"port (default: first free from {_DEFAULT_PORT})")
    parser.add_argument("--host", default="0.0.0.0",
                        help="interface to listen on (default: 0.0.0.0, all interfaces)")
    parser.add_argument("--once", action="store_true",
                        help="print the status in the terminal and exit (no web server)")
    args = parser.parse_args()
    sys.stdout.reconfigure(line_buffering=True)    # messages show up at once also when redirected

    if "/" in args.stem:
        parser.error("the stem is a name, not a path: give the folder with --dir")
    directory = Path(args.dir).expanduser()
    if not directory.is_dir():
        parser.error(f"folder not found: {directory}")
    if args.refresh < 1:
        parser.error("--refresh must be at least 1 s")

    if args.mindist_dt <= 0 or args.rms_dt <= 0:
        parser.error("--mindist-dt and --rms-dt must be positive")

    def existing(option, value):
        path = Path(value).expanduser() if value else None
        if path is not None and not path.is_file():
            parser.error(f"{option}: file not found: {path}")
        return path

    edr = mindist = rms = None
    if (args.edr or args.mindist or args.rms) and not shutil.which(args.gmx):
        print(f"Warning: '{args.gmx}' not found: --edr/--mindist/--rms disabled (source GMXRC, or use "
              f"--gmx).", file=sys.stderr)
    else:
        if args.rms:
            rms = RmsReader(args.gmx, args.rms_group, args.rms_dt * 1000,
                            existing("--rms-index", args.rms_index), existing("--rms-ref", args.rms_ref))
        if args.edr:
            edr = EdrReader(args.gmx, args.edr_interval)
        if args.mindist:
            mindist = MindistReader(args.gmx, args.mindist_group, args.mindist_dt * 1000,
                                    existing("--mindist-index", args.mindist_index))
    run = Run(directory, args.stem, args.end, edr, args.max_points, mindist, rms)
    run.update()

    if args.once:
        if edr:     # wait for the first .edr read, so that --once can report it
            edr._last = 0
            edr.maybe_refresh(run.edr_paths())
            for _ in range(600):
                time.sleep(0.5)
                if not edr._running:
                    break
        for reader, what in ((mindist, "the distance to the periodic image (gmx mindist -pi)"),
                             (rms, "the RMSD (gmx rms)")):
            if not reader:   # wait for the trajectory analyses (a first pass can take minutes)
                continue
            print(f"Computing {what}, one frame per {reader.dt_ps / 1000:g} ns...")
            while True:
                run.update()
                time.sleep(0.5)
                while reader.running:
                    time.sleep(0.5)
                if reader.error or not any(
                        trajectory_due(reader, seg) for seg in run.trajectory_segments()):
                    break
        snap = run.snapshot()
        print(text_status(snap))
        if edr:
            terms = [t for t in EDR_TERMS if t in snap.get("edr", {})]
            print(f"  .edr      : {edr.status}; terms: {', '.join(terms) or 'none'}")
        return

    if not run.order:
        print(f"Note: no mdrun log for '{args.stem}' in {directory} yet; the panel will show it "
              f"when one appears.")
    server = start_server(args.host, args.port, make_handler(run, args.refresh))
    port = server.server_address[1]
    print(f"Monitoring {args.stem} in {directory.resolve()}")
    print("Panel available at:")
    print(f"  http://localhost:{port}/")
    if args.host in ("0.0.0.0", ""):
        print(f"  http://{socket.gethostname()}:{port}/")
    print(f"(refreshes every {args.refresh:g} s; Ctrl-C to stop)")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.")



if __name__ == "__main__":
    main()

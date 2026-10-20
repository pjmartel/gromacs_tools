# gmx_monitor.py

Live web panel for a running GROMACS MD simulation.

## Overview

Follows the mdrun log files of one run and serves a page that refreshes itself: the state of
the run, the progress of the current segment and of the whole run, the speed, the time left and
the finish time, and graphs of temperature, pressure, energies and speed over the whole
simulation. It understands the file layouts of [gmx_continue_grompp.sh](md_continuation.md)
(one log per segment) and [gmx_continue_extend.sh](gmx_continue_extend.md) (one log that is
appended to, or `.partNNNN` logs), and also a single `mdrun -deffnm`.

Start it next to the run, on the machine where it is running, and open the page in a browser
from anywhere on your network (e.g. over Tailscale).

## Features

- ✅ **Status at a glance**: running / finished / between segments / stopped / **stalled?** /
  **error**, also shown in the browser tab title (`md_0 47.3% · 3 h 12 min left`)
- ✅ **Progress and time left**: bars for the segment and the whole run, speed now (last 15 min)
  and on average, time left and finish time for the segment and the whole run
- ✅ **Graphs of the whole run** (segment boundaries marked, current segment shaded):
  temperature with its target, pressure with a running mean and its target, potential energy,
  total and conserved energy, constraint rmsd, speed (ns/day) over the simulation
- ✅ **Conserved-energy drift** of the current segment (kJ/mol/ns per atom), next to the value
  GROMACS reports when a segment finishes; the best early warning of an unstable setup
- ✅ **Problems in the log**: LINCS/SETTLE warnings, atoms that moved too far, NaN energies, and the
  message of a fatal error
- ✅ **Box, volume and density** from the `.edr` files with `--edr` (through `gmx energy`)
- ✅ **Incremental**: only the new part of each log is read on a refresh, so it follows long runs
  at no cost; new segments are picked up as they appear
- ✅ **No dependencies**: Python standard library only (any `python3` ≥ 3.8, no packages, no uv);
  the page draws its graphs itself and needs no internet access
- ✅ **Terminal check**: `--once` prints the status and exits, e.g. over ssh

## Usage

### Basic Syntax

```bash
gmx_monitor.py <stem> [--dir DIR] [--end NS] [--edr] [--refresh S] [--port PORT] [--once]
```

`<stem>` is `<basename>_<replica>`, as given to the continuation scripts (`md 0` → `md_0`). Then
open `http://localhost:8090/` (or the host-name URL it prints) in a browser. Stop it with
`Ctrl+C`; stopping the monitor never affects the simulation.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `<stem>` | required | `<basename>_<replica>` of the run, e.g. `md_0` |
| `-d, --dir <dir>` | `.` | Folder with the run's files |
| `--end <ns>` | detected | Length of the whole run in ns (see [Whole-Run Length](#whole-run-length)) |
| `--refresh <s>` | `10` | Page refresh interval in seconds |
| `--edr` | off | Also plot box size, volume and density from the `.edr` files (see [Box, Volume and Density](#box-volume-and-density---edr)) |
| `--edr-interval <s>` | `60` | Seconds between `.edr` reads |
| `--gmx <exe>` | `gmx` | GROMACS executable for `--edr` |
| `--max-points <n>` | `2000` | Points per graph; longer runs are thinned out evenly |
| `-p, --port <port>` | `8090`, or the next free port | Port of the panel. Without `--port`, a port in use is skipped (8091, 8092, …); an explicitly given port is used, or the monitor stops with a message saying which program holds it |
| `--host <addr>` | `0.0.0.0` | Interface to listen on (`0.0.0.0` = all interfaces) |
| `--once` | off | Print the status in the terminal and exit (no web server) |
| `-h, --help` | - | Show help and exit |

## How It Works

### Which files

| Layout | Files | Made by |
|--------|-------|---------|
| segments | `<stem>_<start>_<end>.log`, in order of `<start>` | `gmx_continue_grompp.sh` |
| parts | `<stem>.log`, `<stem>.part0002.log`, … | `gmx_continue_extend.sh` (append or `-noappend`), `mdrun -deffnm <stem>` |

Only real mdrun logs are read: a file must start with the `GROMACS - gmx mdrun` header. The
combined output log that `gmx_continue_grompp.sh` writes as `<stem>.log` (it starts with
`gmx grompp`) is therefore ignored. With `--edr`, each log's `.edr` has the same name.

### What is read from the log

The log has an energy frame every `nstlog` steps (`Step`/`Time`, then the `Energies (kJ/mol)`
block), the run parameters (`nsteps`, `dt`, `tinit`, `init-step`, `nstlog`, `ref-t`, `ref-p`),
the `Started mdrun` / `Finished mdrun` times, the checkpoint lines (`Writing checkpoint, step …
at …`) and, at the end, the performance and the conserved-energy drift. Parsing is strict: an
energy frame whose names and values don't match is skipped and counted (shown as a problem),
never guessed. A frame cut off by a crash is dropped silently. After a restart from a checkpoint
with `-append`, the frames after the checkpoint appear twice in the log; the newer ones replace
the older ones.

### Speed and time left

The log has no clock time per frame, so the speed is measured:

- **average**: simulation time done since the last `Started mdrun`, divided by the wall time since
  then (up to the last write of the log)
- **now**: over the last 15 minutes, from the checkpoint lines (written every 15 min by default)
  and from what the monitor saw itself on each refresh
- a finished segment shows the speed GROMACS reports (`Performance:`)

Time left = simulation time left / current speed, for the segment and for the whole run (the
latter assumes the next segments run at the same speed). The speed graph uses the checkpoint
lines of every log, so it has a history even if the monitor was started late.

### States

| State | Meaning |
|-------|---------|
| running | the log is being written |
| stalled? | no log output for more than 10 times the usual interval between frames (at least 5 min, at most 1 h): mdrun may have hung, or the machine is overloaded |
| between segments | a segment finished less than 10 min ago and the whole run isn't done: the next segment should start soon |
| stopped | a segment finished more than 10 min ago and no next one started; or mdrun finished before the end of its run (`-maxh`, `-nsteps`) |
| finished | the whole run is done |
| error | the log has a `Fatal error` (its message is shown) |
| starting | mdrun has started, but there is no energy frame yet |

### Whole-run length

The progress of the whole run needs its total length, taken from (first found):

1. `--end <ns>`
2. the invocation recorded in `gmx_continue_grompp_commands*.sh` / `gmx_continue_extend_commands*.sh`
   in the same folder (`# Invocation: … <basename> <replica> <start> <end> <dt> …`, ps if
   `--ps` was used): the newest file whose `<basename>_<replica>` is the stem
3. layout *parts* only: the end of the log's own run (`nsteps`, or `mdrun -nsteps`)

When the length is unknown, the panel says so ("total unknown: use --end").

After an extension with `convert-tpr` and a restart with `-append`, mdrun doesn't print the
parameters again, so the log's `nsteps` is out of date once the run is past it: the end then
comes only from `--end` or the commands file.

### Box, Volume and Density (`--edr`)

These are not in the log. With `--edr`, `gmx energy` reads the `Box-X`, `Box-Y`, `Box-Z`,
`Volume` and `Density` terms of each segment's `.edr`, in the background (the panel never waits
for it), at most every `--edr-interval` seconds, and only for files that changed. `gmx` must be
on `PATH` (source `GMXRC`) or given with `--gmx`. A constant-volume run has none of these terms,
and the panel says so. `gmx energy` runs in a temporary folder, so it leaves no backup files
(`#…#`) next to the run.

## Remote Machines

By default the monitor listens on all interfaces, and on startup it prints both a `localhost`
URL and one using the machine's hostname (e.g. `http://marvin:8090/` over Tailscale).

The panel has no authentication. On a shared or public machine, bind to localhost only and
reach it through an SSH tunnel:

```bash
# On the remote machine
gmx_monitor.py md_0 --host 127.0.0.1 --port 8090

# On your local machine, then open http://localhost:8090/
ssh -L 8090:localhost:8090 user@remote-host
```

To keep it running after you log out, start it in `screen` (or with `nohup … &`).

## Examples

```bash
# Run in the current folder (segments md_0_0_100.log, md_0_100_200.log, ...)
gmx_monitor.py md_0

# Run in another folder, with box, volume and density
gmx_monitor.py MnMT4_apo_0 --dir ~/projects/MnMT4/apo/run0 --edr

# The whole run is 1000 ns (no commands file in the folder)
gmx_monitor.py md_0 --end 1000

# Quick check over ssh, without the web panel
ssh marvin 'gmx_monitor.py md_0 --dir ~/projects/MnMT4/apo/run0 --once'
```

`--once` prints, for example (a finished 2 ns test run):

```
live_0 in /home/martel/projects/TRP_cage/run1  (parts, 1 log file(s))
  state     : finished - finished 2026-10-20 17:51
  time      : 2.000 ns   step 1000000
  segment   : live_0.log  0-2 ns  100.0%  left -  (finish -)
  whole run : end 2 ns (from the log (nsteps))  100.0%  left -  (finish -)
  speed     : 320.6 ns/day now, 320.6 average
  GROMACS 2023.3 on quark, 8339 atoms, dt 2 fs
```

## Troubleshooting

### "no mdrun log for '<stem>'"

Check the stem (`<basename>_<replica>`, without `_<start>_<end>` and without `.log`) and the
folder (`--dir`). The monitor keeps running and shows the run as soon as a log appears.

### "total unknown: use --end"

There is no commands file of the continuation scripts in the folder (or the run was extended with
`-append`, see [Whole-Run Length](#whole-run-length)): give the length with `--end`.

### Speed shows "-"

The speed needs at least 30 s of run (average) or 60 s of observations (now) after a start.

### `--edr` shows no graphs

`gmx` was not found (a warning is printed at startup: source `GMXRC` or use `--gmx`), the run is at
constant volume, or the first read is still running (the status is in the page footer).

### Port already in use

Without `--port`, the next free port is taken. With an explicit `--port`, the monitor stops and
says which program holds the port.

# gmx_continue_extend.sh

Continue GROMACS MD simulations in segments by extending a single TPR file.

## Overview

An alternative to [gmx_continue_grompp.sh](md_continuation.md). Instead of running `grompp`
to build a new `.tpr` for every segment, this script creates **one** `.tpr`, then for each
segment extends it with `gmx convert-tpr -until` and continues from the checkpoint with
`gmx mdrun -cpi`. This is the continuation method recommended by the GROMACS manual. It
gives an exact continuation from the checkpoint and needs no MDP file after the first
segment.

| | `gmx_continue_grompp.sh` | `gmx_continue_extend.sh` |
|---|---|---|
| New `.tpr` per segment | Yes (`grompp`) | No (`convert-tpr -until`) |
| MDP changes between segments | Possible | Not possible (fixed in the `.tpr`) |
| Output file names | `<base>_<rep>_<start>_<end>.*` | `<base>_<rep>.*` (+ `.part000N.*` in noappend mode) |
| Works with [gmx_analysis.sh](gmx_analysis.md) | Yes | No (different file naming) |

## Features

- ✅ **Easy mode**: `gmx_continue_extend.sh 5000` extends the simulation in the current directory by 5000 ps
- ✅ **Segmented runs** from equilibration files, or from an existing `.tpr`/`.cpt`
- ✅ **Append or noappend output**: one continuous set of files, or a new `part000N` file set per segment
- ✅ **Crash recovery**: re-run the same command; segments the checkpoint has already reached are skipped,
  and an interrupted segment is completed without overshooting its end time
- ✅ **STOP file**: create `./STOP` to halt cleanly after the current segment
- ✅ **PLUMED support** via `--plumed`
- ✅ **Command logging** (on by default): every gmx command run, plus the exact invocation used,
  is saved to `gmx_continue_extend_commands.sh` (override with `--commands-file`); the file from an earlier run is kept as `…_commands.1.sh`, `.2.sh`, … (highest number = most recent)

## Usage

### Easy Mode

```bash
gmx_continue_extend.sh <extend_time_ps> [tpr_file]
```

Extends a single simulation in the current directory by `extend_time_ps` **picoseconds**
(an integer):

1. Uses `tpr_file` if given. Otherwise it auto-detects the only `*.tpr` in the current
   directory, and exits with a list of candidates if there is more than one.
2. Runs `gmx convert-tpr -s <tpr> -o <tpr> -extend <time>`, which overwrites the `.tpr` in place.
3. Runs `gmx mdrun -deffnm <base> -cpi <base>.cpt -append`, or without `-cpi` if no checkpoint exists.

Easy mode always appends to the existing output files and takes no options. Each call adds
`extend_time_ps` again, so after a crash resume with `gmx mdrun -deffnm <base> -cpi <base>.cpt -append`
rather than re-running easy mode; for crash-safe multi-segment runs, use normal mode.

### Normal Mode

```bash
gmx_continue_extend.sh <basename> <replica> <start_time> <end_time> <dt> [OPTIONS]
```

| Argument | Description |
|----------|-------------|
| `basename` | Base name for output files (e.g. `md`) |
| `replica` | Replica number (integer, e.g. `0`) |
| `start_time` | Start time, in ns (or ps with `--ps`); integer |
| `end_time` | End time, in ns (or ps with `--ps`); integer |
| `dt` | Segment length, in ns (or ps with `--ps`); integer |

All five positional arguments must come before any option. Output files are named
`<basename>_<replica>.*` (e.g. `md_0.xtc`).

### Options (Normal Mode)

| Option | Default | Description |
|--------|---------|-------------|
| `--template <file>` | `md.mdp` | MDP template used to build the initial `.tpr` (`tinit` and `nsteps` are set automatically) |
| `--initial <basename>` | `npt` | Basename of the equilibration files used as the starting point (`.gro`, `.edr` required; `.cpt` used if present) |
| `--timestep <ps>` | `0.002` | Integration timestep, used to compute `nsteps` per segment |
| `--title <suffix>` | - | Suffix appended to the `title` line of the MDP |
| `--append` | - | Continue writing to the same `.xtc`/`.edr`/`.log` files |
| `--noappend` | *(default)* | Write each continuation segment to new `<base>.part000N.*` files |
| `--tpr <file>` | - | Start from an existing `.tpr` instead of creating one (skips `grompp`) |
| `--cpt <file>` | - | Checkpoint to use with `--tpr` |
| `--topology <file>` | `topol.top` | Topology file (not needed with `--tpr`) |
| `--plumed <file>` | - | PLUMED input file, passed to `mdrun -plumed` |
| `--ps` | off | Interpret `start_time`/`end_time`/`dt` as picoseconds instead of nanoseconds |
| `--commands-file <file>` | `gmx_continue_extend_commands.sh` | Reproducibility script recording every gmx command run |

## Workflow

### First segment (`start_time = 0`, no `--tpr`)

1. Copy `--template` to `<base>_<rep>_init.mdp`, setting `tinit = 0` and `nsteps` to one segment.
2. `gmx grompp -f <base>_<rep>_init.mdp -c <initial>.gro [-t <initial>.cpt] -e <initial>.edr -p <topology> -o <base>_<rep>.tpr`
3. `gmx mdrun -deffnm <base>_<rep>`

### Starting from an existing `.tpr` (`--tpr`)

The base name is taken from the `.tpr` file name (e.g. `--tpr md_0.tpr` → `md_0`), so that
append mode keeps writing to the original files. The checkpoint stores their names and
checksums. `grompp` is not run, and `--template`/`--initial`/`--topology` are ignored.

### Each following segment

1. `gmx convert-tpr -s <base>.tpr -o <base>.tpr -until <segment end in ps>` (`-until` sets an absolute end time, so
   repeating it after a crash does not extend the run twice)
2. `gmx mdrun -deffnm <base> -cpi <base>.cpt -s <base>.tpr -append|-noappend`

The checkpoint is always `<base>.cpt`. The `part000N` suffix only applies to
trajectory/energy/log files.

> [!NOTE]
> With `start_time > 0` and no `--tpr`, the script goes straight to extending
> `<basename>_<replica>.tpr`, so that file (and its `.cpt`) must already exist in the current directory.

## Times and Crash Recovery

`start_time`/`end_time` are **absolute simulation times**, matching the time in the
checkpoint (e.g. to add 500 ns to a run that is at 1000 ns: `1000 1500 100`).

A segment counts as complete when the time stored in `<base>.cpt` (read with
`gmx check -f <base>.cpt`) has reached the segment's end time. This is the same in both output
modes. It reflects how far the simulation actually got, so re-running the same command:

- skips every segment the checkpoint has already reached,
- finishes an interrupted segment from the last checkpoint, up to its original end time (the
  `-until` end time is recomputed, not added again),
- does nothing at all if the whole range is complete.

## Append vs. Noappend

| Mode | Output |
|------|--------|
| `--noappend` (default) | `md_0.xtc` (segment 1), `md_0.part0002.xtc` (segment 2), `md_0.part0003.xtc`, ... (same for `.edr`/`.log`) |
| `--append` | `md_0.xtc`, `md_0.edr`, `md_0.log`, extended continuously |

GROMACS counts the first run as simulation part 1, so the first continuation is
`part0002`. To join noappend output into one trajectory:

```bash
gmx trjcat -f md_0.xtc md_0.part*.xtc -o md_0_complete.xtc -cat
```

## Examples

### Easy Mode

```bash
# Extend the simulation in the current directory by 5 ns
bash bin/gmx_continue_extend.sh 5000

# Several .tpr files present: say which one
bash bin/gmx_continue_extend.sh 5000 md_0.tpr
```

### New Production Run

```bash
# 500 ns in 100 ns segments from npt.gro/cpt/edr (noappend, part files)
bash bin/gmx_continue_extend.sh md 0 0 500 100 --template production.mdp

# Same, single continuous output files
bash bin/gmx_continue_extend.sh md 0 0 500 100 --template production.mdp --append

# All options
bash bin/gmx_continue_extend.sh md 0 0 500 100 --template production.mdp \
    --initial npt --timestep 0.002 --title "MyProtein" --noappend
```

### Continue an Existing Run

```bash
bash bin/gmx_continue_extend.sh md 0 0 500 100 --tpr md_0.tpr --cpt md_0.cpt --append
```

### Times in Picoseconds

```bash
# 10 ns in 2 ns segments
bash bin/gmx_continue_extend.sh md 0 0 10000 2000 --template production.mdp --ps
```

### Crash Recovery and Stopping

```bash
# After a crash: re-run the exact same command; completed segments are skipped
bash bin/gmx_continue_extend.sh md 0 0 500 100 --template production.mdp

# Stop cleanly after the running segment finishes
touch STOP
```

## Troubleshooting

### "Found optional flag ... in required argument position"
Normal mode needs all five positional arguments before any `--option`.

### "Multiple TPR files found in current directory"
In easy mode, pass the `.tpr` explicitly as the second argument.

### "Template MDP file / Initial structure / Initial energy file not found"
The first segment needs `--template` (default `md.mdp`) and `<initial>.gro` +
`<initial>.edr` (default `npt.*`) in the current directory.

### mdrun refuses to append
In append mode the checkpoint must match the existing output files (names and checksums).
Don't rename or edit `<base>.xtc`/`.edr`/`.log` between segments; use `--noappend` if
you need separate files.

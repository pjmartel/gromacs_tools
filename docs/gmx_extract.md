# gmx_extract.sh

Extract, clean up (PBC/center/fit) and convert a GROMACS trajectory for visualization.

## Overview

Wraps a sequence of `gmx trjconv` (and optionally `gmx convert-tpr`) calls into a single
repeatable pipeline that takes a raw trajectory to something ready for viewing in VMD,
PyMOL, or similar tools: whole molecules, centered in the box, fitted to a reference, and
optionally written out as a multi-model PDB.

## Features

- ✅ **Fixed, sensible pipeline order**: extract → pbc → center → fit → pdb
- ✅ **Selectable subset of steps** via `--steps`
- ✅ **Trajectory slicing**: `-b`/`-e`/`--dt`/`--skip`, applied once to the first stage that reads the raw input
- ✅ **Independent group selection** for extraction, centering, and fitting
- ✅ **Non-contiguous group support** (e.g. Protein + Ion) via `--subset-tpr`
- ✅ **Dry-run mode** to preview commands without running them
- ✅ **Script export** (`--save-script`) to save the generated commands as a runnable script

## Pipeline Stages

Stages run in this fixed order; use `--steps` to select a subset.

| Step | Description |
|------|-------------|
| `extract` | Plain group extraction (optional `-b`/`-e`/`--dt`/`--skip` slice), no PBC treatment |
| `pbc` | Group extraction + `-pbc whole` (makes the selected group whole across PBC) |
| `center` | Center the group in the box, keep molecules whole (`-pbc mol -center -ur ...`) |
| `fit` | Least-squares fit of the trajectory onto the reference structure |
| `pdb` | Write the final frame set as a multi-model PDB trajectory |

Default pipeline: `pbc,center,fit,pdb`.

## Usage

### Basic Syntax

```bash
gmx_extract.sh -s <tpr> -f <xtc> [OPTIONS]
```

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `-s, --tpr <file>` | *(required)* | Run input file (.tpr) |
| `-f, --xtc <file>` | *(required)* | Trajectory file (.xtc/.trr/...) |
| `-g, --group <name>` | `Protein` | Group to extract/output |
| `--center-group <name>` | same as `--group` | Group used for centering |
| `--fit-group <name>` | same as `--group` | Group used for the least-squares fit |
| `-b, --start <ps>` | - | Start time (ps) |
| `-e, --end <ps>` | - | End time (ps) |
| `--dt <ps>` | - | Only write frames at this time interval (ps) |
| `--skip <n>` | - | Only write every nth frame |
| `--steps <list>` | `pbc,center,fit,pdb` | Comma-separated subset of: `extract,pbc,center,fit,pdb` |
| `--fit-mode <mode>` | `progressive` | `gmx trjconv -fit` value |
| `--ur <rect\|tric\|compact>` | `compact` | `gmx trjconv -ur` value for the center step |
| `--subset-tpr` | off | Build a reduced `.tpr` for non-contiguous groups (see below) |
| `-o, --output-basename <n>` | derived from `--xtc` | Basename for output files |
| `--outdir <dir>` | `.` | Output directory |
| `-n, --index <file>` | - | Optional index (.ndx) file passed to every trjconv call |
| `--gmx <path>` | `gmx` | gmx binary/command to use |
| `--dry-run` | off | Print commands (and group selections) without running them |
| `--save-script <file>` | - | Also write the commands to a runnable bash script |
| `-h, --help` | - | Show help and exit |

### Available Default Groups

From a typical protein `.tpr`, no index file required:

`System, Protein, Protein-H, C-alpha, Backbone, MainChain, MainChain+Cb, MainChain+H,
SideChain, SideChain-H, Prot-Masses, non-Protein, Water, SOL, non-Water, Ion, Water_and_ions`

## Non-Contiguous Groups (`--subset-tpr`)

If the extract group is non-contiguous (e.g. `Protein_Ion`, a custom group combining
protein and ions), the atom indices in the original `.tpr` no longer match the extracted
`.xtc` once that group is pulled into its own trajectory. Passing `--subset-tpr` builds a
reduced `.tpr` (via `gmx convert-tpr`) containing only the extracted atoms, and switches
subsequent center/fit/pdb steps to use it with output group `System`.

`--center-group`/`--fit-group` default to the extract group, or to `System` once
`--subset-tpr` swaps in the reduced `.tpr` (the extract group itself, e.g. `non-Water`, may
no longer exist as a named group there). Pass `--center-group`/`--fit-group` explicitly
(e.g. `C-alpha`) to use a group that still exists in the reduced `.tpr` instead.

Contiguous selections (e.g. plain `Protein`) don't need `--subset-tpr`, and it is off by
default.

## Examples

### Full Default Pipeline

```bash
# pbc -> center -> fit -> pdb on the Protein group
bash bin/gmx_extract.sh -s md.tpr -f md.xtc
```

### Extract a Time Slice Only

```bash
# Only extract a 0-100 ns slice of the Protein group, no PBC/center/fit/pdb
bash bin/gmx_extract.sh -s md.tpr -f md.xtc --steps extract -b 0 -e 100000
```

### PBC Cleanup Only

```bash
# Extract + remove PBC only (no centering/fit/pdb)
bash bin/gmx_extract.sh -s md.tpr -f md.xtc --steps pbc -b 0 -e 100000
```

### Resampling

```bash
# Resample the whole system every 100 ps, no PBC treatment
bash bin/gmx_extract.sh -s md.tpr -f md.xtc -g System --steps extract --dt 100
```

### Dry Run with Saved Script

```bash
# Preview the full pipeline and save the commands to a runnable script
bash bin/gmx_extract.sh -s md.tpr -f md.xtc --dry-run --save-script run_viz.sh
```

### Non-Contiguous Group

```bash
# Protein + ion group: build a reduced .tpr to keep indices consistent
bash bin/gmx_extract.sh -s md.tpr -f md.xtc -g Protein_Ion --fit-group C-alpha \
    --index index.ndx --subset-tpr
```

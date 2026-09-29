# gmx_analysis.sh

Standard post-processing analysis for a (possibly segmented) GROMACS production run.

## Overview

Takes a production run produced by [gmx_continue_grompp.sh](md_continuation.md) — segment
files named `<basename>_<replica>_<start>_<end>.{xtc,edr,tpr}` — concatenates all segments
into a single trajectory and energy file, and runs a standard set of analyses on it. Each
analysis is written to its own clearly named `.xvg` file in one output directory, ready for
[plot_xvg.py](plot_xvg.md) or the [gmx_panel.py](gmx_panel.md) web dashboard.

## Features

- ✅ **Automatic segment discovery**: finds every `<basename>_<replica>_<start>_<end>.xtc` in
  `--segments-dir`, sorts them numerically by start time, and checks each has a matching `.edr` and `.tpr`
- ✅ **Concatenation** with `gmx trjcat` / `gmx eneconv` (reusable later with `--skip-concat`)
- ✅ **Six analysis categories**: energy terms, RMSD, radius of gyration, SASA, secondary structure (DSSP), RMSF
- ✅ **Selectable analyses** via `--only` / `--skip`
- ✅ **Time slicing** with `--begin`/`--end`/`--dt`, given in any unit via `--tu`
- ✅ **Plot time unit** (`--plot-tu`, default `ns`) independent of the slicing unit, and a separate
  unit for the energy plots (`--energy-tu`, default `ps`)
- ✅ **Custom groups** from an index file (`-n/--index`)
- ✅ **Fault tolerant**: a failing analysis is reported and skipped; the remaining analyses still run
- ✅ **Dry-run mode** to preview commands without running them
- ✅ **Command logging** (on by default): every command run, plus the exact invocation used, is
  saved to `<output-dir>/gmx_analysis_commands.sh` (override with `--save-script`), written even with `--dry-run`

## Pipeline

| Step | Category | Tool | Output |
|------|----------|------|--------|
| 1 | *(always, unless `--skip-concat`)* | `gmx trjcat`, `gmx eneconv` | `<basename>_<replica>_concat.xtc`, `<basename>_<replica>_concat.edr` |
| 2 | `energy` | `gmx energy` | one `.xvg` per term, e.g. `temperature.xvg`, `pressure.xvg`, `potential.xvg`, `total_energy.xvg` |
| 3 | `rmsd` | `gmx rms` | `rmsd.xvg` |
| 4 | `gyration` | `gmx gyrate` | `gyration_radius.xvg` |
| 5 | `sasa` | `gmx sasa` | `sasa.xvg` |
| 6 | `secondary-structure` (alias `dssp`) | `gmx dssp` | `secondary_structure.dat`, `secondary_structure_counts.xvg` |
| 7 | `rmsf` | `gmx rmsf` | `rmsf.xvg` (per residue by default) |

Energy file names are the term name in lower case with `-` replaced by `_`
(`Total-Energy` → `total_energy.xvg`).

`gmx dssp` requires GROMACS ≥ 2023. On older builds the secondary-structure step is skipped
with a message instead of failing.

## Usage

### Basic Syntax

```bash
gmx_analysis.sh <basename> <replica> [OPTIONS]
```

| Argument | Description |
|----------|-------------|
| `basename` | Base name used for the run (e.g. `MnMT4_apo`) |
| `replica` | Replica number (e.g. `0`), so segments are `<basename>_<replica>_<start>_<end>.*` |

Both positional arguments must come before any option.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--segments-dir <dir>` | `.` | Directory containing the segment files |
| `--structure <file>` | `.tpr` of the earliest segment | Reference `.tpr`/`.gro`/`.pdb` used for all analyses (including as the RMSD reference) |
| `--skip-concat` | off | Reuse an existing `<output-dir>/<basename>_<replica>_concat.{xtc,edr}` instead of rebuilding it |
| `--begin <time>` | - | First frame to analyze, in the unit given by `--tu` |
| `--end <time>` | - | Last frame to analyze, in the unit given by `--tu` |
| `--dt <time>` | - | Only use frames spaced by this interval, **always in ps** regardless of `--tu` (not applied to energy terms, see [Time Units](#time-units)) |
| `--tu <ps\|ns\|us\|fs>` | `ps` | Unit for `--begin`/`--end` |
| `--plot-tu <ps\|ns\|us\|fs>` | `ns` | Time unit for the x-axis of the generated `.xvg` files (passed as `-tu` to each tool that accepts it; energy terms use `--energy-tu`) |
| `--energy-tu <ps\|ns\|us\|fs>` | `ps` | Time unit for the x-axis of the energy-term `.xvg` files, applied by rescaling their time column after `gmx energy` runs |
| `-n, --index <file>` | - | Index (`.ndx`) file passed to every tool that accepts one (all except `gmx energy`) |
| `-g, --group <name>` | `Protein` | Group used for RMSD (calculation), gyration, SASA (surface), secondary structure and RMSF |
| `--fit-group <name>` | same as `--group` | Group used for the RMSD least-squares fit |
| `--skip <list>` | - | Comma-separated categories to skip |
| `--only <list>` | - | Comma-separated categories to run exclusively (takes precedence over `--skip`) |
| `--energy-terms <list>` | `Temperature,Pressure,Potential,Total-Energy` | Comma-separated `gmx energy` term names, one `.xvg` per term |
| `--probe-radius <nm>` | `0.14` | SASA solvent probe radius (`gmx sasa -probe`) |
| `--dssp-hmode <gromacs\|dssp>` | `gromacs` | `gmx dssp -hmode` |
| `--per-atom` | off | RMSF per atom instead of the default per residue (`-res`) |
| `--output-dir <dir>` | `<basename>_<replica>_analysis` | Where `.xvg`/`.dat`/log files are written |
| `--xvg-format <xmgrace\|xmgr\|none>` | `xmgrace` | Passed as `-xvg` to every tool that supports it |
| `--gmx <path>` | `gmx` | gmx binary/command to use |
| `--dry-run` | off | Print commands (and stdin group selections) without running them |
| `--save-script <file>` | `<output-dir>/gmx_analysis_commands.sh` | Reproducibility script recording every command run (and the exact invocation of this tool); written even with `--dry-run` |
| `-h, --help` | - | Show help and exit |

Valid category names for `--skip`/`--only`: `energy`, `rmsd`, `gyration`, `sasa`,
`secondary-structure` (or `dssp`), `rmsf`. Unknown names are rejected before anything runs.

### Available Default Groups

From a typical protein `.tpr`, no index file required:

`System, Protein, Protein-H, C-alpha, Backbone, MainChain, MainChain+Cb, MainChain+H,
SideChain, SideChain-H, Prot-Masses, non-Protein, Water, SOL, non-Water, Ion, Water_and_ions`

## Time Units

Two separate options control time:

- **`--tu`** only says how to read `--begin`/`--end`. The values are converted to ps
  internally so that every tool sees the same frame range.
- **`--plot-tu`** sets the x-axis unit of the output plots (`-tu` for each tool). The
  range is converted back to this unit before being passed as `-b`/`-e`/`-dt`, because gmx
  tools read those flags in whatever unit `-tu` specifies.

`--dt` is always in ps.

Not every gmx tool accepts all of these options, and the set differs between GROMACS
versions. Before each analysis, the script reads the tool's `-h` output and passes only the
options it supports:

- Tools without `-tu` get `-b`/`-e`/`-dt` in ps. `gmx rmsf` has no time axis at all; the
  range only selects the frames it averages over.
- **`gmx energy` has neither `-tu` nor `-dt`.** It always gets `-b`/`-e` in ps and writes
  its time axis in ps; energy terms use every frame in the `--begin`/`--end` range. The
  script prints a `Note:` line when an option (such as `-dt`) is left out.

### Energy plot time unit (`--energy-tu`)

Energy plots have their own time unit option, `--energy-tu` (default `ps`), independent of
`--plot-tu`. Since `gmx energy` can't write any other unit, the option is not passed to it.
Instead, after each energy term is calculated, the script:

1. multiplies the first (time) column of the `.xvg` file by the conversion factor
   (e.g. `0.001` for ns), leaving every other column and header line untouched, and
2. changes the x-axis label from `Time (ps)` to the new unit.

Both steps are ordinary logged commands (`awk`, then `mv`), so they appear in `--dry-run` output
and in the reproducibility script.

```bash
# Energy plots in ns, to match the other plots
bash bin/gmx_analysis.sh MnMT4_apo 0 --energy-tu ns
```

```bash
# Analyze 500-1000 ns every 100 ps; plot in ns (the default)
bash bin/gmx_analysis.sh MnMT4_apo 0 --begin 500 --end 1000 --tu ns --dt 100
```

## Output Files

All written to `--output-dir` (default `<basename>_<replica>_analysis/`):

| File | Description |
|------|-------------|
| `<basename>_<replica>_concat.xtc` / `.edr` | Concatenated trajectory and energy files |
| `*.xvg` | One file per analysis (see [Pipeline](#pipeline)) |
| `secondary_structure.dat` | Per-residue DSSP string for every frame |
| `<basename>_<replica>_analysis.log` | Full stdout/stderr of every gmx command |
| `gmx_analysis_commands.sh` | Reproducibility script (see `--save-script`) |

## Examples

### Full Analysis

```bash
# Segments MnMT4_apo_0_0_500.xtc, MnMT4_apo_0_500_1000.xtc, ... in the current directory
bash bin/gmx_analysis.sh MnMT4_apo 0
```

### Time Window

```bash
# Only the last 500 ns, sampled every 100 ps
bash bin/gmx_analysis.sh MnMT4_apo 0 --begin 500 --end 1000 --tu ns --dt 100

# Keep ps on the plot x-axis instead of the ns default
bash bin/gmx_analysis.sh MnMT4_apo 0 --plot-tu ps
```

### Custom Reference and Groups

```bash
# RMSD/RMSF/gyration relative to the equilibrated structure, backbone only
bash bin/gmx_analysis.sh MnMT4_apo 0 --structure equil_npt2.gro --group Backbone

# Fit on C-alpha, compute RMSD for a custom group from an index file
bash bin/gmx_analysis.sh MnMT4_apo 0 -n index.ndx --fit-group C-alpha --group Active_site
```

### Selecting Analyses

```bash
# Only recompute RMSD and RMSF, reusing the already concatenated trajectory
bash bin/gmx_analysis.sh MnMT4_apo 0 --only rmsd,rmsf --skip-concat

# Everything except the (slow) secondary-structure analysis
bash bin/gmx_analysis.sh MnMT4_apo 0 --skip dssp

# Extra energy terms
bash bin/gmx_analysis.sh MnMT4_apo 0 --only energy --energy-terms "Temperature,Density,Volume"
```

### Segments in Another Directory

```bash
bash bin/gmx_analysis.sh MnMT4_apo 0 --segments-dir ../production --output-dir apo_rep0
```

### Dry Run with Saved Script

```bash
bash bin/gmx_analysis.sh MnMT4_apo 0 --dry-run --save-script run_analysis.sh
```

### View the Results

```bash
# Serve all plots as a web dashboard
python bin/gmx_panel.py --dir MnMT4_apo_0_analysis

# Or plot individual files
python bin/plot_xvg.py MnMT4_apo_0_analysis/rmsd.xvg --style lines
```

## Troubleshooting

### "No trajectory segments found"
Segment names must match `<basename>_<replica>_<start>_<end>.xtc` exactly, with integer
start/end values. Check `--segments-dir`, and that `basename`/`replica` match the names
passed to `gmx_continue_grompp.sh`.

### "Missing energy file / run input file for segment"
Every `.xtc` segment needs a matching `.edr` and `.tpr`. Move incomplete segments (e.g. from
a crashed run) out of the directory, or finish them first.

### An analysis failed
Failures are reported as `✗ ... failed (see <log>)` and the next analysis continues. The
gmx error is in `<output-dir>/<basename>_<replica>_analysis.log`. The usual cause is a group
name that doesn't exist in the `.tpr`; pass `-n index.ndx` for custom groups.

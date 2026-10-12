# gmx_analysis.sh

Standard post-processing analysis for a (possibly segmented) GROMACS production run.

## Overview

Takes a production run produced by [gmx_continue_grompp.sh](md_continuation.md) — segment
files named `<basename>_<replica>_<start>_<end>.{xtc,edr,tpr}` — concatenates all segments
into a single trajectory and energy file, and runs a standard set of analyses on it. A
**single trajectory** (`md.xtc`, `md.tpr`, optionally `md.edr`) can be analyzed too: give the
base name without a replica number (see [Single Runs](#single-runs)). Each
analysis is written to its own clearly named `.xvg` file in one output directory, ready for
[plot_xvg.py](plot_xvg.md) or the [gmx_panel.py](gmx_panel.md) web dashboard.

## Features

- ✅ **Automatic segment discovery**: finds every `<basename>_<replica>_<start>_<end>.xtc` in
  `--segments-dir`, sorts them numerically by start time, and checks each has a matching `.edr` and `.tpr`
- ✅ **Concatenation** with `gmx trjcat` / `gmx eneconv` (reusable later with `--skip-concat`), or on its
  own with `--concat-only`
- ✅ **Six analysis categories**: energy terms, RMSD, radius of gyration, SASA, secondary structure (DSSP), RMSF
- ✅ **Distances, angles and dihedrals** between atoms or group centres, defined in a measurements file
  (`--measure`), with every selection checked before anything runs
- ✅ **Selectable analyses** via `--only` / `--skip`
- ✅ **Time slicing** with `--begin`/`--end`/`--dt`, given in any unit via `--tu`
- ✅ **Plot time unit** (`--plot-tu`, default `ns`) independent of the slicing unit, and a separate
  unit for the energy plots (`--energy-tu`, default `ps`)
- ✅ **Custom groups** from an index file (`-n/--index`)
- ✅ **Fault tolerant**: a failing analysis is reported and skipped; the remaining analyses still run
- ✅ **Dry-run mode** to preview commands without running them
- ✅ **Command logging** (on by default): every command run, plus the exact invocation used, is
  saved to `<output-dir>/gmx_analysis_commands.sh` (override with `--save-script`), written even with `--dry-run`; the file from an earlier run is kept as `…_commands.1.sh`, `.2.sh`, … (highest number = most recent)

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
| 8 | `distance`, `angle` (only with `--measure`) | `gmx distance`, `gmx gangle` | one file per measurement: `distance_<name>.xvg`, `angle_<name>.xvg`, `dihedral_<name>.xvg` |

Energy file names are the term name in lower case with `-` replaced by `_`
(`Total-Energy` → `total_energy.xvg`).

`gmx dssp` requires GROMACS ≥ 2023. On older builds the secondary-structure step is skipped
with a message instead of failing.

## Usage

### Basic Syntax

```bash
gmx_analysis.sh <basename> <replica> [OPTIONS]    # segmented run (gmx_continue_grompp.sh)
gmx_analysis.sh <basename> [OPTIONS]              # single run: <basename>.xtc
```

| Argument | Description |
|----------|-------------|
| `basename` | Base name used for the run (e.g. `MnMT4_apo`; for a single run e.g. `md` or `md.xtc`, optionally with a path) |
| `replica` | Replica number (e.g. `0`), so segments are `<basename>_<replica>_<start>_<end>.*`. Without it, the run is analyzed as a single trajectory |

The positional arguments must come before any option.

### Single Runs

Without a replica number, `<basename>.xtc` is analyzed directly, with `<basename>.tpr` as the
reference structure (or `--structure`) and `<basename>.edr` for the energy terms. These are the
files `gmx mdrun -deffnm <basename>` writes. They are looked for in `--segments-dir` (default:
the current directory), unless the base name includes a path.

```bash
gmx_analysis.sh md                          # md.xtc, md.tpr, md.edr -> md_analysis/
gmx_analysis.sh ../run2/md.xtc --only rmsd,rmsf
```

- Nothing is concatenated, so `--skip-concat` and `--concat-only` don't apply (and are
  rejected).
- If there is no `<basename>.edr`, the energy analyses are skipped with a note; everything else
  runs.
- Output goes to `<basename>_analysis/` by default.
- Only one trajectory file is used: `gmx_continue_extend.sh` output in several parts
  (`md_0.part0002.xtc`, ...) needs joining with `gmx trjcat` first.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--segments-dir <dir>` | `.` | Directory containing the segment files, or a single run's files |
| `--structure <file>` | `.tpr` of the earliest segment, or `<basename>.tpr` for a single run | Reference `.tpr`/`.gro`/`.pdb` used for all analyses (including as the RMSD reference) |
| `--skip-concat` | off | Reuse an existing `<output-dir>/<basename>_<replica>_concat.{xtc,edr}` instead of rebuilding it |
| `--concat-only` | off | Only build the concatenated `.xtc`/`.edr` files, then stop without running any analysis. Cannot be combined with `--only`, `--skip` or `--skip-concat` |
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
| `--measure <file>` | - | Distances, angles and dihedrals to measure, see [Distances, Angles and Dihedrals](#distances-angles-and-dihedrals---measure) |
| `--output-dir <dir>` | `<basename>_<replica>_analysis` (single run: `<basename>_analysis`) | Where `.xvg`/`.dat`/log files are written |
| `--xvg-format <xmgrace\|xmgr\|none>` | `xmgrace` | Passed as `-xvg` to every tool that supports it |
| `--gmx <path>` | `gmx` | gmx binary/command to use |
| `--dry-run` | off | Print commands (and stdin group selections) without running them |
| `--save-script <file>` | `<output-dir>/gmx_analysis_commands.sh` | Reproducibility script recording every command run (and the exact invocation of this tool); written even with `--dry-run` |
| `-h, --help` | - | Show help and exit |

Valid category names for `--skip`/`--only`: `energy`, `rmsd`, `gyration`, `sasa`,
`secondary-structure` (or `dssp`), `rmsf`, `distance`, `angle` (the last two need `--measure`).
Unknown names are rejected before anything runs.

### Available Default Groups

From a typical protein `.tpr`, no index file required:

`System, Protein, Protein-H, C-alpha, Backbone, MainChain, MainChain+Cb, MainChain+H,
SideChain, SideChain-H, Prot-Masses, non-Protein, Water, SOL, non-Water, Ion, Water_and_ions`

## Distances, Angles and Dihedrals (`--measure`)

List the measurements in a text file, one per line, and pass it with `--measure`:

```
# type      name        selections, separated by |        ('#' starts a comment line)
distance    SAM_C36     resname SAM and name SD | resid 36 and name SG
distance    SAM_loop    com of resname SAM | com of group "Loop1"
angle       hinge       resid 10 and name CA | com of group "Helix" | resid 90 and name CA
dihedral    chi1_W6     resid 6 and name N | resid 6 and name CA | resid 6 and name CB | resid 6 and name CG
```

```bash
bash bin/gmx_analysis.sh MnMT4_apo 0 --measure measurements.txt -n index.ndx
# Only the measurements, reusing the concatenated trajectory
bash bin/gmx_analysis.sh MnMT4_apo 0 --measure measurements.txt --only distance,angle --skip-concat
```

- **type**: `distance` (2 selections, `gmx distance`), `angle` (3) or `dihedral` (4) (both
  `gmx gangle`). Distances are in nm, angles and dihedrals in degrees.
- **name**: letters, digits, `_ . -`; used in the output file name (`distance_SAM_C36.xvg`) and
  as the plot title. Names must be unique.
- **selections**: [GROMACS selection syntax](https://manual.gromacs.org/current/onlinehelp/selections.html)
  (`gmx help selections`). **Each must give exactly one position**: one atom, or `com of` /
  `cog of` (centre of mass / geometry) a group. Periodic boundaries are taken into account.

Measurements form the categories `distance` and `angle` (dihedrals included) for `--only` /
`--skip`, and they run only when `--measure` is given. They use the same time window and
`--plot-tu` as the other analyses.

### Checks before anything runs

The file is checked line by line (type, name, number of selections, empty selections), then
every selection is evaluated with `gmx select` on the structure file (and index file). Any
problem stops the script, before the concatenation, with all problems listed at once, e.g.:

```
✗ distance whole_res: 'resname TRP' matches 24 atom(s); it needs exactly 1 (one atom, or "com of"/"cog of" a group)
✗ distance typo: gmx select rejected a selection:
      Selection 'resname TRP and name NEE1' never matches any atoms.
✗ angle no_group: gmx select rejected a selection:
      Cannot match 'group "Loop9"', because no such index group can be found.
```

Selections are evaluated on the structure file; ones that depend on positions (e.g. `within`)
could match a different number of atoms in later frames.

### Selection tips

| You want | Selection |
|---|---|
| an atom by residue name (if that residue occurs once) | `resname SAM and name SD` |
| an atom by residue number | `resid 36 and name SG` |
| centre of mass of a residue / of an index group | `com of resname SAM`, `com of group "Loop1"` |
| a group combined with other conditions | `group "Helix" and name CA` (a bare `"Helix" and ...` is a syntax error) |
| logic, ranges, several names | `com of ((resname TRP or resname TYR) and not name N CA C O)`, `resid 1 to 10` |
| wildcards | `name "H*"` (must be quoted) |

- Atom names follow your topology: check them with `grep SAM conf.gro` or try a selection with
  `gmx select -s md.tpr -n index.ndx -select 'resname SAM and name SD'`.
- With `-n/--index`, GROMACS' predefined groups (`Protein`, `Backbone`, ...) are only available
  if the index file contains them (files made with `gmx make_ndx` do).

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

### Energy plot post-processing (titles, `--energy-tu` and statistics)

`gmx energy` titles every plot `GROMACS Energies` and can only write its time axis in ps.
After each energy term is calculated, the script post-processes its `.xvg` file:

1. **Title**: sets the title to the term's name, taken from the legend `gmx energy` writes
   (e.g. `Temperature`, `Total Energy`, `Kinetic En.`). This is always done.
2. **Time unit**: energy plots have their own time unit option, `--energy-tu` (default `ps`),
   independent of `--plot-tu`. It is not passed to `gmx energy`. Instead, if it is not `ps`, the
   first (time) column is multiplied by the conversion factor (e.g. `0.001` for ns) and the
   x-axis label is changed from `Time (ps)` to the new unit.
3. **Statistics**: adds the summary `gmx energy` prints for the term (average, error
   estimate, RMSD and total drift over the analyzed range) as a header comment line:

   ```
   # gmx_analysis stats: Average=299.99 Err.Est.=0.098 RMSD=3.53537 Tot-Drift=0.6368 unit=K
   ```

   Other tools ignore `#` lines, so the file stays a normal `.xvg`. `plot_xvg.py --stats` and
   `gmx_panel.py --stats` show these values in a text box on the plot. To read the table,
   `gmx energy`'s output is saved to a temporary `<term>.xvg.energy.txt` file, which is
   removed afterwards.

Nothing else in the file changes: the data columns and all other header lines are kept as
`gmx energy` wrote them. The post-processing runs as ordinary logged commands (`awk`, then
`mv`), so it appears in `--dry-run` output and in the reproducibility script.

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

### Analyzing One Segment

Each segment of a segmented run is a complete set of files (`.xtc`, `.tpr`, `.edr`), so giving
a segment's full name without a replica number analyzes just that segment, as a
[single run](#single-runs): nothing is concatenated, which makes it fast.

```bash
# Only the 500-1000 ns segment of replica 0 -> MnMT4_apo_0_500_1000_analysis/
bash bin/gmx_analysis.sh MnMT4_apo_0_500_1000 --segments-dir ../prod
```

Plots keep the segment's real times (here 500-1000 ns) on the x-axis. For a time range that
spans several segments, use a time window on the whole run instead (`--begin`/`--end`, see
[Time Window](#time-window)); with `--skip-concat` it reuses an existing concatenated trajectory.

### Custom Reference and Groups

```bash
# RMSD/RMSF/gyration relative to the equilibrated structure, backbone only
bash bin/gmx_analysis.sh MnMT4_apo 0 --structure equil_npt2.gro --group Backbone

# Fit on C-alpha, compute RMSD for a custom group from an index file
bash bin/gmx_analysis.sh MnMT4_apo 0 -n index.ndx --fit-group C-alpha --group Active_site
```

### Selecting Analyses

```bash
# Only concatenate the segments (no analysis), e.g. to inspect the trajectory first
bash bin/gmx_analysis.sh MnMT4_apo 0 --concat-only

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

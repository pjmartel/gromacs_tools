# gmx_pca.py

Automates Principal Component Analysis (PCA) on GROMACS trajectories with automatic plotting of essential results.

## Overview

Performs complete PCA workflow from trajectory preprocessing through eigenvector analysis and visualization. Handles all GROMACS command execution and generates publication-ready plots automatically.

## Features

- ✅ Atom selection with GROMACS selection syntax, for the atoms analysed (`--selection`) and,
  separately, the atoms used for fitting (`--fit-selection`): e.g. fit on the protein backbone
  and analyse a bound ligand, a loop or the N-terminus
- ✅ Choice of frames: `-b`/`-e` (time range, in ps or `--tu`) and `--dt` (spacing)
- ✅ PBC treatment (whole molecules, or no jumps across the box edge), centring and fitting
- ✅ Covariance matrix and eigenvectors
- ✅ Projections on any range of PCs (`--first`/`--last`), one file per PC (ready for 2D and
  3D plots with `plot_xvg.py`) and one with all of them
- ✅ Variance explained per PC: a table in the terminal and `variance.xvg`
- ✅ Automatic plots (scree, cumulative variance, projections over time, scatter of two PCs),
  plus `gmx_pca_plot_commands.sh`: the same plots as `plot_xvg.py` commands, to edit and re-run
- ✅ **Several trajectories** (replicas, apo vs holo): one PCA of all of them, of one group or of
  one trajectory, with every trajectory projected on the same PCs and kept apart in the plots
  and files (colour per trajectory or per group). See [Several Trajectories](#several-trajectories-replicas-apo-vs-holo)
- ✅ **RMSIP** (`--rmsip`): how similar the motions are, from separate PCAs: first vs second half
  of each trajectory (convergence), trajectory vs trajectory, group vs group, within and
  between groups. See [Comparing the Motions](#comparing-the-motions-rmsip---rmsip)
- ✅ **Region split** (`--region`): every projection split into the part of a region (e.g. a floppy
  terminus) and the rest: the projections without the region, and how much of each PC it is.
  See [The Part of a Region](#the-part-of-a-region---region)
- ✅ Safe re-runs: earlier preprocessing is reused only if made with the same files and options
- ✅ **Command logging** (on by default): every GROMACS command run, plus the exact invocation
  used, is saved to `<outdir>/gmx_pca_commands.sh` (override with `--commands-file`); the file from an earlier run is kept as `…_commands.1.sh`, `.2.sh`, … (highest number = most recent)

## Requirements

- GROMACS 2021+ (tested with 2023.3 and 2026.2)
- Python 3.8+
- numpy, matplotlib (installed automatically when run directly with [uv](https://docs.astral.sh/uv/))

## Quick Start

```bash
# Backbone PCA, PC1 and PC2
gmx_pca.py -s topol.tpr -f traj.xtc

# First 5 PCs
gmx_pca.py -s topol.tpr -f traj.xtc --last 5

# A bound ligand, fitted on the protein backbone
gmx_pca.py -s topol.tpr -f traj.xtc --selection "resname SAM" --fit-selection Backbone \
    -o pca_sam --last 3

# Five replicas: one PCA of all of them, each replica in its own colour
gmx_pca.py -s topol.tpr -f rep1.xtc rep2.xtc rep3.xtc rep4.xtc rep5.xtc -o pca_reps
```

## Workflow

1. **Index groups**: the atoms analysed (`--selection`), the atoms fitted (`--fit-selection`,
   default the same), and both together, the only atoms kept from here on
   (`selection.ndx` for the whole system, `pca_groups.ndx` for the reduced structure)
2. **Trajectory preprocessing**: takes the frames chosen with `-b`/`-e`/`--dt` (all by
   default), makes the kept atoms whole and centres them on the fitted atoms
   (PBC treatment, which needs the run's `.tpr`), takes the first frame of that as the reference
   structure `fit.gro`, and fits the trajectory to it on the fitted atoms. The reference is taken
   *after* the PBC treatment: taken from the raw trajectory, it would be split across the box
   whenever the protein crosses the box edge in the first frame, which silently corrupts the
   fit and the PCA
3. **Covariance analysis** (`gmx covar`): covariance matrix of the analysed atoms' coordinates
   (fitting on the fitted atoms), eigenvalues and eigenvectors; the variance per PC is printed
   and written to `variance.xvg`
4. **Projection** (`gmx anaeig`, one pass): the trajectory projected on PCs `--first` to
   `--last`, and the RMSF along them
5. **Plots**

Only the analysed atoms enter the PCA: their 3N coordinates (N atoms) are the dimensions of the
covariance matrix, and their motions are what the PCs describe. Atoms outside both selections
play no part.

## Usage

### Choosing the Atoms

```bash
# Protein backbone (the default)
gmx_pca.py -s topol.tpr -f traj.xtc --selection Backbone

# C-alpha atoms only
gmx_pca.py -s topol.tpr -f traj.xtc --selection "name CA"

# Specific residues
gmx_pca.py -s topol.tpr -f traj.xtc --selection "resid 10 to 50 and name CA"
```

### Fitting on Other Atoms (`--fit-selection`)

By default the same atoms are fitted and analysed. To see how one part moves *relative to* the
rest, fit on the stable part and analyse the moving one:

```bash
# Distribution of a bound ligand in the active site, relative to the protein
gmx_pca.py -s topol.tpr -f traj.xtc --selection "resname SAM" --fit-selection Backbone -o pca_sam

# Conformations of a loop
gmx_pca.py -s topol.tpr -f traj.xtc --selection "resid 45 to 60 and name N CA C" \
    --fit-selection Backbone -o pca_loop

# The N-terminus, fitted on the rest of the backbone
gmx_pca.py -s topol.tpr -f traj.xtc --selection "resid 1 to 15 and name N CA C" \
    --fit-selection "name N CA C and not resid 1 to 15" -o pca_nterm
```

Centring also uses the fitted atoms. With PBC treatment `mol`, every molecule is made whole and
placed in the compact unit cell around them, so a ligand appears in the periodic image closest
to the protein's centre. This is right for ligands **bound** to the protein; a ligand or ion
that diffuses freely through the box mostly shows that wandering, which a PCA can't make much
of.

### Choosing Frames (`-b`, `-e`, `--dt`, `--tu`)

```bash
# 100-500 ns only (e.g. leaving out equilibration)
gmx_pca.py -s topol.tpr -f traj.xtc -b 100 -e 500 --tu ns

# Every 100 ps (fewer frames: faster for long runs)
gmx_pca.py -s topol.tpr -f traj.xtc --dt 100
```

The frames are chosen where the raw trajectory is read (the first `gmx trjconv` step), so the
fitting, the covariance analysis and the projections all use the same frames, and the
reference structure is the first frame of the range. `-b`/`-e` are in the unit of `--tu`
(default ps); `--dt` is always in ps (`trjconv -dt`: frames whose time is a multiple of it, from
the first one). A `-b` beyond the end of the trajectory is reported by GROMACS ("Specified frame
... doesn't exist").

### PCs to Project and Plot (`--first`, `--last`)

```bash
# PC1 to PC5: five projections and time series, the PC1-PC2 scatter, RMSF along PC1-PC5
gmx_pca.py -s topol.tpr -f traj.xtc --last 5

# Only PC3 and PC4
gmx_pca.py -s topol.tpr -f traj.xtc --first 3 --last 4
```

The number of PCs is at most 3 × (number of analysed atoms), and smaller than the number of
frames; a `--last` beyond that is reported with the actual number.

### PBC and Fitting Options

```bash
# Molecules whole, then no jumps across the box edge between frames (long runs)
gmx_pca.py -s topol.tpr -f traj.xtc --pbc nojump

# No PBC treatment (trajectory already whole and centred)
gmx_pca.py -s topol.tpr -f traj.xtc --pbc no

# Progressive fitting (each frame to the previous one)
gmx_pca.py -s topol.tpr -f traj.xtc --fit progressive

# No fitting (pre-aligned trajectories)
gmx_pca.py -s topol.tpr -f traj.xtc --fit none --no-center
```

### Output Directory

```bash
gmx_pca.py -s topol.tpr -f traj.xtc -o pca_results

# Make the preprocessing again even when it could be reused
gmx_pca.py -s topol.tpr -f traj.xtc -o pca_out --overwrite
```

### Custom GROMACS Binary

```bash
gmx_pca.py -s topol.tpr -f traj.xtc --gmx-bin gmx_2024
gmx_pca.py -s topol.tpr -f traj.xtc --gmx-bin /usr/local/gromacs/2024/bin/gmx
```

### Plotting Control

```bash
# Skip plotting (only run GROMACS steps)
gmx_pca.py -s topol.tpr -f traj.xtc --no-plots

# Then use plot_xvg.py for custom plots: all projections over time
plot_xvg.py pca_out/projections.xvg
```

## Several Trajectories (Replicas, Apo vs Holo)

Each PCA has its own principal components: PC1 of one replica is in general another direction
than PC1 of another, so the projections of separate PCAs can't be compared. To compare
trajectories, they must all be projected on the **same** PCs. With several trajectories,
`gmx_pca.py` does that, and keeps each trajectory's frames apart, so the plots and files show
which points belong to which trajectory:

```bash
# Five holo replicas: one PCA of all of them, each replica in its own colour
gmx_pca.py -s topol.tpr -f holo1.xtc holo2.xtc holo3.xtc holo4.xtc holo5.xtc -o pca_holo

# Apo vs holo, from a file listing the ten trajectories (below):
# one PCA of all ten, apo in one colour and holo in another
gmx_pca.py --trajectories replicas.txt --color-by group -o pca_apo_holo

# The PCA of the apo replicas only, with the holo ones projected on it
gmx_pca.py --trajectories replicas.txt --pca-from apo --color-by group -o pca_from_apo
```

### Giving the Trajectories

**On the command line** (simple cases): several files after `-f`, and after `-s` one `.tpr`
for all of them or one per trajectory (in the same order). Each trajectory gets a label, by
default its file name (`holo1.xtc` → `holo1`; if those are all the same, as in `rep1/md.xtc`,
`rep2/md.xtc`, the folder too: `rep1_md`), or give them with `--labels`. `--groups` gives each
a group:

```bash
gmx_pca.py -s apo.tpr apo.tpr holo.tpr holo.tpr -f apo_r1.xtc apo_r2.xtc holo_r1.xtc holo_r2.xtc \
    --labels apo1 apo2 holo1 holo2 --groups apo apo holo holo
```

**In a file** (`--trajectories FILE`), one trajectory per line: label, group, trajectory and
structure (`.tpr`). `#` starts a comment; a group `-` means none; without a structure (or with
`-`) the one given with `-s` is used. Relative paths are relative to the file's folder; put
paths with spaces in quotes.

```
# label  group  trajectory              structure
apo1     apo    apo/rep1/md_all.xtc     apo/rep1/md.tpr
apo2     apo    apo/rep2/md_all.xtc     apo/rep2/md.tpr
holo1    holo   holo/rep1/md_all.xtc    holo/rep1/md.tpr
holo2    holo   holo/rep2/md_all.xtc    holo/rep2/md.tpr
```

Labels and groups name files, so they may only contain letters, digits and `_ . + -`; a label
can't also be a group name, and neither can be `all`. Each trajectory is one file: concatenate
a replica's segments first (e.g. with `gmx_analysis.sh` or `gmx trjcat`).

### What the PCA Is Made Of (`--pca-from`)

- `--pca-from all` (default): one PCA of all trajectories together. Symmetric: no trajectory is
  the reference, and the PCs describe the motions over all of them, including the differences
  between them (e.g. apo vs holo, if those are large). Each trajectory weighs by its number of
  frames: give replicas of similar length, or use `-b`/`-e`/`--dt` to even them out
- `--pca-from GROUP` (e.g. `apo`): the PCA of that group's trajectories, with all
  trajectories projected on it. Asks "how do the holo replicas look in the apo motions?":
  holo motions that apo doesn't have are not shown
- `--pca-from LABEL` (e.g. `apo1`): the PCA of one trajectory, all projected on it

For apo vs holo, both views are worth having: e.g. `-o pca_all` with the default and
`-o pca_apo` with `--pca-from apo`.

### Same Options, Same Atoms, One Reference

`--selection`, `--fit-selection`, `--pbc`, centring, `--fit` and `-b`/`-e`/`--dt` apply to every
trajectory. Each trajectory's own `.tpr` is used for its index groups and PBC treatment, so
the systems may differ (e.g. holo has a ligand and apo doesn't), but the **selected atoms must
be the same** in all of them: same number, same atom names, same order (the PCA compares the
coordinates atom by atom). `gmx_pca.py` checks that and stops with the first difference, e.g.
when `--selection Protein` includes a ligand present in only some trajectories, or when the
proteins have other residues. Residue names may differ (e.g. `HIS`/`HSD` protonation), with a
note.

All trajectories are fitted to **one** reference structure, `fit.gro`: the first frame of the
first trajectory of the PCA (first in the list for `all`), with the PBC treatment. Fitting each
trajectory to its own first frame would put them in different orientations, and the
differences would show up as false motions. For the same reason, `--fit progressive` is not
accepted with several trajectories.

### How It Runs

1. Per trajectory, in `trajectories/<label>/`: the index groups, from its own `.tpr`
   (`selection.ndx`, `pca_groups.ndx`), the PBC treatment and centring of the selected atoms
   (`pbc.xtc`), and its first frame treated the same way (`first.gro`). This is the slow part,
   and is **reused** in later runs in the same output folder as long as the trajectory, its
   structure and the options it depends on stay the same (`preprocessing.json` in that folder;
   `--overwrite` makes it again). Changing `--pca-from`, `--last` or `--color-by` reuses it
2. The check that all trajectories give the same atoms; the common reference `fit.gro`; each
   trajectory fitted to it (`trajectories/<label>/fit.xtc`)
3. `gmx covar` on the trajectories of the PCA, put one after the other (`gmx trjcat -cat`, a
   temporary file); the variance table and `variance.xvg` as for one trajectory
4. `gmx anaeig` on each trajectory separately, so each keeps its identity; the RMSF along the
   PCs once (it depends only on the PCs)
5. The mean and standard deviation of each trajectory's (and each group's) projections, in the
   terminal and in `projection_stats.txt`; the plots

### Outputs

As for one trajectory, except:

- `projection_pc<N>_<label>.xvg`, `projections_<label>.xvg` - Each trajectory's projections
- `projection_pc<N>_group_<group>.xvg` - A group's projections, its trajectories one after the
  other (so `plot_xvg.py` can draw a group in one colour; times repeat)
- `projection_stats.txt` - Mean and standard deviation on each PC, per trajectory and per group:
  where each replica sits in the PC space, and how much of it it explores
- `trajectories/<label>/` - Each trajectory's index groups, PBC-treated and fitted trajectories,
  first frame and logs. For `gmx_pca_movie.py`, use `fit.gro` and `eigenvectors.trr` of the
  output folder with any trajectory's `trajectories/<label>/fit.xtc`
- `plot_pc<N>_timeseries` - One panel per trajectory, on the same scale
- `plot_pc<A>_pc<B>_scatter` - Every frame of every trajectory, one colour per trajectory
  (`--color-by trajectory`, default) or per group (`--color-by group`)
- `gmx_pca_plot_commands.sh` - The scatter plot (and a 3D one with three or more PCs) as
  `plot_xvg.py --xy-correlation --multi` commands, in the same colours; the projections over
  time as one plot with all trajectories (`plot_pc<N>_timeseries_overlay.png`), since
  `plot_xvg.py` doesn't make panels

The projections and plots of an earlier run in the same folder are removed first, so files of
another PCA, other PCs or other trajectories are not left behind among the new ones.

## Comparing the Motions: RMSIP (`--rmsip`)

Projecting everything on one set of PCs shows **where** the trajectories go. A different
question is whether they make the **same motions**. For that, each trajectory (or group) gets
its own PCA, and their eigenvectors are compared. The projections of separate PCAs can't be
compared, but the eigenvectors themselves can: they are directions in the same space of atom
coordinates, as long as the atoms and the reference structure are the same. `gmx_pca.py`
makes sure of both: every separate PCA uses the main PCA's atoms and reference, `fit.gro`.

```bash
# Convergence of one trajectory: its first half vs its second half
gmx_pca.py -s topol.tpr -f traj.xtc --rmsip

# Two force fields, five replicas each: all of the comparisons below
gmx_pca.py --trajectories ff_replicas.txt --color-by group --rmsip -o pca_ff
```

The inner product of two eigenvectors (both of unit length) is 1 for the same motion and 0 for
unrelated ones (its sign means nothing: an eigenvector can point either way). The **RMSIP**
(root mean square inner product) of the first N PCs (`--rmsip-pcs`, default 10) of two PCAs,

RMSIP = √( (1/N) Σᵢ Σⱼ (vᵢ · wⱼ)² ),  i, j = 1 … N

is 1 when the two sets of PCs span the same motions and 0 when they share none. Values above
about 0.7 are usually read as similar.

What `--rmsip` computes:

- **Convergence**, for every trajectory: its first half vs its second half (split halfway
  between the two middle frames). A low value means that the trajectory hasn't sampled its
  main motions consistently: it is still changing, or too short
- With several trajectories, **trajectory vs trajectory**: every pair. With groups, the
  terminal summarises them **within** each group (e.g. replica vs replica of ff14sb) and
  **between** groups (an ff14sb replica vs an ff19sb replica)
- With groups, **group vs group**: a PCA of each group (its trajectories one after the other),
  their RMSIP and the overlap matrix of their PCs (|vᵢ · wⱼ| for each pair of PCs), which shows
  which motion matches which: e.g. PC1 of one group may be PC2 of the other

**Read the values against the within-group ones.** The RMSIP is never 1, even for two replicas
of the same system: each samples only part of the motions, and shorter samples agree less. A
value between groups means something only next to the values within them:

```
Trajectory vs trajectory, within and between groups
  within ff14sb: mean 0.62  (0.55-0.70, 10 pairs)
  within ff19sb: mean 0.60  (0.52-0.68, 10 pairs)
  between ff14sb and ff19sb: mean 0.58  (0.49-0.66, 25 pairs)
```

If the values between groups are like those within them, the groups make the same motions,
as far as the replicas can tell; if they are clearly lower, the motions differ. The group PCAs
(more frames each) usually give higher values than single trajectories.

Each PCA needs more frames than PCs: with `--rmsip-pcs 10`, each half of each trajectory needs
at least 11 frames (more is much better; `gmx_pca.py` stops with a message otherwise). The
separate PCAs reuse the fitted trajectories and only write their first N eigenvectors: with C-
alphas they take seconds; with the backbone of a large protein, each one takes as long as the
main `gmx covar`.

Outputs (in the output folder):

- `rmsip_summary.txt` - The summary printed in the terminal
- `rmsip_halves.dat` - First vs second half of each trajectory, with the numbers of frames
- `rmsip_trajectories.dat` - The matrix of all trajectory pairs (several trajectories)
- `rmsip_groups.dat` - The matrix of all group pairs (two groups or more)
- `overlap_<groupA>_<groupB>.dat` - |inner product| of each PC of group A (rows) with each PC of
  group B (columns)
- `plot_rmsip_trajectories` - Heatmap of the trajectory matrix, with lines between the groups
- `plot_overlap_<groupA>_<groupB>` - Heatmap of the overlap matrix
- `rmsip/` - The separate PCAs (eigenvectors, eigenvalues, averages, logs)

The heatmaps are made by `gmx_pca.py` only (`plot_xvg.py` doesn't make heatmaps), so they are
not in `gmx_pca_plot_commands.sh`. The RMSIP is computed from the eigenvectors read with
`gmx dump`; for one pair, `gmx anaeig -v A.trr -v2 B.trr -first 1 -last N -over overlap.xvg`
gives the same number as the square root of the last value in `overlap.xvg`.

## The Part of a Region (`--region`)

How much of the PCA is due to one part of the protein, e.g. a floppy N-terminus? A projection is a
sum over the analysed atoms (each atom's displacement from the average structure times its part of
the eigenvector), so it splits exactly into the part of the region's atoms and the part of the
rest:

  projection(t) = region's part(t) + rest's part(t)

The rest's part is the trajectory **without the region**, projected on the same PCs: the
difference from the full projection is the effect of the region.

```bash
# Backbone PCA of 5 replicas: how much is the N-terminus? (fitted without it, see below)
gmx_pca.py --trajectories replicas.txt --selection Backbone \
    --fit-selection 'group "Backbone" and not resnr 1 to 10' \
    --region 'resnr 1 to 10' --region-name Nterm --last 3
```

Only the region's **analysed** atoms count: with `--selection 'name CA'`, `--region 'resnr 1 to
10'` means their C-alphas. It works with one trajectory or several (each split separately, and
each group's frames together), with `--pca-from` and `--rmsip`, and needs `--fit rot+trans`
(the default).

In the terminal and `region_split.txt`:

```
Weight of Nterm in each PC (share of the eigenvector's squared length on its atoms):
  PC1  76.6%  PC2  54.0%  PC3  15.1%

Variance of the projections, split: Nterm's part, the rest's part, and together
(twice their covariance: + they move together along the PC, - they cancel out), as % of the variance:
                  PC1: Nterm   rest together      PC2: Nterm   rest together
  apo1                   72%     6%      22%             44%    27%      29%
```

- **Weight**: how much of the PC's eigenvector is on the region (the same for every
  trajectory). Here PC1 is mostly the N-terminus
- **Variance split**, per trajectory and group: the region's part, the rest's part and their
  covariance term add up to 100 %. A negative "together" means that along this PC the region and
  the rest move against each other (their parts partly cancel); the rest's part can then be over
  100 %

**Leave the region out of the fit** (`--fit-selection`). If a floppy region is fitted, the fit
rotates the rest of the protein to follow it, and the region's motion is smeared over all atoms:
in the example above, with the N-terminus in the fit, PC1 is 17 % N-terminus instead of 77 %. A
warning is printed when the region's atoms are also fitted.

Outputs:

- `projection_pc<N>_parts[_<label>].xvg` - Time, the projection with all atoms, the region's
  part and the rest's part (without the region)
- `region_split.txt`, `region_weights.dat` - The tables above
- `plot_pc<N>_<name>_parts` - The three over time, a panel per trajectory
- `plot_pc<A>_pc<B>_without_<name>` - The scatter of the first two PCs with all atoms (grey)
  and without the region (colour), a panel per trajectory
- `gmx_pca_plot_commands.sh` - Also the time series (`plot_xvg.py` makes no panels nor this
  scatter)

This is the full PCA seen without the region's atoms: the PCs stay those of the full PCA. A
PCA made without the region has PCs of its own (run it with a `--selection` that leaves it out).

How it is checked: the split is computed in Python from the fitted trajectories and the
eigenvectors (read with `gmx dump`). The fitted frames must reproduce gmx covar's average
structure (to 0.002 nm, apart from one translation), and all atoms together must reproduce gmx
anaeig's projections. These agree to rounding (about 0.0005 nm) when the fitted and analysed
atoms are the same; when they differ, anaeig fits each frame itself, slightly differently from
covar (about 1 % of the projections' range; the split follows covar).

## Command-Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `-s, --structure` | Required | The run input file (`.tpr`): needed to make molecules whole. A GRO/PDB is accepted only with `--pbc no` (trajectory already whole and centred). Several trajectories: one for all, or one per trajectory |
| `-f, --trajectory` | Required | Trajectory file (XTC/TRR); several for a PCA of several trajectories |
| `--trajectories` | | A file listing several trajectories: `label group trajectory [structure]` per line (see [Giving the Trajectories](#giving-the-trajectories)) |
| `--labels` | file names | With several `-f`: a label per trajectory |
| `--groups` | none | With several `-f`: a group per trajectory (`-` for none) |
| `--pca-from` | `all` | Several trajectories: the PCA of `all`, of a group or of one trajectory (label); all are projected on it |
| `--color-by` | `trajectory` | Several trajectories: one colour per `trajectory` or per `group` in the plots |
| `--rmsip` | `False` | Also compare the motions with separate PCAs (RMSIP): first vs second half of each trajectory and, with several, trajectory vs trajectory and group vs group (see [Comparing the Motions](#comparing-the-motions-rmsip---rmsip)) |
| `--rmsip-pcs` | `10` | Number of PCs compared with `--rmsip` |
| `--region` | | Split every projection into the part of these analysed atoms and the rest (see [The Part of a Region](#the-part-of-a-region---region)) |
| `--region-name` | `region` | Name of `--region` in files, tables and plots, e.g. `Nterm` |
| `-o, --outdir` | `pca_out` | Output directory |
| `--selection` | `Backbone` | GROMACS selection of the atoms analysed: the PCA is of their coordinates |
| `--fit-selection` | `--selection` | GROMACS selection of the atoms used for fitting and centring |
| `-b, --begin` | first frame | First frame used, in the unit of `--tu` |
| `-e, --end` | last frame | Last frame used, in the unit of `--tu` |
| `--dt` | every frame | Only frames spaced by this interval, always in ps |
| `--tu` | `ps` | Unit of `-b`/`-e`: `fs`, `ps`, `ns` or `us` |
| `--first` | `1` | First PC projected and plotted |
| `--last` | `2` | Last PC projected and plotted |
| `--pbc` | `mol` | PBC treatment: `mol`, `res`, `atom`, `nojump` (whole, then no jumps between frames) or `no` |
| `--no-center` | `False` | Disable centring |
| `--fit` | `rot+trans` | Fitting mode (rot+trans/progressive/none) |
| `--gmx-bin` | `gmx` | GROMACS binary name/path |
| `--no-plots` | `False` | Skip matplotlib plotting |
| `--overwrite` | `False` | Make the index, fitted trajectory and reference again, even when the ones in `--outdir` could be reused (see [Re-running](#re-running-in-the-same-output-folder)) |

## Output Files

### GROMACS Outputs

All in `--outdir` (GROMACS runs there, so its own logs stay there too, and it makes no
`#file.N#` backups: re-running replaces the files).

- `selection.ndx` - Index groups of the whole system: all kept atoms, analysed, fitted
  (`select_out.ndx` is what `gmx select` wrote)
- `pca_groups.ndx` - The same groups in the reduced structure (`fit.gro`): fit, analysis, all.
  `gmx_pca_movie.py` uses it
- `fit.xtc` - Whole, centred and fitted trajectory of the kept atoms
- `fit.gro` - Reference structure: the first frame, whole and centred like `fit.xtc`
- `preprocessing.json` - The input files and options the files above were made with
- `eigenvalues.xvg` - Eigenvalues (nm²) from the covariance analysis
- `variance.xvg` - Per PC: eigenvalue, % of the total variance, cumulative %
- `eigenvectors.trr` - Eigenvectors (for projections/movies)
- `average.pdb` - Average structure of the analysed atoms
- `projection_pc<N>.xvg` - The trajectory projected on PC N (time, value), for PCs `--first` to `--last`
- `projections.xvg` - The same, all PCs as columns (time, PC first, ..., PC last)
- `rmsf.xvg` - RMSF per atom along each eigenvector `--first` to `--last` (one data set per PC)
- `covar.log` - Log of `gmx covar`; `logs/` - output of every GROMACS step (and, when one
  fails, its full output as `failed_gmx_<tool>.log`; the error message itself is printed)
- `gmx_pca_commands.sh` - The GROMACS commands run (see Features)
- `gmx_pca_plot_commands.sh` - `plot_xvg.py` commands for the plots (see [Tweaking the Plots](#tweaking-the-plots))

### Generated Plots

Each as `.png` and `.pdf`:

- `plot_eigenvalues` - Scree plot: variance (%) of each PC (the first 25 PCs)
- `plot_cumulative_variance` - Cumulative variance (%) (the first 25 PCs)
- `plot_pc<N>_timeseries` - Projection on PC N over time, for each PC `--first` to `--last`
- `plot_pc<A>_pc<B>_scatter` - The first two projected PCs against each other, coloured by time

### Tweaking the Plots

`gmx_pca_plot_commands.sh` holds one `plot_xvg.py` command per plot, with the same file names,
titles and labels (and, with three or more PCs, a 3D plot of the first three, which
`gmx_pca.py` itself doesn't make). Edit it — titles, axis labels, `--xlim`/`--ylim`, `--figsize`,
`--fontsize`, `--colormap`, `--output plot.pdf`, ... (`plot_xvg.py --help`) — and run it:

```bash
bash pca_out/gmx_pca_plot_commands.sh      # works from any folder
```

It replaces the PNG plots in the output folder with the new versions, which look much like
`gmx_pca.py`'s own (same sizes, no legends, the scatter coloured by time with a colour bar in ps:
`--figsize`, `--no-legend`, `--color-by time`). Running `gmx_pca.py` again writes a new script and keeps
the earlier one, maybe edited, as `gmx_pca_plot_commands.1.sh`, `.2.sh`, ...

It uses `plot_xvg.py` from `PATH`, or else the one next to `gmx_pca.py`.

### Variance per PC

The terminal shows the first 10 PCs (or up to `--last`), e.g.

```
    PC  eigenvalue (nm^2)   variance  cumulative
     1            0.11607     35.08%      35.08%
     2           0.045379     13.71%      48.80%
     3           0.038716     11.70%      60.50%
```

and `variance.xvg` has all of them. The percentages are each eigenvalue over the sum of all
eigenvalues; that sum is checked against the trace of the covariance matrix that `gmx covar`
reports (a warning is printed if they differ). They are relative to the analysed atoms and the
fit: the same simulation gives other numbers for `Backbone` than for `name CA`.

### Re-running in the Same Output Folder

The index, fitted trajectory and reference structure of an earlier run are reused only if
they were made from the same structure and trajectory files (same path, size and modification
time) with the same `--selection`, `--fit-selection`, `--pbc`, `--fit`, centring and frames
(`-b`/`-e`/`--dt`), as recorded in `preprocessing.json`. Otherwise they are made again, with a note saying what changed (e.g.
"changed since: selection"). `--overwrite` makes them again in any case. The covariance
analysis, projections and plots are always made again.

## Generated Plots Explained

### Scree Plot
Variance of each PC in descending order. Helps identify how many PCs are significant ("elbow method").

### Cumulative Variance
Cumulative percentage of variance explained by the first PCs.

### Scatter of Two PCs
2D projection showing conformational sampling. Clusters indicate metastable states.

### Time Series
How the conformations evolve along each PC over time. Useful for identifying transitions.

## Example Workflows

### Basic Backbone PCA

```bash
gmx_pca.py -s topol.tpr -f traj.xtc --last 3

# Outputs in pca_out/:
#    - variance.xvg and the table in the terminal: how much each PC explains
#    - projection_pc1/2/3.xvg, projections.xvg: the projections
#    - fit.gro, eigenvectors.trr: for movies (gmx_pca_movie.py)
#    - plot_*.png: automatic visualizations
```

### 3D Plot of PC1, PC2 and PC3

```bash
gmx_pca.py -s topol.tpr -f traj.xtc --last 3
cd pca_out
plot_xvg.py projection_pc1.xvg projection_pc2.xvg projection_pc3.xvg --xy-correlation
```

### PC1 vs PC2 with plot_xvg.py

```bash
plot_xvg.py pca_out/projections.xvg --xcol 1 --columns 2 --scatter
```

### Comparing Replicas

Give all replicas at once: one PCA, every replica projected on the same PCs and shown in its
own colour (see [Several Trajectories](#several-trajectories-replicas-apo-vs-holo)):

```bash
gmx_pca.py -s topol.tpr -f rep1.xtc rep2.xtc rep3.xtc --last 3 -o pca_reps
```

### Pre-processed Trajectory

```bash
# If the trajectory is already whole, centred and aligned
gmx_pca.py -s aligned.gro -f aligned.xtc --fit none --no-center --pbc no
```

## GROMACS Selection Syntax

Common selections:
- `Backbone` - Protein backbone atoms (N, CA, C)
- `name CA` or `group "C-alpha"` - Only C-alpha atoms
- `Protein` - All protein atoms
- `resid 10 to 50` - Residues 10-50
- `resname TRP TYR PHE` - Aromatic residues
- `group "SideChain-H"` - Heavy sidechain atoms
- `resid 45 and within 0.5 of resname LIG` - Atoms near ligand

Combine with `and`, `or`, `not`:
```bash
--selection "resid 10 to 50 and name CA"
--selection "(resid 10 to 30 or resid 60 to 80) and name N CA C"
```

Default group names (`Backbone`, `C-alpha`, `SideChain-H`, ...) are safest written as
`group "Backbone"`: a plain name works only on its own and without a hyphen (`Backbone` does;
`C-alpha` is read as "C minus alpha"), and inside an expression the `group` form is required
(`group "Backbone" and resid 1 to 15`). Or name the atoms: `name N CA C` for the backbone.
`Backbone and resid 1 to 15` and `resid 1 to 15 and backbone` are syntax errors.

## Tips & Best Practices

1. **Atom Selection**:
   - C-alpha only: Faster, captures global motion
   - Backbone: More detailed, standard choice
   - All-atom: Most detailed but slow and noisy

2. **Number of PCs**:
   - Check the variance table and the cumulative variance plot
   - Analyze more PCs (5-10) for complex systems

3. **PBC Handling**:
   - `mol`: Good default for proteins
   - `nojump`: molecules made whole, then kept continuous: for long runs where something
     crosses the box edge
   - `no`: Only for pre-processed trajectories

4. **Fitting**:
   - `rot+trans`: Standard
   - `progressive`: each frame to the previous one
   - `none`: For pre-aligned data

5. **Output Organization**:
   - Use descriptive output directory names
   - Keep different selections in separate directories

6. **Performance**:
   - C-alpha selection is much faster than all-atom
   - Consider reducing trajectory size with `gmx trjconv -skip`

## Troubleshooting

### "Selection is empty" / "--selection selects no atoms"
- Check your selection string syntax
- Verify atoms exist in your structure: `gmx select -s topol.tpr -select "your_selection"`

### "--pbc mol needs a .tpr"
- Give the run's `.tpr` as `-s`: making molecules whole needs the topology. With a `.gro` or
  `.pdb`, GROMACS can only warn and leaves molecules broken, which would corrupt the PCA.
- If the trajectory is already whole and centred, use `--pbc no` (see Pre-processed Trajectory)

### "--last N, but there are only M PCs"
- There are at most 3 × (analysed atoms) PCs, and fewer than the number of frames

### "X gives N atoms ..., Y M" / "atom K of the selection is ... in X but ... in Y"
- With several trajectories, the selections must give the same atoms in each. Make
  `--selection` (and `--fit-selection`) leave out what only some systems have, e.g. `Backbone`
  or `Protein` rather than a group that includes a ligand, and check that the proteins have
  the same residues (e.g. the same termini)

### GROMACS not found
- Ensure `gmx` is in PATH: `which gmx`
- Or specify path: `--gmx-bin /path/to/gmx`

### Memory issues
- Reduce selection (e.g., use C-alpha instead of all-atom)
- Reduce trajectory size with `gmx trjconv -skip`

## Further Analysis

### PC Movies
```bash
# Generate movie along PC1
gmx_pca_movie.py -s pca_out/fit.gro -f pca_out/fit.xtc \
    -v pca_out/eigenvectors.trr -o pc1_movie.pdb --pc 1 --nframes 100
```

## References

For more on PCA in MD:
- Amadei et al. (1993) Proteins 17:412-425
- Garcia (1992) Phys Rev Lett 68:2696-2699
- Ichiye & Karplus (1991) Proteins 11:205-217

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
- ✅ Automatic plots (scree, cumulative variance, projections over time, scatter of two PCs)
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

## Command-Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `-s, --structure` | Required | The run input file (`.tpr`): needed to make molecules whole. A GRO/PDB is accepted only with `--pbc no` (trajectory already whole and centred) |
| `-f, --trajectory` | Required | Trajectory file (XTC/TRR) |
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

### Generated Plots

Each as `.png` and `.pdf`:

- `plot_eigenvalues` - Scree plot: variance (%) of each PC (the first 25 PCs)
- `plot_cumulative_variance` - Cumulative variance (%) (the first 25 PCs)
- `plot_pc<N>_timeseries` - Projection on PC N over time, for each PC `--first` to `--last`
- `plot_pc<A>_pc<B>_scatter` - The first two projected PCs against each other, coloured by time

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

Each PCA has its own principal components: PC1 of one replica is in general another direction
than PC1 of another, so the projections of two separate PCAs can't be compared directly. To
compare replicas, do **one** PCA of all of them together, so that all frames are projected on
the same PCs:

```bash
# Concatenate the replicas (same system, same atoms); -settime/-cat options as needed
gmx trjcat -f rep1.xtc rep2.xtc -o reps.xtc -cat
gmx_pca.py -s topol.tpr -f reps.xtc -o pca_reps --last 3
# The projections then contain rep1's frames followed by rep2's
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

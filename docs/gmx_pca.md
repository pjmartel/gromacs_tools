# gmx_pca.py

Automates Principal Component Analysis (PCA) on GROMACS trajectories with automatic plotting of essential results.

## Overview

Performs complete PCA workflow from trajectory preprocessing through eigenvector analysis and visualization. Handles all GROMACS command execution and generates publication-ready plots automatically.

## Features

- ✅ Automatic index file generation from selection strings
- ✅ Trajectory centering and fitting
- ✅ Covariance matrix computation
- ✅ Eigenvector analysis and projection
- ✅ Automatic plot generation (scree, cumulative variance, PC projections, timeseries)
- ✅ Flexible atom selection using GROMACS selection syntax
- ✅ Customizable PBC and fitting options
- ✅ Batch processing support
- ✅ **Command logging** (on by default): every GROMACS command run, plus the exact invocation
  used, is saved to `<outdir>/gmx_pca_commands.sh` (override with `--commands-file`); the file from an earlier run is kept as `…_commands.1.sh`, `.2.sh`, … (highest number = most recent)

## Requirements

- GROMACS 2021+
- Python 3.7+
- numpy, matplotlib

## Quick Start

```bash
# Basic PCA on backbone atoms
python bin/gmx_pca.py -s topol.tpr -f traj.xtc --selection "Backbone"

# Analyze first 5 principal components
python bin/gmx_pca.py -s topol.tpr -f traj.xtc \
    --selection "Backbone" --first 1 --last 5
```

## Workflow

The script automatically performs these steps:

1. **Index Generation**: Creates atom selection from GROMACS selection string
2. **Trajectory Preprocessing**: makes the selected atoms whole and centres them (PBC
   treatment, which needs the run's `.tpr`), takes the first frame of that as the reference
   structure `fit.gro`, and fits the trajectory to it (removes rotation/translation). The
   reference is taken *after* the PBC treatment: taken from the raw trajectory, it would be
   split across the box whenever the protein crosses the box edge in the first frame, which
   silently corrupts the fit and the PCA
3. **Covariance Analysis**: Computes covariance matrix and eigenvectors
4. **Projection**: Projects trajectory onto principal components
5. **Visualization**: Generates plots for analysis

## Usage

### Basic Usage

```bash
python bin/gmx_pca.py -s topol.tpr -f traj.xtc
```

### Custom Selection

```bash
# Protein backbone
python bin/gmx_pca.py -s topol.tpr -f traj.xtc --selection "Backbone"

# C-alpha atoms only
python bin/gmx_pca.py -s topol.tpr -f traj.xtc --selection "name CA"

# Specific residues
python bin/gmx_pca.py -s topol.tpr -f traj.xtc \
    --selection "resid 10 to 50 and name CA"

# Binding site region
python bin/gmx_pca.py -s topol.tpr -f traj.xtc \
    --selection "resid 15 23 45 67 and sidechain"
```

### Output Directory

```bash
# Custom output directory
python bin/gmx_pca.py -s topol.tpr -f traj.xtc -o pca_results

# Overwrite existing output
python bin/gmx_pca.py -s topol.tpr -f traj.xtc -o pca_out --overwrite
```

### Eigenvector Range

```bash
# Analyze first 10 PCs
python bin/gmx_pca.py -s topol.tpr -f traj.xtc --first 1 --last 10

# Analyze only PC1
python bin/gmx_pca.py -s topol.tpr -f traj.xtc --first 1 --last 1
```

### PBC and Fitting Options

```bash
# Different PBC handling
python bin/gmx_pca.py -s topol.tpr -f traj.xtc --pbc nojump

# No PBC correction
python bin/gmx_pca.py -s topol.tpr -f traj.xtc --pbc no

# Progressive fitting (for long trajectories)
python bin/gmx_pca.py -s topol.tpr -f traj.xtc --fit progressive

# No fitting (pre-aligned trajectories)
python bin/gmx_pca.py -s topol.tpr -f traj.xtc --fit none --no-center
```

### Custom GROMACS Binary

```bash
# Use specific GROMACS version
python bin/gmx_pca.py -s topol.tpr -f traj.xtc --gmx-bin gmx_2024

# Full path to GROMACS
python bin/gmx_pca.py -s topol.tpr -f traj.xtc \
    --gmx-bin /usr/local/gromacs/2024/bin/gmx
```

### Plotting Control

```bash
# Skip plotting (only run GROMACS steps)
python bin/gmx_pca.py -s topol.tpr -f traj.xtc --no-plots

# Then use plot_xvg.py for custom plots
python bin/plot_xvg.py pca_out/projection.xvg --scatter --colormap viridis
```

## Command-Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `-s, --structure` | Required | The run input file (`.tpr`): needed to make molecules whole. A GRO/PDB is accepted only with `--pbc no` (trajectory already whole and centred) |
| `-f, --trajectory` | Required | Trajectory file (XTC/TRR) |
| `-o, --outdir` | `pca_out` | Output directory |
| `--selection` | `Backbone` | GROMACS selection string |
| `--first` | `1` | First eigenvector to analyze |
| `--last` | `2` | Last eigenvector to analyze |
| `--pbc` | `mol` | PBC handling (no/mol/res/atom/nojump) |
| `--no-center` | `False` | Disable trajectory centering |
| `--fit` | `rot+trans` | Fitting mode (rot+trans/progressive/none) |
| `--gmx-bin` | `gmx` | GROMACS binary name/path |
| `--no-plots` | `False` | Skip matplotlib plotting |
| `--overwrite` | `False` | Make the index, fitted trajectory and reference again, even when the ones in `--outdir` could be reused (see [Re-running](#re-running-in-the-same-output-folder)) |

## Output Files

### GROMACS Outputs

All in `--outdir` (GROMACS runs there, so its own logs stay there too, and it makes no
`#file.N#` backups: re-running replaces the files).

- `selection.ndx` - Index file with selected atoms
- `fit.xtc` - Whole, centred and fitted trajectory of the selected atoms
- `fit.gro` - Reference structure: the first frame, whole and centred like `fit.xtc`
- `preprocessing.json` - The input files and options the three files above were made with
- `eigenvalues.xvg` - Eigenvalues from covariance analysis
- `eigenvectors.trr` - Eigenvectors (for projections/movies)
- `average.pdb` - Average structure
- `projection_pc1.xvg`, `projection_pc2.xvg` - Trajectory projected onto PC1 and PC2
- `rmsf.xvg` - RMSF per atom along eigenvectors `--first` to `--last`
- `covar.log` - Log of `gmx covar`; `logs/` - output of every GROMACS step
- `gmx_pca_commands.sh` - The GROMACS commands run (see Features)

### Generated Plots

Each as `.png` and `.pdf`:

- `plot_eigenvalues` - Scree plot (eigenvalues ranked)
- `plot_cumulative_variance` - Cumulative variance explained
- `plot_pc1_pc2_scatter` - PC1 vs PC2 scatter plot, coloured by time
- `plot_pc1_timeseries`, `plot_pc2_timeseries` - PC1 and PC2 projections over time

### Re-running in the Same Output Folder

The index, fitted trajectory and reference structure of an earlier run are reused only if
they were made from the same structure and trajectory files (same path, size and modification
time) with the same `--selection`, `--pbc`, `--fit` and centring, as recorded in
`preprocessing.json`. Otherwise they are made again, with a note saying what changed (e.g.
"changed since: selection"). `--overwrite` makes them again in any case. The covariance
analysis, projections and plots are always made again.

## Generated Plots Explained

### Scree Plot
Shows eigenvalues in descending order. Helps identify how many PCs are significant ("elbow method").

### Cumulative Variance
Shows cumulative percentage of variance explained. Typically 80-90% captured in first few PCs.

### PC1 vs PC2 Scatter
2D projection showing conformational sampling. Clusters indicate metastable states.

### PC1 Timeseries
Shows how conformations evolve along PC1 over time. Useful for identifying transitions.

## Example Workflows

### Basic Backbone PCA

```bash
# 1. Run PCA
python bin/gmx_pca.py -s topol.tpr -f traj.xtc

# 2. Outputs in pca_out/:
#    - fit.gro, eigenvectors.trr (for movies)
#    - projection.xvg (for custom plots)
#    - plot_*.png (automatic visualizations)
```

### Focused Region Analysis

```bash
# Analyze specific binding pocket
python bin/gmx_pca.py -s topol.tpr -f traj.xtc \
    --selection "resid 45 46 47 89 90 91 and sidechain" \
    -o binding_site_pca --first 1 --last 5
```

### Multiple Replicas Comparison

```bash
# Replica 1
python bin/gmx_pca.py -s topol.tpr -f rep1.xtc -o pca_rep1

# Replica 2
python bin/gmx_pca.py -s topol.tpr -f rep2.xtc -o pca_rep2

# Compare projections
python bin/plot_xvg.py pca_rep1/projection.xvg pca_rep2/projection.xvg \
    --xy-correlation --multi --legends "Replica 1" "Replica 2" \
    --title "PC1 vs PC2: Replicas Comparison" --aspect equal
```

### Pre-processed Trajectory

```bash
# If trajectory already aligned
python bin/gmx_pca.py -s aligned.gro -f aligned.xtc \
    --fit none --no-center --pbc no
```

## GROMACS Selection Syntax

Common selections:
- `Backbone` - Protein backbone atoms (N, CA, C, O)
- `C-alpha` - Only C-alpha atoms
- `Protein` - All protein atoms
- `name CA` - Atoms named CA
- `resid 10 to 50` - Residues 10-50
- `resname TRP TYR PHE` - Aromatic residues
- `sidechain and not hydrogen` - Heavy sidechain atoms
- `resid 45 and within 0.5 of resname LIG` - Atoms near ligand

Combine with `and`, `or`, `not`:
```bash
--selection "resid 10 to 50 and name CA"
--selection "(resid 10 to 30 or resid 60 to 80) and backbone"
```

## Tips & Best Practices

1. **Atom Selection**:
   - C-alpha only: Faster, captures global motion
   - Backbone: More detailed, standard choice
   - All-atom: Most detailed but slow and noisy

2. **Eigenvector Range**:
   - First 2-3 PCs usually capture >70% variance
   - Analyze more (5-10) for complex systems
   - Check cumulative variance plot

3. **PBC Handling**:
   - `mol`: Good default for proteins
   - `nojump`: Better for membrane systems
   - `no`: Only for pre-processed trajectories

4. **Fitting**:
   - `rot+trans`: Standard for equilibrated trajectories
   - `progressive`: Better for very long trajectories
   - `none`: For pre-aligned data

5. **Output Organization**:
   - Use descriptive output directory names
   - Keep different selections in separate directories
   - Save `projection.xvg` for custom analysis

6. **Performance**:
   - C-alpha selection is much faster than all-atom
   - Large trajectories benefit from compressed XTC format
   - Consider reducing trajectory size with `gmx trjconv -skip`

## Troubleshooting

### "Selection is empty"
- Check your selection string syntax
- Verify atoms exist in your structure: `gmx select -s topol.tpr -select "your_selection"`

### "--pbc mol needs a .tpr"
- Give the run's `.tpr` as `-s`: making molecules whole needs the topology. With a `.gro` or
  `.pdb`, GROMACS can only warn and leaves molecules broken, which would corrupt the PCA.
- If the trajectory is already whole and centred, use `--pbc no` (see Pre-processed Trajectory)

### GROMACS not found
- Ensure `gmx` is in PATH: `which gmx`
- Or specify path: `--gmx-bin /path/to/gmx`

### Trajectory not aligned
- Don't use `--no-center` or `--fit none` unless pre-aligned
- Default settings handle alignment automatically

### Memory issues
- Reduce selection (e.g., use C-alpha instead of all-atom)
- Reduce trajectory size with `gmx trjconv -skip`

## Further Analysis

After running PCA, use these for additional analysis:

### Custom Projections
```bash
# More detailed PC1 vs PC2 plot
python bin/plot_xvg.py pca_out/projection.xvg \
    --scatter --colormap plasma --aspect equal \
    --title "Conformational Space" --markersize 2
```

### PC Movies
```bash
# Generate movie along PC1
python bin/gmx_pca_movie.py -s pca_out/fit.gro \
    -v pca_out/eigenvectors.trr -o pc1_movie.pdb --pc 1 --nframes 100
```

### Multiple PC Correlation
```bash
# Extract individual PC columns from projection.xvg
# Then correlate PC1 from different systems
python bin/plot_xvg.py system1_proj.xvg system2_proj.xvg \
    --xy-correlation --legends "System 1" "System 2"
```

## References

For more on PCA in MD:
- Amadei et al. (1993) Proteins 17:412-425
- Garcia (1992) Phys Rev Lett 68:2696-2699
- Ichiye & Karplus (1991) Proteins 11:205-217

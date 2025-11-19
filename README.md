# GROMACS Tools

A collection of Python tools for GROMACS molecular dynamics simulations.

## Scripts

### gromacs_prepare.py

A comprehensive preparation script for GROMACS molecular dynamics simulations. This script automates the entire setup workflow from initial structure to equilibrated system.

#### Features

- Automatic protein/ligand structure preparation
- Topology generation
- System solvation and ionization
- Energy minimization
- NVT and NPT equilibration
- Production MD setup
- Comprehensive error handling and logging

#### Requirements

- Python 3.x
- GROMACS (tested with version 2021 or later)
- Required Python packages: (add your dependencies here)

#### Usage

```bash
python gromacs_prepare.py [options]
```

For detailed options and configuration, see the script's help documentation.

### plot_xvg.py

A versatile command-line tool for plotting GROMACS XVG files with Matplotlib. Supports single or multiple XVG files, various plot modes (xy plots, scatter plots, histograms, correlations), moving averages, and extensive customization options.

#### Features

- **Multiple plot modes:**
  - Standard XY plots with customizable styles (dots, lines, lines+dots)
  - Scatter plots with time/frame-ordered color gradient (ideal for PCA projections)
  - Histogram mode for distribution analysis
  - XY correlation mode to plot two datasets against each other
- Plot single or multiple XVG files on the same axes
- Apply moving average smoothing to time series data
- Custom legends, titles, axis labels, and colormaps
- Backend selection for different display environments (Qt5Agg, TkAgg, etc.)
- Save plots to various formats (PNG, PDF, SVG, EPS, etc.)
- Automatically parses XVG metadata (titles, labels, legends)

#### Requirements

- Python 3.x with `numpy` and `matplotlib`

#### Usage

**Basic XY plots:**
```bash
# Plot a single XVG file (default: dots)
python plot_xvg.py energy.xvg

# Plot with line style
python plot_xvg.py energy.xvg --style lines

# Plot with lines and dots
python plot_xvg.py rmsd.xvg --style lines+dots
```

**Moving average:**
```bash
# Apply moving average smoothing
python plot_xvg.py energy.xvg --moving-avg --window 50 --style lines
```

**Scatter plots (colored by time/frame order):**
```bash
# Scatter plot colored by frame/time order (like PCA plots)
python plot_xvg.py pc_projection.xvg --scatter --colormap viridis

# Scatter with custom colormap
python plot_xvg.py pc1_vs_pc2.xvg --scatter --colormap plasma \\
  --title "PC1 vs PC2" --output pca.pdf
```

**Histogram mode:**
```bash
# Plot histogram of second column (distribution analysis)
python plot_xvg.py rmsd.xvg --histogram --bins 50

# Histogram with custom labels
python plot_xvg.py rmsd.xvg --histogram --bins 100 \\
  --title "RMSD Distribution" --xlabel "RMSD (nm)" --output rmsd_dist.png

# Compare distributions from multiple files
python plot_xvg.py run1.xvg run2.xvg run3.xvg --multi --histogram \\
  --bins 40 --legends "Wild Type" "Mutant A" "Mutant B"
```

**XY correlation (plot two datasets against each other):**
```bash
# 2D: Plot second column of file1 vs second column of file2
python plot_xvg.py rmsd_protein.xvg rmsd_ligand.xvg --xy-correlation

# 2D correlation with scatter coloring by time
python plot_xvg.py pc1.xvg pc2.xvg --xy-correlation --scatter \\
  --title "PC1 vs PC2" --colormap coolwarm --output pc_correlation.png

# Square aspect ratio for proper correlation visualization
python plot_xvg.py pc1.xvg pc2.xvg --xy-correlation --aspect equal

# 3D: Plot three datasets as x, y, z coordinates
python plot_xvg.py pc1.xvg pc2.xvg pc3.xvg --xy-correlation --scatter \\
  --title "PC1 vs PC2 vs PC3" --colormap viridis --output pc_3d.png

# 3D trajectory path
python plot_xvg.py x.xvg y.xvg z.xvg --xy-correlation --style lines \\
  --aspect equal --title "3D Trajectory"

# Multiple 2D correlations on same plot (different color per dataset)
python plot_xvg.py rmsd1a.xvg rmsd2a.xvg rmsd1b.xvg rmsd2b.xvg \\
  --xy-correlation --multi --legends "System A" "System B" \\
  --title "RMSD Correlations"

# Multiple 3D correlations overlaid
python plot_xvg.py pc1a.xvg pc2a.xvg pc3a.xvg pc1b.xvg pc2b.xvg pc3b.xvg \\
  --xy-correlation --multi --style lines \\
  --legends "Wild Type" "Mutant" --title "PC Space Comparison"

# Will error if files have different number of rows (validation included)
```

**Marker size and aspect ratio control:**
```bash
# Adjust marker/dot size for dense data
python plot_xvg.py trajectory.xvg --markersize 1 --style dots

# Larger markers for visibility
python plot_xvg.py sparse_data.xvg --markersize 6

# Square plot for correlation analysis
python plot_xvg.py data.xvg --aspect equal

# Custom aspect ratio (width/height)
python plot_xvg.py data.xvg --aspect 1.5
```

**Multiple files on same axes:**
```bash
# Plot multiple XVG files overlaid
python plot_xvg.py file1.xvg file2.xvg file3.xvg --multi

# With custom legends and style
python plot_xvg.py rmsd1.xvg rmsd2.xvg --multi \\
  --legends "System A" "System B" --style lines+dots
```

**Output and customization:**
```bash
# Save plot to file with high DPI
python plot_xvg.py energy.xvg --output energy_plot.png --dpi 300

# Custom title and labels
python plot_xvg.py rmsd.xvg --title "Backbone RMSD" \\
  --xlabel "Time (ns)" --ylabel "RMSD (nm)"

# Specify matplotlib backend (useful for headless/remote systems)
python plot_xvg.py --backend Qt5Agg energy.xvg
python plot_xvg.py --backend Agg energy.xvg --output energy.png
```

**Trajectory slicing (view specific portions):**
```bash
# Plot only frames 100-500 (exclude equilibration and late trajectory)
python plot_xvg.py rmsd.xvg --start 100 --end 500 --style lines

# Plot from frame 1000 onwards (skip equilibration period)
python plot_xvg.py energy.xvg --start 1000 --style lines

# Plot up to frame 2000 (early trajectory analysis)
python plot_xvg.py trajectory.xvg --end 2000 --scatter

# Combine slicing with other options
python plot_xvg.py rmsd.xvg --start 500 --end 2000 --moving-avg --window 50
```

#### Key Options

**Input and mode:**
- `files`: One or more XVG files to plot (positional argument)
- `--style, -s`: Plot style - `dots` (default), `lines`, or `lines+dots`
- `--scatter`: Enable scatter mode with color gradient by frame/time order
- `--histogram, --hist`: Plot histogram of second column (y-values)
- `--xy-correlation, --xycorr`: Plot y-values from files as coordinates
  - 2 files: 2D correlation plot (x from file1, y from file2)
  - 3 files: 3D correlation plot (x, y, z from file1, file2, file3)
  - With `--multi`: overlay multiple correlations (pairs or triplets) with different colors

**Customization:**
- `--colormap, --cmap`: Colormap for scatter mode (default: `viridis`). Options: `viridis`, `plasma`, `inferno`, `magma`, `coolwarm`, `rainbow`
- `--bins`: Number of bins for histogram mode (default: 50)
- `--markersize, --ms`: Size of markers/dots (default: 3.0). Smaller for dense data (1-2), larger for visibility (4-6)
- `--aspect`: Aspect ratio of plot axes
  - `equal`: 1:1 square aspect (perfect for correlations)
  - `auto`: Default matplotlib behavior
  - Numeric (e.g., `1.5`): Custom width/height ratio
- `--moving-avg, --ma`: Apply moving average filter (only for single dataset in xy mode)
- `--window, -w`: Window size for moving average (default: 10)
- `--title, -t`: Custom plot title (overrides XVG title)
- `--xlabel`: Custom x-axis label (overrides XVG label)
- `--ylabel`: Custom y-axis label (overrides XVG label)

**Data slicing:**
- `--start`: First data row to plot (0-indexed). Allows viewing specific trajectory portions by skipping initial frames
- `--end`: Last data row to plot (exclusive, 0-indexed). Allows viewing specific trajectory portions by limiting to early frames
- These options work with all plot modes (single, multi, correlation, histogram, etc.)
- Useful for excluding equilibration, focusing on specific time windows, or comparing trajectory segments

**Multiple files:**
- `--multi, -m`: Plot multiple XVG files on the same axes
  - Standard mode: each file is plotted as a separate series
  - With `--xy-correlation`: expects pairs (2D) or triplets (3D) of files
- `--legends, -l`: Custom legend labels
  - Standard/histogram mode: must match number of files
  - XY correlation with `--multi`: must match number of correlation sets (pairs/triplets)

**Output:**
- `--output, -o`: Save plot to file instead of displaying interactively
- `--figsize WIDTH HEIGHT`: Figure size in inches (default: 10 6)
- `--dpi`: Resolution for saved figures (default: 300)
- `--backend, --mpl-backend`: Matplotlib backend (e.g., Qt5Agg, TkAgg, Agg)

#### Outputs

- Interactive plot window (if `--output` not specified)
- Saved plot file in specified format (PNG, PDF, SVG, EPS, etc.)

#### Use Cases

- **Time series analysis**: Plot energy, RMSD, RMSF, distances over time
- **PCA visualization**: 
  - 2D: Scatter plots with time-ordered coloring for PC1 vs PC2
  - 3D: Interactive 3D plots for PC1 vs PC2 vs PC3 conformational space
- **Distribution analysis**: Histograms of RMSD, energy, angles, distances
- **Single correlation studies**: 
  - 2D: Compare two observables (e.g., protein vs ligand RMSD)
  - 3D: Explore three observables simultaneously (e.g., three distances, PC1-PC2-PC3)
- **Multiple correlation studies**: Overlay correlation plots from different systems
  - Compare WT vs mutant in PC space
  - Multiple replicas or conditions on same correlation plot
  - Different color for each dataset
- **3D trajectories**: Visualize center-of-mass motion or particle tracking in 3D space
- **Multi-system comparison**: Overlay multiple datasets with custom legends
- **Trajectory slicing**: Analyze specific simulation segments
  - Exclude equilibration period using `--start`
  - Focus on early/late trajectory portions using `--start` and `--end`
  - Compare different time windows from same trajectory
  - Production run analysis (skip equilibration frames)
- **Publication figures**: Control marker size, aspect ratio, DPI for journal-quality plots

#### Notes

- **Default style is `dots`** which works well for most GROMACS data
- **Scatter mode** (`--scatter`) shows temporal evolution via color gradient - ideal for PCA projections
- **Histogram mode** analyzes distributions of the second column (y-values)
- **XY correlation mode** supports both 2D (2 files) and 3D (3 files) plots:
  - Single correlation: 2 or 3 files for one correlation plot
  - Multiple correlations with `--multi`: 4, 6, 8... files (pairs) or 6, 9, 12... files (triplets)
  - Each correlation set gets a different color (no time-based coloring in multi mode)
  - Validates that all files have matching row counts (errors if different)
  - 3D plots are interactive: rotate, zoom, and pan to explore data
  - Use `--aspect equal` for proper geometric scaling in correlation/PCA plots
- **Marker size** (`--markersize`):
  - Dense data with many points: 1-2
  - Default/medium: 3 (current default)
  - Sparse data or presentations: 4-6
- **Aspect ratio** (`--aspect`):
  - Use `equal` for correlation plots where x and y have the same units
  - Critical for PC1 vs PC2 plots to see true geometric relationships
  - For 3D plots with `equal`, creates cube-like scaling
- The script automatically parses XVG metadata (titles, axis labels, legends)
- **Data slicing** (`--start` and `--end`):
  - Row indices are 0-based (first data row is 0)
  - `--end` is exclusive (like Python slicing)
  - Works with all modes: single, multi, correlation, histogram, scatter
  - Applied after parsing, before plotting (preserves file integrity)
  - Useful for equilibration removal, time window analysis, or reducing plot density
- Moving average only applies to single datasets in standard xy plot mode
- When using `--legends`:
  - Standard mode: number of legends must match number of files
  - Multi-correlation mode: number of legends must match number of correlation sets
- Use `--backend` to control display behavior (useful for remote/headless systems)
- Scatter mode adds a colorbar showing frame/time progression (single correlations only)

### gromacs_pca.py

Automates Principal Component Analysis (PCA) on GROMACS trajectories and produces common plots (scree, cumulative variance, PC1–PC2 scatter, PC1 timeseries).

#### Requirements

- GROMACS 2021+ in PATH (binary `gmx`) or specify with `--gmx-bin`
- Python 3.x with `numpy` and `matplotlib` (see `requirements.txt`)

#### Typical workflow performed

1. Build an atom selection index using GROMACS selection syntax (default: `Backbone`).
2. Center and fit the trajectory on the selection (`gmx trjconv`).
3. Compute covariance matrix and eigenvectors (`gmx covar`).
4. Analyze eigenvectors and project the trajectory (`gmx anaeig`).
5. Parse `.xvg` outputs and generate plots.

#### Usage

```bash
python gromacs_pca.py -s topol.tpr -f traj.xtc -o pca_out \
	--selection "Backbone" --first 1 --last 2
```

Key options:

- `-s, --structure`: Structure/topology file (TPR recommended; GRO/PDB also works)
- `-f, --trajectory`: Trajectory file (XTC/TRR)
- `-o, --outdir`: Output directory (default: `pca_out`)
- `--selection`: GROMACS selection string (default: `Backbone`)
- `--first/--last`: Eigenvector range to analyze (default: 1–2)
- `--pbc`: PBC handling for `trjconv` (`no|mol|res|atom`; default: `mol`)
- `--no-center`: Disable centering in `trjconv`
- `--fit`: Fitting mode for `trjconv` (`rot+trans|progressive|none`; default: `rot+trans`)
- `--gmx-bin`: GROMACS binary to use (e.g., `gmx` or full path)
- `--no-plots`: Skip Matplotlib plotting (still runs GROMACS steps)
- `--overwrite`: Overwrite existing outputs

Outputs (in the chosen `outdir`):

- `selection.ndx` — index generated from your selection
- `fit.xtc` — centered/fitted trajectory
- `eigenvalues.xvg`, `eigenvectors.trr`, `average.pdb` — from `gmx covar`
- `projection.xvg`, `eig.xvg`, `rmsf.xvg`, `comp.xvg` — from `gmx anaeig`
- `plot_*.png/.pdf` — Matplotlib plots (scree, cumulative, PC1–PC2, PC1 timeseries)

Notes:

- The script feeds group selections programmatically using the generated `selection.ndx`, which contains a single group. If you want different atoms, change `--selection` (GROMACS selection language).
- Ensure GROMACS is installed and callable (try `gmx --version`).

### gromacs_pca_movie.py

Creates trajectory movies showing conformational changes along principal components. This generates interpolated structures that can be animated in molecular visualization software.

#### Requirements

- GROMACS 2021+ in PATH (binary `gmx`) or specify with `--gmx-bin`
- Outputs from `gromacs_pca.py` (or equivalent PCA workflow): `fit.gro` and `eigenvectors.trr`

#### What it does

Uses `gmx anaeig -extr` to interpolate structures along a principal component, creating a trajectory that shows the motion captured by that PC. This is useful for:
- Visualizing what conformational change a PC represents
- Creating animations for presentations
- Understanding the dominant motions in your system

#### Usage

```bash
# Generate a 50-frame movie along PC1
python gromacs_pca_movie.py -s pca_out/fit.gro -v pca_out/eigenvectors.trr \\
  -o pc1_movie.pdb --pc 1 --nframes 50 --extreme 2.0

# Generate only the extreme structures (min and max)
python gromacs_pca_movie.py -s pca_out/fit.gro -v pca_out/eigenvectors.trr \\
  -o extremes.pdb --pc 1 --extremes-only --extreme 2.0
```

Key options:

- `-s, --structure`: Structure file matching PCA atoms (use `fit.gro` from `gromacs_pca.py`)
- `-v, --eigenvec`: Eigenvectors file (use `eigenvectors.trr` from `gromacs_pca.py`)
- `-o, --output`: Output trajectory file (PDB/XTC/TRR)
- `--pc`: Principal component number (default: 1)
- `--nframes`: Number of interpolated frames (default: 50)
- `--extreme`: Extent along eigenvector in nm (default: 2.0)
- `--extremes-only`: Generate only min/max structures instead of full movie
- `--outdir`: Output directory for logs
- `--gmx-bin`: GROMACS binary to use

Output:
- A trajectory file showing interpolated motion along the selected PC
- Can be opened in VMD, PyMOL, ChimeraX, etc., and played as an animation

Workflow example:
```bash
# 1. Run PCA analysis
python gromacs_pca.py -s topol.tpr -f traj.xtc -o pca_out

# 2. Generate movie for PC1
python gromacs_pca_movie.py -s pca_out/fit.gro -v pca_out/eigenvectors.trr \\
  -o pc1_movie.pdb --pc 1 --nframes 100

# 3. Visualize in VMD
vmd pc1_movie.pdb

# 4. In VMD: Graphics → Representations, then play the trajectory
```

## Installation

Clone this repository:

```bash
git clone https://github.com/yourusername/gromacs_tools.git
cd gromacs_tools
```

## License

[Add your license here]

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Author

[Paulo Martel]

## Acknowledgments

Built for streamlining GROMACS MD simulation workflows.

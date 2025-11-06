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

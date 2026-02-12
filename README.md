# GROMACS Tools

A collection of Python and Bash tools for automating GROMACS molecular dynamics simulations, analysis, and visualization.

## 🚀 Quick Start

```bash
git clone https://github.com/yourusername/gromacs_tools.git
cd gromacs_tools
pip install -r requirements.txt

# Add to PATH (optional)
export PATH="$PATH:$(pwd)/bin"
```

## 📦 Tools

| Tool | Purpose | Documentation |
|------|---------|---------------|
| **gromacs_prepare.py** | Automate system preparation pipeline (pdb2gmx, solvate, ions) | [docs/gromacs_prepare.md](docs/gromacs_prepare.md) |
| **gromacs_pca.py** | Principal component analysis workflow with automatic plotting | [docs/gromacs_pca.md](docs/gromacs_pca.md) |
| **gromacs_pca_movie.py** | Generate trajectory movies along principal components | [docs/gromacs_pca_movie.md](docs/gromacs_pca_movie.md) |
| **plot_xvg.py** | Versatile XVG plotting tool with extensive customization | [docs/plot_xvg.md](docs/plot_xvg.md) |
| **gmx_continue_grompp.sh** | Intelligent MD simulation continuation with crash recovery | [docs/md_continuation.md](docs/md_continuation.md) |

## 💡 Quick Usage Examples

### System Preparation
```bash
# Preview preparation workflow
python bin/gromacs_prepare.py protein.pdb --forcefield amber14sb --dry-run

# Execute preparation
python bin/gromacs_prepare.py protein.pdb --forcefield charmm36
```

### PCA Analysis
```bash
# Full PCA workflow with automatic plotting
python bin/gromacs_pca.py -s topol.tpr -f traj.xtc --selection "Backbone"

# Generate PC motion movie
python bin/gromacs_pca_movie.py fit.gro eigenvec.trr --pc 1 --frames 50
```

### Plotting XVG Files
```bash
# Simple time series
python bin/plot_xvg.py rmsd.xvg --style lines

# Multiple files with legends
python bin/plot_xvg.py file1.xvg file2.xvg --multi --legends "WT" "Mutant"

# Scatter plot with moving average
python bin/plot_xvg.py energy.xvg --scatter --moving-avg --window 50

# Plot specific columns with custom legends
python bin/plot_xvg.py multicolumn.xvg --columns 1 3 --legends "RMSD" "RMSF"
```

### MD Continuation
```bash
# Run 500 ns in 100 ns segments
bash bin/gmx_continue_grompp.sh md 0 0 500 100 production.mdp npt

# With custom title
bash bin/gmx_continue_grompp.sh md 0 0 500 100 production.mdp npt 0.002 "TRP_cage"

# Automatic crash recovery: just re-run the same command
```

## 📋 Requirements

- **Python**: 3.7+
- **GROMACS**: 2021+ (tested with 2021-2025)
- **Python packages**: numpy, matplotlib (see `requirements.txt`)

## 📂 Repository Structure

```
gromacs_tools/
├── bin/               # Executable scripts
├── docs/              # Detailed documentation per tool
├── mdp_templates/     # Force field-specific MDP parameter files
│   ├── charmm36/     # CHARMM36 templates
│   ├── amber14sb/    # AMBER14SB templates
│   └── plumed/       # PLUMED configuration files
└── examples/          # Usage examples and workflows
```

## 🎯 Key Features

- **Automation**: Reduce manual GROMACS command execution
- **Reproducibility**: Automatic command logging to `commands.sh`
- **Safety**: Dry-run modes, crash recovery, validation checks
- **Flexibility**: Extensive customization options
- **Documentation**: Comprehensive guides and examples

## 📖 MDP Templates

Production-ready MDP parameter files optimized for specific force fields:
- **CHARMM36**: 1.2 nm cutoffs, force-switch, all-bonds constraints
- **AMBER14SB**: 1.0 nm cutoffs, potential-shift, h-bonds constraints
- See [mdp_templates/README.md](mdp_templates/README.md) for details

## 🤝 Contributing

Contributions welcome! Feel free to open issues or submit pull requests.

## 📜 License

MIT License - see [LICENSE](LICENSE) file for details.

## 🔗 Resources

- [GROMACS Documentation](https://manual.gromacs.org/)
- [GROMACS Tutorials](http://www.mdtutorials.com/gmx/)

# Reference Molecular Systems

This directory contains pre-equilibrated molecular systems with topology files and starting configurations for testing and demonstration purposes.

## Directory Structure

Each system subdirectory should contain:

### Essential Files (tracked in git)
- `*.top` - GROMACS topology file
- `*.itp` - Include topology files (force field parameters, molecule definitions)
- `*.gro` or `*.pdb` - Structure/coordinate files
- `*.mdp` - MDP parameter files for different simulation stages
  - `em.mdp` - Energy minimization
  - `nvt.mdp` - NVT equilibration
  - `npt.mdp` - NPT equilibration
  - `md.mdp` - Production MD
- `*.tpr` - Pre-compiled run input files (small, useful for quick restarts)
- `*.cpt` - Checkpoint files from equilibration (for seamless continuation)
- `plumed.dat` - PLUMED input files (if using enhanced sampling)
- `README.md` - System-specific documentation

### Excluded Files (not tracked)
- `*.xtc`, `*.trr` - Trajectory files (too large)
- `*.edr` - Energy files (large binary files)
- `*.log` - Log files (can be regenerated)
- `#*#`, `*.bak` - Backup files

## Available Systems

### Small Molecules
- **cyclohexane/** - Simple cyclic hydrocarbon, useful for testing liquid simulations
- **ethane/** - Minimal organic molecule, fast test system
- **aspirin/** - Small drug molecule with rings

### Peptides
- **kyotorphin/** - Dipeptide (Tyr-Arg), useful for peptide simulation testing
- **trp_cage/** - Small fast-folding protein (20 residues), protein simulation benchmark

## Usage Examples

### Continue from equilibrated state
```bash
cd systems/trp_cage/
gmx_continue_extend.sh 10000  # Extend by 10 ns using auto-detected files
```

### Start new production run from checkpoint
```bash
cd systems/aspirin/
gmx_continue_grompp.sh md 0 0 100 10 --template md.mdp --initial npt
```

### Test with PLUMED
```bash
cd systems/kyotorphin/
gmx_continue_extend.sh 5000 --plumed plumed.dat
```

## Adding New Systems

When adding a new system:

1. Create subdirectory: `systems/system_name/`
2. Add topology and structure files
3. Include MDP files for all simulation stages
4. Run equilibration (EM → NVT → NPT)
5. Save final checkpoint and TPR from NPT
6. Add system-specific README.md with:
   - System description
   - Force field used
   - Box size and composition
   - Equilibration protocol
   - Any special considerations

## Notes

- All systems assume GROMACS 2021 or later
- Checkpoint files enable seamless continuation
- TPR files allow immediate restart without rerunning grompp
- Keep files compressed if they're large but still needed (`.gro.gz`, `.tpr.gz`)

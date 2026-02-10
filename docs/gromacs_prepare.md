# gromacs_prepare.py

Comprehensive automation script for GROMACS system preparation workflow.

## Overview

Automates the complete setup pipeline from initial PDB structure to a fully prepared, solvated, and neutralized system ready for equilibration.

## Features

- ✅ Automatic topology generation (`pdb2gmx`)
- ✅ Box definition (`editconf`)
- ✅ System solvation (`solvate`)
- ✅ Ion addition for neutralization (`genion`)
- ✅ Comprehensive error handling
- ✅ Command logging to `commands.sh` for reproducibility
- ✅ **Dry-run mode**: Preview workflow without execution

## Requirements

- Python 3.7+
- GROMACS 2021+
- Input: PDB structure file

## Usage

### Basic Usage

```bash
# Full preparation pipeline
python bin/gromacs_prepare.py protein.pdb
```

### Force Field Selection

```bash
# AMBER14SB force field
python bin/gromacs_prepare.py protein.pdb --forcefield amber14sb

# CHARMM36 force field
python bin/gromacs_prepare.py protein.pdb --forcefield charmm36

# List available force fields
gmx pdb2gmx -h  # See force field options
```

### Water Models

```bash
# TIP3P water (default)
python bin/gromacs_prepare.py protein.pdb --water-model tip3p

# TIP4P water
python bin/gromacs_prepare.py protein.pdb --water-model tip4p

# SPC/E water
python bin/gromacs_prepare.py protein.pdb --water-model spce
```

### Box Configuration

```bash
# Custom box distance (1.5 nm from protein)
python bin/gromacs_prepare.py protein.pdb --box-distance 1.5

# Dodecahedral box (fewer water molecules)
python bin/gromacs_prepare.py protein.pdb --box-type dodecahedron

# Cubic box (default)
python bin/gromacs_prepare.py protein.pdb --box-type cubic
```

### Ion Selection

```bash
# Use potassium and chloride ions
python bin/gromacs_prepare.py protein.pdb --cation K --anion CL

# Use sodium and chloride (default)
python bin/gromacs_prepare.py protein.pdb --cation NA --anion CL
```

### Dry-Run Mode

```bash
# Generate commands.sh without executing
python bin/gromacs_prepare.py protein.pdb --dry-run

# Review the generated commands
cat commands.sh

# Execute when ready
bash commands.sh
```

### Complete Example

```bash
python bin/gromacs_prepare.py protein.pdb \
    --forcefield amber14sb \
    --water-model tip3p \
    --box-distance 1.2 \
    --box-type dodecahedron \
    --cation K \
    --anion CL \
    --dry-run
```

## Command-Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `input_pdb` | Required | Input PDB structure file |
| `--forcefield` | `amber14sb` | Force field selection |
| `--water-model` | `tip3p` | Water model selection |
| `--box-type` | `cubic` | Simulation box type |
| `--box-distance` | `1.0` | Distance from protein to box edge (nm) |
| `--cation` | `NA` | Cation type for neutralization |
| `--anion` | `CL` | Anion type for neutralization |
| `--dry-run` | `False` | Generate commands without executing |

## Output Files

- `topol.top` - System topology
- `processed.gro` - Processed structure
- `newbox.gro` - Structure with defined box
- `solvated.gro` - Solvated system
- `ionized.gro` - Final system with ions
- `commands.sh` - Log of all executed commands

## Typical Workflow

```bash
# 1. Prepare system
python bin/gromacs_prepare.py protein.pdb --forcefield amber14sb

# 2. Energy minimization
gmx grompp -f mdp_templates/amber14sb/minimization/em_steep.mdp \
           -c ionized.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em

# 3. NVT equilibration
gmx grompp -f mdp_templates/amber14sb/equilibration/nvt_equil.mdp \
           -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

# 4. NPT equilibration
gmx grompp -f mdp_templates/amber14sb/equilibration/npt_equil.mdp \
           -c nvt.gro -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

# 5. Production run
bash bin/gmx_continue_grompp.sh md 0 0 500 100 \
     mdp_templates/amber14sb/production/md_npt.mdp npt
```

## Error Handling

The script includes comprehensive error checking:
- Validates GROMACS installation
- Checks input file existence
- Monitors command execution
- Reports clear error messages
- Preserves intermediate files on failure

## Tips

- Always use `--dry-run` first to preview the workflow
- Check `commands.sh` for reproducibility
- Review force field compatibility with your system
- Adjust box distance based on simulation length
- Consider dodecahedral boxes for better efficiency

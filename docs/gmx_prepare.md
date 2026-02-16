# gmx_prepare.py

Comprehensive automation script for GROMACS system preparation workflow.

## Overview

Automates the complete setup pipeline from initial PDB structure to a fully prepared, solvated, and neutralized system ready for equilibration.

## Features

- ✅ Automatic topology generation (`pdb2gmx`)
- ✅ **Skip pdb2gmx**: Use existing topology/structure files to start from step 2
- ✅ Box definition (`editconf`)
- ✅ System solvation (`solvate`)
- ✅ Ion addition for neutralization (`genion`)
- ✅ **Input validation**: PDB file structure and force field/water compatibility
- ✅ **Progress indicators**: Clear [1/5], [2/5]... step tracking (or [1/4], [2/4]... when skipping pdb2gmx)
- ✅ **Professional logging**: Python logging module with INFO/WARNING/ERROR levels
- ✅ **System summary**: Post-pipeline statistics (atoms, box volume, composition)
- ✅ **Enhanced error hints**: Step-specific troubleshooting guidance
- ✅ Command logging to `commands.sh` for reproducibility
- ✅ **Dry-run mode**: Preview workflow with configuration summary
- ✅ **Verbose mode**: Detailed debug output with `--verbose` flag

## Requirements

- Python 3.7+
- GROMACS 2021+
- Input: PDB structure file

## Usage

### Basic Usage

```bash
# Full preparation pipeline
python bin/gmx_prepare.py protein.pdb
```

### Force Field Selection

```bash
# AMBER14SB force field
python bin/gmx_prepare.py protein.pdb --forcefield amber14sb

# CHARMM36 force field
python bin/gmx_prepare.py protein.pdb --forcefield charmm36

# List available force fields
gmx pdb2gmx -h  # See force field options
```

### Water Models

```bash
# TIP3P water (default)
python bin/gmx_prepare.py protein.pdb --water-model tip3p

# TIP4P water
python bin/gmx_prepare.py protein.pdb --water-model tip4p

# SPC/E water
python bin/gmx_prepare.py protein.pdb --water-model spce
```

### Box Configuration

```bash
# Custom box distance (1.5 nm from protein)
python bin/gmx_prepare.py protein.pdb --box-distance 1.5

# Dodecahedral box (fewer water molecules)
python bin/gmx_prepare.py protein.pdb --box-type dodecahedron

# Cubic box (default)
python bin/gmx_prepare.py protein.pdb --box-type cubic
```

### Ion Selection

```bash
# Use potassium and chloride ions
python bin/gmx_prepare.py protein.pdb --cation K --anion CL

# Use sodium and chloride (default)
python bin/gmx_prepare.py protein.pdb --cation NA --anion CL
```

### Handling PDB Files with Hydrogens

```bash
# If PDB already contains hydrogen atoms, use --ignore-hydrogens
# This prevents pdb2gmx from crashing when hydrogens are present
python bin/gmx_prepare.py protein.pdb --ignore-hydrogens

# Combine with other options
python bin/gmx_prepare.py protein.pdb --ignore-hydrogens --forcefield amber99sb-ildn
```

### Skip pdb2gmx (Use Existing Topology/Structure)

```bash
# If you already have a topology and structure file, skip pdb2gmx
python bin/gmx_prepare.py --topology protein.top --structure protein.gro

# This starts the pipeline from step 2 (editconf)
# Useful when:
#  - You've manually created/edited topology files
#  - You want to re-process an existing system with different box/ion parameters
#  - You're working with non-standard topologies

# Combine with other options
python bin/gmx_prepare.py --topology protein.top --structure protein.gro \
    --box-distance 1.5 --box-type dodecahedron
```

### Dry-Run Mode

```bash
# Generate commands.sh without executing
python bin/gmx_prepare.py protein.pdb --dry-run

# Review the generated commands
cat commands.sh

# Execute when ready
bash commands.sh
```

### Verbose Mode

```bash
# Enable detailed debug output
python bin/gmx_prepare.py protein.pdb --verbose

# Combine with dry-run for detailed preview
python bin/gmx_prepare.py protein.pdb --dry-run --verbose
```

### Complete Example

```bash
python bin/gmx_prepare.py protein.pdb \
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
| `pdb_file` | Optional | Input PDB structure file (or use --topology/--structure) |
| `-t, --topology` | None | Use existing topology file (with --structure, skips pdb2gmx) |
| `-s, --structure` | None | Use existing structure/GRO file (with --topology, skips pdb2gmx) |
| `-ff, --forcefield` | `charmm27` | Force field selection (only with pdb_file) |
| `-wm, --water-model` | `tip3p` | Water model selection (only with pdb_file) |
| `-bt, --box-type` | `cubic` | Simulation box type |
| `-d, --box-distance` | `1.0` | Distance from protein to box edge (nm) |
| `-cs, --solvent-file` | `spc216.gro` | Solvent model file |
| `-pname, --cation` | `NA` | Cation type for neutralization |
| `-nname, --anion` | `CL` | Anion type for neutralization |
| `--no-neutral` | `False` | Do not add neutralizing ions |
| `--ignore-hydrogens` | `False` | Ignore H atoms in PDB (pass -ignh to pdb2gmx) |
| `--dry-run` | `False` | Generate commands without executing (preview mode) |
| `-v, --verbose` | `False` | Verbose output (show detailed progress) |

## Output Files

Output files use the input PDB basename (e.g., for `protein.pdb`):

- `protein.top` - System topology
- `protein.gro` - Processed structure with hydrogens
- `protein_box.gro` - Structure with defined box
- `protein_solvent.gro` - Solvated system
- `protein_ions.gro` - Final system with ions
- `protein_ions.tpr` - Input file for genion
- `minimal.mdp` - Generated MDP file for preprocessing
- `commands.sh` - Log of all executed commands

## Pipeline Output

The script provides detailed progress tracking:

```
============================================================
GROMACS SYSTEM PREPARATION PIPELINE
============================================================

📋 Configuration:
  Input PDB:      protein.pdb
  Force field:    charmm27
  Water model:    tip3p
  Box type:       cubic
  Box distance:   1.0 nm
  Ions:           NA/CL
  Neutralize:     Yes

✓ PDB file validation passed: protein.pdb

============================================================
[1/5] Generate topology (pdb2gmx)
============================================================
Running: gmx pdb2gmx -f protein.pdb -o protein.gro ...
✅ Step completed successfully

[2/5] Define simulation box (editconf)
...

============================================================
📊 SYSTEM SUMMARY
============================================================
Total atoms:     23456
Box dimensions:  6.2 x 6.2 x 6.2 nm
Box volume:      238.3 nm³

Molecule composition:
  Protein_chain_A    1
  SOL                7542
  NA                 12
  CL                 8
============================================================
```

## Typical Workflow

```bash
# 1. Prepare system
python bin/gmx_prepare.py protein.pdb --forcefield amber14sb

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

## Validation & Error Handling

### Pre-flight Validation

The script performs validation before execution:

- **PDB File Validation**:
  - Checks file exists and is readable
  - Verifies file contains ATOM or HETATM records
  - Validates file size is reasonable (>0 bytes, <1GB)
  - **Note**: If your PDB has hydrogen atoms, use `--ignore-hydrogens` to prevent pdb2gmx crashes

- **Force Field / Water Model Compatibility**:
  - Validates 18+ common force fields
  - Warns if water model is suboptimal for chosen force field
  - Example: Using TIP4P with CHARMM27 triggers warning
  
  ```
  ⚠️  Water model 'tip4p' may not be optimal with force field 'charmm27'
      Recommended: tip3p, spc, spce
  ```

### Error Hints

Each GROMACS tool has specific troubleshooting guidance:

- **pdb2gmx errors**: Missing atoms, unrecognized residues, terminal issues
- **editconf errors**: Invalid coordinates, box dimension problems
- **solvate errors**: Box too small, solvent file issues
- **grompp errors**: MDP syntax, topology issues, atom mismatches
- **genion errors**: Insufficient solvent molecules, group selection

Example error output:
```
❌ Command failed with exit code 1

💡 Troubleshooting hints for pdb2gmx:
  • Check for missing atoms in PDB file
  • Verify all residues are recognized by force field
  • Use pdb2gmx interactively to see detailed errors
  • Check terminal groups (NH3+/COO-)
```

### Error Recovery

- Validates GROMACS installation
- Checks input file existence before starting
- Monitors command execution at each step
- Reports clear error messages with context
- Preserves intermediate files on failure
- Provides step-specific troubleshooting hints

## Tips

- **Always** use `--dry-run` first to preview the workflow and validate configuration
- Use `--verbose` for detailed debugging output
- Check `commands.sh` for reproducibility and manual execution
- Review force field compatibility warnings—they help avoid simulation issues
- Adjust box distance based on simulation length (larger for longer simulations)
- Consider dodecahedral boxes for ~25% fewer water molecules (better efficiency)
- The system summary at the end provides useful statistics for documentation
- Force field defaults: AMBER variants use TIP3P, CHARMM uses TIP3P, GROMOS uses SPC/E

## Compatibility Matrix

Common force field and water model combinations:

| Force Field | Recommended Water Models |
|-------------|-------------------------|
| AMBER (all variants) | TIP3P, TIP4P-Ew, SPC/E |
| CHARMM27 | TIP3P, SPC, SPC/E |
| CHARMM36 | TIP3P only |
| GROMOS (all variants) | SPC, SPC/E |
| OPLS-AA | TIP3P, TIP4P, TIP5P, SPC/E |

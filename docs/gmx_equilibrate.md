# gmx_equilibrate.py

Automated equilibration pipeline for GROMACS molecular dynamics simulations.

## Overview

Automates the complete equilibration workflow from energy minimization through NVT and NPT equilibration with position restraints, preparing systems for production MD runs.

## Features

- ✅ **4-stage pipeline**: Energy minimization → NVT → NPT (restrained) → NPT (unrestrained)
- ✅ **Built-in MDP templates**: Sensible defaults for all simulation parameters
- ✅ **Automatic position restraints**: Generate restraints for heavy atoms with `gmx genrestr`
- ✅ **Proper continuation**: Handles velocity generation and checkpoint files correctly
- ✅ **Easy mode**: Simple CLI flags for common parameters (temperature, pressure, time)
- ✅ **Advanced mode**: Fine-tune coupling parameters, cutoffs, constraints, timestep
- ✅ **Professional logging**: Python logging module with INFO/WARNING/ERROR levels
- ✅ **Progress tracking**: Clear [1/4], [2/4], [3/4], [4/4] stage indicators
- ✅ **Dry-run mode**: Preview all commands before execution
- ✅ **Command logging**: Reproducible pipeline saved to `equilibrate_commands.sh`
- ✅ **Error hints**: Context-specific troubleshooting guidance

## Requirements

- Python 3.7+
- GROMACS 2021+
- Input: Structure file (.gro) and topology (.top) from gmx_prepare.py

## Usage

### Basic Usage

```bash
# Complete equilibration with defaults
python bin/gmx_equilibrate.py protein_ions.gro protein.top
```

**Default pipeline:**
- Energy minimization: 50,000 steps
- NVT (restrained): 100 ps at 300 K
- NPT (restrained): 200 ps at 300 K, 1 bar
- NPT (unrestrained): 500 ps at 300 K, 1 bar

### Custom Temperature and Pressure

```bash
# Human body temperature and pressure
python bin/gmx_equilibrate.py protein_ions.gro protein.top --temp 310 --pressure 1.0

# Custom temperature with extended final equilibration
python bin/gmx_equilibrate.py protein_ions.gro protein.top --temp 310 --time-npt2 1000
```

### Custom Run Lengths

```bash
# Thorough equilibration protocol
python bin/gmx_equilibrate.py protein_ions.gro protein.top \
    --time-em 100000 \
    --time-nvt 200 \
    --time-npt1 500 \
    --time-npt2 2000
```

### Temperature and Pressure Coupling

```bash
# Use Nose-Hoover thermostat
python bin/gmx_equilibrate.py protein_ions.gro protein.top --tcoupl nose-hoover

# Use Parrinello-Rahman barostat (for final NPT)
python bin/gmx_equilibrate.py protein_ions.gro protein.top --pcoupl Parrinello-Rahman

# Combine both
python bin/gmx_equilibrate.py protein_ions.gro protein.top \
    --tcoupl nose-hoover \
    --pcoupl Parrinello-Rahman
```

### Custom Output Prefix

```bash
# Use custom prefix for output files
python bin/gmx_equilibrate.py protein_ions.gro protein.top --prefix myprotein
# Produces: myprotein_em.gro, myprotein_nvt.gro, myprotein_npt1.gro, myprotein_npt2.gro
```

### Dry-Run Mode

```bash
# Preview pipeline without executing
python bin/gmx_equilibrate.py protein_ions.gro protein.top --dry-run

# Review generated commands
cat equilibrate_commands.sh

# Execute manually if desired
bash equilibrate_commands.sh
```

### Verbose Mode

```bash
# Enable detailed debug output
python bin/gmx_equilibrate.py protein_ions.gro protein.top --verbose

# Combine with dry-run for detailed preview
python bin/gmx_equilibrate.py protein_ions.gro protein.top --dry-run --verbose
```

### Advanced Fine-Tuning

```bash
# Custom coupling time constants
python bin/gmx_equilibrate.py protein_ions.gro protein.top \
    --tau-t 0.5 \
    --tau-p 5.0

# Custom cutoffs
python bin/gmx_equilibrate.py protein_ions.gro protein.top \
    --rcoulomb 1.2 \
    --rvdw 1.2

# Custom timestep and constraints
python bin/gmx_equilibrate.py protein_ions.gro protein.top \
    --dt 0.001 \
    --constraints all-bonds \
    --lincs-order 6

# Custom energy minimization tolerance
python bin/gmx_equilibrate.py protein_ions.gro protein.top --emtol 500
```

### Complete Example

```bash
python bin/gmx_equilibrate.py protein_ions.gro protein.top \
    --temp 310 \
    --pressure 1.0 \
    --time-em 100000 \
    --time-nvt 200 \
    --time-npt1 500 \
    --time-npt2 1000 \
    --tcoupl V-rescale \
    --pcoupl Berendsen \
    --prefix equil_310K \
    --dry-run
```

## Command-Line Options

### Positional Arguments

| Argument | Description |
|----------|-------------|
| `structure` | Input structure file (.gro) |
| `topology` | Topology file (.top) |

### Basic Options

| Option | Default | Description |
|--------|---------|-------------|
| `--prefix` | `equil` | Output file prefix |
| `--temp` | `300` | Temperature in K |
| `--pressure` | `1.0` | Pressure in bar |
| `--time-em` | `50000` | Energy minimization max steps |
| `--time-nvt` | `100` | NVT equilibration time (ps) |
| `--time-npt1` | `200` | NPT restrained equilibration time (ps) |
| `--time-npt2` | `500` | NPT unrestrained equilibration time (ps) |
| `--tcoupl` | `V-rescale` | Temperature coupling: V-rescale, berendsen, nose-hoover |
| `--pcoupl` | `Berendsen` | Pressure coupling: Berendsen, Parrinello-Rahman, C-rescale |
| `--dry-run` | `False` | Generate commands without executing (preview mode) |
| `-v, --verbose` | `False` | Verbose output (detailed progress) |
| `--commands-file` | `equilibrate_commands.sh` | Command log output file |

### Advanced Options (Fine-Tuning)

| Option | Default | Description |
|--------|---------|-------------|
| `--dt` | `0.002` | Timestep in ps (2 fs) |
| `--tau-t` | `0.1` | Temperature coupling time constant (ps) |
| `--tau-p` | `2.0` | Pressure coupling time constant (ps) |
| `--rcoulomb` | `1.0` | Coulomb cutoff in nm |
| `--rvdw` | `1.0` | Van der Waals cutoff in nm |
| `--lincs-order` | `4` | LINCS constraint order |
| `--constraints` | `h-bonds` | Constraint type: h-bonds, all-bonds, none |
| `--emtol` | `1000` | Energy minimization tolerance (kJ/mol/nm) |

## Pipeline Stages

### Stage 1: Energy Minimization
- **Algorithm**: Steepest descent
- **Purpose**: Remove steric clashes and bad contacts
- **Output**: `{prefix}_em.gro`
- **Typical duration**: Variable (until convergence or max steps)

**MDP parameters:**
- Integrator: steep
- Tolerance: 1000 kJ/mol/nm (adjustable with `--emtol`)
- Max steps: 50,000 (adjustable with `--time-em`)
- Constraints: h-bonds

### Stage 2: NVT Equilibration (Restrained)
- **Ensemble**: Canonical (constant N, V, T)
- **Purpose**: Stabilize temperature with position restraints
- **Restraints**: Heavy atoms (1000 kJ/mol/nm²)
- **Velocity generation**: Yes (random at target temperature)
- **Output**: `{prefix}_nvt.gro`, `{prefix}_nvt.cpt`
- **Default duration**: 100 ps

**MDP parameters:**
- Integrator: md (leap-frog)
- Timestep: 2 fs
- Temperature coupling: V-rescale (τ = 0.1 ps)
- Pressure coupling: None
- Constraints: h-bonds (LINCS)
- Position restraints: `-DPOSRES` (heavy atoms)

### Stage 3: NPT Equilibration (Restrained)
- **Ensemble**: Isothermal-isobaric (constant N, P, T)
- **Purpose**: Stabilize pressure and density with restraints
- **Restraints**: Heavy atoms (1000 kJ/mol/nm²)
- **Velocity generation**: No (continuation from NVT)
- **Output**: `{prefix}_npt1.gro`, `{prefix}_npt1.cpt`
- **Default duration**: 200 ps

**MDP parameters:**
- Integrator: md
- Timestep: 2 fs
- Temperature coupling: V-rescale (τ = 0.1 ps)
- Pressure coupling: Berendsen (τ = 2.0 ps, isotropic)
- Reference pressure: 1.0 bar
- Constraints: h-bonds (LINCS)
- Position restraints: `-DPOSRES` (heavy atoms)

### Stage 4: NPT Equilibration (Unrestrained)
- **Ensemble**: Isothermal-isobaric (constant N, P, T)
- **Purpose**: Final equilibration without restraints
- **Restraints**: None
- **Velocity generation**: No (continuation from NPT1)
- **Output**: `{prefix}_npt2.gro`, `{prefix}_npt2.cpt` ⭐
- **Default duration**: 500 ps

**MDP parameters:**
- Integrator: md
- Timestep: 2 fs
- Temperature coupling: V-rescale (τ = 0.1 ps)
- Pressure coupling: User-specified (default Berendsen, can use Parrinello-Rahman)
- Reference pressure: 1.0 bar
- Constraints: h-bonds (LINCS)
- Position restraints: None

## Output Files

All files use the specified prefix (default: `equil`):

### Generated MDP Files
- `{prefix}_em.mdp` - Energy minimization parameters
- `{prefix}_nvt.mdp` - NVT equilibration parameters
- `{prefix}_npt1.mdp` - NPT restrained parameters
- `{prefix}_npt2.mdp` - NPT unrestrained parameters

### Trajectory and Structure Files
- `{prefix}_em.gro` - Minimized structure
- `{prefix}_nvt.gro` - NVT equilibrated structure
- `{prefix}_npt1.gro` - NPT restrained structure
- `{prefix}_npt2.gro` - **Final equilibrated structure** ⭐

### Binary Files
- `{prefix}_*.tpr` - Portable binary run input files
- `{prefix}_*.cpt` - Checkpoint files (for continuation)
- `{prefix}_*.edr` - Energy files (for analysis)
- `{prefix}_*.xtc` - Compressed trajectory files

### Log Files
- `{prefix}_*_grompp.log` - Preprocessing logs
- `{prefix}_*_mdrun.log` - Simulation run logs

### Additional Files
- `posre_heavy.itp` - Position restraint topology
- `equilibrate_commands.sh` - Reproducibility script

## Pipeline Output

The script provides detailed progress tracking:

```
============================================================
GROMACS EQUILIBRATION PIPELINE
============================================================

🔍 DRY RUN MODE: Commands will be logged but not executed

📋 Configuration:
  Structure:       protein_ions.gro
  Topology:        protein.top
  Output prefix:   equil
  Temperature:     300 K
  Pressure:        1.0 bar
  T-coupling:      V-rescale
  P-coupling:      Berendsen

📊 Pipeline stages:
  [1/4] Energy minimization:    50000 steps
  [2/4] NVT (restrained):       100 ps
  [3/4] NPT (restrained):       200 ps
  [4/4] NPT (unrestrained):     500 ps

All commands logged to: equilibrate_commands.sh

[1/4] Energy Minimization
------------------------------------------------------------
✓ Created equil_em.mdp
✓ Success! (grompp and mdrun)

[2/4] NVT Equilibration (restrained)
------------------------------------------------------------
✓ Created posre_heavy.itp
✓ Created equil_nvt.mdp
✓ Success! (grompp and mdrun)

[3/4] NPT Equilibration (restrained)
------------------------------------------------------------
✓ Created equil_npt1.mdp
✓ Success! (grompp and mdrun)

[4/4] NPT Equilibration (unrestrained)
------------------------------------------------------------
✓ Created equil_npt2.mdp
✓ Success! (grompp and mdrun)

============================================================
✅ EQUILIBRATION COMPLETED SUCCESSFULLY!
============================================================

📁 Output files:
  Final structure:      equil_npt2.gro ⭐
  Checkpoint:           equil_npt2.cpt
  Energy files:         equil_*.edr
  Log files:            equil_*_mdrun.log
```

## Error Handling

The script includes comprehensive error detection and hints:

### Pre-flight Checks
- Verifies structure file exists and is readable
- Verifies topology file exists and is readable
- Checks GROMACS is installed and accessible

### Runtime Error Hints

**grompp errors:**
- Check MDP file syntax and parameters
- Verify topology file is complete and valid
- Ensure structure and topology match (atom counts)
- Check for missing position restraint files if using restraints

**mdrun errors:**
- Check if system is stable (minimize first)
- Verify temperature/pressure parameters are reasonable
- Check if previous step completed successfully
- Ensure sufficient disk space for output files

**genrestr errors:**
- Verify structure file exists and is valid
- Check if protein/heavy atoms are present
- Ensure index file has correct groups if using custom selection

## Analysis After Equilibration

After equilibration completes, analyze convergence:

```bash
# Temperature convergence (should stabilize around target)
gmx energy -f equil_nvt.edr -o temperature.xvg
# Select: Temperature

# Pressure convergence (NPT stages)
gmx energy -f equil_npt1.edr -o pressure.xvg
# Select: Pressure

# Density convergence
gmx energy -f equil_npt2.edr -o density.xvg
# Select: Density

# Potential energy
gmx energy -f equil_npt2.edr -o potential.xvg
# Select: Potential

# RMSD from starting structure
gmx rms -s equil_em.tpr -f equil_npt2.gro -o rmsd.xvg
# Select: Backbone, Backbone
```

## Tips and Best Practices

### Temperature Control
- **V-rescale** is recommended for most equilibrations (good balance of speed and accuracy)
- **Berendsen** heats/cools quickly but can cause artifacts (use for initial heating only)
- **Nose-Hoover** is more rigorous but requires longer equilibration

### Pressure Control
- **Berendsen** is excellent for equilibration (fast density relaxation)
- **Parrinello-Rahman** is better for production (correct ensemble) but can be unstable during equilibration
- Consider using Berendsen for restrained NPT, then switch to Parrinello-Rahman for final stage

### Equilibration Length
- **Minimum**: EM → 100 ps NVT → 200 ps NPT → 500 ps NPT (default)
- **Standard**: EM → 200 ps NVT → 500 ps NPT → 1000 ps NPT
- **Thorough**: EM → 500 ps NVT → 1 ns NPT → 2 ns NPT
- Monitor temperature, pressure, and density to confirm convergence

### Position Restraints
- Heavy atom restraints (1000 kJ/mol/nm²) are standard for equilibration
- Allows solvent and ions to relax while keeping protein structure stable
- Final unrestrained stage ensures system is fully equilibrated

### Common Workflows

**Fast equilibration** (testing):
```bash
gmx_equilibrate.py protein_ions.gro protein.top \
    --time-em 10000 --time-nvt 50 --time-npt1 100 --time-npt2 200
```

**Standard equilibration** (most systems):
```bash
gmx_equilibrate.py protein_ions.gro protein.top \
    --time-nvt 200 --time-npt1 500 --time-npt2 1000
```

**Conservative equilibration** (challenging systems):
```bash
gmx_equilibrate.py protein_ions.gro protein.top \
    --time-em 100000 --time-nvt 500 --time-npt1 1000 --time-npt2 2000 \
    --tcoupl berendsen  # Use Berendsen for better stability
```

## Integration with Other Tools

### Complete Workflow Example

```bash
# 1. System preparation
python bin/gmx_prepare.py protein.pdb --forcefield charmm36

# 2. Equilibration (this tool)
python bin/gmx_equilibrate.py protein_ions.gro protein.top --temp 310

# 3. Production MD
bash bin/gmx_continue_grompp.sh md 0 0 500 100 \
     mdp_templates/charmm36/production/md_npt.mdp equil_npt2

# 4. Analysis
python bin/gromacs_pca.py -s md_0.tpr -f md_0.xtc
```

### Resuming Failed Equilibration

If a stage fails, you can manually continue from the last successful checkpoint:

```bash
# If NVT completed but NPT1 failed:
gmx grompp -f equil_npt1.mdp -c equil_nvt.gro -r equil_em.gro \
           -t equil_nvt.cpt -p protein.top -o equil_npt1.tpr
gmx mdrun -deffnm equil_npt1
```

## Troubleshooting

### "System blowing up" during MD
- Energy minimization may need more steps: `--time-em 100000`
- Try smaller timestep: `--dt 0.001`
- Use stricter constraints: `--constraints all-bonds`

### Temperature not stabilizing
- Increase NVT time: `--time-nvt 200` or longer
- Check thermostat parameters: `--tau-t 0.1` (standard)
- Use Berendsen for initial heating if V-rescale struggles

### Pressure oscillating wildly
- Increase pressure coupling time: `--tau-p 5.0`
- Use Berendsen instead of Parrinello-Rahman during equilibration
- Ensure NVT completed successfully before starting NPT

### "Restraints not applied" warning
- Check that posre_heavy.itp was generated successfully
- Verify topology includes position restraint directive
- Check MDP has `-DPOSRES` in define line

## See Also

- [gmx_prepare.md](gmx_prepare.md) - System preparation (upstream step)
- [md_continuation.md](md_continuation.md) - Production MD continuation
- [GROMACS documentation](http://manual.gromacs.org/) - Official reference

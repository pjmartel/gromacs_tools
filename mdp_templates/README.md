# GROMACS MDP Templates

This directory contains well-tested MDP parameter files organized by force field and simulation type.

## Directory Structure

```
mdp_templates/
├── charmm36/          # CHARMM36 force field templates
├── amber14sb/         # AMBER14SB force field templates
└── plumed/            # PLUMED configuration files for enhanced sampling
```

## Force Field Specific Settings

### CHARMM36
- **Cutoffs**: 1.2 nm (recommended)
- **van der Waals**: Force-switch from 1.0 to 1.2 nm
- **Electrostatics**: PME with 0.12 nm grid spacing
- **Constraints**: All bonds (lincs-order 4)

### AMBER14SB
- **Cutoffs**: 1.0 nm (standard)
- **van der Waals**: Potential-shift-verlet
- **Electrostatics**: PME with 0.12 nm grid spacing
- **Constraints**: H-bonds (lincs-order 4)

## Template Organization

Each force field directory contains:
- **minimization/**: Energy minimization protocols
- **equilibration/**: NVT/NPT equilibration with various position restraints
- **production/**: Production MD settings

## Usage with gmx_continue_grompp.sh

All production templates have `nsteps = 0` and `tinit = 0` by default.
The continuation script will automatically set these values based on segment duration.

## Adding Custom Templates

When creating new templates:
1. Place in appropriate force field directory
2. Use descriptive filenames
3. Include header comments explaining purpose and settings
4. Set `nsteps = 0` and `tinit = 0` for production runs
5. Use title format: `<forcefield>,<purpose>`

## PLUMED Directory

Contains `.dat` files for enhanced sampling methods:
- Steered MD (SMD)
- Metadynamics
- Umbrella sampling
- Collective variable definitions (native contacts, helicity, etc.)

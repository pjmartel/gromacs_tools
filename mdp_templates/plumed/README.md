# PLUMED Configuration Files

This directory contains PLUMED `.dat` files for enhanced sampling and collective variable analysis.

## Usage

To use PLUMED with GROMACS:
```bash
gmx mdrun -deffnm production -plumed plumed.dat
```

## File Organization

Add your PLUMED files here organized by method:
- **Steered MD**: Force-based pulling simulations
- **Metadynamics**: Free energy surface exploration
- **Umbrella Sampling**: Enhanced sampling along reaction coordinates
- **Collective Variables**: Analysis tools (native contacts, helicity, RMSD, etc.)

## Example Structure (to be added)

```
plumed/
├── smd_pull_z.dat              # Pull along z-axis
├── metad_distance.dat          # Metadynamics on distance CV
├── umbrella_windows.dat        # Umbrella sampling template
├── cv_native_contacts.dat      # Native contacts definition
├── cv_helicity.dat             # Alpha helix content
└── cv_gyration.dat             # Radius of gyration
```

## Resources

- PLUMED Documentation: https://www.plumed.org/
- PLUMED Masterclass: https://www.plumed.org/masterclass
- CV Library: https://www.plumed.org/doc-master/user-doc/html/_colvar.html

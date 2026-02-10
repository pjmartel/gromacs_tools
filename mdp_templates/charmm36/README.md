# CHARMM36 Force Field Templates

MDP parameter files optimized for the CHARMM36 force field.

## Key Settings for CHARMM36

- **Cutoff scheme**: Verlet
- **vdW cutoff**: 1.2 nm
- **vdW modifier**: Force-switch (from 1.0 to 1.2 nm)
- **Coulomb cutoff**: 1.2 nm
- **Coulomb type**: PME
- **Constraints**: all-bonds (allows 2 fs timestep)
- **LINCS order**: 4 (recommended for all-bonds)

## Reference

MacKerell et al., J. Phys. Chem. B 1998, 102, 3586-3616
Best et al., J. Chem. Theory Comput. 2012, 8, 3257-3273

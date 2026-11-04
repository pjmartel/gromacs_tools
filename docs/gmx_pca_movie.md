# gmx_pca_movie.py

Generate trajectory movies showing conformational changes along principal components.

## Overview

Creates animated trajectories that visualize the motion captured by specific principal components. Uses GROMACS `gmx anaeig -extr` to interpolate structures along eigenvectors, producing movies that can be visualized in VMD, PyMOL, ChimeraX, etc.

## Features

- ✅ Generate smooth interpolated trajectories along any PC
- ✅ Takes everything from a `gmx_pca.py` output folder (`--pca`), also of a PCA of several
  trajectories (replicas, apo vs holo): the movie then spans all of them
- ✅ Full movies or extreme structures only
- ✅ Extent from the trajectories (default) or chosen (`--extreme`), and frame count
- ✅ Compatible with all molecular viewers
- ✅ Multiple output formats (PDB, XTC, TRR)
- ✅ **Command logging** (on by default): every GROMACS command run, plus the exact invocation
  used, is saved to `<outdir>/gmx_pca_movie_commands.sh` (override with `--commands-file`); the file from an earlier run is kept as `…_commands.1.sh`, `.2.sh`, … (highest number = most recent)

## Requirements

- GROMACS 2021+
- Outputs from PCA analysis: the output folder of `gmx_pca.py` (`fit.gro`, `eigenvectors.trr`,
  `pca_groups.ndx` and the fitted trajectories)
- Molecular viewer (VMD, PyMOL, ChimeraX, etc.)

## Quick Start

```bash
# Generate 50-frame movie along PC1
python bin/gmx_pca_movie.py --pca pca_out -o pc1_movie.pdb --pc 1 --nframes 50
```

## Usage

### Basic Movie Generation

```bash
# PC1 movie (default: 50 frames)
python bin/gmx_pca_movie.py --pca pca_out -o pc1_movie.pdb --pc 1

# PC2 movie with 100 frames
python bin/gmx_pca_movie.py --pca pca_out -o pc2_movie.pdb --pc 2 --nframes 100
```

`--pca pca_out` stands for `-s pca_out/fit.gro -f pca_out/fit.xtc -v pca_out/eigenvectors.trr
-n pca_out/pca_groups.ndx`; any of these given explicitly takes precedence. Without `--pca`,
`-s`, `-f` and `-v` are needed (and `-n` defaults to `pca_groups.ndx` next to `-s`).

### Several Trajectories (Replicas, Apo vs Holo)

A PCA of several trajectories (`gmx_pca.py -f t1 t2 ...` or `--trajectories`) has no `fit.xtc`:
each trajectory is fitted in `trajectories/<label>/`, and `fitted_trajectories.txt` lists them.
`--pca` uses all of them, joined, so the movie spans the whole range the trajectories explore on
that PC, as in the PCA's plots:

```bash
python bin/gmx_pca.py --trajectories replicas.txt -o pca_reps
python bin/gmx_pca_movie.py --pca pca_reps -o pc1_movie.pdb --pc 1
```

The extremes of only some trajectories: give them with `-f` (several are joined):

```bash
# PC1 of the common PCA, over the range of the apo replicas only
python bin/gmx_pca_movie.py --pca pca_reps -o pc1_apo.pdb --pc 1 \
    -f pca_reps/trajectories/apo1/fit.xtc pca_reps/trajectories/apo2/fit.xtc
```

A folder made by an older `gmx_pca.py` (before `fitted_trajectories.txt`): run `gmx_pca.py`
again with the same command (it reuses the preprocessed trajectories, so it is quick), or give
the trajectories with `-f`.

### Extreme Structures Only

```bash
# Generate only min and max structures (2 frames total)
python bin/gmx_pca_movie.py --pca pca_out -o pc1_extremes.pdb --pc 1 --extremes-only
```

### Custom Extent

By default the movie goes from the lowest to the highest projection of the trajectories on the
PC. `--extreme X` makes it go from -X to +X nm instead, symmetric around the average structure:

```bash
# From -1.5 to +1.5 nm along the eigenvector
python bin/gmx_pca_movie.py --pca pca_out -o pc1_wide.pdb --pc 1 --extreme 1.5
```

The projections are in nm summed over all analysed atoms, so their typical size depends on the
system and the selection: look at the range of `projection_pc<N>*.xvg` (or, for several
trajectories, the standard deviations in `projection_stats.txt`) before choosing X. Beyond the
range the trajectories reach, the motion is extrapolated and exaggerated.

### Different Output Formats

```bash
# XTC format (compressed, smaller files)
python bin/gmx_pca_movie.py --pca pca_out -o pc1_movie.xtc --pc 1

# TRR format (full precision)
python bin/gmx_pca_movie.py --pca pca_out -o pc1_movie.trr --pc 1
```

### Custom Output Directory

```bash
python bin/gmx_pca_movie.py --pca pca_out -o movies/pc1.pdb --pc 1 --outdir movies
```

## Command-Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `--pca DIR` | - | Output folder of `gmx_pca.py`: the defaults of `-s`, `-f`, `-v` and `-n` (below) |
| `-s, --structure` | `DIR/fit.gro` | Structure file (`fit.gro` from PCA) |
| `-f, --trajectory` | `DIR/fit.xtc`, or the trajectories in `DIR/fitted_trajectories.txt` | The fitted trajectory or trajectories (several are joined); the movie spans their lowest to highest projection |
| `-v, --eigenvec` | `DIR/eigenvectors.trr` | Eigenvectors file (`eigenvectors.trr` from PCA) |
| `-n, --index` | `pca_groups.ndx` next to `-s`, if present | Index with the fit and analysis groups, as written by `gmx_pca.py`. Needed when the PCA was run with `--fit-selection`: `fit.gro` then holds the fitted and the analysed atoms, the eigenvectors only the analysed ones (the movie shows those). Without it (older PCA folders) all atoms of `-s` are used |
| `-o, --output` | Required | Output trajectory (PDB/XTC/TRR) |
| `--pc` | `1` | Principal component number |
| `--nframes` | `50` | Number of interpolated frames |
| `--extreme NM` | from the trajectories | Go from -NM to +NM nm along the eigenvector (`gmx anaeig -max`) |
| `--extremes-only` | `False` | Generate only min/max structures, as `-o` |
| `--outdir` | `.` | Output directory for logs and the commands file |
| `--gmx-bin` | `gmx` | GROMACS binary name/path |
| `--commands-file` | `gmx_pca_movie_commands.sh` | Commands file, relative to `--outdir` |

## Understanding the Movie

### What Does It Show?

The frames are evenly spaced along the eigenvector, from the average structure + (minimum ×
eigenvector) to the average structure + (maximum × eigenvector):

- **Frame 1**: the minimum (most negative projection)
- **Last frame**: the maximum (most positive projection)
- **In between**: straight-line interpolation; the average structure is where the projection
  is 0 (the middle frame with `--extreme`)

Only the analysed atoms (the PCA's `--selection`) are in the movie. To see the motion back and
forth, play it as a rock/bounce loop in the viewer.

### Interpreting Motion

- **Large displacement**: Important conformational change
- **Collective motion**: Multiple regions moving together
- **Hinge regions**: Areas with minimal motion between moving domains
- **Breathing motions**: Opening/closing of cavities or pockets

## Complete Workflow

```bash
# 1. Run PCA analysis
python bin/gmx_pca.py -s topol.tpr -f traj.xtc \
    --selection "Backbone" -o pca_out

# 2. Generate movies for first 3 PCs
for pc in 1 2 3; do
    python bin/gmx_pca_movie.py --pca pca_out -o pc${pc}_movie.pdb --pc $pc --nframes 100
done

# 3. Visualize in VMD
vmd pc1_movie.pdb

# 4. In VMD: Graphics → Representations, then play trajectory
```

## Visualization Tips

The movie's B-factors are all zero: colour by residue, structure or chain instead.

### VMD

```tcl
# Load movie
vmd pc1_movie.pdb

# In VMD Tk Console:
# Nice cartoon representation
mol modstyle 0 0 NewCartoon
mol modcolor 0 0 Structure

# Smooth playback, back and forth
animate goto 0
animate style rock
animate speed 0.5
```

### PyMOL

```python
# Load movie
load pc1_movie.pdb

# Cartoon representation
hide everything
show cartoon
spectrum resi

# Play
mplay
```

### ChimeraX

```bash
# Open movie
open pc1_movie.pdb

# Style
cartoon style width 2.0 thick 0.5
rainbow

# Play
coordset slider #1
```

## Use Cases

### 1. Understanding Dominant Motions

```bash
# Visualize top 3 PCs to understand main conformational changes
for pc in 1 2 3; do
    python bin/gmx_pca_movie.py --pca pca_out -o pc${pc}_movie.pdb --pc $pc
done
```

### 2. Comparing WT vs Mutant

```bash
# One PCA of both: the same PCs, each movie over that system's range (a selection with the
# same atoms in both, e.g. C-alpha: residue names may differ, atom names not)
python bin/gmx_pca.py -s wt.tpr mut.tpr -f wt.xtc mut.xtc --labels wt mut \
    --selection 'name CA' -o pca_both
python bin/gmx_pca_movie.py --pca pca_both -o wt_pc1.pdb --pc 1 -f pca_both/trajectories/wt/fit.xtc
python bin/gmx_pca_movie.py --pca pca_both -o mut_pc1.pdb --pc 1 -f pca_both/trajectories/mut/fit.xtc

# Compare in VMD: load both and play side-by-side
```

With separate PCAs (`pca_wt/`, `pca_mut/`), PC1 of one need not be the same motion as PC1 of
the other: compare them with `gmx_pca.py --rmsip` first.

### 3. Presentation/Publication

```bash
# High quality, smooth 200-frame movie
python bin/gmx_pca_movie.py --pca pca_out -o presentation_pc1.pdb --pc 1 --nframes 200
```

### 4. Structural Analysis

```bash
# Extract extreme conformations for docking
python bin/gmx_pca_movie.py --pca pca_out -o pc1_extremes.pdb --pc 1 --extremes-only

# Result: 2-frame trajectory (the two ends of the motion, e.g. open and closed)
```

## Advanced Usage

### Batch Generate Multiple PCs

```bash
#!/bin/bash
# generate_all_pc_movies.sh

for pc in {1..5}; do
    echo "Generating PC${pc} movie..."
    python bin/gmx_pca_movie.py --pca pca_out \
        -o movies/pc${pc}_movie.pdb \
        --pc $pc \
        --nframes 100 \
        --outdir movies
done
```

### Multiple Extents

```bash
# Small, medium, large extent movies
for extent in 0.5 1.0 2.0; do
    python bin/gmx_pca_movie.py --pca pca_out \
        -o pc1_extent${extent}.pdb --pc 1 --extreme $extent
done
```

## Tips & Best Practices

1. **Frame Count**:
   - 50 frames: Quick preview
   - 100 frames: Standard smooth animation
   - 200+ frames: Very smooth, publication quality

2. **Extent Selection**:
   - Default (the trajectories' range): shows what the simulations actually reach
   - `--extreme` below the range: the central, best-sampled part of the motion
   - `--extreme` beyond the range: exaggerated, may show unphysical structures

3. **File Format**:
   - PDB: Universal compatibility, larger files
   - XTC: Compressed, smaller, good for long movies
   - TRR: Full precision, largest files

4. **Which PCs to Visualize**:
   - Always visualize PC1-2 (most variance)
   - Check scree plot to see if PC3+ are significant
   - Usually first 3-5 PCs are interpretable

5. **Visualization**:
   - Use cartoon for proteins (clearer motion)
   - Add rock/bounce loop animation for continuous playback

6. **Combining with Analysis**:
   - Use extremes-only for docking studies
   - Extract representative structures from movies
   - Overlay with ligand binding sites to see accessibility changes

## Troubleshooting

### "Error: eigenvector file not found" / "-s, -f, -v needed"
- Ensure you've run `gmx_pca.py` first
- Give its output folder with `--pca`, or check the paths

### "... has neither fit.xtc nor fitted_trajectories.txt"
- The `--pca` folder is not a `gmx_pca.py` output folder, or comes from an older
  `gmx_pca.py` run of several trajectories: run `gmx_pca.py` again, or give `-f`

### GROMACS error about the number of atoms
- Use the structure and index of the same PCA (`fit.gro`, `pca_groups.ndx`): `--pca` does
- Don't use original `topol.tpr` (may have different atom selection)

### Movie looks unrealistic
- Use a smaller `--extreme` (below the trajectories' range)
- Some PCs capture noise rather than real motion

### No visible motion
- Check if you selected the right PC
- Verify PC explains significant variance (check scree plot)
- A larger `--extreme` exaggerates the motion

### GROMACS error during generation
- Verify input files are correct
- Check that eigenvector file has enough components (`--pc` value exists)
- The anaeig log in `{outdir}/logs/` has the details

## Output Files

- `{output}` - Trajectory movie (specified with `-o`)
- `{outdir}/logs/anaeig_pc<N>_movie.log` (`anaeig_pc<N>_extremes.log` with `--extremes-only`) -
  GROMACS log; `trjcat.log` when several trajectories are joined
- `{outdir}/gmx_pca_movie_commands.sh` - The commands, to re-run

## See Also

- [gmx_pca.md](gmx_pca.md) - For running PCA analysis
- [plot_xvg.md](plot_xvg.md) - For plotting PC projections
- GROMACS `gmx anaeig` documentation - For understanding interpolation

## References

- GROMACS manual: `gmx anaeig -extr`
- Essential dynamics: Amadei et al. (1993) Proteins 17:412-425

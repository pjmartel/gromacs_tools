# gromacs_pca_movie.py

Generate trajectory movies showing conformational changes along principal components.

## Overview

Creates animated trajectories that visualize the motion captured by specific principal components. Uses GROMACS `gmx anaeig -extr` to interpolate structures along eigenvectors, producing movies that can be visualized in VMD, PyMOL, ChimeraX, etc.

## Features

- ✅ Generate smooth interpolated trajectories along any PC
- ✅ Full movies or extreme structures only
- ✅ Customizable extent and frame count
- ✅ Compatible with all molecular viewers
- ✅ Multiple output formats (PDB, XTC, TRR)

## Requirements

- GROMACS 2021+
- Outputs from PCA analysis: `fit.gro` and `eigenvectors.trr`
- Molecular viewer (VMD, PyMOL, ChimeraX, etc.)

## Quick Start

```bash
# Generate 50-frame movie along PC1
python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro \
    -v pca_out/eigenvectors.trr \
    -o pc1_movie.pdb --pc 1 --nframes 50
```

## Usage

### Basic Movie Generation

```bash
# PC1 movie (default: 50 frames)
python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro \
    -v pca_out/eigenvectors.trr \
    -o pc1_movie.pdb --pc 1

# PC2 movie with 100 frames
python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro \
    -v pca_out/eigenvectors.trr \
    -o pc2_movie.pdb --pc 2 --nframes 100
```

### Extreme Structures Only

```bash
# Generate only min and max structures (2 frames total)
python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro \
    -v pca_out/eigenvectors.trr \
    -o pc1_extremes.pdb --pc 1 --extremes-only
```

### Custom Extent

```bash
# Larger amplitude (3.0 nm along eigenvector)
python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro \
    -v pca_out/eigenvectors.trr \
    -o pc1_large.pdb --pc 1 --extreme 3.0

# Smaller amplitude (1.0 nm)
python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro \
    -v pca_out/eigenvectors.trr \
    -o pc1_small.pdb --pc 1 --extreme 1.0
```

### Different Output Formats

```bash
# XTC format (compressed, smaller files)
python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro \
    -v pca_out/eigenvectors.trr \
    -o pc1_movie.xtc --pc 1

# TRR format (full precision)
python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro \
    -v pca_out/eigenvectors.trr \
    -o pc1_movie.trr --pc 1
```

### Custom Output Directory

```bash
python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro \
    -v pca_out/eigenvectors.trr \
    -o movies/pc1.pdb --pc 1 --outdir movies
```

## Command-Line Options

| Option | Default | Description |
|--------|---------|-------------|
| `-s, --structure` | Required | Structure file (use `fit.gro` from PCA) |
| `-v, --eigenvec` | Required | Eigenvectors file (`eigenvectors.trr` from PCA) |
| `-o, --output` | Required | Output trajectory (PDB/XTC/TRR) |
| `--pc` | `1` | Principal component number |
| `--nframes` | `50` | Number of interpolated frames |
| `--extreme` | `2.0` | Extent along eigenvector (nm) |
| `--extremes-only` | `False` | Generate only min/max structures |
| `--outdir` | `.` | Output directory for logs |
| `--gmx-bin` | `gmx` | GROMACS binary name/path |

## Understanding the Movie

### What Does It Show?

The movie interpolates structures from:
- **Negative extreme**: Average structure - (extreme × eigenvector)
- **Average structure**: The mean conformation
- **Positive extreme**: Average structure + (extreme × eigenvector)

The animation goes: min → average → max → average → min (full cycle)

### Interpreting Motion

- **Large displacement**: Important conformational change
- **Collective motion**: Multiple regions moving together
- **Hinge regions**: Areas with minimal motion between moving domains
- **Breathing motions**: Opening/closing of cavities or pockets

## Complete Workflow

```bash
# 1. Run PCA analysis
python bin/gromacs_pca.py -s topol.tpr -f traj.xtc \
    --selection "Backbone" -o pca_out

# 2. Generate movies for first 3 PCs
python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro -v pca_out/eigenvectors.trr \
    -o pc1_movie.pdb --pc 1 --nframes 100

python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro -v pca_out/eigenvectors.trr \
    -o pc2_movie.pdb --pc 2 --nframes 100

python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro -v pca_out/eigenvectors.trr \
    -o pc3_movie.pdb --pc 3 --nframes 100

# 3. Visualize in VMD
vmd pc1_movie.pdb

# 4. In VMD: Graphics → Representations, then play trajectory
```

## Visualization Tips

### VMD

```tcl
# Load movie
vmd pc1_movie.pdb

# In VMD Tk Console:
# Nice cartoon representation
mol modstyle 0 0 NewCartoon
mol modcolor 0 0 Structure

# Color by B-factor (RMSF)
mol modcolor 0 0 Beta

# Smooth playback
animate goto 0
animate style loop
animate speed 0.5
```

### PyMOL

```python
# Load movie
load pc1_movie.pdb

# Cartoon representation
hide everything
show cartoon

# Color by B-factor
spectrum b, blue_white_red, minimum=0, maximum=5

# Play
mplay
```

### ChimeraX

```bash
# Open movie
open pc1_movie.pdb

# Style
cartoon style width 2.0 thick 0.5
color byattribute bfactor palette cyan:yellow:red

# Play
coordset slider #1
```

## Use Cases

### 1. Understanding Dominant Motions

```bash
# Visualize top 3 PCs to understand main conformational changes
for pc in 1 2 3; do
    python bin/gromacs_pca_movie.py \
        -s pca_out/fit.gro -v pca_out/eigenvectors.trr \
        -o pc${pc}_movie.pdb --pc $pc
done
```

### 2. Comparing WT vs Mutant

```bash
# Wild type PC1
python bin/gromacs_pca_movie.py \
    -s pca_wt/fit.gro -v pca_wt/eigenvectors.trr \
    -o wt_pc1.pdb --pc 1

# Mutant PC1
python bin/gromacs_pca_movie.py \
    -s pca_mut/fit.gro -v pca_mut/eigenvectors.trr \
    -o mut_pc1.pdb --pc 1

# Compare in VMD: load both and play side-by-side
```

### 3. Presentation/Publication

```bash
# High quality, smooth 200-frame movie
python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro -v pca_out/eigenvectors.trr \
    -o presentation_pc1.pdb --pc 1 --nframes 200 --extreme 2.5
```

### 4. Structural Analysis

```bash
# Extract extreme conformations for docking
python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro -v pca_out/eigenvectors.trr \
    -o pc1_extremes.pdb --pc 1 --extremes-only

# Result: 2-frame trajectory (open and closed conformations)
# Extract frames with gmx trjconv or viewer
```

## Advanced Usage

### Batch Generate Multiple PCs

```bash
#!/bin/bash
# generate_all_pc_movies.sh

for pc in {1..5}; do
    echo "Generating PC${pc} movie..."
    python bin/gromacs_pca_movie.py \
        -s pca_out/fit.gro \
        -v pca_out/eigenvectors.trr \
        -o movies/pc${pc}_movie.pdb \
        --pc $pc \
        --nframes 100 \
        --outdir movies
done
```

### Extract Frames for Further Analysis

```bash
# Generate movie
python bin/gromacs_pca_movie.py \
    -s pca_out/fit.gro -v pca_out/eigenvectors.trr \
    -o pc1.xtc --pc 1 --nframes 100

# Extract specific frame (e.g., frame 50 = max extreme)
gmx trjconv -f pc1.xtc -s pca_out/fit.gro -o pc1_max.pdb -dump 50
```

### Multiple Extents

```bash
# Small, medium, large extent movies
for extent in 1.0 2.0 3.0; do
    python bin/gromacs_pca_movie.py \
        -s pca_out/fit.gro -v pca_out/eigenvectors.trr \
        -o pc1_extent${extent}.pdb --pc 1 --extreme $extent
done
```

## Tips & Best Practices

1. **Frame Count**:
   - 50 frames: Quick preview
   - 100 frames: Standard smooth animation
   - 200+ frames: Very smooth, publication quality

2. **Extent Selection**:
   - 2.0 nm: Good default, covers typical motion
   - 1.0 nm: Conservative, for subtle motions
   - 3.0+ nm: May show unphysical structures

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
   - Color by B-factor (RMSF) to show flexible regions
   - Add loop animation for continuous playback

6. **Combining with Analysis**:
   - Use extremes-only for docking studies
   - Extract representative structures from movies
   - Overlay with ligand binding sites to see accessibility changes

## Troubleshooting

### "Error: eigenvectors file not found"
- Ensure you've run `gromacs_pca.py` first
- Check path to `eigenvectors.trr`

### "Error: atoms mismatch"
- Use the same structure file as PCA analysis (`fit.gro`)
- Don't use original `topol.tpr` (may have different atom selection)

### Movie looks unrealistic
- Reduce `--extreme` value (try 1.5 or 1.0)
- Some PCs capture noise rather than real motion

### No visible motion
- Increase `--extreme` value (try 3.0 or 4.0)
- Check if you selected the right PC
- Verify PC explains significant variance (check scree plot)

### GROMACS error during generation
- Verify input files are correct
- Check that eigenvector file has enough components (`--pc` value exists)

## Output Files

- `{output}` - Trajectory movie (specified with `-o`)
- `{outdir}/extr.log` - GROMACS log file

## See Also

- [gromacs_pca.md](gromacs_pca.md) - For running PCA analysis
- [plot_xvg.md](plot_xvg.md) - For plotting PC projections
- GROMACS `gmx anaeig` documentation - For understanding interpolation

## References

- GROMACS manual: `gmx anaeig -extr`
- Essential dynamics: Amadei et al. (1993) Proteins 17:412-425

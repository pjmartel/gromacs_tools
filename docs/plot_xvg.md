# plot_xvg.py

Versatile command-line tool for plotting GROMACS XVG files with extensive customization options.

## Overview

A powerful plotting tool that handles everything from simple time series to complex multi-dimensional correlations, histograms, and scatter plots with color gradients.

## Features

- ✅ **Multiple plot types**: Line plots, scatter plots, histograms, correlations
- ✅ **Column selection**: Plot specific columns from multi-column files
- ✅ **Custom legends**: Essential for files without XVG headers
- ✅ **Moving averages**: Smooth noisy data
- ✅ **Trajectory slicing**: Analyze specific simulation portions
- ✅ **Multi-file overlay**: Compare multiple datasets
- ✅ **2D/3D correlations**: Correlation analysis
- ✅ **Extensive customization**: Styles, colors, labels, sizes
- ✅ **Publication-ready**: High DPI, multiple formats (PNG, PDF, SVG)

## Quick Examples

```bash
# Simple line plot
python bin/plot_xvg.py rmsd.xvg --style lines

# Scatter plot with time-colored points
python bin/plot_xvg.py pca_projection.xvg --scatter --colormap viridis

# Moving average
python bin/plot_xvg.py energy.xvg --moving-avg --window 50 --style lines

# Multiple files overlay
python bin/plot_xvg.py wt.xvg mutant.xvg --multi --legends "WT" "Mutant"

# Histogram
python bin/plot_xvg.py rmsd.xvg --histogram --bins 100

# Specific columns with custom legends
python bin/plot_xvg.py data.xvg --columns 1 3 --legends "RMSD" "RMSF"
```

## Plot Types

### 1. Time Series (Default)

```bash
# Dots (default)
python bin/plot_xvg.py energy.xvg

# Lines
python bin/plot_xvg.py energy.xvg --style lines

# Lines + dots
python bin/plot_xvg.py energy.xvg --style lines+dots
```

### 2. Scatter Plots

```bash
# Color by time/frame order
python bin/plot_xvg.py projection.xvg --scatter --colormap viridis

# Custom marker size
python bin/plot_xvg.py data.xvg --scatter --markersize 5
```

### 3. Histograms

```bash
# Distribution of values
python bin/plot_xvg.py rmsd.xvg --histogram --bins 50

# Multiple datasets
python bin/plot_xvg.py rmsd1.xvg rmsd2.xvg --multi --histogram --bins 100
```

### 4. Correlations

```bash
# 2D correlation (y-values from 2 files)
python bin/plot_xvg.py rmsd1.xvg rmsd2.xvg --xy-correlation

# 3D correlation (y-values from 3 files)
python bin/plot_xvg.py pc1.xvg pc2.xvg pc3.xvg --xy-correlation --scatter

# Multiple correlations overlaid
python bin/plot_xvg.py r1a.xvg r2a.xvg r1b.xvg r2b.xvg \
    --xy-correlation --multi --legends "Set A" "Set B"
```

## Column Selection & Legends

### Working with Multi-Column Files

```bash
# Plot only columns 1, 3, and 5 (1-indexed)
python bin/plot_xvg.py multicolumn.xvg --columns 1 3 5

# Add custom legends
python bin/plot_xvg.py data.xvg --columns 1 2 3 \
    --legends "RMSD" "RMSF" "Gyration"

# Essential for files without XVG headers
python bin/plot_xvg.py raw_data.dat \
    --columns 1 2 --legends "Temperature" "Pressure" --style lines
```

## Trajectory Slicing

```bash
# Plot frames 100-500 (0-indexed)
python bin/plot_xvg.py trajectory.xvg --start 100 --end 500

# Skip equilibration (first 1000 frames)
python bin/plot_xvg.py rmsd.xvg --start 1000 --style lines

# Analyze early trajectory only
python bin/plot_xvg.py energy.xvg --end 2000
```

## Customization Options

### Labels and Titles

```bash
python bin/plot_xvg.py data.xvg \
    --title "RMSD vs Time" \
    --xlabel "Time (ns)" \
    --ylabel "RMSD (nm)"
```

### Figure Size and DPI

```bash
# Custom figure size (width height in inches)
python bin/plot_xvg.py data.xvg --figsize 12 8

# High DPI for publications
python bin/plot_xvg.py data.xvg --output figure.png --dpi 600
```

### Fonts

```bash
# Set base font size
python bin/plot_xvg.py data.xvg --fontsize 14
```

### Axis Limits and Aspect

```bash
# Set axis limits
python bin/plot_xvg.py data.xvg --xlim 0 100 --ylim 0 5

# Equal aspect ratio (1:1)
python bin/plot_xvg.py correlation.xvg --aspect equal

# Custom aspect ratio
python bin/plot_xvg.py data.xvg --aspect 1.5
```

### Styles and Colors

```bash
# Matplotlib style
python bin/plot_xvg.py data.xvg --plot-style seaborn

# List available styles
python bin/plot_xvg.py --plot-style available

# Custom colormap for scatter
python bin/plot_xvg.py pca.xvg --scatter --colormap plasma
```

### Marker Size

```bash
# Larger markers
python bin/plot_xvg.py data.xvg --markersize 5

# Smaller for dense data
python bin/plot_xvg.py dense_data.xvg --scatter --markersize 1
```

## Saving Plots

```bash
# PNG (default)
python bin/plot_xvg.py data.xvg --output plot.png

# PDF (vector, publication-ready)
python bin/plot_xvg.py data.xvg --output plot.pdf --dpi 300

# SVG (editable in Inkscape/Illustrator)
python bin/plot_xvg.py data.xvg --output plot.svg
```

## Complete Examples

### Publication Figure

```bash
python bin/plot_xvg.py rmsd.xvg \
    --style lines \
    --moving-avg --window 100 \
    --title "Backbone RMSD" \
    --xlabel "Time (ns)" \
    --ylabel "RMSD (nm)" \
    --figsize 8 6 \
    --fontsize 12 \
    --output rmsd_figure.pdf \
    --dpi 300 \
    --plot-style seaborn
```

### Multi-System Comparison

```bash
python bin/plot_xvg.py wt_rmsd.xvg mutant_rmsd.xvg drug_rmsd.xvg \
    --multi \
    --legends "Wild Type" "F6A Mutant" "+ Drug" \
    --style lines \
    --title "RMSD Comparison" \
    --xlabel "Time (ns)" \
    --ylabel "Backbone RMSD (nm)" \
    --output comparison.png
```

### PCA Projection

```bash
python bin/plot_xvg.py pc1_vs_pc2.xvg \
    --scatter \
    --colormap viridis \
    --title "PC1 vs PC2" \
    --xlabel "PC1 (nm)" \
    --ylabel "PC2 (nm)" \
    --aspect equal \
    --markersize 2 \
    --output pca_projection.pdf
```

### Raw Data File (No Headers)

```bash
python bin/plot_xvg.py raw_multicolumn.dat \
    --columns 1 2 3 \
    --legends "Temperature" "Pressure" "Volume" \
    --style lines \
    --title "Thermodynamic Properties" \
    --output thermodynamics.png
```

## Command-Line Options Summary

| Option | Description |
|--------|-------------|
| `--style` | Plot style: dots, lines, lines+dots |
| `--scatter` | Scatter plot with time-colored points |
| `--histogram` | Plot distribution histogram |
| `--xy-correlation` | 2D/3D correlation plot |
| `--moving-avg` | Apply moving average filter |
| `--window` | Moving average window size |
| `--columns` | Select specific columns (1-indexed) |
| `--legends` | Custom legend labels |
| `--multi` | Plot multiple files on same axes |
| `--start` / `--end` | Trajectory slicing (row indices) |
| `--title` / `--xlabel` / `--ylabel` | Custom labels |
| `--figsize` | Figure size (width height) |
| `--aspect` | Aspect ratio (equal, auto, or number) |
| `--xlim` / `--ylim` | Axis limits |
| `--markersize` | Marker/dot size |
| `--colormap` | Colormap for scatter plots |
| `--plot-style` | Matplotlib style theme |
| `--fontsize` | Base font size |
| `--output` | Save to file (PNG, PDF, SVG) |
| `--dpi` | Output resolution |

## Tips

- Use `--scatter --colormap viridis` for PCA projections
- Add `--moving-avg` for noisy energy plots
- Use `--columns` with `--legends` for files without headers
- Set `--aspect equal` for correlation/PCA plots
- Save as PDF for publications (vector graphics)
- Use `--start` to skip equilibration periods

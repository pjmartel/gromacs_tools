# Example Workflows

This directory contains example scripts and workflows for common GROMACS analysis tasks.

## Basic Workflows

### 1. Complete MD Simulation Workflow

```bash
#!/bin/bash
# complete_md_workflow.sh - Full simulation pipeline

# Step 1: System preparation
python bin/gromacs_prepare.py protein.pdb \
    --forcefield amber14sb \
    --water-model tip3p \
    --box-distance 1.2 \
    --box-type dodecahedron

# Step 2: Energy minimization
gmx grompp -f mdp_templates/amber14sb/minimization/em_steep.mdp \
           -c ionized.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em

# Step 3: NVT equilibration (100 ps)
gmx grompp -f mdp_templates/amber14sb/equilibration/nvt_equil.mdp \
           -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

# Step 4: NPT equilibration (1 ns)
gmx grompp -f mdp_templates/amber14sb/equilibration/npt_equil.mdp \
           -c nvt.gro -p topol.top -o npt.tpr \
           -t nvt.cpt -maxwarn 1
gmx mdrun -deffnm npt

# Step 5: Production MD (500 ns in 100 ns segments)
bash bin/gmx_continue_grompp.sh md 0 0 500 100 \
     mdp_templates/amber14sb/production/md_npt.mdp npt 0.002 "MyProtein"
```

### 2. Post-Simulation Analysis

```bash
#!/bin/bash
# analysis_workflow.sh - Common analysis tasks

TRAJ="md_complete.xtc"
TPR="npt.tpr"

# Concatenate trajectory segments first
gmx trjcat -f md_0_*_*.xtc -o $TRAJ -cat

# 1. RMSD calculation
echo "Backbone Backbone" | gmx rms -s $TPR -f $TRAJ -o rmsd.xvg

# 2. Plot RMSD
python bin/plot_xvg.py rmsd.xvg --style lines \
    --title "Backbone RMSD" --xlabel "Time (ns)" --ylabel "RMSD (nm)" \
    --output rmsd_plot.png

# 3. Radius of gyration
echo "Protein" | gmx gyrate -s $TPR -f $TRAJ -o gyrate.xvg

# 4. RMSF calculation
echo "Backbone" | gmx rmsf -s $TPR -f $TRAJ -o rmsf.xvg -res

# 5. Secondary structure
echo "Protein" | gmx do_dssp -s $TPR -f $TRAJ -o ss.xpm -sc scount.xvg

# 6. Hydrogen bonds
echo "Protein Protein" | gmx hbond -s $TPR -f $TRAJ -num hbnum.xvg
```

### 3. PCA Analysis Workflow

```bash
#!/bin/bash
# pca_workflow.sh - Principal component analysis

TRAJ="md_complete.xtc"
TPR="npt.tpr"

# 1. Run PCA on backbone
python bin/gromacs_pca.py -s $TPR -f $TRAJ \
    --selection "Backbone" \
    --first 1 --last 5 \
    -o pca_backbone

# 2. Generate movies for top 3 PCs
for pc in 1 2 3; do
    python bin/gromacs_pca_movie.py \
        -s pca_backbone/fit.gro \
        -v pca_backbone/eigenvectors.trr \
        -o pca_backbone/pc${pc}_movie.pdb \
        --pc $pc --nframes 100
done

# 3. Custom PC projection plot
python bin/plot_xvg.py pca_backbone/projection.xvg \
    --scatter --colormap viridis --aspect equal \
    --title "PC1 vs PC2" --xlabel "PC1 (nm)" --ylabel "PC2 (nm)" \
    --output pc1_pc2_scatter.pdf --dpi 300

# 4. PC timeseries
python bin/plot_xvg.py pca_backbone/projection.xvg --columns 1 \
    --legends "PC1" --style lines \
    --title "PC1 Evolution" --xlabel "Frame" --ylabel "PC1 (nm)" \
    --output pc1_timeseries.png
```

### 4. Multi-System Comparison

```bash
#!/bin/bash
# compare_systems.sh - Compare WT vs mutant

# System names
SYSTEMS=("wildtype" "mutant_F6A" "mutant_L10A")

# Run simulations for each system
for sys in "${SYSTEMS[@]}"; do
    # Preparation
    python bin/gromacs_prepare.py ${sys}.pdb \
        --forcefield amber14sb -o prep_${sys}
    
    cd prep_${sys}
    
    # Minimization + equilibration
    gmx grompp -f ../mdp_templates/amber14sb/minimization/em_steep.mdp \
               -c ionized.gro -p topol.top -o em.tpr
    gmx mdrun -deffnm em
    
    gmx grompp -f ../mdp_templates/amber14sb/equilibration/nvt_equil.mdp \
               -c em.gro -p topol.top -o nvt.tpr
    gmx mdrun -deffnm nvt
    
    gmx grompp -f ../mdp_templates/amber14sb/equilibration/npt_equil.mdp \
               -c nvt.gro -p topol.top -o npt.tpr -t nvt.cpt
    gmx mdrun -deffnm npt
    
    # Production
    bash ../bin/gmx_continue_grompp.sh md 0 0 200 50 \
         ../mdp_templates/amber14sb/production/md_npt.mdp npt 0.002 "$sys"
    
    cd ..
done

# Analyze RMSD for all systems
for sys in "${SYSTEMS[@]}"; do
    cd prep_${sys}
    gmx trjcat -f md_0_*_*.xtc -o md_complete.xtc -cat
    echo "Backbone Backbone" | gmx rms -s npt.tpr -f md_complete.xtc -o rmsd_${sys}.xvg
    cd ..
done

# Compare RMSD plots
python bin/plot_xvg.py \
    prep_wildtype/rmsd_wildtype.xvg \
    prep_mutant_F6A/rmsd_mutant_F6A.xvg \
    prep_mutant_L10A/rmsd_mutant_L10A.xvg \
    --multi --legends "WT" "F6A" "L10A" \
    --style lines --title "Backbone RMSD Comparison" \
    --output rmsd_comparison.png
```

### 5. Multiple Replicas Workflow

```bash
#!/bin/bash
# replicas_workflow.sh - Run multiple replicas

N_REPLICAS=3
BASENAME="md"

# After equilibration (npt.gro, npt.cpt exist)

# Run replicas in parallel
for rep in $(seq 0 $(($N_REPLICAS - 1))); do
    (
        bash bin/gmx_continue_grompp.sh $BASENAME $rep 0 200 50 \
             mdp_templates/amber14sb/production/md_npt.mdp npt 0.002 "rep${rep}"
    ) &
done

wait

echo "All replicas completed"

# Analyze each replica
for rep in $(seq 0 $(($N_REPLICAS - 1))); do
    gmx trjcat -f ${BASENAME}_${rep}_*_*.xtc -o replica${rep}.xtc -cat
    echo "Backbone Backbone" | gmx rms -s npt.tpr -f replica${rep}.xtc -o rmsd_rep${rep}.xvg
done

# Compare replicas
python bin/plot_xvg.py rmsd_rep*.xvg --multi \
    --legends "Replica 0" "Replica 1" "Replica 2" \
    --style lines --title "RMSD Across Replicas" \
    --output replicas_rmsd.png
```

### 6. Trajectory Preprocessing: PBC Removal and Centering

```bash
#!/bin/bash
# trajectory_preprocessing.sh - Remove PBC artifacts and center

TRAJ="md_complete.xtc"
TPR="npt.tpr"

# 1. Remove PBC artifacts (make molecules whole)
echo "Protein System" | gmx trjconv -s $TPR -f $TRAJ \
    -o traj_nojump.xtc -pbc nojump

# 2. Center protein in box and remove PBC
echo "Protein System" | gmx trjconv -s $TPR -f $TRAJ \
    -o traj_centered.xtc -pbc mol -center

# 3. Center and make compact (protein in center, minimum image)
echo "Protein System" | gmx trjconv -s $TPR -f $TRAJ \
    -o traj_compact.xtc -pbc mol -center -ur compact

# 4. Full preprocessing: whole molecules + center + compact
echo "Protein Protein" | gmx trjconv -s $TPR -f $TRAJ \
    -o traj_processed.xtc -pbc mol -center -ur compact

# Note: First selection = group to center, Second = output group
```

**PBC options explained:**
- `nojump` - Make molecules whole (no broken molecules across boundaries)
- `mol` - Put molecules back in box
- `res` - Put residues back in box
- `atom` - Put atoms back in box

**Unit cell representation (-ur):**
- `rect` - Rectangular box (default)
- `compact` - Compact representation (minimum image)
- `tric` - Triclinic box

### 7. Trajectory Fitting and Alignment

```bash
#!/bin/bash
# trajectory_fitting.sh - Fit trajectory to remove rotation/translation

TRAJ="md_complete.xtc"
TPR="npt.tpr"

# 1. Fit to initial structure (standard approach)
echo "Backbone Protein" | gmx trjconv -s $TPR -f $TRAJ \
    -o traj_fit.xtc -fit rot+trans

# 2. Progressive fitting (each frame to previous frame)
# Better for large conformational changes
echo "Backbone Protein" | gmx trjconv -s $TPR -f $TRAJ \
    -o traj_progressive.xtc -fit progressive

# 3. Fit to average structure
# First: create average structure
echo "Protein" | gmx trjconv -s $TPR -f $TRAJ -o average.pdb -dump 0

# Then: fit all frames to average
echo "Backbone Protein" | gmx trjconv -s average.pdb -f $TRAJ \
    -o traj_fit_avg.xtc -fit rot+trans

# 4. Rotation only (no translation)
echo "Backbone Protein" | gmx trjconv -s $TPR -f $TRAJ \
    -o traj_rot.xtc -fit rotxy+transxy

# 5. Complete preprocessing + fitting workflow
echo "Backbone Protein" | gmx trjconv -s $TPR -f $TRAJ \
    -o traj_final.xtc -pbc mol -center -ur compact -fit rot+trans

# Note: First selection = fitting group, Second = output group
```

**Fitting modes:**
- `rot+trans` - Remove rotation and translation (standard)
- `rotxy+transxy` - Remove rotation/translation in XY plane only
- `progressive` - Fit each frame to previous frame
- `none` - No fitting

**Common use cases:**
- **Standard analysis**: `rot+trans` with backbone fitting
- **Large conformational changes**: `progressive` fitting
- **Membrane systems**: `rotxy+transxy` (keep Z-dimension)
- **Visualization**: Combine with `-pbc mol -center`

### 8. Trajectory Format Conversion

```bash
#!/bin/bash
# trajectory_conversion.sh - Convert between trajectory formats

TRAJ="md_complete.xtc"
TPR="npt.tpr"

# 1. Convert XTC to PDB (entire trajectory)
echo "Protein" | gmx trjconv -s $TPR -f $TRAJ -o trajectory.pdb

# 2. Convert to PDB with specific time range
echo "Protein" | gmx trjconv -s $TPR -f $TRAJ -o traj_100-200ns.pdb \
    -b 100000 -e 200000  # Time in ps

# 3. Convert to TRR (full precision)
echo "System" | gmx trjconv -s $TPR -f $TRAJ -o trajectory.trr

# 4. Extract specific frames
# Single frame (e.g., frame at 50 ns)
echo "Protein" | gmx trjconv -s $TPR -f $TRAJ -o frame_50ns.pdb -dump 50000

# Every 10th frame
echo "Protein" | gmx trjconv -s $TPR -f $TRAJ -o traj_skip10.xtc -skip 10

# Every 100 ps
echo "Protein" | gmx trjconv -s $TPR -f $TRAJ -o traj_dt100.xtc -dt 100

# 5. Convert with selection (e.g., only protein without water)
echo "Protein" | gmx trjconv -s $TPR -f $TRAJ -o protein_only.pdb

# 6. High-quality PDB for VMD/PyMOL
echo "Protein" | gmx trjconv -s $TPR -f $TRAJ -o movie.pdb \
    -pbc mol -center -ur compact -fit rot+trans -skip 5

# 7. Convert GRO to PDB
gmx trjconv -s structure.gro -f structure.gro -o structure.pdb

# 8. Create multi-model PDB (for NMR-style visualization)
echo "Protein" | gmx trjconv -s $TPR -f $TRAJ -o multimodel.pdb -sep
```

**Format options:**
- `.xtc` - Compressed, lossy (typical for MD)
- `.trr` - Full precision, larger files
- `.pdb` - Visualization, single/multi-model
- `.gro` - GROMACS structure format
- `.dcd` - CHARMM/NAMD format

**Useful conversion flags:**
- `-b` / `-e` - Begin/end time (ps)
- `-dt` - Output every X ps
- `-skip` - Output every Nth frame
- `-dump` - Extract single frame at time X
- `-sep` - Separate file for each frame
- `-pbc mol -center -ur compact` - Clean up for visualization

**Example: Prepare trajectory for publication figure**
```bash
# 1. Extract key transition (100-120 ns)
# 2. Fit and center
# 3. Output only protein
# 4. Reduce to 1 frame per ns
echo "Backbone Protein" | gmx trjconv -s npt.tpr -f md_complete.xtc \
    -o figure_trajectory.pdb \
    -b 100000 -e 120000 \
    -dt 1000 \
    -pbc mol -center -ur compact -fit rot+trans
```

## Quick Reference Commands

### System Preparation
```bash
# Preview workflow
python bin/gromacs_prepare.py protein.pdb --dry-run

# Execute with CHARMM36
python bin/gromacs_prepare.py protein.pdb --forcefield charmm36
```

### Plotting
```bash
# Basic line plot
python bin/plot_xvg.py data.xvg --style lines

# Multiple files comparison
python bin/plot_xvg.py file1.xvg file2.xvg --multi --legends "A" "B"

# Scatter with color gradient
python bin/plot_xvg.py data.xvg --scatter --colormap viridis

# Specific columns
python bin/plot_xvg.py data.xvg --columns 1 3 --legends "Col1" "Col3"
```

### MD Continuation
```bash
# Basic usage
bash bin/gmx_continue_grompp.sh md 0 0 500 100 production.mdp npt

# With custom title
bash bin/gmx_continue_grompp.sh md 0 0 500 100 production.mdp npt 0.002 "MySystem"
```

### PCA
```bash
# Full PCA workflow
python bin/gromacs_pca.py -s topol.tpr -f traj.xtc --selection "Backbone"

# Generate PC movie
python bin/gromacs_pca_movie.py -s pca_out/fit.gro -v pca_out/eigenvectors.trr -o pc1.pdb --pc 1
```

## Tips

1. **Always use dry-run first** to preview commands
2. **Keep organized** - separate directories for different systems
3. **Monitor disk space** - trajectories can be large
4. **Use descriptive names** - include system/replica info in titles
5. **Save intermediate files** - checkpoints allow recovery
6. **Concatenate carefully** - use `gmx trjcat` with `-cat` flag
7. **Check convergence** - plot RMSD, energy before analysis

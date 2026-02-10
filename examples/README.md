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

# gmx_continue_grompp.sh

Intelligent bash script for running long MD simulations in sequential segments with automatic crash recovery.

## Overview

Automates the continuation of molecular dynamics simulations by:
- Running simulations in manageable segments
- Automatic checkpoint-based continuation
- Intelligent crash recovery (just re-run the same command)
- Proper time management (tinit, nsteps)
- Optional custom titles for tracking

## Features

- ✅ **Segmented simulations**: Break long runs into manageable chunks
- ✅ **Automatic continuation**: Uses previous segment's checkpoint
- ✅ **Crash recovery**: Detects incomplete segments and resumes
- ✅ **Checkpoint handling**: Works with or without checkpoint files
- ✅ **Time management**: Automatically sets tinit and nsteps
- ✅ **Template-based**: Uses MDP templates from `mdp_templates/`
- ✅ **Flexible**: Optional parameters for customization
- ✅ **Idempotent**: Re-running same command is safe
- ✅ **Command logging** (on by default): every gmx command run, plus the exact invocation used,
  is saved to `gmx_continue_grompp_commands.sh` / `gmx_continue_extend_commands.sh`; the file from an earlier run is kept as `…_commands.1.sh`, `.2.sh`, … (highest number = most recent)
  (override with `--commands-file`)

## Quick Start

```bash
# Run 500 ns in 100 ns segments, starting from npt.gro/cpt/edr
bash bin/gmx_continue_grompp.sh md 0 0 500 100 \
    --template mdp_templates/amber14sb/production/md_npt.mdp --initial npt

# After crash: just re-run the exact same command
bash bin/gmx_continue_grompp.sh md 0 0 500 100 \
    --template mdp_templates/amber14sb/production/md_npt.mdp --initial npt
```

## Usage

### Basic Syntax

```bash
gmx_continue_grompp.sh <basename> <replica> <start> <end> <dt> [OPTIONS]
```

The five positional arguments come first; everything else is an option starting with `--`
(any other extra argument is an error).

### Positional Arguments

| Argument | Description |
|----------|-------------|
| `basename` | Base name for output files (e.g., "md") |
| `replica` | Replica number (e.g., 0, 1, 2) |
| `start` | Start time in ns (0 for a new run, >0 to extend) |
| `end` | End time in ns (e.g., 500) |
| `dt` | Segment length in ns (e.g., 100) |

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--template FILE` | `md.mdp` | Production MDP template (see [Which MDP is used](#which-mdp-is-used)) |
| `--initial NAME` | `npt` | Starting files of the first segment: `NAME.gro`, `NAME.edr`, `NAME.cpt` (optional) |
| `--timestep PS` | 0.002 | Integration timestep in ps, used to compute `nsteps` |
| `--title TEXT` | - | Replaces the MDP's `title` |
| `--topology FILE` | `topol.top` | Topology file |
| `--posres-ref FILE` | - | Reference structure for position restraints (`grompp -r FILE`); needed when the topology/MDP define position restraints (e.g. `-DPOSRES`) |
| `--plumed FILE` | - | PLUMED input file, passed to `mdrun -plumed` |
| `--ps` | off | Times are in ps instead of ns |
| `--ignore-initial-edr` | off | Don't pass `NAME.edr` to grompp (e.g. when the barostat changes, Berendsen → Parrinello-Rahman) |
| `--commands-file FILE` | `gmx_continue_grompp_commands.sh` | Reproducibility script of every gmx command run |

### Which MDP is used

- **First segment** (`start` = 0): `--template` if given; otherwise a copy of `NAME.mdp` (the
  equilibration's own MDP, `--initial NAME`) if it exists; otherwise `md.mdp`. **Give
  `--template` when starting from equilibration**, or the run will use the equilibration's
  settings (position restraints, barostat, ...).
- **Later segments**: a copy of the previous segment's MDP (so extending needs no template), or
  `--template` if given.

In every case `tinit` and `nsteps` are set for the segment, and `--title` replaces the title.

### File Naming Convention

Segments are named: `{basename}_{replica}_{start}_{end}`

Example with `basename=md`, `replica=0`, 100 ns segments:
- `md_0_0_100.gro/cpt/edr` (0-100 ns)
- `md_0_100_200.gro/cpt/edr` (100-200 ns)
- `md_0_200_300.gro/cpt/edr` (200-300 ns)
- etc.

The output of grompp and mdrun for all segments is appended to `{basename}_{replica}.log`.

## Examples

### Basic Usage

```bash
# 500 ns in 100 ns segments: md.mdp (or npt.mdp, if present) and npt.gro/cpt/edr
bash bin/gmx_continue_grompp.sh md 0 0 500 100
```

### With Template MDP

```bash
# AMBER14SB production template
bash bin/gmx_continue_grompp.sh md 0 0 500 100 \
    --template mdp_templates/amber14sb/production/md_npt.mdp

# CHARMM36 production template
bash bin/gmx_continue_grompp.sh md 0 0 500 100 \
    --template mdp_templates/charmm36/production/md_npt.mdp
```

### Custom Initial Files

```bash
# Start from nvt.gro/cpt/edr instead of npt.gro/cpt/edr
bash bin/gmx_continue_grompp.sh md 0 0 500 100 \
    --template mdp_templates/amber14sb/production/md_npt.mdp --initial nvt
```

### Custom Timestep

```bash
# 4 fs timestep (0.004 ps), e.g. with hydrogen mass repartitioning
bash bin/gmx_continue_grompp.sh md 0 0 500 100 \
    --template mdp_templates/amber14sb/production/md_npt.mdp --timestep 0.004
```

The timestep only sets `nsteps`; the MDP's `dt` must match it.

### Topology, Position Restraints and PLUMED

```bash
bash bin/gmx_continue_grompp.sh md 0 0 500 100 \
    --template prod.mdp --topology system.top \
    --posres-ref npt.gro --plumed plumed.dat
```

### With Custom Title

```bash
# The MDP's title becomes "TRP_cage_WT"
bash bin/gmx_continue_grompp.sh md 0 0 500 100 \
    --template mdp_templates/amber14sb/production/md_npt.mdp --title "TRP_cage_WT"
```

### Multiple Replicas

```bash
for rep in 0 1 2; do
    bash bin/gmx_continue_grompp.sh md $rep 0 500 100 \
        --template prod.mdp --initial npt --title "replica_$rep" &
done
wait
```

Each replica writes its own files (`md_0_*`, `md_1_*`, ...) and log, but they share the command
file; use `--commands-file md_${rep}_commands.sh` to keep one per replica. Running several in
the same folder also shares the `STOP` file.

### Extending from a Later Time

```bash
# Extend from 10000 to 20000 ns: copies md_0_9500_10000.mdp
bash bin/gmx_continue_grompp.sh md 0 10000 20000 500
```

## Crash Recovery

### How It Works

The script automatically detects:
1. **Completed segments**: Has `.gro`, `.edr`, `.log` with "Finished mdrun"
2. **Incomplete segments**: Has `.tpr` but segment didn't complete
3. **Not started**: No `.tpr` file exists yet

### Recovery Process

```bash
# Initial run (crashes at segment 2)
bash bin/gmx_continue_grompp.sh md 0 0 500 100 --template prod.mdp --initial npt

# Segments completed: 0-100, 100-200
# Segment crashed: 200-300

# Recovery: Just re-run the EXACT SAME COMMAND
bash bin/gmx_continue_grompp.sh md 0 0 500 100 --template prod.mdp --initial npt
```

**What happens:**
```
✓ Segment md_0_0_100: Already complete, skipping
✓ Segment md_0_100_200: Already complete, skipping
⚠ Segment md_0_200_300: Incomplete, resuming from checkpoint
→ gmx mdrun -deffnm md_0_200_300 -cpi md_0_200_300.cpt
✓ Segment md_0_300_400: Running normally
✓ Segment md_0_400_500: Running normally
```

### No Checkpoint Available

If a checkpoint file is missing (rare), the script:
- Warns you
- Continues using only `.gro` and `.edr` files
- Takes the velocities from the `.gro` file (if it has them), or generates them if the MDP has `gen_vel = yes`

## Extending Simulations

To extend an existing simulation:

```bash
# Original: 0-500 ns
bash bin/gmx_continue_grompp.sh md 0 0 500 100 --template prod.mdp --initial npt

# Later: Extend to 1000 ns (same command, just change end time)
bash bin/gmx_continue_grompp.sh md 0 0 1000 100 --template prod.mdp --initial npt
```

The script automatically:
- Skips completed segments (0-500)
- Only runs new segments (500-1000)

## Required Files

### Before First Segment

- `topol.top` - System topology
- `npt.gro` - Equilibrated structure (or `NAME.gro` with `--initial NAME`)
- `npt.edr` - Energy file from equilibration
- `npt.cpt` - Checkpoint (optional but recommended)
- Production MDP template (`--template`)

### During Continuation

The script automatically uses outputs from previous segments.

## MDP Template Requirements

Production MDP templates should have:
```
nsteps = 0    ; Will be set automatically
tinit = 0     ; Will be set automatically
```

See `mdp_templates/` for examples.

## Output Files Per Segment

Each segment produces:
- `{name}.mdp` - Modified MDP (tinit, nsteps set)
- `{name}.tpr` - Run input file
- `{name}.gro` - Final coordinates
- `{name}.cpt` - Checkpoint for continuation
- `{name}.edr` - Energy file
- `{name}.log` - Run log
- `{name}.xtc` - Compressed trajectory
- `{name}.trr` - Full precision trajectory (if enabled)

## Workflow Integration

### Complete Workflow

```bash
# 1. System preparation
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

# 5. Production (500 ns in 100 ns segments)
bash bin/gmx_continue_grompp.sh md 0 0 500 100 \
     --template mdp_templates/amber14sb/production/md_npt.mdp --initial npt
```

### Concatenating Trajectories

After completion, concatenate segment trajectories:

```bash
# Concatenate XTC files
gmx trjcat -f md_0_*_*.xtc -o md_complete.xtc -cat

# Concatenate EDR files
gmx eneconv -f md_0_*_*.edr -o md_complete.edr
```

## Tips & Best Practices

1. **Segment Length**: 
   - 100 ns segments: Good for most simulations
   - 50 ns segments: Better crash recovery, more overhead
   - 200 ns segments: Less overhead, longer recovery time

2. **Checkpoint Files**:
   - Ensure `nstxout-compressed` in MDP is reasonable
   - Checkpoints are written at end of run + periodically

3. **Crash Recovery**:
   - Just re-run the exact same command
   - No need to calculate where it stopped

4. **Multiple Replicas**:
   - Use different replica numbers (0, 1, 2, etc.)
   - Add custom titles for easier tracking

5. **Disk Space**:
   - Monitor space, especially for long simulations
   - Consider deleting intermediate `.trr` files if not needed

6. **Time Management**:
   - The script handles `tinit` and `nsteps` automatically
   - No manual MDP editing needed

## Troubleshooting

### "Error: Template MDP not found"
- Check path to template MDP file
- Use absolute or relative path from working directory

### "Error: Missing initial files"
- Ensure `npt.gro` and `npt.edr` exist (or `NAME.gro`/`NAME.edr` with `--initial NAME`)
- Check that NPT equilibration completed successfully

### "Error: Missing output files from previous segment"
- Previous segment may have failed
- Check log files for errors
- May need to restart from earlier segment

### Simulation runs but no checkpoint
- Check `nsteps` calculation (should be > 0)
- Verify simulation actually ran (check log file)
- Very short segments may not write checkpoint during run

## Advanced: Custom STOP File

Create a `STOP` file in the working directory to gracefully halt:

```bash
# While simulation is running
touch STOP

# Script will finish current segment and exit
```

# Testing gmx_continue_extend.sh

This directory contains test scripts for `gmx_continue_extend.sh`. Test data files are **not committed** to the repository (see `.gitignore`).

## Setup Test Data

You'll need a small, fast system for testing. TRP-cage is ideal (58 residues, ~1000 atoms).

### Option 1: Use Your Existing TRP-cage Files

Copy your TRP-cage simulation files here:
```bash
mkdir -p tests/trp_cage
cp /path/to/your/trp_cage/{md_0.tpr,md_0.cpt,md_0.xtc,md_0.edr,md_0.log,md_0.gro} tests/trp_cage/
```

### Option 2: Download TRP-cage from PDB

```bash
mkdir -p tests/trp_cage
cd tests/trp_cage

# Download structure
wget https://files.rcsb.org/download/1L2Y.pdb

# Set up system (requires GROMACS)
gmx pdb2gmx -f 1L2Y.pdb -o conf.gro -water tip3p -ff amber99sb-ildn
gmx editconf -f conf.gro -o box.gro -c -d 1.0 -bt cubic
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top
gmx grompp -f em.mdp -c solvated.gro -p topol.top -o em.tpr
# ... (run minimization, equilibration, short MD to get TPR/CPT)
```

## Test Scenarios

### Test 1: Easy Mode (auto-detect single TPR)
```bash
cd tests/trp_cage
../../bin/gmx_continue_extend.sh 1000  # Extend by 1 ns
```

**Expected:**
- Auto-detects `md_0.tpr` and `md_0.cpt`
- Extends TPR by 1000 ps
- Runs mdrun in append mode

---

### Test 2: Easy Mode (explicit TPR, multiple files)
```bash
cd tests/trp_cage
# Assume both md_0.tpr and npt.tpr exist
../../bin/gmx_continue_extend.sh 1000 md_0.tpr
```

**Expected:**
- Uses specified `md_0.tpr`
- Extends and continues as above

---

### Test 3: Normal Mode - Noappend with Segmentation
```bash
../../bin/gmx_continue_extend.sh md 0 0 5 1 --tpr md_0.tpr --noappend
```

**Expected:**
- Uses existing `md_0.tpr`
- Runs 5 segments of 1 ns each
- Creates: `md_0.xtc`, `md_0.part0001.xtc`, `md_0.part0002.xtc`, ...
- Checkpoint remains: `md_0.cpt` (updated in place)

---

### Test 4: Normal Mode - Append with Segmentation
```bash
../../bin/gmx_continue_extend.sh md 0 0 5 1 --tpr md_0.tpr --append
```

**Expected:**
- Uses existing `md_0.tpr`
- Runs 5 segments of 1 ns each
- All output appended to same files: `md_0.xtc`, `md_0.edr`, `md_0.log`
- Useful for HPC scheduling control

---

### Test 5: Crash Recovery (Easy Mode)
```bash
cd tests/trp_cage
../../bin/gmx_continue_extend.sh 2000 md_0.tpr

# Kill with Ctrl+C during mdrun

# Re-run same command - should resume from checkpoint
../../bin/gmx_continue_extend.sh 2000 md_0.tpr
```

**Expected:**
- Second run detects existing checkpoint
- Resumes with `-cpi` flag

---

### Test 6: Crash Recovery (Segmented Mode)
```bash
../../bin/gmx_continue_extend.sh md 0 0 10 2 --tpr md_0.tpr --noappend

# Kill during segment 3

# Re-run same command
../../bin/gmx_continue_extend.sh md 0 0 10 2 --tpr md_0.tpr --noappend
```

**Expected:**
- Script detects completed segments (checks log files)
- Resumes from incomplete segment
- Uses checkpoint for crash recovery

---

## Automated Test Script

Run `./test_extend_modes.sh` to verify all modes work correctly (coming soon).

## Test Results

Document your test results:

| Test | Mode | Status | Notes |
|------|------|--------|-------|
| Easy mode (1 TPR) | Append | ✅ | |
| Easy mode (multiple TPR) | Append | ✅ | |
| Segmented noappend | Part files | ✅ | |
| Segmented append | Same file | ✅ | |
| Crash recovery easy | Append | ✅ | |
| Crash recovery segment | Part files | ✅ | |

## Cleanup

```bash
# Remove all test outputs, keep original TPR
cd tests/trp_cage
rm -f md_0.part*.* *.backup md_0_prev.cpt
```

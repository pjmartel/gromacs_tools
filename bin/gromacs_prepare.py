#!/usr/bin/env python3
"""
GROMACS Automation Script with command-line arguments and logging
Usage: python gromacs_prepare.py <pdb_file> [options]
"""

import argparse
import subprocess
import sys
import os
import logging
import re
from pathlib import Path
from datetime import datetime

# Global variable to store the command log script path
_command_script_path = None

# Error hints for common failures
ERROR_HINTS = {
    'pdb2gmx': [
        "- Check if PDB has missing atoms/residues",
        "- Verify force field supports all residues",
        "- Try different force field if current one fails",
        "- Check for non-standard residues or HETATM entries",
        "- If PDB contains hydrogens, use --ignore-hydrogens flag"
    ],
    'editconf': [
        "- Verify input structure file exists",
        "- Check box distance isn't too small",
        "- Try different box type if geometry issues occur"
    ],
    'solvate': [
        "- Check if box is large enough for solvent",
        "- Verify solvent file exists and is accessible",
        "- Check topology file wasn't corrupted in previous step"
    ],
    'grompp': [
        "- Check MDP file syntax",
        "- Verify all files referenced exist",
        "- Check topology is complete and valid"
    ],
    'genion': [
        "- Verify SOL group was detected correctly",
        "- Check if system has enough water for ion replacement",
        "- Try running 'gmx make_ndx' manually to inspect groups",
        "- Verify ion names match force field (NA vs NA+, etc.)"
    ]
}

# Force field and water model compatibility
FF_WATER_COMPATIBILITY = {
    'amber03': ['tip3p', 'tip4p', 'tip4pew', 'tip5p', 'spc', 'spce'],
    'amber94': ['tip3p', 'tip4p', 'spc', 'spce'],
    'amber96': ['tip3p', 'tip4p', 'spc', 'spce'],
    'amber99': ['tip3p', 'tip4p', 'tip4pew', 'spc', 'spce'],
    'amber99sb': ['tip3p', 'tip4p', 'tip4pew', 'spc', 'spce'],
    'amber99sb-ildn': ['tip3p', 'tip4pew', 'spce'],
    'amberGS': ['tip3p', 'tip4pew', 'spce'],
    'charmm27': ['tip3p', 'spc', 'spce'],
    'charmm36': ['tip3p'],
    'charmm36-feb2021': ['tip3p'],
    'charmm36-jul2021': ['tip3p'],
    'charmm36m': ['tip3p'],
    'gromos43a1': ['spc', 'spce'],
    'gromos43a2': ['spc', 'spce'],
    'gromos45a3': ['spc', 'spce'],
    'gromos53a5': ['spc', 'spce'],
    'gromos53a6': ['spc', 'spce'],
    'gromos54a7': ['spc', 'spce'],
    'oplsaa': ['tip3p', 'tip4p', 'tip4pew', 'tip5p', 'spc', 'spce']
}

def setup_logging(verbose=False):
    """Configure logging for the script.
    
    Parameters
    ----------
    verbose : bool
        If True, set logging level to DEBUG for detailed output
    """
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%H:%M:%S'
    )

def validate_forcefield_water_compatibility(forcefield, water_model):
    """Validate force field and water model compatibility.
    
    Parameters
    ----------
    forcefield : str
        Force field name
    water_model : str
        Water model name
        
    Returns
    -------
    bool
        True if compatible or unknown, False with warning if incompatible
    """
    logger = logging.getLogger(__name__)
    
    if forcefield in FF_WATER_COMPATIBILITY:
        if water_model not in FF_WATER_COMPATIBILITY[forcefield]:
            logger.warning(f"⚠️  Water model '{water_model}' may not be optimal with force field '{forcefield}'")
            logger.warning(f"    Recommended: {', '.join(FF_WATER_COMPATIBILITY[forcefield])}")
            return False
    else:
        logger.debug(f"Force field '{forcefield}' not in compatibility database, skipping check")
    
    return True

def validate_pdb_file(pdb_file):
    """Validate PDB file structure and contents.
    
    Parameters
    ----------
    pdb_file : str
        Path to PDB file
        
    Raises
    ------
    ValueError
        If PDB file is invalid or malformed
    """
    logger = logging.getLogger(__name__)
    
    # Check file exists
    if not os.path.exists(pdb_file):
        raise ValueError(f"PDB file not found: {pdb_file}")
    
    # Check file size
    file_size = os.path.getsize(pdb_file)
    if file_size < 100:
        raise ValueError(f"PDB file is suspiciously small ({file_size} bytes)")
    
    # Check for ATOM or HETATM records
    has_atoms = False
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                has_atoms = True
                break
    
    if not has_atoms:
        raise ValueError(f"PDB file doesn't contain ATOM or HETATM records")
    
    logger.info(f"✓ PDB file validation passed: {pdb_file}")

def create_minimal_mdp(filename="minimal.mdp"):
    """Create minimal MDP file for pre-processing.
    
    This MDP is only used for genion and doesn't need to be complex.
    It's not used for actual MD simulation.
    
    Parameters
    ----------
    filename : str
        Output filename for MDP file
    """
    content = """; Minimal MDP for pre-processing (genion)
; This is NOT used for actual MD simulation
; Only used to create TPR file for adding ions

integrator  = steep      ; Steepest descent minimizer
emtol       = 1000.0     ; Stop when max force < 1000 kJ/mol/nm
emstep      = 0.01       ; Energy step size (nm)
nsteps      = 1          ; Just need structure, not simulation
"""
    with open(filename, 'w') as f:
        f.write(content)

def parse_gro_file(gro_file):
    """Parse GRO file to extract system information.
    
    Parameters
    ----------
    gro_file : str
        Path to GRO file
        
    Returns
    -------
    dict
        Dictionary with atom count, box dimensions, etc.
    """
    if not os.path.exists(gro_file):
        return None
    
    try:
        with open(gro_file, 'r') as f:
            lines = f.readlines()
            
        if len(lines) < 3:
            return None
            
        n_atoms = int(lines[1].strip())
        box_line = lines[-1].strip().split()
        
        info = {
            'n_atoms': n_atoms,
            'box_x': float(box_line[0]) if len(box_line) > 0 else 0.0,
            'box_y': float(box_line[1]) if len(box_line) > 1 else 0.0,
            'box_z': float(box_line[2]) if len(box_line) > 2 else 0.0,
        }
        
        # Calculate volume
        info['volume'] = info['box_x'] * info['box_y'] * info['box_z']
        
        return info
    except (ValueError, IndexError):
        return None

def parse_topology_file(top_file):
    """Parse topology file to extract molecule counts.
    
    Parameters
    ----------
    top_file : str
        Path to topology file
        
    Returns
    -------
    dict
        Dictionary with molecule types and counts
    """
    if not os.path.exists(top_file):
        return {}
    
    molecules = {}
    in_molecules_section = False
    
    try:
        with open(top_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Skip comments and empty lines
                if not line or line.startswith(';'):
                    continue
                
                # Check for molecules section
                if line.startswith('[ molecules ]'):
                    in_molecules_section = True
                    continue
                
                # Exit molecules section on new section
                if line.startswith('[') and in_molecules_section:
                    break
                
                # Parse molecule counts
                if in_molecules_section:
                    parts = line.split()
                    if len(parts) >= 2:
                        mol_name = parts[0]
                        mol_count = int(parts[1])
                        molecules[mol_name] = mol_count
    except (ValueError, IndexError):
        pass
    
    return molecules

def print_system_summary(final_gro, topology_file):
    """Print comprehensive system summary after pipeline completion.
    
    Parameters
    ----------
    final_gro : str
        Path to final GRO structure file
    topology_file : str
        Path to topology file
    """
    logger = logging.getLogger(__name__)
    
    logger.info("")
    logger.info("="*60)
    logger.info("📊 SYSTEM SUMMARY")
    logger.info("="*60)
    
    # Parse GRO file
    gro_info = parse_gro_file(final_gro)
    if gro_info:
        logger.info(f"Total atoms:      {gro_info['n_atoms']:,}")
        logger.info(f"Box dimensions:   {gro_info['box_x']:.3f} × {gro_info['box_y']:.3f} × {gro_info['box_z']:.3f} nm")
        logger.info(f"Box volume:       {gro_info['volume']:.3f} nm³")
    
    # Parse topology
    molecules = parse_topology_file(topology_file)
    if molecules:
        logger.info("")
        logger.info("Molecule composition:")
        for mol_name, count in molecules.items():
            logger.info(f"  {mol_name:15s} {count:,} molecules")
    
    logger.info("="*60)

def set_command_script_path(path):
    """Set the global path for logging commands."""
    global _command_script_path
    _command_script_path = path
    # Initialize the script with a header
    with open(path, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write("# GROMACS preparation pipeline commands\n")
        f.write("# Generated by gromacs_prepare.py\n")
        f.write("#\n")
        f.write("# This script reproduces all GROMACS commands executed by the Python script.\n")
        f.write("# You can run this independently to get the same results.\n")
        f.write("#\n\n")
        f.write("set -e  # Exit on error\n\n")

def log_command_to_script(command, comment=None):
    """Append a command to the reproducibility script."""
    global _command_script_path
    if _command_script_path is None:
        return
    
    with open(_command_script_path, 'a') as f:
        if comment:
            f.write(f"\n# {comment}\n")
        f.write(f"{command}\n")

def run_gromacs_command(command, description, log_filename, dry_run=False, step_info=None):
    """Run a GROMACS command, capture output to log file, and handle errors
    
    Parameters
    ----------
    command : str
        The command to execute
    description : str
        Human-readable description of the command
    log_filename : str
        Path to log file for output
    dry_run : bool
        If True, only log the command without executing it
    step_info : str, optional
        Step progress information (e.g., "[1/5]")
    """
    logger = logging.getLogger(__name__)
    
    step_prefix = f"{step_info} " if step_info else ""
    
    logger.info("")
    logger.info("="*60)
    if dry_run:
        logger.info(f"{step_prefix}[DRY RUN] {description}")
    else:
        logger.info(f"{step_prefix}{description}")
    logger.debug(f"Command: {command}")
    if not dry_run:
        logger.debug(f"Log file: {log_filename}")
    logger.info("="*60)
    
    # Log to reproducibility script
    log_command_to_script(command, description)
    
    # In dry-run mode, skip execution
    if dry_run:
        logger.info("[DRY RUN] Command logged but not executed")
        return True
    
    # Write command and timestamp to log file
    with open(log_filename, 'w') as log_file:
        log_file.write(f"GROMACS Automation Log\n")
        log_file.write(f"Tool: {description}\n")
        log_file.write(f"Timestamp: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        log_file.write(f"Command: {command}\n")
        log_file.write(f"{'='*60}\n\n")
    
    try:
        # Run command and capture both stdout and stderr
        result = subprocess.run(command, shell=True, check=True, 
                              capture_output=True, text=True)
        
        # Append output to log file
        with open(log_filename, 'a') as log_file:
            if result.stdout:
                log_file.write("STDOUT:\n")
                log_file.write(result.stdout)
                log_file.write("\n")
            if result.stderr:
                log_file.write("STDERR:\n")
                log_file.write(result.stderr)
                log_file.write("\n")
            log_file.write(f"\nReturn code: {result.returncode}\n")
            log_file.write("✓ Command completed successfully\n")
        
        logger.info("✓ Success!")
        return True
        
    except subprocess.CalledProcessError as e:
        # Append error information to log file
        with open(log_filename, 'a') as log_file:
            if e.stdout:
                log_file.write("STDOUT:\n")
                log_file.write(e.stdout)
                log_file.write("\n")
            if e.stderr:
                log_file.write("STDERR:\n")
                log_file.write(e.stderr)
                log_file.write("\n")
            log_file.write(f"\nReturn code: {e.returncode}\n")
            log_file.write("✗ Command failed!\n")
        
        logger.error(f"❌ {description} failed!")
        logger.error(f"📋 Check log file: {log_filename}")
        logger.error(f"Return code: {e.returncode}")
        
        # Print troubleshooting hints
        if description in ERROR_HINTS:
            logger.error("")
            logger.error("💡 Troubleshooting hints:")
            for hint in ERROR_HINTS[description]:
                logger.error(f"   {hint}")
        
        return False

def get_solvent_group_number(tpr_file):
    """Get the group number for solvent (SOL) from the TPR using gmx make_ndx.

    Returns the numeric group id as a string (e.g., '13') if found, else None.
    """
    logger = logging.getLogger(__name__)
    
    try:
        result = subprocess.run(f"gmx make_ndx -f {tpr_file} -o index.ndx", 
                                shell=True, check=True, 
                                input="q\n", capture_output=True, text=True)
        output = result.stdout
        for line in output.splitlines():
            if "SOL" in line:
                parts = line.split()
                if parts and parts[0].isdigit():
                    logger.debug(f"Found SOL group: {parts[0]}")
                    return parts[0]
    except subprocess.CalledProcessError as e:
        logger.error("Error determining solvent group number.")
    return None

def automate_gromacs(args):
    """Automate GROMACS preparation steps for MD simulations.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments containing at least:
        pdb_file, forcefield, water_model, box_type, box_distance,
        solvent_file, cation, anion, dry_run.
    """
    logger = logging.getLogger(__name__)
    
    # Initialize command logging script
    commands_file = Path.cwd() / "commands.sh"
    set_command_script_path(commands_file)
    
    # Validate inputs
    logger.info("")
    logger.info("="*60)
    logger.info("GROMACS SYSTEM PREPARATION PIPELINE")
    logger.info("="*60)
    
    if args.dry_run:
        logger.info("")
        logger.info("🔍 DRY RUN MODE: Commands will be logged but not executed")
        logger.info("")
        logger.info("📋 Configuration:")
        logger.info(f"  Input PDB:      {args.pdb_file}")
        logger.info(f"  Force field:    {args.forcefield}")
        logger.info(f"  Water model:    {args.water_model}")
        logger.info(f"  Box type:       {args.box_type}")
        logger.info(f"  Box distance:   {args.box_distance} nm")
        logger.info(f"  Ions:           {args.cation}/{args.anion}")
        logger.info(f"  Neutralize:     {'No' if args.no_neutral else 'Yes'}")
    else:
        logger.info(f"Input: {args.pdb_file}")
        logger.info(f"Force field: {args.forcefield}, Water: {args.water_model}")
    
    logger.info(f"")
    logger.info(f"All commands will be logged to: {commands_file}")
    
    # Validate force field and water model compatibility
    validate_forcefield_water_compatibility(args.forcefield, args.water_model)
    
    # Validate PDB file
    try:
        validate_pdb_file(args.pdb_file)
    except ValueError as e:
        logger.error(f"❌ PDB validation failed: {e}")
        return False

    # Derive stem from input pdb file
    input_pdb_path = Path(args.pdb_file)
    stem = input_pdb_path.stem

    # Total pipeline steps
    total_steps = 5
    
    # Step 1: pdb2gmx
    topology_file = f"{stem}.top"
    gro_output = f"{stem}.gro"
    pdb2gmx_log = f"{stem}_pdb2gmx.log"
    
    # Build pdb2gmx command with optional -ignh flag
    ignh_flag = " -ignh" if getattr(args, "ignore_hydrogens", False) else ""
    pdb2gmx_cmd = (
        f"gmx pdb2gmx -f {args.pdb_file} -o {gro_output} -p {topology_file} "
        f"-ff {args.forcefield} -water {args.water_model}{ignh_flag}"
    )
    if not run_gromacs_command(pdb2gmx_cmd, "Generate topology (pdb2gmx)", pdb2gmx_log, 
                               dry_run=args.dry_run, step_info=f"[1/{total_steps}]"):
        return False

    # Step 2: editconf (define box)
    box_output = f"{stem}_box.gro"
    editconf_log = f"{stem}_editconf.log"
    editconf_cmd = (
        f"gmx editconf -f {gro_output} -o {box_output} -c -d {args.box_distance} "
        f"-bt {args.box_type}"
    )
    if not run_gromacs_command(editconf_cmd, "Define simulation box (editconf)", editconf_log, 
                               dry_run=args.dry_run, step_info=f"[2/{total_steps}]"):
        return False

    # Step 3: solvate
    solvent_output = f"{stem}_solvent.gro"
    solvate_log = f"{stem}_solvate.log"
    solvate_cmd = (
        f"gmx solvate -cp {box_output} -cs {args.solvent_file} -o {solvent_output} "
        f"-p {topology_file}"
    )
    if not run_gromacs_command(solvate_cmd, "Add solvent molecules (solvate)", solvate_log, 
                               dry_run=args.dry_run, step_info=f"[3/{total_steps}]"):
        return False

    # Step 4: prepare for ion addition (minimal mdp + grompp)
    create_minimal_mdp("minimal.mdp")

    grompp_log = f"{stem}_grompp.log"
    tpr_file = f"{stem}_ions.tpr"
    grompp_cmd = (
        f"gmx grompp -f minimal.mdp -c {solvent_output} -p {topology_file} -o {tpr_file}"
    )
    if not run_gromacs_command(grompp_cmd, "Prepare for ion addition (grompp)", grompp_log, 
                               dry_run=args.dry_run, step_info=f"[4/{total_steps}]"):
        return False

    # Step 5: genion (determine solvent group dynamically)
    ions_output = f"{stem}_ions.gro"
    genion_log = f"{stem}_genion.log"
    
    # In dry-run mode, use placeholder for solvent group number
    if args.dry_run:
        solvent_group_number = "SOL_GROUP"  # Placeholder for dry-run
        logger.info("[DRY RUN] Using placeholder 'SOL_GROUP' for solvent group number")
        log_command_to_script("", "NOTE: Replace SOL_GROUP below with actual solvent group number from 'gmx make_ndx'")
    else:
        solvent_group_number = get_solvent_group_number(tpr_file)
        if not solvent_group_number:
            logger.error("❌ Could not detect SOL group in index (from make_ndx). Aborting.")
            return False

    # Build genion command with optional neutralization
    neutral_flag = "" if getattr(args, "no_neutral", False) else " -neutral"
    genion_cmd = (
        f"echo {solvent_group_number} | gmx genion -s {tpr_file} -o {ions_output} -p {topology_file} "
        f"-pname {args.cation} -nname {args.anion}{neutral_flag}"
    )
    if not run_gromacs_command(genion_cmd, "Add ions and neutralize (genion)", genion_log, 
                               dry_run=args.dry_run, step_info=f"[5/{total_steps}]"):
        return False

    # Cleanup temporary files (skip in dry-run mode)
    if not args.dry_run:
        for temp_file in ["minimal.mdp", "mdout.mdp", tpr_file]:
            if os.path.exists(temp_file):
                os.remove(temp_file)

    # Summary
    logger.info("")
    logger.info("="*60)
    if args.dry_run:
        logger.info("✅ DRY RUN COMPLETED SUCCESSFULLY!")
        logger.info("="*60)
        logger.info("")
        logger.info(f"📝 Commands saved to: {commands_file}")
        logger.info("")
        logger.info("📋 Pipeline will produce:")
        logger.info(f"  1. Topology:        {topology_file}")
        logger.info(f"  2. Protein + H:     {gro_output}")
        logger.info(f"  3. In box:          {box_output}")
        logger.info(f"  4. Solvated:        {solvent_output}")
        logger.info(f"  5. Final (+ ions):  {ions_output}")
        logger.info("")
        logger.info("🚀 To execute the pipeline:")
        logger.info(f"   bash {commands_file}")
        logger.info("   OR")
        logger.info(f"   python {sys.argv[0]} {args.pdb_file} [options]")
        logger.info("="*60)
    else:
        logger.info("✅ PIPELINE COMPLETED SUCCESSFULLY!")
        logger.info("="*60)
        logger.info("")
        logger.info("📁 Output files:")
        logger.info(f"  Topology:           {topology_file}")
        logger.info(f"  Final structure:    {ions_output} ⭐")
        logger.info(f"  Intermediates:      {box_output}, {solvent_output}")
        logger.info("")
        logger.info("📋 Log files:")
        logger.info(f"  {pdb2gmx_log}, {editconf_log}, {solvate_log}")
        logger.info(f"  {grompp_log}, {genion_log}")
        logger.info(f"")
        logger.info(f"📝 Reproducibility:  {commands_file}")
        
        # Print system summary
        print_system_summary(ions_output, topology_file)
        
        logger.info("")
        logger.info("🎯 Next steps:")
        logger.info("  1. Energy minimization")
        logger.info("  2. NVT equilibration")
        logger.info("  3. NPT equilibration")
        logger.info("  4. Production MD")
        logger.info("="*60)
    return True

def main():
    parser = argparse.ArgumentParser(
        description="Automate GROMACS preprocessing pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with defaults
  %(prog)s protein.pdb
  
  # Custom force field and water model
  %(prog)s protein.pdb --forcefield amber99sb-ildn --water-model tip3p
  
  # Larger box with dodecahedral shape
  %(prog)s protein.pdb --box-distance 1.5 --box-type dodecahedron
  
  # Preview commands without executing (dry-run)
  %(prog)s protein.pdb --dry-run
  
  # Verbose output for debugging
  %(prog)s protein.pdb --verbose
        """)
    
    parser.add_argument("pdb_file", help="Input PDB file")
    parser.add_argument("-ff", "--forcefield", default="charmm27", 
                       help="Force field (default: charmm27)")
    parser.add_argument("-d", "--box-distance", type=float, default=1.0, 
                       help="Distance to box wall in nm (default: 1.0)")
    parser.add_argument("-bt", "--box-type", default="cubic", 
                       help="Box type: cubic, dodecahedron, etc. (default: cubic)")
    parser.add_argument("-cs", "--solvent-file", default="spc216.gro", 
                       help="Solvent model file (default: spc216.gro)")
    parser.add_argument("-wm", "--water-model", default="tip3p", 
                       help="Water model for pdb2gmx (default: tip3p)")
    parser.add_argument("-pname", "--cation", default="NA", 
                       help="Cation type (default: NA)")
    parser.add_argument("-nname", "--anion", default="CL", 
                       help="Anion type (default: CL)")
    parser.add_argument("--no-neutral", action="store_true", 
                       help="Do not add neutralizing ions (omit -neutral in genion)")
    parser.add_argument("--ignore-hydrogens", action="store_true",
                       help="Ignore hydrogen atoms in input PDB (pass -ignh to pdb2gmx)")
    parser.add_argument("--dry-run", action="store_true", 
                       help="Generate commands.sh without executing (preview mode)")
    parser.add_argument("-v", "--verbose", action="store_true",
                       help="Verbose output (show detailed progress)")
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(verbose=args.verbose)
    logger = logging.getLogger(__name__)
    
    # Check PDB file exists
    if not os.path.exists(args.pdb_file):
        logger.error(f"❌ PDB file not found: {args.pdb_file}")
        sys.exit(1)
    
    # Run pipeline
    success = automate_gromacs(args)
    
    if not success:
        logger.error("")
        logger.error("💥 PIPELINE FAILED!")
        logger.error("Check the log files above for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()


#!/usr/bin/env python3
"""
GROMACS Automation Script with command-line arguments and logging
Usage: python gromacs_prepare.py <pdb_file> [options]
"""

import argparse
import subprocess
import sys
import os
from pathlib import Path
from datetime import datetime

def run_gromacs_command(command, description, log_filename):
    """Run a GROMACS command, capture output to log file, and handle errors"""
    print(f"\n{'='*50}")
    print(f"Running: {description}")
    print(f"Command: {command}")
    print(f"Log file: {log_filename}")
    print(f"{'='*50}")
    
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
        
        print("✓ Success!")
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
        
        print(f"✗ Error running {description}:")
        print(f"Check log file: {log_filename}")
        print(f"Return code: {e.returncode}")
        return False

def get_solvent_group_number(tpr_file):
    """Get the group number for solvent from the tpr file using gmx make_ndx"""
    try:
        result = subprocess.run(f"gmx make_ndx -f {tpr_file} -o index.ndx", 
                                shell=True, check=True, 
                                input="q\n", capture_output=True, text=True)
        output = result.stdout
        for line in output.splitlines():
            if "SOL" in line: # or "Water" in line: don't think I need "Water"
                parts = line.split()
                group_number = parts[0]
                return group_number
    except subprocess.CalledProcessError as e:
        print("Error determining solvent group number.")
        return None

def automate_gromacs(args):
    """
    Automate GROMACS preprocessing pipeline
    """
    stem = Path(args.pdb_file).stem
    
    print(f"Starting GROMACS automation for: {args.pdb_file}")
    print(f"File stem: {stem}")
    print(f"Log files will be saved as: {stem}_<tool>.log")
    
    # Step 1: pdb2gmx
    topology_file = f"{stem}.top"
    gro_output = f"{stem}.gro"
    pdb2gmx_log = f"{stem}_pdb2gmx.log"
    
    pdb2gmx_cmd = f"gmx pdb2gmx -f {args.pdb_file} -o {gro_output} -p {topology_file} -ff {args.forcefield} -water {args.water_model}"
    
    if not run_gromacs_command(pdb2gmx_cmd, "pdb2gmx", pdb2gmx_log):
        return False
    
    # Step 2: editconf
    box_output = f"{stem}_box.gro"
    editconf_log = f"{stem}_editconf.log"
    
    editconf_cmd = f"gmx editconf -f {gro_output} -o {box_output} -c -d {args.box_distance} -bt {args.box_type}"
    
    if not run_gromacs_command(editconf_cmd, "editconf", editconf_log):
        return False
    
    # Step 3: solvate
    solvent_output = f"{stem}_solvent.gro"
    solvate_log = f"{stem}_solvate.log"
    
    solvate_cmd = f"gmx solvate -cp {box_output} -cs {args.solvent_file} -o {solvent_output} -p {topology_file}"
    
    if not run_gromacs_command(solvate_cmd, "solvate", solvate_log):
        return False
    
    # Step 4: genion
    ions_output = f"{stem}_ions.gro"
    genion_log = f"{stem}_genion.log"
    
    # Create minimal mdp file
    with open("minimal.mdp", "w") as f:
        f.write(
"""
integrator  = md
dt          = 0.001
nsteps      = 1
""")
    
    # Grompp for genion (with its own log file)
    grompp_log = f"{stem}_grompp.log"
    tpr_file = f"{stem}_ions.tpr"
    grompp_cmd = f"gmx grompp -f minimal.mdp -c {solvent_output} -p {topology_file} -o {tpr_file}"
    
    if not run_gromacs_command(grompp_cmd, "grompp", grompp_log):
        return False

    solvent_group_number = get_solvent_group_number(tpr_file)
    # Genion
    genion_cmd = f"echo {solvent_group_number} | gmx genion -s {tpr_file} -o {ions_output} -p {topology_file} -pname {args.cation} -nname {args.anion} -neutral"

    if not run_gromacs_command(genion_cmd, "genion", genion_log):
        return False
    
    # Cleanup
    temp_files = ["minimal.mdp", "mdout.mdp", tpr_file]
    for temp_file in temp_files:
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    print(f"\n{'='*60}")
    print("✓ GROMACS AUTOMATION COMPLETED SUCCESSFULLY!")
    print(f"{'='*60}")
    print(f"Final output files:")
    print(f"- Topology: {topology_file}")
    print(f"- Final structure: {ions_output}")
    print(f"- Intermediate structures: {stem}_box.gro, {stem}_solvent.gro")
    print(f"\nLog files:")
    print(f"- pdb2gmx: {stem}_pdb2gmx.log")
    print(f"- editconf: {stem}_editconf.log")
    print(f"- solvate: {stem}_solvate.log")
    print(f"- grompp: {stem}_grompp.log")
    print(f"- genion: {stem}_genion.log")
    
    return True

def main():
    parser = argparse.ArgumentParser(description="Automate GROMACS preprocessing pipeline")
    parser.add_argument("pdb_file", help="Input PDB file")
    parser.add_argument("-ff", "--forcefield", default="charmm27", help="Force field (default: charmm27)")
    parser.add_argument("-d", "--box-distance", type=float, default=1.0, help="Distance to box wall in nm (default: 1.0)")
    parser.add_argument("-bt", "--box-type", default="cubic", help="Box type (default: cubic)")
    parser.add_argument("-cs", "--solvent-file", default="spc216.gro", help="Solvent model file (default: spc216.gro)")
    parser.add_argument("-wm", "--water-model", default="tip3p", help="Water model for pdb2gmx (default: tip3p)")
    parser.add_argument("-pname", "--cation", default="NA", help="Cation type (default: NA)")
    parser.add_argument("-nname", "--anion", default="CL", help="Anion type (default: CL)")
    
    args = parser.parse_args()
    
    if not os.path.exists(args.pdb_file):
        print(f"Error: PDB file '{args.pdb_file}' not found!")
        sys.exit(1)
    
    success = automate_gromacs(args)
    
    if not success:
        print("\n💥 Pipeline failed! Check the log files above for details.")
        sys.exit(1)

if __name__ == "__main__":
    main()


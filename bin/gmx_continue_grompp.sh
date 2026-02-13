#!/bin/bash -e
# Script for continuing a GROMACS molecular dynamics simulation
# It uses grompp to prepare the input files for continuation
# Required files: topology (.top), structure (.gro), checkpoint (.cpt), and energy (.edr)
#
# Usage: ./md_continue_grompp.sh <basename> <replica> <start_time> <end_time> <dt> [template_mdp] [initial_basename] [timestep] [title_suffix]
#   basename:         Base name for output files (e.g., "md")
#   replica:          Replica number (e.g., 0, 1, 2)
#   start_time:       Start time in ns (e.g., 0)
#   end_time:         End time in ns (e.g., 500)
#   dt:               Time step per segment in ns (e.g., 100)
#   template_mdp:     (Optional) MDP template for production runs, default: "md.mdp"
#   initial_basename: (Optional) Basename for initial equilibration files, default: "npt"
#   timestep:         (Optional) Integration timestep in ps, default: 0.002 (2 fs)
#   title_suffix:     (Optional) System-specific title to append (e.g., "TRP_cage_replica_1")
#
# Example: ./md_continue_grompp.sh md 0 0 500 100 md_production.mdp npt 0.002 "TRP_cage"
#   Uses npt.gro/cpt/edr as starting point, md_production.mdp as template
#   Appends "TRP_cage" to MDP title
#   Then runs: md_0_0_100, md_0_100_200, md_0_200_300, etc.
#
# CRASH RECOVERY:
#   If a simulation crashes or is interrupted, simply re-run the script with 
#   the EXACT SAME ARGUMENTS. The script will:
#   - Automatically detect which segments completed successfully
#   - Skip already-completed segments
#   - Resume interrupted segments from their checkpoint files
#   
#   Example after crash:
#     ./md_continue_grompp.sh md 0 0 500 100 md_production.mdp npt
#     (Same command - the script figures out where to continue)

# Check arguments
if [[ $# -lt 5 ]] || [[ $# -gt 9 ]]; then
    echo "Error: Incorrect number of arguments"
    echo "Usage: $0 <basename> <replica> <start_time> <end_time> <dt> [template_mdp] [initial_basename] [timestep] [title_suffix]"
    echo "  basename:         Base name for output files (e.g., 'md')"
    echo "  replica:          Replica number (e.g., 0)"
    echo "  start_time:       Start time in ns (e.g., 0)"
    echo "  end_time:         End time in ns (e.g., 500)"
    echo "  dt:               Time step per segment in ns (e.g., 100)"
    echo "  template_mdp:     (Optional) MDP template for production runs (default: 'md.mdp')"
    echo "  initial_basename: (Optional) Basename for initial equilibration files (default: 'npt')"
    echo "  timestep:         (Optional) Integration timestep in ps (default: 0.002)"
    echo "  title_suffix:     (Optional) System-specific title to append (e.g., 'TRP_cage_replica_1')"
    echo "  --plumed <file>:  (Optional) PLUMED input file for enhanced sampling/analysis"
    echo ""
    echo "Example: $0 md 0 0 500 100 md_production.mdp npt 0.002 'TRP_cage'"
    echo "  Uses npt.gro/cpt/edr as starting point, md_production.mdp as template"
    echo "  Appends 'TRP_cage' to MDP title, then runs md_0_0_100, md_0_100_200, etc."
    exit 1
fi

# Parse arguments
basename_arg="$1"
replica="$2"
tstart="$3"
tend="$4"
dt="$5"
template_mdp="${6:-md.mdp}"      # Default to "md.mdp" if not provided
initial_basename="${7:-npt}"     # Default to "npt" if not provided
timestep_ps="${8:-0.002}"        # Default to 0.002 ps (2 fs) if not provided
title_suffix="${9:-}"            # Optional title suffix
topology="topol.top"
plumed_file=""
base_name="${basename_arg}_${replica}"

# Parse optional flags (--plumed)
shift 9
while [[ $# -gt 0 ]]; do
    case "$1" in
        --plumed)
            plumed_file="$2"
            shift 2
            ;;
        *)
            echo "Warning: Unknown option '$1' ignored"
            shift
            ;;
    esac
done

# Setup PLUMED flag if provided
if [[ -n "${plumed_file}" ]]; then
    if [[ ! -f "${plumed_file}" ]]; then
        echo "Error: PLUMED file '${plumed_file}' not found"
        exit 1
    fi
    plumed_flag="-plumed ${plumed_file}"
    echo "Using PLUMED file: ${plumed_file}"
else
    plumed_flag=""
fi

# Calculate nsteps per segment: (dt_ns * 1000 ps/ns) / timestep_ps
nsteps_per_segment=$(awk -v dt="$dt" -v ts="$timestep_ps" 'BEGIN {printf "%.0f\n", (dt * 1000) / ts}')

echo "=== MD Continuation Script ==="
echo "Base name:       ${basename_arg}"
echo "Replica:         ${replica}"
echo "Time range:      ${tstart} -> ${tend} ns"
echo "Segment dt:      ${dt} ns (${nsteps_per_segment} steps)"
echo "Timestep:        ${timestep_ps} ps"
echo "Template MDP:    ${template_mdp}"
echo "Initial files:   ${initial_basename}.gro/cpt/edr"
echo "=============================="

# Source GROMACS if not already available
if ! command -v gmx &> /dev/null; then
    if [[ -f /programs/gromacs-2025.2/bin/GMXRC.bash ]]; then
        . /programs/gromacs-2025.2/bin/GMXRC.bash
    else
        echo "Error: GROMACS not found. Please source GMXRC.bash manually or add gmx to PATH."
        exit 1
    fi
fi

# Check for topology file
if [[ ! -f ${topology} ]]; then
    echo "Error: Topology file '${topology}' not found"
    exit 1
fi

# Check for template MDP file
if [[ ! -f ${template_mdp} ]]; then
    echo "Error: Template MDP file '${template_mdp}' not found"
    echo "This should be your production run MDP file (NOT the NPT/NVT equilibration MDP)."
    exit 1
fi

# Function to check if a segment completed successfully
# Returns 0 (success) if segment is complete, 1 (failure) otherwise
check_segment_complete() {
    local segment_name="$1"
    
    # Check if essential output files exist (.cpt is optional)
    if [[ ! -f ${segment_name}.gro ]] || [[ ! -f ${segment_name}.edr ]] || \
       [[ ! -f ${segment_name}.log ]]; then
        return 1
    fi
    
    # Check if log file indicates successful completion
    if grep -q "Finished mdrun" ${segment_name}.log 2>/dev/null; then
        return 0
    else
        return 1
    fi
}

# Determine where to start based on existing files (crash recovery)
actual_start=${tstart}

# Check if we're restarting from a crash
for ((time=${tstart} ; time<${tend} ; time+=${dt})) ; do
    check_prev="${base_name}_$((time))_$((time+dt))"
    
    # If TPR exists but segment didn't complete, this segment was interrupted
    if [[ -f ${check_prev}.tpr ]]; then
        if ! check_segment_complete "${check_prev}"; then
            echo "Found incomplete segment: ${check_prev}"
            echo "Will resume from this segment..."
            actual_start=$time
            break
        fi
    fi
    
    # If segment completed successfully, move to next
    if check_segment_complete "${check_prev}"; then
        actual_start=$((time + dt))
    else
        # Missing files or incomplete, start from here
        break
    fi
done

if [[ ${actual_start} -ge ${tend} ]]; then
    echo "All segments already completed (0 -> ${tend} ns). Nothing to do."
    exit 0
fi

if [[ ${actual_start} -gt ${tstart} ]]; then
    echo "Resuming from ${actual_start} ns (previous segments already completed)"
fi

# Initial segment (time 0 to dt) - only if starting from beginning
if [[ ${actual_start} -eq ${tstart} ]] && [[ ${tstart} -eq 0 ]]; then
    initial_cur="${base_name}_0_${dt}"
    echo "Setting up initial segment: 0 -> ${dt} ns (using ${initial_basename} as starting point)"

    # Check for initial equilibration files
    if [[ ! -f ${initial_basename}.gro ]]; then
        echo "Error: Missing initial structure file (${initial_basename}.gro)"
        echo "This file should come from your NPT/NVT equilibration step."
        exit 1
    fi
    
    if [[ ! -f ${initial_basename}.edr ]]; then
        echo "Error: Missing initial energy file (${initial_basename}.edr)"
        echo "This file should come from your NPT/NVT equilibration step."
        exit 1
    fi

    # Check if TPR exists but segment didn't complete (crashed during this segment)
    if [[ -f ${initial_cur}.tpr ]]; then
        if ! check_segment_complete "${initial_cur}"; then
            echo "Resuming interrupted initial segment..."
            echo "Note: Existing incomplete output files will be backed up with .bak extension"
            gmx mdrun -deffnm ${initial_cur} -cpi ${initial_cur}.cpt ${plumed_flag}
        else
            echo "Initial segment already completed successfully."
        fi
    else
        echo "Running grompp for initial segment (using ${template_mdp})..."
        # Create MDP with tinit=0 and nsteps for first segment
        sed -e "s/\(tinit\s*=\s*\)[0-9]\+/\10/" \
            -e "s/\(nsteps\s*=\s*\)[0-9]\+/\1${nsteps_per_segment}/" \
            ${template_mdp} > ${initial_cur}.mdp
        
        # Append title suffix if provided
        if [[ -n "${title_suffix}" ]]; then
            sed -i "s/\(title\s*=\s*.*\)/\1,${title_suffix}/" ${initial_cur}.mdp
        fi
        
        # Use checkpoint if available, otherwise continue without it
        if [[ -f ${initial_basename}.cpt ]]; then
            echo "Using checkpoint from ${initial_basename}"
            gmx grompp -f ${initial_cur}.mdp -c ${initial_basename}.gro -t ${initial_basename}.cpt \
                       -e ${initial_basename}.edr -p ${topology} -o ${initial_cur}.tpr
        else
            echo "Warning: No checkpoint file found (${initial_basename}.cpt)"
            echo "Continuing without checkpoint. Velocities will be regenerated."
            gmx grompp -f ${initial_cur}.mdp -c ${initial_basename}.gro \
                       -e ${initial_basename}.edr -p ${topology} -o ${initial_cur}.tpr
        fi

        echo "Running mdrun for initial segment..."
        gmx mdrun -deffnm ${initial_cur} ${plumed_flag}
    fi

    if [[ -f ./STOP ]]; then
        echo "Found STOP file. Stopping."
        exit 0
    fi
fi

# Continuation segments
start_time=$(( actual_start > dt ? actual_start : dt ))
for ((time=${start_time} ; time<${tend} ; time+=${dt})) ; do
    prev="${base_name}_$((time-dt))_${time}"
    cur="${base_name}_${time}_$((time+dt))"
    timeps=$((time * 1000))  # Convert ns to ps for tinit parameter
    
    echo "Setting up segment: ${time} -> $((time+dt)) ns"

    # Check if this segment was interrupted (TPR exists but segment incomplete)
    if [[ -f ${cur}.tpr ]]; then
        if ! check_segment_complete "${cur}"; then
            echo "Resuming interrupted segment..."
            echo "Note: Existing incomplete output files will be backed up with .bak extension"
            gmx mdrun -deffnm ${cur} -cpi ${cur}.cpt ${plumed_flag}
        else
            echo "Segment already completed successfully. Skipping..."
            continue
        fi
    else
        # Normal workflow: prepare and run
        if [[ ! -f ${prev}.gro ]]; then
            echo "Error: Missing structure file from previous segment (${prev}.gro). Aborting."
            exit 1
        fi
        
        if [[ ! -f ${prev}.edr ]]; then
            echo "Error: Missing energy file from previous segment (${prev}.edr). Aborting."
            exit 1
        fi

        # Update tinit and nsteps in template MDP file for current segment
        sed -e "s/\(tinit\s*=\s*\)[0-9]\+/\1${timeps}/" \
            -e "s/\(nsteps\s*=\s*\)[0-9]\+/\1${nsteps_per_segment}/" \
            ${template_mdp} > ${cur}.mdp
        
        # Append title suffix if provided
        if [[ -n "${title_suffix}" ]]; then
            sed -i "s/\(title\s*=\s*.*\)/\1,${title_suffix}/" ${cur}.mdp
        fi

        echo "Running grompp (using ${template_mdp} with tinit=${timeps} ps, nsteps=${nsteps_per_segment})..."
        
        # Check if checkpoint exists (optional for continuation)
        if [[ -f ${prev}.cpt ]]; then
            echo "Using checkpoint from previous segment"
            gmx grompp -f ${cur}.mdp -c ${prev}.gro -t ${prev}.cpt \
                       -e ${prev}.edr -p ${topology} -o ${cur}.tpr
        else
            echo "Warning: No checkpoint file found (${prev}.cpt), continuing without it"
            echo "Velocities will be regenerated. This is fine but less seamless."
            gmx grompp -f ${cur}.mdp -c ${prev}.gro \
                       -e ${prev}.edr -p ${topology} -o ${cur}.tpr
        fi
        
        echo "Running mdrun..."
        gmx mdrun -deffnm ${cur} ${plumed_flag}
    fi
    
    if [[ -f ./STOP ]]; then
        echo "Found STOP file. Stopping."
        exit 0
    fi
done

echo "All segments completed successfully (0 -> ${tend} ns)"

#!/bin/bash -e
# Script for continuing a GROMACS molecular dynamics simulation
# It uses grompp to prepare the input files for continuation
# Required files: topology (.top), structure (.gro), checkpoint (.cpt), and energy (.edr)
#
# Usage: ./md_continue_grompp.sh <basename> <replica> <start_time> <end_time> <dt>
#   basename:   Base name for output files (e.g., "md")
#   replica:    Replica number (e.g., 0, 1, 2)
#   start_time: Start time in ns (e.g., 0)
#   end_time:   End time in ns (e.g., 500)
#   dt:         Time step per segment in ns (e.g., 100)
#
# Example: ./md_continue_grompp.sh md 0 0 500 100
#   This will run: md_0_0_100, md_0_100_200, md_0_200_300, etc.

# Check arguments
if [[ $# -ne 5 ]]; then
    echo "Error: Incorrect number of arguments"
    echo "Usage: $0 <basename> <replica> <start_time> <end_time> <dt>"
    echo "  basename:   Base name for output files (e.g., 'md')"
    echo "  replica:    Replica number (e.g., 0)"
    echo "  start_time: Start time in ns (e.g., 0)"
    echo "  end_time:   End time in ns (e.g., 500)"
    echo "  dt:         Time step per segment in ns (e.g., 100)"
    echo ""
    echo "Example: $0 md 0 0 500 100"
    exit 1
fi

# Parse arguments
basename_arg="$1"
replica="$2"
tstart="$3"
tend="$4"
dt="$5"
topology="topol.top"
base_name="${basename_arg}_${replica}"

echo "=== MD Continuation Script ==="
echo "Base name:  ${basename_arg}"
echo "Replica:    ${replica}"
echo "Time range: ${tstart} -> ${tend} ns"
echo "Segment dt: ${dt} ns"
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

# Determine where to start based on existing files (crash recovery)
actual_start=${tstart}

# Check if we're restarting from a crash
for ((time=${tstart} ; time<${tend} ; time+=${dt})) ; do
    check_prev="${base_name}_$((time))_$((time+dt))"
    
    # If TPR exists but output files don't, this segment was interrupted
    if [[ -f ${check_prev}.tpr ]] && [[ ! -f ${check_prev}.gro ]]; then
        echo "Found incomplete segment: ${check_prev}"
        echo "Restarting from this segment..."
        actual_start=$time
        break
    fi
    
    # If all output files exist, this segment completed
    if [[ -f ${check_prev}.gro ]] && [[ -f ${check_prev}.cpt ]] && [[ -f ${check_prev}.edr ]]; then
        actual_start=$((time + dt))
    else
        # Missing files, start from here
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
    initial_prev="${base_name}_0_0"
    initial_cur="${base_name}_0_${dt}"
    echo "Setting up initial segment: 0 -> ${dt} ns"

    if [[ ! -f ${initial_prev}.gro ]] || [[ ! -f ${initial_prev}.cpt ]] || [[ ! -f ${initial_prev}.edr ]]; then
        echo "Error: Missing initial files (${initial_prev}.gro, .cpt, or .edr)"
        exit 1
    fi

    if [[ ! -f ${initial_cur}.mdp ]]; then
        echo "Error: ${initial_cur}.mdp not found"
        exit 1
    fi

    # Check if TPR exists but outputs don't (crashed during this segment)
    if [[ -f ${initial_cur}.tpr ]] && [[ ! -f ${initial_cur}.gro ]]; then
        echo "Resuming interrupted initial segment..."
        gmx mdrun -deffnm ${initial_cur} -cpi ${initial_cur}.cpt
    else
        echo "Running grompp for initial segment..."
        gmx grompp -f ${initial_cur}.mdp -c ${initial_prev}.gro -t ${initial_prev}.cpt \
                   -e ${initial_prev}.edr -p ${topology} -o ${initial_cur}.tpr

        echo "Running mdrun for initial segment..."
        gmx mdrun -deffnm ${initial_cur}
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

    # Check if this segment was interrupted (TPR exists but outputs don't)
    if [[ -f ${cur}.tpr ]] && [[ ! -f ${cur}.gro ]]; then
        echo "Resuming interrupted segment..."
        gmx mdrun -deffnm ${cur} -cpi ${cur}.cpt
    else
        # Normal workflow: prepare and run
        if [[ ! -f ${prev}.mdp ]]; then
            echo "Error: ${prev}.mdp not found. Aborting."
            exit 1
        fi
        
        if [[ ! -f ${prev}.gro ]] || [[ ! -f ${prev}.cpt ]] || [[ ! -f ${prev}.edr ]]; then
            echo "Error: Missing output files from previous segment (${prev}.gro, .cpt, or .edr). Aborting."
            exit 1
        fi

        # Update tinit in MDP file for current segment
        sed -e "s/\(tinit\s*=\s*\)[0-9]\+/\1${timeps}/" ${prev}.mdp > ${cur}.mdp

        echo "Running grompp..."
        gmx grompp -f ${cur}.mdp -c ${prev}.gro -t ${prev}.cpt \
                   -e ${prev}.edr -p ${topology} -o ${cur}.tpr
        
        echo "Running mdrun..."
        gmx mdrun -deffnm ${cur}
    fi
    
    if [[ -f ./STOP ]]; then
        echo "Found STOP file. Stopping."
        exit 0
    fi
done

echo "All segments completed successfully (0 -> ${tend} ns)"

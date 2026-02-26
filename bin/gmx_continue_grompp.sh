#!/bin/bash -e
# Script for continuing/extending GROMACS molecular dynamics simulations
# It uses grompp to prepare the input files for continuation
# Required files: topology (.top), and either:
#   - Initial run (start_time=0): initial equilibration files (.gro/.cpt/.edr/.mdp)
#   - Continuation (start_time>0): previous segment files - MDP copied from previous segment
#
# Usage: ./gmx_continue_grompp.sh <basename> <replica> <start_time> <end_time> <dt> [OPTIONS]
#   basename:         Base name for output files (e.g., "md" or "trp_fold")
#   replica:          Replica number (e.g., 0, 1, 2)
#   start_time:       Start time in ns (e.g., 0 for new run, 10000 for continuation)
#   end_time:         End time in ns (e.g., 500 or 20000)
#   dt:               Time step per segment in ns (e.g., 100 or 500)
#
# CONTINUATION MODE (start_time > 0):
#   The script automatically copies the MDP file from the previous segment.
#   No template MDP needed! Simply specify the time range you want to extend.
#
#   Example - extend from 10000 to 20000 ns:
#     ./gmx_continue_grompp.sh trp_fold 0 10000 20000 500
#     Copies trp_fold_0_9500_10000.mdp (previous segment)
#     Continues: trp_fold_0_10000_10500, trp_fold_0_10500_11000, etc.
#
# INITIAL RUN (start_time = 0):
#   Starting a new production run from equilibration files.
#   MDP copied from equilibration MDP, or uses template if provided.
#
#   Example - start new 500 ns run:
#     ./gmx_continue_grompp.sh md 0 0 500 100 --initial npt_final
#     Copies npt_final.mdp and uses npt_final.gro/cpt/edr as starting point
#
# CRASH RECOVERY:
#   If a simulation crashes, simply re-run with the EXACT SAME ARGUMENTS.
#   The script will automatically detect and resume from the interruption point.

# Check minimum arguments
if [[ $# -lt 5 ]]; then
    echo "Error: Insufficient arguments"
    echo "Usage: $0 <basename> <replica> <start_time> <end_time> <dt> [OPTIONS]"
    echo "  basename:         Base name for output files (e.g., 'md' or 'trp_fold')"
    echo "  replica:          Replica number (e.g., 0)"
    echo "  start_time:       Start time in ns (0 for new run, >0 for continuation)"
    echo "  end_time:         End time in ns (e.g., 500 or 20000)"
    echo "  dt:               Time step per segment in ns (e.g., 100 or 500)"
    echo ""
    echo "Optional arguments:"
    echo "  --template <file>      MDP template (only needed if no previous .mdp to copy from)"
    echo "  --initial <basename>   Basename for initial equilibration files (default: 'npt')"
    echo "  --timestep <ps>        Integration timestep in ps (default: 0.002)"
    echo "  --title <suffix>       System-specific title to append (e.g., 'TRP_cage_replica_1')"
    echo "  --topology <file>      Topology file (default: 'topol.top')"
    echo "  --plumed <file>        PLUMED input file for enhanced sampling/analysis"
    echo "  --ps                   Interpret times as picoseconds (default: nanoseconds)"
    echo ""
    echo "Examples:"
    echo "  # Extend existing run from 10000 to 20000 ns (MDP auto-copied)"
    echo "  $0 trp_fold 0 10000 20000 500"
    echo ""
    echo "  # Start new run from 0 to 500 ns (copies npt_final.mdp)"
    echo "  $0 md 0 0 500 100 --initial npt_final"
    exit 1
fi

# Validate that first 5 arguments are not flags (before shifting)
for i in 1 2 3 4 5; do
    eval "arg=\${$i}"
    if [[ "$arg" == --* ]]; then
        echo "Error: Argument $i appears to be a flag ('$arg')"
        echo "The first 5 arguments must be: <basename> <replica> <start_time> <end_time> <dt>"
        echo ""
        echo "Correct format:"
        echo "  $0 md 0 0 500 100 [OPTIONS]"
        echo ""
        echo "Your command had a flag in position $i. Did you:"
        echo "  - Forget a required positional argument?"
        echo "  - Misspell a flag name?"
        exit 1
    fi
done

# Parse required positional arguments
basename_arg="$1"
replica="$2"
tstart="$3"
tend="$4"
dt="$5"
shift 5

# Set default values for optional arguments
template_mdp="md.mdp"
template_provided=false  # Track if --template was explicitly provided
initial_basename="npt"
timestep_ps="0.002"
title_suffix=""
plumed_file=""
topology="topol.top"
time_unit="ns"  # Default to nanoseconds

# Parse optional arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --template)
            if [[ -z "$2" ]] || [[ "$2" == --* ]]; then
                echo "Error: --template requires a value"
                exit 1
            fi
            template_mdp="$2"
            template_provided=true
            shift 2
            ;;
        --initial)
            if [[ -z "$2" ]] || [[ "$2" == --* ]]; then
                echo "Error: --initial requires a value"
                exit 1
            fi
            initial_basename="$2"
            shift 2
            ;;
        --timestep)
            if [[ -z "$2" ]] || [[ "$2" == --* ]]; then
                echo "Error: --timestep requires a value"
                exit 1
            fi
            timestep_ps="$2"
            shift 2
            ;;
        --title)
            if [[ -z "$2" ]] || [[ "$2" == --* ]]; then
                echo "Error: --title requires a value"
                exit 1
            fi
            title_suffix="$2"
            shift 2
            ;;
        --topology)
            if [[ -z "$2" ]] || [[ "$2" == --* ]]; then
                echo "Error: --topology requires a value"
                exit 1
            fi
            topology="$2"
            shift 2
            ;;
        --plumed)
            if [[ -z "$2" ]] || [[ "$2" == --* ]]; then
                echo "Error: --plumed requires a value"
                exit 1
            fi
            plumed_file="$2"
            shift 2
            ;;
        --ps)
            time_unit="ps"
            shift
            ;;
        --*)
            echo "Error: Unknown option '$1'"
            echo ""
            echo "Valid options: --template, --initial, --timestep, --title, --topology, --plumed, --ps"
            echo ""
            echo "Did you misspell an option? Common typos:"
            echo "  --intitial  → should be --initial"
            echo "  --templte   → should be --template"
            exit 1
            ;;
        *)
            echo "Error: Unexpected argument '$1'"
            echo "All optional arguments must start with --"
            exit 1
            ;;
    esac
done

base_name="${basename_arg}_${replica}"

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

# Calculate nsteps per segment
if [[ "${time_unit}" == "ps" ]]; then
    # Times already in ps
    nsteps_per_segment=$(awk -v dt="$dt" -v ts="$timestep_ps" 'BEGIN {printf "%.0f\n", dt / ts}')
else
    # Times in ns, convert to ps
    nsteps_per_segment=$(awk -v dt="$dt" -v ts="$timestep_ps" 'BEGIN {printf "%.0f\n", (dt * 1000) / ts}')
fi

echo "=== MD Continuation Script ==="
echo "Base name:       ${basename_arg}"
echo "Replica:         ${replica}"
echo "Time range:      ${tstart} -> ${tend} ${time_unit}"
echo "Segment dt:      ${dt} ${time_unit} (${nsteps_per_segment} steps)"
echo "Timestep:        ${timestep_ps} ps"
if [[ ${tstart} -eq 0 ]]; then
    echo "Initial files:   ${initial_basename}.gro/cpt/edr/mdp"
else
    echo "Continuing from: ${tstart} ${time_unit} (MDP copied from previous segment)"
fi
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
    echo ""
    echo "The topology file is required for gmx grompp to create new TPR files."
    echo "Please either:"
    echo "  1. Create/copy a topol.top file in this directory, or"
    echo "  2. Use --topology <file> to specify a different topology file"
    echo ""
    echo "Example: $0 md 0 0 5000 1000 --topology /path/to/topol.top --initial npt_free"
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

# Check if we actually need to start from segment 0 (template required)
if [[ ${actual_start} -eq 0 ]] && [[ ! -f ${template_mdp} ]]; then
    echo "Error: Template MDP file '${template_mdp}' not found"
    echo ""
    echo "Starting from time 0 requires a template MDP file."
    echo "Please either:"
    echo "  1. Create/copy '${template_mdp}' in this directory, or"
    echo "  2. Use --template <file> to specify a different MDP file"
    echo ""
    echo "Example: $0 ${basename_arg} ${replica} 0 ${tend} ${dt} --template md_production.mdp"
    exit 1
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
        # Priority 1: If --template explicitly provided, use it
        if [[ ${template_provided} == true ]]; then
            echo "Using template MDP: ${template_mdp}"
            sed -e "s/\(tinit\s*=\s*\)[0-9]\+/\10/" \
                -e "s/\(nsteps\s*=\s*\)[0-9]\+/\1${nsteps_per_segment}/" \
                ${template_mdp} > ${initial_cur}.mdp
        # Priority 2: Try to copy MDP from initial equilibration
        elif [[ -f ${initial_basename}.mdp ]]; then
            echo "Copying MDP from ${initial_basename}.mdp"
            cp "${initial_basename}.mdp" "${initial_cur}.mdp"
            # Update tinit=0 and nsteps for first segment
            sed -i -e "s/\(tinit\s*=\s*\)[0-9]\+/\10/" \
                   -e "s/\(nsteps\s*=\s*\)[0-9]\+/\1${nsteps_per_segment}/" \
                   ${initial_cur}.mdp
        # Priority 3: Use default template (already checked to exist)
        else
            echo "Using template MDP: ${template_mdp}"
            sed -e "s/\(tinit\s*=\s*\)[0-9]\+/\10/" \
                -e "s/\(nsteps\s*=\s*\)[0-9]\+/\1${nsteps_per_segment}/" \
                ${template_mdp} > ${initial_cur}.mdp
        fi
        
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
    # Convert to ps for tinit parameter if needed
    if [[ "${time_unit}" == "ps" ]]; then
        timeps=${time}  # Already in ps
    else
        timeps=$((time * 1000))  # Convert ns to ps
    fi
    
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

        # Priority 1: If --template explicitly provided, always use it
        if [[ ${template_provided} == true ]]; then
            echo "Using template MDP: ${template_mdp}"
            sed -e "s/\(tinit\s*=\s*\)[0-9]\+/\1${timeps}/" \
                -e "s/\(nsteps\s*=\s*\)[0-9]\+/\1${nsteps_per_segment}/" \
                ${template_mdp} > ${cur}.mdp
        # Priority 2: Copy MDP from previous segment
        elif [[ -f ${prev}.mdp ]]; then
            echo "Copying MDP from previous segment: ${prev}.mdp"
            cp "${prev}.mdp" "${cur}.mdp"
            # Update tinit and nsteps for current segment
            sed -i -e "s/\(tinit\s*=\s*\)[0-9]\+/\1${timeps}/" \
                   -e "s/\(nsteps\s*=\s*\)[0-9]\+/\1${nsteps_per_segment}/" \
                   ${cur}.mdp
        # Priority 3: No previous MDP and no explicit template
        else
            echo "Error: Previous MDP file '${prev}.mdp' not found and no --template provided"
            echo "Either the previous segment is missing, or you need to specify --template"
            echo "Example: $0 ${basename_arg} ${replica} ${tstart} ${tend} ${dt} --template md_production.mdp"
            exit 1
        fi
        
        # Append title suffix if provided
        if [[ -n "${title_suffix}" ]]; then
            sed -i "s/\(title\s*=\s*.*\)/\1,${title_suffix}/" ${cur}.mdp
        fi

        echo "Running grompp (tinit=${timeps} ps, nsteps=${nsteps_per_segment})..."
        
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

echo "All segments completed successfully (0 -> ${tend} ${time_unit})"

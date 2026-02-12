#!/bin/bash
# Script for continuing GROMACS MD simulations using TPR extension method
# Unlike gmx_continue_grompp.sh which creates new TPR for each segment,
# this script extends the existing TPR and uses mdrun continuation
#
# Usage: 
#   ./gmx_continue_extend.sh <extend_time>  # EASY MODE: auto-detect files
#   ./gmx_continue_extend.sh <basename> <replica> <start_time> <end_time> <dt> [OPTIONS]
#
# EASY MODE (single argument):
#   extend_time:      Time in ps to extend (auto-detects TPR/CPT in current directory)
#
# NORMAL MODE (5+ arguments):
#   basename:         Base name for output files (e.g., "md")
#   replica:          Replica number (e.g., 0, 1, 2)
#   start_time:       Start time in ps (e.g., 0)
#   end_time:         End time in ps (e.g., 500)
#   dt:               Time step per segment in ps (e.g., 100)
#
# Optional arguments (order-independent after required args):
#   --template <file>      MDP template for initial TPR creation (default: md.mdp)
#   --initial <basename>   Basename for initial equilibration files (default: npt)
#   --timestep <ps>        Integration timestep in ps (default: 0.002)
#   --title <suffix>       System-specific title to append
#   --append               Use append mode (single output file, default: noappend)
#   --noappend             Use noappend mode (creates part000X files)
#   --tpr <file>           Start from existing TPR instead of creating new one
#   --cpt <file>           Checkpoint file to use with --tpr (optional)
#   --topology <file>      Topology file (default: topol.top)
#
# MODES:
#   --append:   Extends same TPR, appends to same .xtc/.edr/.log files
#               Output: md_0.xtc, md_0.edr, md_0.log (continuously extended)
#   
#   --noappend: (DEFAULT) Extends TPR, creates new files per segment with part000X suffix
#               Output: md_0.xtc, md_0.part0001.xtc, md_0.part0002.xtc, etc.
#
# WORKFLOW:
#   1st segment: Creates TPR from template MDP (or uses --tpr if provided)
#   2nd+ segments: Extends TPR with gmx convert-tpr, continues with mdrun -cpi
#
# Examples:
#   # Easy mode: extend by 5 ns (5000 ps)
#   ./gmx_continue_extend.sh 5000
#
#   # Default (noappend mode, creates part files)
#   ./gmx_continue_extend.sh md 0 0 500 100 --template production.mdp
#
#   # Append mode (single output file)
#   ./gmx_continue_extend.sh md 0 0 500 100 --template production.mdp --append
#
#   # Start from existing TPR
#   ./gmx_continue_extend.sh md 0 0 500 100 --tpr md_0.tpr --append
#
#   # Full options
#   ./gmx_continue_extend.sh md 0 0 500 100 --template production.mdp \
#       --initial npt --timestep 0.002 --title "MyProtein" --noappend

# =============================================================================
# EASY MODE: Single argument = extend time (auto-detect files)
# =============================================================================
if [[ $# -eq 1 ]] && [[ "$1" =~ ^[0-9]+$ ]]; then
    extend_time="$1"
    
    echo "=== Easy Mode: Extending simulation by ${extend_time} ps ==="
    echo ""
    
    # Find TPR file(s) in current directory
    tpr_files=(*.tpr)
    
    if [[ ${#tpr_files[@]} -eq 0 ]] || [[ ! -f "${tpr_files[0]}" ]]; then
        echo "Error: No TPR file found in current directory"
        echo "Please run this command in the directory containing your simulation files"
        exit 1
    fi
    
    if [[ ${#tpr_files[@]} -gt 1 ]]; then
        echo "Error: Multiple TPR files found in current directory:"
        for tpr in "${tpr_files[@]}"; do
            echo "  - $tpr"
        done
        echo ""
        echo "Please specify which one to use with full syntax:"
        echo "  $0 <basename> <replica> <start_time> <end_time> <dt> --tpr <file>"
        exit 1
    fi
    
    # Use the found TPR file
    found_tpr="${tpr_files[0]}"
    found_base="${found_tpr%.tpr}"
    found_cpt="${found_base}.cpt"
    
    echo "Auto-detected TPR file: ${found_tpr}"
    
    # Check for checkpoint
    if [[ -f "${found_cpt}" ]]; then
        echo "Auto-detected checkpoint: ${found_cpt}"
        existing_cpt="${found_cpt}"
    else
        echo "Note: No checkpoint found (${found_cpt})"
        echo "Will start from TPR end time"
        existing_cpt=""
    fi
    
    # Get current end time from TPR
    tpr_info=$(gmx check -s "${found_tpr}" 2>&1)
    current_end=$(echo "$tpr_info" | grep "Last frame" | awk '{print $NF}')
    
    if [[ -z "${current_end}" ]]; then
        echo "Error: Could not determine current end time from TPR"
        exit 1
    fi
    
    # Calculate new end time
    new_end=$(awk "BEGIN {printf \"%.0f\", ${current_end} + ${extend_time}}")
    
    echo ""
    echo "Current end time: ${current_end} ps"
    echo "Extension: ${extend_time} ps"
    echo "New end time: ${new_end} ps"
    echo "Mode: append (default for easy mode)"
    echo ""
    
    # Set up arguments for main script logic
    basename_arg="${found_base}"
    replica=0
    tstart=0
    tend="${new_end}"
    dt="${extend_time}"
    existing_tpr="${found_tpr}"
    append_mode="yes"
    
    # Set default values
    template_mdp="md.mdp"
    initial_basename="npt"
    timestep_ps="0.002"
    title_suffix=""
    topology="topol.top"

# =============================================================================
# NORMAL MODE: 5+ arguments
# =============================================================================
elif [[ $# -lt 5 ]]; then
    echo "Error: Invalid arguments"
    echo ""
    echo "Usage: $0 <extend_time>  # Easy mode"
    echo "   or: $0 <basename> <replica> <start_time> <end_time> <dt> [OPTIONS]"
    echo ""
    echo "Easy mode:"
    echo "  extend_time       Time in ps to extend (auto-detects files)"
    echo ""
    echo "Normal mode:"
    echo "  basename          Base name for output files"
    echo "  replica           Replica number"
    echo "  start_time        Start time in ps"
    echo "  end_time          End time in ps"
    echo "  dt                Segment length in ps"
    echo ""
    echo "Optional:"
    echo "  --template <file>     MDP template (default: md.mdp)"
    echo "  --initial <name>      Initial files basename (default: npt)"
    echo "  --timestep <ps>       Timestep in ps (default: 0.002)"
    echo "  --title <suffix>      Title suffix for MDP"
    echo "  --append              Single output file (extends continuously)"
    echo "  --noappend            Multiple part files (default)"
    echo "  --tpr <file>          Start from existing TPR"
    echo "  --cpt <file>          Checkpoint file for --tpr (optional)"
    echo "  --topology <file>     Topology file (default: topol.top)"
    echo ""
    echo "Examples:"
    echo "  $0 5000                                          # Easy: extend by 5 ns"
    echo "  $0 md 0 0 500 100 --template production.mdp     # Normal mode"
    echo "  $0 md 0 0 500 100 --tpr md_0.tpr --append       # Continue existing"
    exit 1
else
    echo "Error: Missing required arguments"
    echo "Usage: $0 <basename> <replica> <start_time> <end_time> <dt> [OPTIONS]"
    echo ""
    echo "Required:"
    echo "  basename          Base name for output files"
    echo "  replica           Replica number"
    echo "  start_time        Start time in ns"
    echo "  end_time          End time in ns"
    echo "  dt                Segment length in ns"
    echo ""
    echo "Optional:"
    echo "  --template <file>     MDP template (default: md.mdp)"
    echo "  --initial <name>      Initial files basename (default: npt)"
    echo "  --timestep <ps>       Timestep in ps (default: 0.002)"
    echo "  --title <suffix>      Title suffix for MDP"
    echo "  --append              Single output file (extends continuously)"
    echo "  --noappend            Multiple part files (default)"
    echo "  --tpr <file>          Start from existing TPR"
    echo "  --cpt <file>          Checkpoint file for --tpr (optional)"
    echo "  --topology <file>     Topology file (default: topol.top)"
    echo ""
    echo "Examples:"
    echo "  $0 md 0 0 500 100 --template production.mdp --noappend"
    echo "  $0 md 0 0 500 100 --tpr md_0.tpr --cpt md_0.cpt --append"
else
    # Check if first 5 args look like optional flags (common mistake)
    for arg in "$1" "$2" "$3" "$4" "$5"; do
        if [[ "$arg" == --* ]]; then
            echo "Error: Found optional flag '$arg' in required argument position"
            echo ""
            echo "You must provide 5 required arguments BEFORE any optional flags:"
            echo "  $0 <basename> <replica> <start_time> <end_time> <dt> [OPTIONS]"
            echo ""
            echo "Example:"
            echo "  $0 md 0 0 500 100 --template production.mdp"
            echo "     ^^ ^  ^   ^   ^^  (5 required args first)"
            exit 1
        fi
    done

    basename_arg="$1"
    replica="$2"
    tstart="$3"
    tend="$4"
    dt="$5"
    shift 5

    # Default values
    template_mdp="md.mdp"
    initial_basename="npt"
    timestep_ps="0.002"
    title_suffix=""
    append_mode="no"
    existing_tpr=""
    existing_cpt=""
    topology="topol.top"

    # Parse optional arguments
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --template)
                template_mdp="$2"
                shift 2
                ;;
            --initial)
                initial_basename="$2"
                shift 2
                ;;
            --timestep)
                timestep_ps="$2"
                shift 2
                ;;
            --title)
                title_suffix="$2"
                shift 2
                ;;
            --append)
                append_mode="yes"
                shift
                ;;
            --noappend)
                append_mode="no"
                shift
                ;;
            --tpr)
                existing_tpr="$2"
                shift 2
                ;;
            --cpt)
                existing_cpt="$2"
                shift 2
                ;;
            --topology)
                topology="$2"
                shift 2
                ;;
            *)
                echo "Error: Unknown option '$1'"
                exit 1
                ;;
        esac
    done
fi

# =============================================================================
# Common setup for both modes
# =============================================================================

base_name="${basename_arg}_${replica}"

# Validate numeric arguments
if ! [[ "${tstart}" =~ ^[0-9]+$ ]] || ! [[ "${tend}" =~ ^[0-9]+$ ]] || ! [[ "${dt}" =~ ^[0-9]+$ ]]; then
    echo "Error: start_time, end_time, and dt must be integers"
    echo "Got: start=${tstart}, end=${tend}, dt=${dt}"
    exit 1
fi

if ! [[ "${replica}" =~ ^[0-9]+$ ]]; then
    echo "Error: replica must be an integer, got: ${replica}"
    exit 1
fi

# Calculate nsteps per segment
nsteps_per_segment=$(awk -v dt="$dt" -v ts="$timestep_ps" 'BEGIN {printf "%.0f\n", (dt * 1000) / ts}')

echo "=== MD Extension Script ==="
echo "Base name:       ${basename_arg}"
echo "Replica:         ${replica}"
echo "Time range:      ${tstart} -> ${tend} ns"
echo "Segment dt:      ${dt} ns (${nsteps_per_segment} steps)"
echo "Timestep:        ${timestep_ps} ps"
echo "Append mode:     ${append_mode}"
if [[ -n "${existing_tpr}" ]]; then
    echo "Starting TPR:    ${existing_tpr}"
    if [[ -n "${existing_cpt}" ]]; then
        echo "Starting CPT:    ${existing_cpt}"
    fi
else
    echo "Template MDP:    ${template_mdp}"
    echo "Initial files:   ${initial_basename}.gro/cpt/edr"
fi
echo "Topology:        ${topology}"
echo "=============================="

# Source GROMACS if not already available
if ! command -v gmx &> /dev/null; then
    if [[ -f /programs/gromacs-2025.2/bin/GMXRC.bash ]]; then
        . /programs/gromacs-2025.2/bin/GMXRC.bash
    else
        echo "Error: GROMACS not found. Please source GMXRC.bash or add gmx to PATH."
        exit 1
    fi
fi

# Check topology file (not needed if starting from existing TPR)
if [[ -z "${existing_tpr}" ]] && [[ ! -f ${topology} ]]; then
    echo "Error: Topology file '${topology}' not found"
    exit 1
fi

# Function to check if a segment completed successfully
# Args: segment_number, append_mode, target_time_ns
check_segment_complete() {
    local seg_num="$1"
    local append="$2"
    local target_time_ns="$3"
    local log_file
    
    if [[ "${append}" == "yes" ]]; then
        # In append mode, check if we've reached the target time
        log_file="${base_name}.log"
        
        if [[ ! -f ${log_file} ]]; then
            return 1
        fi
        
        # Check if log indicates completion
        if ! grep -q "Finished mdrun" ${log_file} 2>/dev/null; then
            return 1
        fi
        
        # Check the TPR file to see what time it's set to run until
        local tpr_end_time=$(gmx check -f ${base_name}.tpr 2>&1 | grep "Last frame" | tail -1 | awk '{print $NF}')
        
        if [[ -n "${tpr_end_time}" ]]; then
            # tpr_end_time is in ps, target_time_ns is in ns
            local tpr_end_ns=$(awk -v t="${tpr_end_time}" 'BEGIN {printf "%.0f", t/1000}')
            
            # If TPR is set to at least target time, segment is complete
            if [[ ${tpr_end_ns} -ge ${target_time_ns} ]]; then
                return 0
            fi
        fi
        
        return 1
    else
        # In noappend mode, first segment is main log, rest are partXXXX
        if [[ ${seg_num} -eq 0 ]]; then
            log_file="${base_name}.log"
        else
            log_file="${base_name}.part$(printf '%04d' ${seg_num}).log"
        fi
        
        if [[ ! -f ${log_file} ]]; then
            return 1
        fi
        
        if grep -q "Finished mdrun" ${log_file} 2>/dev/null; then
            return 0
        else
            return 1
        fi
    fi
}

# Determine output file pattern based on append mode
if [[ "${append_mode}" == "yes" ]]; then
    main_tpr="${base_name}.tpr"
    main_log="${base_name}.log"
    main_cpt="${base_name}.cpt"
    mdrun_append_flag="-append"
else
    main_tpr="${base_name}.tpr"
    main_log="${base_name}.log"
    main_cpt="${base_name}.cpt"
    mdrun_append_flag="-noappend"
fi

# Initialize or use existing TPR for first segment
segment_num=0
current_time=${tstart}
total_segments=$(( (tend - tstart) / dt ))

if [[ -n "${existing_tpr}" ]]; then
    # Starting from existing TPR
    if [[ ! -f ${existing_tpr} ]]; then
        echo "Error: Specified TPR file '${existing_tpr}' not found"
        exit 1
    fi
    
    echo "Using existing TPR: ${existing_tpr}"
    
    echo "Using existing TPR: ${existing_tpr}"
    
    # When using existing TPR, derive base_name from TPR filename
    # This ensures we work with the original files (critical for append mode
    # since checkpoint contains filenames and checksums)
    base_name="${existing_tpr%.tpr}"
    main_tpr="${base_name}.tpr"
    main_cpt="${base_name}.cpt"
    
    echo "Base name set to: ${base_name}"
    
    # Copy TPR only if path is different (e.g., if in subdirectory)
    if [[ "${existing_tpr}" != "${main_tpr}" ]]; then
        cp "${existing_tpr}" "${main_tpr}"
        echo "Copied to ${main_tpr}"
    fi
    
    # Handle checkpoint file if provided
    if [[ -n "${existing_cpt}" ]]; then
        if [[ ! -f ${existing_cpt} ]]; then
            echo "Error: Specified checkpoint file '${existing_cpt}' not found"
            exit 1
        fi
        
        echo "Using existing checkpoint: ${existing_cpt}"
        
        # Copy to expected location if different
        if [[ "${existing_cpt}" != "${main_cpt}" ]]; then
            cp "${existing_cpt}" "${main_cpt}"
            echo "Copied to ${main_cpt}"
        fi
    else
        echo "No checkpoint specified. Will run without -cpi if none exists."
    fi
    
    # Check if first segment already completed
    if check_segment_complete 0 "${append_mode}" ${dt}; then
        echo "First segment already completed. Will continue from segment 2..."
        segment_num=1
        current_time=$((tstart + dt))
    fi
    
elif [[ ${tstart} -eq 0 ]]; then
    # Create initial TPR from template and equilibration files
    
    # Check if already completed
    if [[ -f ${main_tpr} ]] && check_segment_complete 0 "${append_mode}" ${dt}; then
        echo "First segment already completed. Will continue from segment 2..."
        segment_num=1
        current_time=${dt}
    else
        echo "Creating initial TPR for segment 1 (0 -> ${dt} ns)..."
        
        # Check for required input files
        if [[ ! -f ${template_mdp} ]]; then
            echo "Error: Template MDP file '${template_mdp}' not found"
            exit 1
        fi
        
        if [[ ! -f ${initial_basename}.gro ]]; then
            echo "Error: Initial structure '${initial_basename}.gro' not found"
            exit 1
        fi
        
        if [[ ! -f ${initial_basename}.edr ]]; then
            echo "Error: Initial energy file '${initial_basename}.edr' not found"
            exit 1
        fi
        
        # Create MDP with correct parameters
        first_mdp="${base_name}_init.mdp"
        sed -e "s/\(tinit\s*=\s*\)[0-9]\+/\10/" \
            -e "s/\(nsteps\s*=\s*\)[0-9]\+/\1${nsteps_per_segment}/" \
            ${template_mdp} > ${first_mdp}
        
        if [[ -n "${title_suffix}" ]]; then
            sed -i "s/\(title\s*=\s*.*\)/\1,${title_suffix}/" ${first_mdp}
        fi
        
        echo "Running grompp for initial segment..."
        if [[ -f ${initial_basename}.cpt ]]; then
            gmx grompp -f ${first_mdp} -c ${initial_basename}.gro \
                       -t ${initial_basename}.cpt -e ${initial_basename}.edr \
                       -p ${topology} -o ${main_tpr}
        else
            echo "Warning: No checkpoint file (${initial_basename}.cpt), continuing without"
            gmx grompp -f ${first_mdp} -c ${initial_basename}.gro \
                       -e ${initial_basename}.edr -p ${topology} -o ${main_tpr}
        fi
        
        echo "Running mdrun for segment 1..."
        if [[ -f ${base_name}.cpt ]]; then
            # Checkpoint exists, use it for crash recovery
            gmx mdrun -deffnm ${base_name} -cpi ${base_name}.cpt
        else
            # No checkpoint yet, first run
            gmx mdrun -deffnm ${base_name}
        fi
        
        segment_num=1
        current_time=${dt}
    fi
fi

# Continue with remaining segments using TPR extension
for ((seg=segment_num; seg<total_segments; seg++)); do
    segment_time=$((tstart + (seg + 1) * dt))
    segment_time_ps=$((segment_time * 1000))
    
    echo ""
    echo "=== Segment $((seg + 1))/${total_segments}: ${current_time} -> ${segment_time} ns ==="
    
    # Check if segment already completed
    if check_segment_complete $((seg + 1)) "${append_mode}" ${segment_time}; then
        echo "Segment $((seg + 1)) already completed. Skipping..."
        current_time=${segment_time}
        continue
    fi
    
    # Extend TPR by fixed increment (dt in ps)
    # Note: -extend is INCREMENTAL, not absolute time
    extend_by=$((dt * 1000))
    echo "Extending TPR by ${dt} ns (${extend_by} ps)..."
    extended_tpr="${main_tpr}"
    
    gmx convert-tpr -s ${main_tpr} -o ${extended_tpr} -extend ${extend_by}
    
    # Determine checkpoint file for this segment
    # IMPORTANT: In both append and noappend modes, checkpoint is always ${base_name}.cpt
    # The part000X suffix only applies to trajectory/energy/log files, NOT checkpoint!
    current_cpt="${base_name}.cpt"
    
    # Run mdrun with checkpoint continuation
    echo "Running mdrun (segment $((seg + 1)))..."
    
    if [[ "${append_mode}" == "yes" ]]; then
        # Append mode: continue to same files
        if [[ -f ${current_cpt} ]]; then
            gmx mdrun -deffnm ${base_name} -cpi ${current_cpt} -s ${extended_tpr} -append
        else
            echo "Warning: Checkpoint ${current_cpt} not found, starting from structure"
            gmx mdrun -deffnm ${base_name} -s ${extended_tpr} -append
        fi
    else
        # Noappend mode: creates part000X files
        if [[ -f ${current_cpt} ]]; then
            gmx mdrun -deffnm ${base_name} -cpi ${current_cpt} -s ${extended_tpr} -noappend
        else
            echo "Warning: Checkpoint ${current_cpt} not found, starting from structure"
            gmx mdrun -deffnm ${base_name} -s ${extended_tpr} -noappend
        fi
    fi
    
    current_time=${segment_time}
    
    # Check for STOP file
    if [[ -f ./STOP ]]; then
        echo "Found STOP file. Stopping."
        exit 0
    fi
done

echo ""
echo "=== All segments completed successfully (${tstart} -> ${tend} ns) ==="

if [[ "${append_mode}" == "yes" ]]; then
    echo "Output files: ${base_name}.xtc, ${base_name}.edr, ${base_name}.log"
else
    echo "Output files: ${base_name}.xtc, ${base_name}.part0001.xtc, ... ${base_name}.part$(printf '%04d' $((total_segments - 1))).xtc"
    echo ""
    echo "To concatenate trajectories:"
    echo "  gmx trjcat -f ${base_name}.xtc ${base_name}.part*.xtc -o ${base_name}_complete.xtc -cat"
fi

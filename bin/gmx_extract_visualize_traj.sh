#!/bin/bash -e
set -o pipefail
# Extract, clean up (PBC/center/fit) and convert a GROMACS trajectory for visualization.
#
# Pipeline stages (run in this fixed order, subset selectable via --steps):
#   extract  Plain group extraction (optional -b/-e/--dt/--skip slice), no PBC treatment
#   pbc      Group extraction + "-pbc whole" (makes the selected group whole across PBC)
#   center   Center the group in the box, keep molecules whole ("-pbc mol -center -ur ...")
#   fit      Least-squares fit of the trajectory onto the reference structure
#   pdb      Write the final frame set as a multi-model PDB trajectory
#
# Default pipeline (matches the reference 4-command workflow): pbc,center,fit,pdb
#
# Usage: ./gmx_extract_visualize_traj.sh -s <tpr> -f <xtc> [OPTIONS]
#
# Required:
#   -s, --tpr <file>          Run input file (.tpr)
#   -f, --xtc <file>          Trajectory file (.xtc/.trr/...)
#
# Group selection:
#   -g, --group <name>        Group to extract/output (default: "Protein")
#   --center-group <name>     Group used for centering (default: same as --group)
#   --fit-group <name>        Group used for the least-squares fit (default: same as --group)
#
# Trajectory slicing (applied only to the first stage that reads the raw input):
#   -b, --start <ps>          Start time (ps)
#   -e, --end <ps>            End time (ps)
#   --dt <ps>                 Only write frames at this time interval (ps)
#   --skip <n>                Only write every nth frame
#
# Pipeline control:
#   --steps <list>            Comma-separated subset of: extract,pbc,center,fit,pdb
#                             (default: pbc,center,fit,pdb)
#   --fit-mode <mode>         gmx trjconv -fit value (default: progressive)
#   --ur <rect|tric|compact>  gmx trjconv -ur value for the center step (default: compact)
#
# Non-contiguous atom groups (e.g. Protein + Ion):
#   --subset-tpr              After extraction, build a reduced .tpr (via "gmx convert-tpr")
#                             containing only the extracted atoms, and use it (with output
#                             group "System") for the center/fit/pdb steps. Needed because
#                             trjconv/tpr atom indices no longer match once a non-contiguous
#                             group has been extracted into its own .xtc. Off by default,
#                             since contiguous selections (e.g. plain "Protein") don't need it.
#                             --center-group/--fit-group default to the extract group, or to
#                             "System" once --subset-tpr swaps in the reduced .tpr (since the
#                             extract group itself, e.g. "non-Water", may no longer exist as a
#                             named group there). Pass --center-group/--fit-group explicitly
#                             (e.g. "C-alpha") to use a group that still exists in the reduced
#                             .tpr instead.
#
# Output:
#   -o, --output-basename <n> Basename for output files (default: derived from --xtc)
#   --outdir <dir>            Output directory (default: current directory)
#   -n, --index <file>        Optional index (.ndx) file passed to every trjconv call
#
# Misc:
#   --gmx <path>              gmx binary/command to use (default: "gmx")
#   --dry-run                 Print commands (and group selections) without running them
#   --save-script <file>      Also write the commands to a runnable bash script
#   -h, --help                Show this help and exit
#
# Available default groups (from a typical protein .tpr, no index file required):
#   System, Protein, Protein-H, C-alpha, Backbone, MainChain, MainChain+Cb,
#   MainChain+H, SideChain, SideChain-H, Prot-Masses, non-Protein, Water, SOL,
#   non-Water, Ion, Water_and_ions
#
# Examples:
#   # Full default pipeline (pbc -> center -> fit -> pdb) on the Protein group
#   ./gmx_extract_visualize_traj.sh -s md.tpr -f md.xtc
#
#   # Only extract a 0-100 ns slice of the Protein group, no PBC/center/fit/pdb
#   ./gmx_extract_visualize_traj.sh -s md.tpr -f md.xtc --steps extract -b 0 -e 100000
#
#   # Extract + remove PBC only (no centering/fit/pdb)
#   ./gmx_extract_visualize_traj.sh -s md.tpr -f md.xtc --steps pbc -b 0 -e 100000
#
#   # Resample the whole system every 100 ps, no PBC treatment
#   ./gmx_extract_visualize_traj.sh -s md.tpr -f md.xtc -g System --steps extract --dt 100
#
#   # Dry run of the full pipeline, saving the commands to a script
#   ./gmx_extract_visualize_traj.sh -s md.tpr -f md.xtc --dry-run --save-script run_viz.sh
#
#   # Non-contiguous group (protein + ion): build a reduced .tpr to keep indices consistent
#   ./gmx_extract_visualize_traj.sh -s md.tpr -f md.xtc -g Protein_Ion --fit-group C-alpha \
#       --index index.ndx --subset-tpr

print_usage() {
    sed -n '2,60p' "$0" | sed -e 's/^# \{0,1\}//'
}

if [[ $# -eq 0 ]]; then
    print_usage
    exit 1
fi

# Defaults
tpr_file=""
xtc_file=""
group="Protein"
center_group=""
fit_group=""
start_time=""
end_time=""
dt=""
skip=""
steps_arg="pbc,center,fit,pdb"
fit_mode="progressive"
ur_mode="compact"
output_basename=""
outdir="."
index_file=""
gmx_bin="gmx"
dry_run=false
save_script=""
subset_tpr=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -s|--tpr)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            tpr_file="$2"; shift 2 ;;
        -f|--xtc)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            xtc_file="$2"; shift 2 ;;
        -g|--group)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            group="$2"; shift 2 ;;
        --center-group)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            center_group="$2"; shift 2 ;;
        --fit-group)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            fit_group="$2"; shift 2 ;;
        -b|--start)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            start_time="$2"; shift 2 ;;
        -e|--end)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            end_time="$2"; shift 2 ;;
        --dt)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            dt="$2"; shift 2 ;;
        --skip)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            skip="$2"; shift 2 ;;
        --steps)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            steps_arg="$2"; shift 2 ;;
        --fit-mode)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            fit_mode="$2"; shift 2 ;;
        --ur)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            ur_mode="$2"; shift 2 ;;
        -o|--output-basename)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            output_basename="$2"; shift 2 ;;
        --outdir)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            outdir="$2"; shift 2 ;;
        -n|--index)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            index_file="$2"; shift 2 ;;
        --gmx)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            gmx_bin="$2"; shift 2 ;;
        --dry-run)
            dry_run=true; shift ;;
        --save-script)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            save_script="$2"; shift 2 ;;
        --subset-tpr)
            subset_tpr=true; shift ;;
        -h|--help)
            print_usage; exit 0 ;;
        --*)
            echo "Error: Unknown option '$1'"
            echo "Run with --help to see valid options"
            exit 1 ;;
        *)
            echo "Error: Unexpected argument '$1'"
            echo "All options must start with - or --"
            exit 1 ;;
    esac
done

# Validate required arguments
if [[ -z "${tpr_file}" ]]; then
    echo "Error: --tpr <file> is required"
    exit 1
fi
if [[ -z "${xtc_file}" ]]; then
    echo "Error: --xtc <file> is required"
    exit 1
fi
if [[ ! -f "${tpr_file}" ]]; then
    echo "Error: TPR file '${tpr_file}' not found"
    exit 1
fi
if [[ ! -f "${xtc_file}" ]]; then
    echo "Error: Trajectory file '${xtc_file}' not found"
    exit 1
fi
if [[ -n "${index_file}" ]] && [[ ! -f "${index_file}" ]]; then
    echo "Error: Index file '${index_file}' not found"
    exit 1
fi

# Validate and expand requested steps against the canonical pipeline order
IFS=',' read -ra requested_steps <<< "${steps_arg}"
canonical_steps=(extract pbc center fit pdb)
declare -A step_requested
for s in "${requested_steps[@]}"; do
    s="$(echo "$s" | xargs)"  # trim whitespace
    [[ -z "$s" ]] && continue
    valid=false
    for c in "${canonical_steps[@]}"; do
        [[ "$s" == "$c" ]] && valid=true
    done
    if [[ "${valid}" == false ]]; then
        echo "Error: Unknown step '${s}'"
        echo "Valid steps: ${canonical_steps[*]}"
        exit 1
    fi
    step_requested[$s]=true
done
if [[ ${#step_requested[@]} -eq 0 ]]; then
    echo "Error: --steps produced an empty pipeline"
    exit 1
fi

mkdir -p "${outdir}"

# Derive output basename from the xtc filename if not provided
if [[ -z "${output_basename}" ]]; then
    output_basename="$(basename "${xtc_file}")"
    output_basename="${output_basename%.*}"
fi

# Time-slice suffix, applied only to the first stage that reads the raw input
slice_suffix=""
if [[ -n "${start_time}" ]] || [[ -n "${end_time}" ]]; then
    slice_suffix="${slice_suffix}_${start_time:-0}_${end_time:-end}ps"
fi
if [[ -n "${dt}" ]]; then
    slice_suffix="${slice_suffix}_dt${dt}ps"
fi
if [[ -n "${skip}" ]]; then
    slice_suffix="${slice_suffix}_skip${skip}"
fi

# Source GROMACS if not already available (skip check in dry-run mode)
if [[ "${dry_run}" == false ]] && ! command -v "${gmx_bin%% *}" &> /dev/null; then
    if [[ -f /programs/gromacs-2025.2/bin/GMXRC.bash ]]; then
        . /programs/gromacs-2025.2/bin/GMXRC.bash
    else
        echo "Error: GROMACS not found. Please source GMXRC.bash manually or add gmx to PATH."
        exit 1
    fi
fi

if [[ -n "${save_script}" ]]; then
    {
        echo "#!/bin/bash -e"
        echo "# Generated by gmx_extract_visualize_traj.sh on $(date)"
    } > "${save_script}"
    chmod +x "${save_script}"
fi

index_flag=""
[[ -n "${index_file}" ]] && index_flag="-n ${index_file}"

echo "=== Trajectory Extraction/Visualization Prep ==="
echo "TPR:             ${tpr_file}"
echo "XTC:             ${xtc_file}"
echo "Extract group:   ${group}"
echo "Center group:    ${center_group:-<default: extract group, or System after --subset-tpr>}"
echo "Fit group:       ${fit_group:-<default: extract group, or System after --subset-tpr>}"
echo "Steps:           ${steps_arg}"
[[ -n "${start_time}" ]] && echo "Start:           ${start_time} ps"
[[ -n "${end_time}" ]]   && echo "End:             ${end_time} ps"
[[ -n "${dt}" ]]         && echo "dt:              ${dt} ps"
[[ -n "${skip}" ]]       && echo "skip:            ${skip}"
echo "Output basename: ${output_basename}"
echo "=================================================="

current_input="${xtc_file}"
current_tag=""
raw_input_consumed=false
# active_tpr/active_group switch to the reduced subset .tpr/"System" once
# run_convert_tpr fires (only relevant for center/fit/pdb, and only if --subset-tpr)
active_tpr="${tpr_file}"
active_group="${group}"

# Runs one gmx trjconv call, feeding the given group selections via stdin.
# Args: out_file, group1 [group2], extra gmx trjconv flags...
run_trjconv() {
    local out_file="$1"; shift
    local groups_csv="$1"; shift
    local extra_flags=("$@")

    local cmd=("${gmx_bin}" trjconv -s "${active_tpr}" -f "${current_input}" -o "${out_file}")
    [[ -n "${index_flag}" ]] && cmd+=(${index_flag})
    cmd+=("${extra_flags[@]}")

    # Apply slice flags only to the first stage that reads the raw input
    if [[ "${raw_input_consumed}" == false ]]; then
        [[ -n "${start_time}" ]] && cmd+=(-b "${start_time}")
        [[ -n "${end_time}" ]]   && cmd+=(-e "${end_time}")
        [[ -n "${dt}" ]]         && cmd+=(-dt "${dt}")
        [[ -n "${skip}" ]]       && cmd+=(-skip "${skip}")
        raw_input_consumed=true
    fi

    IFS=',' read -ra groups <<< "${groups_csv}"

    echo ""
    echo "--- Command: ${cmd[*]}"
    echo "--- Group selection(s): ${groups[*]}"

    if [[ -n "${save_script}" ]]; then
        {
            echo ""
            echo "printf '%s\\n' \\"
            for g in "${groups[@]}"; do
                echo "    '${g}' \\"
            done
            echo "  | ${cmd[*]}"
        } >> "${save_script}"
    fi

    if [[ "${dry_run}" == true ]]; then
        echo "(dry run: not executed)"
        current_input="${out_file}"
        return
    fi

    printf '%s\n' "${groups[@]}" | "${cmd[@]}"
    current_input="${out_file}"
}

# Builds a reduced .tpr (via "gmx convert-tpr") containing only the extracted atoms,
# so tpr/xtc indices stay aligned for non-contiguous groups. Switches active_tpr and
# active_group so subsequent steps read it with output group "System".
run_convert_tpr() {
    local subset_tpr_out="${outdir}/${output_basename}${slice_suffix}_subset.tpr"
    local cmd=("${gmx_bin}" convert-tpr -s "${tpr_file}" -o "${subset_tpr_out}")
    [[ -n "${index_flag}" ]] && cmd+=(${index_flag})

    echo ""
    echo "--- Command: ${cmd[*]}"
    echo "--- Group selection(s): ${group}   (builds a reduced .tpr matching the extracted atoms)"

    if [[ -n "${save_script}" ]]; then
        echo "" >> "${save_script}"
        echo "printf '%s\n' '${group}' | ${cmd[*]}" >> "${save_script}"
    fi

    if [[ "${dry_run}" == true ]]; then
        echo "(dry run: not executed)"
    else
        printf '%s\n' "${group}" | "${cmd[@]}"
    fi

    active_tpr="${subset_tpr_out}"
    active_group="System"
    echo "Using reduced TPR for subsequent steps: ${active_tpr} (output group now: ${active_group})"
}

# 1. extract - plain group extraction, no PBC treatment
if [[ "${step_requested[extract]}" == true ]]; then
    current_tag="extract${slice_suffix}"
    out_file="${outdir}/${output_basename}${slice_suffix}_extract.xtc"
    run_trjconv "${out_file}" "${group}"
fi

# 2. pbc - group extraction + make whole across PBC
if [[ "${step_requested[pbc]}" == true ]]; then
    out_file="${outdir}/${output_basename}${slice_suffix}_whole.xtc"
    run_trjconv "${out_file}" "${group}" -pbc whole
    current_tag="whole"
fi

# Build the reduced .tpr before center/fit/pdb, if requested and needed
if [[ "${subset_tpr}" == true ]] && { [[ "${step_requested[center]}" == true ]] || [[ "${step_requested[fit]}" == true ]] || [[ "${step_requested[pdb]}" == true ]]; }; then
    run_convert_tpr
fi

# 3. center - center group in box, keep molecules whole
if [[ "${step_requested[center]}" == true ]]; then
    out_file="${outdir}/${output_basename}${slice_suffix}_center.xtc"
    run_trjconv "${out_file}" "${center_group:-${active_group}},${active_group}" -pbc mol -center -ur "${ur_mode}"
    current_tag="center"
fi

# 4. fit - least-squares fit onto the reference structure
if [[ "${step_requested[fit]}" == true ]]; then
    out_file="${outdir}/${output_basename}${slice_suffix}_fit.xtc"
    run_trjconv "${out_file}" "${fit_group:-${active_group}},${active_group}" -fit "${fit_mode}"
    current_tag="fit"
fi

# 5. pdb - write the (final) trajectory as a PDB trajectory
if [[ "${step_requested[pdb]}" == true ]]; then
    [[ -z "${current_tag}" ]] && current_tag="extract${slice_suffix}"
    out_file="${outdir}/${output_basename}${slice_suffix}_${current_tag}.pdb"
    # avoid duplicating the tag if it was already appended to the xtc name above
    case "${current_tag}" in
        extract*) out_file="${outdir}/${output_basename}${slice_suffix}_extract.pdb" ;;
        *)        out_file="${outdir}/${output_basename}${slice_suffix}_${current_tag}.pdb" ;;
    esac
    run_trjconv "${out_file}" "${active_group}"
fi

echo ""
if [[ "${dry_run}" == true ]]; then
    echo "Dry run complete. No files were written."
else
    echo "Pipeline complete. Final output: ${current_input}"
fi
[[ -n "${save_script}" ]] && echo "Commands saved to: ${save_script}"

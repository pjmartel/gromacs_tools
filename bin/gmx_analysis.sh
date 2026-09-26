#!/bin/bash -e
set -o pipefail
# Standard post-processing analysis for a (possibly segmented) GROMACS production run.
#
# Given a run produced by gmx_continue_grompp.sh (files named
# <basename>_<replica>_<start>_<end>.{xtc,edr,tpr,...}), this script:
#   1. Finds every trajectory segment for <basename>/<replica> in --segments-dir
#   2. Concatenates the segments (gmx trjcat / gmx eneconv) into one trajectory+energy file
#   3. Runs the standard analyses on the (optionally time-sliced) concatenated data:
#        - energy terms (Temperature, Pressure, Potential, Total-Energy by default)
#        - RMSD to a reference structure (gmx rms)
#        - radius of gyration (gmx gyrate)
#        - solvent accessible surface area (gmx sasa)
#        - secondary structure over time (gmx dssp, GROMACS >= 2023)
#        - per-residue RMS fluctuations (gmx rmsf)
#   Each analysis is written to its own clearly-named .xvg file in --output-dir.
#
# Usage: ./gmx_analysis.sh <basename> <replica> [OPTIONS]
#   basename:  Base name used for the run (e.g. "MnMT4_apo")
#   replica:   Replica number (e.g. 0), so segments are "<basename>_<replica>_<start>_<end>"
#
# File discovery:
#   --segments-dir <dir>       Directory containing the segment files (default: .)
#   --structure <file>         Reference .tpr/.gro/.pdb used for all analyses
#                              (default: the .tpr of the earliest/first segment found)
#   --skip-concat              Reuse an existing "<output-dir>/<basename>_<replica>_concat.{xtc,edr}"
#                              instead of rebuilding it from the segments
#
# Time range / sampling:
#   --begin <time>             First frame to analyze (-b), in the unit given by --tu
#   --end <time>                Last frame to analyze (-e), in the unit given by --tu
#   --dt <time>                 Only use frames spaced by this interval (-dt), always in ps
#                              regardless of --tu
#   --tu <ps|ns|us|fs>          Unit for --begin/--end (default: ps). Values are
#                              converted to ps internally before being passed to any tool,
#                              so results stay consistent across every analysis.
#   --plot-tu <ps|ns|us|fs>     Time unit for the x-axis of every generated .xvg plot
#                              (passed as -tu to each gmx tool). Independent of --tu,
#                              which only controls how --begin/--end are interpreted
#                              (default: ns)
#
# Group/selection:
#   -n, --index <file>         Optional index (.ndx) file passed to every tool that accepts one
#   -g, --group <name>         Group/selection used for RMSD (calc), gyration, SASA (surface),
#                              secondary structure and RMSF (default: Protein)
#   --fit-group <name>         Group used for the RMSD least-squares fit (default: same as --group)
#
# Analysis selection:
#   --skip <list>               Comma-separated categories to skip: energy,rmsd,gyration,sasa,
#                              secondary-structure,rmsf ("dssp" is an alias for secondary-structure)
#   --only <list>                Comma-separated categories to run exclusively (same names)
#   --energy-terms <list>       Comma-separated gmx energy term names, one .xvg per term
#                              (default: Temperature,Pressure,Potential,Total-Energy)
#   --probe-radius <nm>          SASA solvent probe radius, gmx sasa -probe (default: 0.14)
#   --dssp-hmode <gromacs|dssp>  gmx dssp -hmode (default: gromacs)
#   --per-atom                   RMSF per-atom instead of the default per-residue (-res)
#
# Output:
#   --output-dir <dir>          Where .xvg/.dat/log files are written
#                              (default: "<basename>_<replica>_analysis")
#   --xvg-format <xmgrace|xmgr|none>  Passed as -xvg to every tool that supports it (default: xmgrace)
#
# Misc:
#   --gmx <path>                gmx binary/command to use (default: "gmx")
#   --dry-run                    Print commands (and stdin group selections) without running them
#   --save-script <file>         Reproducibility script recording every command run (and the
#                              exact invocation of this tool), written even with --dry-run
#                              (default: "<output-dir>/gmx_analysis_commands.sh")
#   -h, --help                   Show this help and exit
#
# Available default groups (from a typical protein .tpr, no index file required):
#   System, Protein, Protein-H, C-alpha, Backbone, MainChain, MainChain+Cb,
#   MainChain+H, SideChain, SideChain-H, Prot-Masses, non-Protein, Water, SOL,
#   non-Water, Ion, Water_and_ions
#
# Examples:
#   # Full analysis of MnMT4_apo replica 0 (segments MnMT4_apo_0_0_500.xtc, _500_1000.xtc, ...)
#   ./gmx_analysis.sh MnMT4_apo 0
#
#   # Only the last 500 ns, sampled every 100 ps
#   ./gmx_analysis.sh MnMT4_apo 0 --begin 500 --end 1000 --tu ns --dt 100
#
#   # Full 1 us trajectory, plots in ns (the default) instead of ps
#   ./gmx_analysis.sh MnMT4_apo 0 --plot-tu ns
#
#   # Keep the old ps-based plots instead of the new ns default
#   ./gmx_analysis.sh MnMT4_apo 0 --plot-tu ps
#
#   # RMSD/RMSF/gyration relative to the equilibrated structure, backbone-only RMSD
#   ./gmx_analysis.sh MnMT4_apo 0 --structure equil_npt2.gro --group Backbone
#
#   # Skip the (slow) secondary-structure analysis, only recompute RMSD and RMSF
#   ./gmx_analysis.sh MnMT4_apo 0 --only rmsd,rmsf
#
#   # Dry run, saving the exact commands to a script for later/offline execution
#   ./gmx_analysis.sh MnMT4_apo 0 --dry-run --save-script run_analysis.sh

print_usage() {
    sed -n '3,84p' "$0" | sed -e 's/^# \{0,1\}//'
}

original_invocation="$0 $*"

if [[ $# -eq 0 ]]; then
    print_usage
    exit 1
fi

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    print_usage
    exit 0
fi

# Validate that first 2 arguments are not flags (before shifting)
for i in 1 2; do
    eval "arg=\${$i}"
    if [[ "$arg" == --* ]]; then
        echo "Error: Argument $i appears to be a flag ('$arg')"
        echo "The first 2 arguments must be: <basename> <replica>"
        exit 1
    fi
done

basename_arg="$1"
replica="$2"
shift 2

# Defaults
segments_dir="."
output_dir=""
structure_file=""
index_file=""
group="Protein"
fit_group=""
begin_time=""
end_time=""
dt_time=""
time_unit="ps"
plot_time_unit="ns"
skip_list=""
only_list=""
energy_terms="Temperature,Pressure,Potential,Total-Energy"
probe_radius="0.14"
dssp_hmode="gromacs"
per_residue=true
skip_concat=false
xvg_format="xmgrace"
gmx_bin="gmx"
dry_run=false
save_script=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --segments-dir)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            segments_dir="$2"; shift 2 ;;
        --output-dir)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            output_dir="$2"; shift 2 ;;
        --structure)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            structure_file="$2"; shift 2 ;;
        -n|--index)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            index_file="$2"; shift 2 ;;
        -g|--group)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            group="$2"; shift 2 ;;
        --fit-group)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            fit_group="$2"; shift 2 ;;
        --begin)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            begin_time="$2"; shift 2 ;;
        --end)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            end_time="$2"; shift 2 ;;
        --dt)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            dt_time="$2"; shift 2 ;;
        --tu)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            time_unit="$2"; shift 2 ;;
        --plot-tu)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            plot_time_unit="$2"; shift 2 ;;
        --skip)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            skip_list="$2"; shift 2 ;;
        --only)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            only_list="$2"; shift 2 ;;
        --energy-terms)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            energy_terms="$2"; shift 2 ;;
        --probe-radius)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            probe_radius="$2"; shift 2 ;;
        --dssp-hmode)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            dssp_hmode="$2"; shift 2 ;;
        --per-atom)
            per_residue=false; shift ;;
        --skip-concat)
            skip_concat=true; shift ;;
        --xvg-format)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            xvg_format="$2"; shift 2 ;;
        --gmx)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            gmx_bin="$2"; shift 2 ;;
        --dry-run)
            dry_run=true; shift ;;
        --save-script)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            save_script="$2"; shift 2 ;;
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

fit_group="${fit_group:-${group}}"
base_name="${basename_arg}_${replica}"
output_dir="${output_dir:-${base_name}_analysis}"

# Default reproducibility script name, unless overridden with --save-script
save_script="${save_script:-${output_dir}/gmx_analysis_commands.sh}"

if [[ ! -d "${segments_dir}" ]]; then
    echo "Error: Segments directory '${segments_dir}' not found"
    exit 1
fi
if [[ -n "${index_file}" ]] && [[ ! -f "${index_file}" ]]; then
    echo "Error: Index file '${index_file}' not found"
    exit 1
fi

# Validate and normalize --skip/--only categories (one per line via mapfile)
valid_categories=(energy rmsd gyration sasa secondary-structure rmsf)
normalize_categories() {
    local raw="$1"
    IFS=',' read -ra parts <<< "${raw}"
    for p in "${parts[@]}"; do
        p="$(echo "${p}" | xargs)"
        [[ -z "${p}" ]] && continue
        [[ "${p}" == "dssp" ]] && p="secondary-structure"
        local ok=false
        for c in "${valid_categories[@]}"; do [[ "${p}" == "${c}" ]] && ok=true; done
        if [[ "${ok}" == false ]]; then
            echo "Error: Unknown analysis category '${p}'" >&2
            echo "Valid categories: ${valid_categories[*]} (dssp is an alias for secondary-structure)" >&2
            exit 1
        fi
        echo "${p}"
    done
}
skip_categories=()
only_categories=()
[[ -n "${skip_list}" ]] && mapfile -t skip_categories < <(normalize_categories "${skip_list}")
[[ -n "${only_list}" ]] && mapfile -t only_categories < <(normalize_categories "${only_list}")

should_run() {
    local name="$1"
    if [[ ${#only_categories[@]} -gt 0 ]]; then
        for c in "${only_categories[@]}"; do [[ "${c}" == "${name}" ]] && return 0; done
        return 1
    fi
    for c in "${skip_categories[@]}"; do [[ "${c}" == "${name}" ]] && return 1; done
    return 0
}

# Discover trajectory segments: <segments_dir>/<base_name>_<start>_<end>.xtc
shopt -s nullglob
candidates=("${segments_dir}/${base_name}_"[0-9]*_[0-9]*.xtc)
if [[ ${#candidates[@]} -eq 0 ]]; then
    echo "Error: No trajectory segments found matching '${base_name}_<start>_<end>.xtc' in '${segments_dir}'"
    exit 1
fi

seg_xtc=()
seg_edr=()
seg_tpr=()
while IFS=$'\t' read -r start end xtc; do
    edr="${segments_dir}/${base_name}_${start}_${end}.edr"
    tpr="${segments_dir}/${base_name}_${start}_${end}.tpr"
    if [[ ! -f "${edr}" ]]; then
        echo "Error: Missing energy file for segment ${start}-${end}: '${edr}'"
        exit 1
    fi
    if [[ ! -f "${tpr}" ]]; then
        echo "Error: Missing run input file for segment ${start}-${end}: '${tpr}'"
        exit 1
    fi
    seg_xtc+=("${xtc}")
    seg_edr+=("${edr}")
    seg_tpr+=("${tpr}")
done < <(
    for f in "${candidates[@]}"; do
        fname="$(basename "${f}")"
        if [[ "${fname}" =~ ^${base_name}_([0-9]+)_([0-9]+)\.xtc$ ]]; then
            printf '%s\t%s\t%s\n' "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}" "${f}"
        fi
    done | sort -n -k1,1
)

structure_file="${structure_file:-${seg_tpr[0]}}"
if [[ ! -f "${structure_file}" ]]; then
    echo "Error: Reference structure '${structure_file}' not found"
    exit 1
fi

# Convert --begin/--end/--dt to ps so every tool (even ones without -tu) agrees
convert_to_ps() {
    local value="$1" unit="$2"
    case "${unit}" in
        ps) echo "${value}" ;;
        ns) awk -v v="${value}" 'BEGIN{printf "%.6f", v*1000}' ;;
        us) awk -v v="${value}" 'BEGIN{printf "%.6f", v*1000000}' ;;
        fs) awk -v v="${value}" 'BEGIN{printf "%.6f", v/1000}' ;;
        *)
            echo "Error: Unsupported --tu unit '${unit}' (use ps, ns, us, or fs)" >&2
            exit 1 ;;
    esac
}
begin_ps=""; end_ps=""; dt_ps=""
[[ -n "${begin_time}" ]] && begin_ps="$(convert_to_ps "${begin_time}" "${time_unit}")"
[[ -n "${end_time}" ]]   && end_ps="$(convert_to_ps "${end_time}" "${time_unit}")"
# --dt is always in ps, regardless of --tu (which only affects --begin/--end)
[[ -n "${dt_time}" ]]    && dt_ps="$(convert_to_ps "${dt_time}" "ps")"

# Convert back from ps to --plot-tu units, since gmx tools interpret -b/-e/-dt
# according to whatever unit is passed via -tu
convert_from_ps() {
    local value_ps="$1" unit="$2"
    case "${unit}" in
        ps) echo "${value_ps}" ;;
        ns) awk -v v="${value_ps}" 'BEGIN{printf "%.6f", v/1000}' ;;
        us) awk -v v="${value_ps}" 'BEGIN{printf "%.6f", v/1000000}' ;;
        fs) awk -v v="${value_ps}" 'BEGIN{printf "%.6f", v*1000}' ;;
        *)
            echo "Error: Unsupported --plot-tu unit '${unit}' (use ps, ns, us, or fs)" >&2
            exit 1 ;;
    esac
}
begin_plot=""; end_plot=""; dt_plot=""
[[ -n "${begin_ps}" ]] && begin_plot="$(convert_from_ps "${begin_ps}" "${plot_time_unit}")"
[[ -n "${end_ps}" ]]   && end_plot="$(convert_from_ps "${end_ps}" "${plot_time_unit}")"
[[ -n "${dt_ps}" ]]    && dt_plot="$(convert_from_ps "${dt_ps}" "${plot_time_unit}")"
time_flags=(-tu "${plot_time_unit}")
[[ -n "${begin_plot}" ]] && time_flags+=(-b "${begin_plot}")
[[ -n "${end_plot}" ]]   && time_flags+=(-e "${end_plot}")
[[ -n "${dt_plot}" ]]    && time_flags+=(-dt "${dt_plot}")

# rmsf's x-axis is residue/atom index, not time, so it has no -tu flag;
# it still needs -b/-e/-dt (in ps) to select the trajectory range to average over.
traj_range_flags=()
[[ -n "${begin_ps}" ]] && traj_range_flags+=(-b "${begin_ps}")
[[ -n "${end_ps}" ]]   && traj_range_flags+=(-e "${end_ps}")
[[ -n "${dt_ps}" ]]    && traj_range_flags+=(-dt "${dt_ps}")

mkdir -p "${output_dir}"
log_file="${output_dir}/${base_name}_analysis.log"
{
    echo "GROMACS analysis log for ${base_name}"
    echo "Generated: $(date)"
    echo "=========================================="
} > "${log_file}"

if [[ -n "${save_script}" ]]; then
    {
        echo "#!/bin/bash -e"
        echo "# Commands executed by gmx_analysis.sh"
        echo "# Invocation: ${original_invocation}"
        echo "# Generated: $(date)"
    } > "${save_script}"
    chmod +x "${save_script}"
fi

echo "=== GROMACS Trajectory Analysis ==="
echo "Basename:        ${basename_arg}  (replica ${replica})"
echo "Segments dir:    ${segments_dir}"
echo "Segments found:  ${#seg_xtc[@]}"
for i in "${!seg_xtc[@]}"; do
    echo "  - $(basename "${seg_xtc[$i]}")"
done
echo "Structure:       ${structure_file}"
echo "Group:           ${group} (fit group for RMSD: ${fit_group})"
[[ -n "${index_file}" ]] && echo "Index:           ${index_file}"
[[ -n "${begin_ps}" ]] && echo "Begin:           ${begin_ps} ps"
[[ -n "${end_ps}" ]]   && echo "End:             ${end_ps} ps"
[[ -n "${dt_ps}" ]]    && echo "dt:              ${dt_ps} ps"
echo "Plot time unit:  ${plot_time_unit}"
echo "Energy terms:    ${energy_terms}"
echo "Output dir:      ${output_dir}"
echo "===================================="

index_flag=()
[[ -n "${index_file}" ]] && index_flag=(-n "${index_file}")

# Runs a plain gmx command (array), honoring --dry-run/--save-script, logging to log_file.
run_cmd() {
    local desc="$1"; shift
    echo ""
    echo "--- ${desc} ---"
    echo "Command: $*"
    if [[ -n "${save_script}" ]]; then
        printf '%q ' "$@" >> "${save_script}"
        echo "" >> "${save_script}"
    fi
    if [[ "${dry_run}" == true ]]; then
        echo "(dry run: not executed)"
        return 0
    fi
    if "$@" >> "${log_file}" 2>&1; then
        echo "✓ ${desc} completed"
        return 0
    else
        echo "✗ ${desc} failed (see ${log_file})"
        return 1
    fi
}

# Same as run_cmd, but feeds $stdin_data (already containing literal newlines) to the tool.
run_piped_cmd() {
    local desc="$1"; shift
    local stdin_data="$1"; shift
    echo ""
    echo "--- ${desc} ---"
    echo "Command: $*"
    echo "Group selection(s): $(echo "${stdin_data}" | tr '\n' ',' | sed 's/,$//')"
    if [[ -n "${save_script}" ]]; then
        {
            echo "printf '%s\\n' \\"
            while IFS= read -r line; do echo "    '${line}' \\"; done <<< "${stdin_data}"
            printf '  | '
            printf '%q ' "$@"
            echo ""
        } >> "${save_script}"
    fi
    if [[ "${dry_run}" == true ]]; then
        echo "(dry run: not executed)"
        return 0
    fi
    if printf '%s\n' "${stdin_data}" | "$@" >> "${log_file}" 2>&1; then
        echo "✓ ${desc} completed"
        return 0
    else
        echo "✗ ${desc} failed (see ${log_file})"
        return 1
    fi
}

has_gmx_dssp() {
    "${gmx_bin}" dssp -h >/dev/null 2>&1
}

# Source GROMACS if not already available (skip check in dry-run mode)
if [[ "${dry_run}" == false ]] && ! command -v "${gmx_bin%% *}" &> /dev/null; then
    if [[ -f /programs/gromacs-2025.2/bin/GMXRC.bash ]]; then
        . /programs/gromacs-2025.2/bin/GMXRC.bash
    else
        echo "Error: GROMACS not found. Please source GMXRC.bash manually or add gmx to PATH."
        exit 1
    fi
fi

# --- 1. Concatenate segments -------------------------------------------------
concat_xtc="${output_dir}/${base_name}_concat.xtc"
concat_edr="${output_dir}/${base_name}_concat.edr"

if [[ "${skip_concat}" == true ]]; then
    if [[ "${dry_run}" == false ]] && { [[ ! -f "${concat_xtc}" ]] || [[ ! -f "${concat_edr}" ]]; }; then
        echo "Error: --skip-concat given but concatenated files not found:"
        echo "  ${concat_xtc}"
        echo "  ${concat_edr}"
        exit 1
    fi
    echo ""
    echo "Skipping concatenation, reusing: ${concat_xtc}, ${concat_edr}"
else
    cmd=("${gmx_bin}" trjcat -f "${seg_xtc[@]}" -o "${concat_xtc}")
    run_cmd "Concatenate ${#seg_xtc[@]} trajectory segment(s)" "${cmd[@]}"

    cmd=("${gmx_bin}" eneconv -f "${seg_edr[@]}" -o "${concat_edr}")
    run_cmd "Concatenate ${#seg_edr[@]} energy segment(s)" "${cmd[@]}"
fi

# --- 2. Energy terms (temperature, pressure, potential, total energy, ...) --
if should_run energy; then
    IFS=',' read -ra terms <<< "${energy_terms}"
    for term in "${terms[@]}"; do
        term="$(echo "${term}" | xargs)"
        [[ -z "${term}" ]] && continue
        slug="$(echo "${term}" | tr '[:upper:]' '[:lower:]' | tr '-' '_')"
        out="${output_dir}/${slug}.xvg"
        cmd=("${gmx_bin}" energy -f "${concat_edr}" -s "${structure_file}" -o "${out}" -xvg "${xvg_format}" "${time_flags[@]}")
        run_piped_cmd "Energy term: ${term}" "${term}" "${cmd[@]}" || true
    done
fi

# --- 3. RMSD to reference structure -----------------------------------------
if should_run rmsd; then
    out="${output_dir}/rmsd.xvg"
    cmd=("${gmx_bin}" rms -s "${structure_file}" -f "${concat_xtc}" -o "${out}" -xvg "${xvg_format}" "${index_flag[@]}" "${time_flags[@]}")
    stdin_data="${fit_group}"$'\n'"${group}"
    run_piped_cmd "RMSD (fit=${fit_group}, calc=${group})" "${stdin_data}" "${cmd[@]}" || true
fi

# --- 4. Radius of gyration ---------------------------------------------------
if should_run gyration; then
    out="${output_dir}/gyration_radius.xvg"
    cmd=("${gmx_bin}" gyrate -s "${structure_file}" -f "${concat_xtc}" -o "${out}" -sel "${group}" -xvg "${xvg_format}" "${index_flag[@]}" "${time_flags[@]}")
    run_cmd "Radius of gyration (${group})" "${cmd[@]}" || true
fi

# --- 5. Solvent accessible surface area --------------------------------------
if should_run sasa; then
    out="${output_dir}/sasa.xvg"
    cmd=("${gmx_bin}" sasa -s "${structure_file}" -f "${concat_xtc}" -o "${out}" -surface "${group}" -probe "${probe_radius}" -xvg "${xvg_format}" "${index_flag[@]}" "${time_flags[@]}")
    run_cmd "Solvent accessible surface area (${group})" "${cmd[@]}" || true
fi

# --- 6. Secondary structure over time (gmx dssp, GROMACS >= 2023) ----------
if should_run secondary-structure; then
    if [[ "${dry_run}" == true ]] || has_gmx_dssp; then
        out_dat="${output_dir}/secondary_structure.dat"
        out_num="${output_dir}/secondary_structure_counts.xvg"
        cmd=("${gmx_bin}" dssp -s "${structure_file}" -f "${concat_xtc}" -sel "${group}" \
             -o "${out_dat}" -num "${out_num}" -hmode "${dssp_hmode}" -xvg "${xvg_format}" \
             "${index_flag[@]}" "${time_flags[@]}")
        run_cmd "Secondary structure (${group})" "${cmd[@]}" || true
    else
        echo ""
        echo "--- Secondary structure ---"
        echo "'gmx dssp' is not available in this GROMACS build (requires GROMACS >= 2023)."
        echo "Skipping secondary structure analysis."
    fi
fi

# --- 7. RMS fluctuations ------------------------------------------------------
if should_run rmsf; then
    out="${output_dir}/rmsf.xvg"
    cmd=("${gmx_bin}" rmsf -s "${structure_file}" -f "${concat_xtc}" -o "${out}" -xvg "${xvg_format}" "${index_flag[@]}" "${traj_range_flags[@]}")
    [[ "${per_residue}" == true ]] && cmd+=(-res)
    run_piped_cmd "RMS fluctuations (${group})" "${group}" "${cmd[@]}" || true
fi

echo ""
echo "=== Analysis complete ==="
if [[ "${dry_run}" == true ]]; then
    echo "Dry run: no files were written."
else
    echo "Output files in ${output_dir}:"
    for f in "${output_dir}"/*.xvg "${output_dir}"/*.dat; do
        [[ -f "${f}" ]] && echo "  - $(basename "${f}")"
    done
    echo "Full command log: ${log_file}"
fi

#!/bin/bash -e
set -o pipefail
# Standard post-processing analysis for a GROMACS production run: either a segmented run
# produced by gmx_continue_grompp.sh, or a single trajectory.
#
# Segmented run (files named <basename>_<replica>_<start>_<end>.{xtc,edr,tpr,...}):
#   1. Finds every trajectory segment for <basename>/<replica> in --segments-dir
#   2. Concatenates the segments (gmx trjcat / gmx eneconv) into one trajectory+energy file
# Single run (files <basename>.xtc, <basename>.tpr and optionally <basename>.edr, as written
# e.g. by gmx mdrun -deffnm): uses these files directly, nothing is concatenated.
# Then, in both cases:
#   3. Runs the standard analyses on the (optionally time-sliced) data:
#        - energy terms (Temperature, Pressure, Potential, Total-Energy by default)
#        - RMSD to a reference structure (gmx rms)
#        - radius of gyration (gmx gyrate)
#        - solvent accessible surface area (gmx sasa)
#        - secondary structure over time (gmx dssp, GROMACS >= 2023)
#        - per-residue RMS fluctuations (gmx rmsf)
#   Each analysis is written to its own clearly-named .xvg file in --output-dir.
#
# Usage: ./gmx_analysis.sh <basename> <replica> [OPTIONS]    segmented run
#        ./gmx_analysis.sh <basename> [OPTIONS]              single run
#   basename:  Base name used for the run (e.g. "MnMT4_apo"; for a single run e.g. "md" or
#              "md.xtc", optionally with a path)
#   replica:   Replica number (e.g. 0), so segments are "<basename>_<replica>_<start>_<end>".
#              Without it, <basename>.xtc is analyzed as a single run.
#
# File discovery:
#   --segments-dir <dir>       Directory containing the segment files, or a single run's files
#                              (default: .)
#   --structure <file>         Reference .tpr/.gro/.pdb used for all analyses (default: the .tpr
#                              of the earliest segment, or <basename>.tpr for a single run)
#   --skip-concat              Reuse an existing "<output-dir>/<basename>_<replica>_concat.{xtc,edr}"
#                              instead of rebuilding it from the segments (segmented runs only)
#   --concat-only              Only build the concatenated trajectory and energy files, then
#                              stop without running any analysis (cannot be combined with
#                              --only, --skip or --skip-concat)
#
# Time range / sampling:
#   --begin <time>             First frame to analyze (-b), in the unit given by --tu
#   --end <time>                Last frame to analyze (-e), in the unit given by --tu
#   --dt <time>                 Only use frames spaced by this interval (-dt), always in ps
#                              regardless of --tu (not applied to gmx energy, which has no -dt)
#   --tu <ps|ns|us|fs>          Unit for --begin/--end (default: ps). Values are
#                              converted to ps internally before being passed to any tool,
#                              so results stay consistent across every analysis.
#   --plot-tu <ps|ns|us|fs>     Time unit for the x-axis of every generated .xvg plot
#                              (as -tu, for tools that accept it; energy terms use --energy-tu). Independent of --tu,
#                              which only controls how --begin/--end are interpreted
#                              (default: ns)
#   --energy-tu <ps|ns|us|fs>   Time unit for the x-axis of the energy-term plots (default: ps).
#                              gmx energy has no -tu option, so the time column of its .xvg
#                              output is rescaled after the calculation
#
# Group/selection:
#   -n, --index <file>         Optional index (.ndx) file passed to every tool that accepts one
#   -g, --group <name>         Group/selection used for RMSD (calc), gyration, SASA (surface),
#                              secondary structure and RMSF (default: Protein)
#   --fit-group <name>         Group used for the RMSD least-squares fit (default: same as --group)
#
# Analysis selection:
#   --skip <list>               Comma-separated categories to skip: energy,rmsd,gyration,sasa,
#                              secondary-structure,rmsf,distance,angle ("dssp" is an alias for
#                              secondary-structure; distance and angle need --measure)
#   --only <list>                Comma-separated categories to run exclusively (same names)
#   --energy-terms <list>       Comma-separated gmx energy term names, one .xvg per term
#                              (default: Temperature,Pressure,Potential,Total-Energy)
#   --probe-radius <nm>          SASA solvent probe radius, gmx sasa -probe (default: 0.14)
#   --dssp-hmode <gromacs|dssp>  gmx dssp -hmode (default: gromacs)
#   --per-atom                   RMSF per-atom instead of the default per-residue (-res)
#
# Distances, angles and dihedrals:
#   --measure <file>            Measure the distances, angles and dihedrals listed in <file>, one
#                              .xvg per measurement (categories "distance" and "angle", the latter
#                              including dihedrals). Format, one measurement per line ('#' starts a
#                              comment line):
#                                <type>  <name>  <selection> | <selection> [| ...]
#                              e.g.
#                                distance  SAM_C36  resname SAM and name SD | resid 36 and name SG
#                                angle     hinge    resid 10 and name CA | com of group "Helix" | resid 90 and name CA
#                                dihedral  chi1_W6  resid 6 and name N | resid 6 and name CA | resid 6 and name CB | resid 6 and name CG
#                              type: distance (2 selections), angle (3) or dihedral (4); name: letters,
#                              digits, _ . - (used in the file name, e.g. distance_SAM_C36.xvg).
#                              Selections use GROMACS selection syntax (see 'gmx help selections'),
#                              with groups from -n/--index; each must match exactly one atom or be a
#                              "com of"/"cog of" a group. All are checked with gmx select first.
#
# Output:
#   --output-dir <dir>          Where .xvg/.dat/log files are written
#                              (default: "<basename>_<replica>_analysis", or "<basename>_analysis"
#                              for a single run)
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
#   # A single run: md.xtc, md.tpr (and md.edr, if present) in the current directory
#   ./gmx_analysis.sh md
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

# Prints the header comment block (from line 3 up to the first non-comment line)
print_usage() {
    awk 'NR >= 3 { if (!/^#/) exit; sub(/^# ?/, ""); print }' "$0"
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

# First argument: <basename>. Second (optional): <replica>. Without a replica (no second
# argument, or it is an option), <basename>.xtc is analyzed as a single run.
if [[ "$1" == -* ]]; then
    echo "Error: The first argument must be <basename>, got the option '$1'"
    echo "Usage: $0 <basename> <replica> [OPTIONS]   (segmented run)"
    echo "       $0 <basename> [OPTIONS]             (single run)"
    exit 1
fi
basename_arg="$1"
shift
if [[ $# -gt 0 && "$1" != -* ]]; then
    replica="$1"
    shift
    single_run=false
else
    replica=""
    single_run=true
    basename_arg="${basename_arg%.xtc}"     # accept "md.xtc" as well as "md"
fi

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
energy_time_unit="ps"
skip_list=""
only_list=""
energy_terms="Temperature,Pressure,Potential,Total-Energy"
probe_radius="0.14"
dssp_hmode="gromacs"
per_residue=true
skip_concat=false
concat_only=false
measure_file=""
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
        --energy-tu)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            energy_time_unit="$2"; shift 2 ;;
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
        --concat-only)
            concat_only=true; shift ;;
        --measure)
            [[ -z "$2" || "$2" == --* ]] && { echo "Error: $1 requires a value"; exit 1; }
            measure_file="$2"; shift 2 ;;
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
if [[ "${single_run}" == true ]]; then
    base_name="${basename_arg##*/}"
else
    base_name="${basename_arg}_${replica}"
fi
output_dir="${output_dir:-${base_name}_analysis}"

# Default reproducibility script name, unless overridden with --save-script
save_script="${save_script:-${output_dir}/gmx_analysis_commands.sh}"

# --concat-only stops before the analyses, so options selecting analyses (or skipping the
# concatenation itself) would silently do nothing: reject the combination instead
if [[ "${concat_only}" == true ]]; then
    conflict=""
    if [[ -n "${only_list}" ]]; then conflict="--only"; fi
    if [[ -n "${skip_list}" ]]; then conflict="--skip"; fi
    if [[ "${skip_concat}" == true ]]; then conflict="--skip-concat"; fi
    if [[ -n "${conflict}" ]]; then
        echo "Error: --concat-only cannot be combined with ${conflict}"
        exit 1
    fi
fi

# A single run has nothing to concatenate
if [[ "${single_run}" == true && "${concat_only}" == true ]]; then
    echo "Error: --concat-only only applies to segmented runs (a single run is used as is)"
    exit 1
fi
if [[ "${single_run}" == true && "${skip_concat}" == true ]]; then
    echo "Error: --skip-concat only applies to segmented runs (a single run is used as is)"
    exit 1
fi

if [[ ! -d "${segments_dir}" ]]; then
    echo "Error: Segments directory '${segments_dir}' not found"
    exit 1
fi
if [[ -n "${index_file}" ]] && [[ ! -f "${index_file}" ]]; then
    echo "Error: Index file '${index_file}' not found"
    exit 1
fi

# Validate and normalize --skip/--only categories (one per line via mapfile)
valid_categories=(energy rmsd gyration sasa secondary-structure rmsf distance angle)
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
# Use a command substitution (not a process substitution, whose exit status is lost)
# so that an unknown category aborts the script
skip_categories=()
only_categories=()
if [[ -n "${skip_list}" ]]; then
    normalized="$(normalize_categories "${skip_list}")" || exit 1
    if [[ -n "${normalized}" ]]; then
        mapfile -t skip_categories <<< "${normalized}"
    fi
fi
if [[ -n "${only_list}" ]]; then
    normalized="$(normalize_categories "${only_list}")" || exit 1
    if [[ -n "${normalized}" ]]; then
        mapfile -t only_categories <<< "${normalized}"
    fi
fi

# Read the --measure file into measure_types/_names/_sels (selections of one measurement
# joined by $'\x1f'). Strict: any malformed line stops the script with its line number.
measure_types=(); measure_names=(); measure_sels=()
if [[ -n "${measure_file}" ]]; then
    if [[ ! -f "${measure_file}" ]]; then
        echo "Error: Measurements file '${measure_file}' not found"
        exit 1
    fi
    measure_error() { echo "Error: ${measure_file}, line ${lineno}: $*" >&2; exit 1; }
    lineno=0
    while IFS= read -r line || [[ -n "${line}" ]]; do
        lineno=$((lineno + 1))
        line="${line%$'\r'}"
        [[ "${line}" =~ ^[[:space:]]*(#|$) ]] && continue
        read -r m_type m_name m_rest <<< "${line}"
        case "${m_type}" in
            distance) m_nsel=2 ;;
            angle)    m_nsel=3 ;;
            dihedral) m_nsel=4 ;;
            *) measure_error "unknown type '${m_type}' (use distance, angle or dihedral)" ;;
        esac
        [[ "${m_name}" =~ ^[A-Za-z0-9_.-]+$ ]] \
            || measure_error "invalid name '${m_name}' (use letters, digits, _ . - only)"
        for n in "${measure_names[@]}"; do
            [[ "${n}" == "${m_name}" ]] && measure_error "name '${m_name}' is used twice"
        done
        # Empty parts ("a | | b", or a leading/trailing |) would be dropped silently by read below
        [[ "|${m_rest}|" =~ \|[[:space:]]*\| ]] \
            && measure_error "empty selection (check the | separators)"
        IFS='|' read -ra m_parts <<< "${m_rest}"
        m_joined=""
        for i in "${!m_parts[@]}"; do
            m_part="${m_parts[i]#"${m_parts[i]%%[![:space:]]*}"}"   # trim leading spaces
            m_part="${m_part%"${m_part##*[![:space:]]}"}"           # trim trailing spaces
            [[ -n "${m_part}" ]] || measure_error "empty selection (check the | separators)"
            m_joined+="${m_joined:+$'\x1f'}${m_part}"
        done
        [[ ${#m_parts[@]} -eq ${m_nsel} ]] \
            || measure_error "$( [[ ${m_type} == angle ]] && echo an || echo a ) ${m_type} needs ${m_nsel} selections separated by |, found ${#m_parts[@]}"
        measure_types+=("${m_type}"); measure_names+=("${m_name}"); measure_sels+=("${m_joined}")
    done < "${measure_file}"
    [[ ${#measure_types[@]} -gt 0 ]] || { echo "Error: No measurements found in '${measure_file}'"; exit 1; }
    unset -f measure_error
fi
for c in "${only_categories[@]}"; do
    if [[ ( "${c}" == distance || "${c}" == angle ) && -z "${measure_file}" ]]; then
        echo "Error: --only ${c} needs --measure <file> with the measurements to make"
        exit 1
    fi
done

should_run() {
    local name="$1"
    if [[ ${#only_categories[@]} -gt 0 ]]; then
        for c in "${only_categories[@]}"; do [[ "${c}" == "${name}" ]] && return 0; done
        return 1
    fi
    for c in "${skip_categories[@]}"; do [[ "${c}" == "${name}" ]] && return 1; done
    return 0
}

if [[ "${single_run}" == true ]]; then
    # Single run: <basename>.xtc (+ .tpr, optional .edr), in --segments-dir unless a path is given
    if [[ "${basename_arg}" == */* ]]; then
        single_prefix="${basename_arg}"
else
        single_prefix="${segments_dir}/${basename_arg}"
    fi
    single_xtc="${single_prefix}.xtc"
    single_edr="${single_prefix}.edr"
    if [[ ! -f "${single_xtc}" ]]; then
        echo "Error: Trajectory '${single_xtc}' not found"
        echo "(For a segmented run from gmx_continue_grompp.sh, also give the replica: $0 ${basename_arg} <replica>)"
        exit 1
    fi
    [[ -f "${single_edr}" ]] || single_edr=""
    default_structure="${single_prefix}.tpr"
else
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

    default_structure="${seg_tpr[0]}"
fi

structure_file="${structure_file:-${default_structure}}"
if [[ ! -f "${structure_file}" ]]; then
    echo "Error: Reference structure '${structure_file}' not found"
    [[ "${single_run}" == true ]] && echo "Give one with --structure (a .tpr is best)"
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
convert_from_ps "0" "${plot_time_unit}" > /dev/null  # validate --plot-tu early
case "${energy_time_unit}" in
    ps|ns|us|fs) ;;
    *)
        echo "Error: Unsupported --energy-tu unit '${energy_time_unit}' (use ps, ns, us, or fs)" >&2
        exit 1 ;;
esac
# The per-tool -tu/-b/-e/-dt flags are built by set_time_flags (below), since not
# every gmx tool accepts all of them

mkdir -p "${output_dir}"
log_file="${output_dir}/${base_name}_analysis.log"
{
    echo "GROMACS analysis log for ${base_name}"
    echo "Generated: $(date)"
    echo "=========================================="
} > "${log_file}"

# If <file> exists, rename it to the next numbered backup before it is overwritten:
# name.sh -> name.1.sh, name.2.sh, ... (highest number = most recent), so the commands
# of earlier runs are kept. Files not ending in .sh become name.1, name.2, ...
backup_existing_file() {
    local file="$1"
    [[ -e "${file}" ]] || return 0
    local base="${file##*/}"
    local prefix="${file%"${base}"}"
    local stem="${base}" ext=""
    if [[ "${base}" == *.sh ]]; then
        stem="${base%.sh}"; ext=".sh"
    fi
    local max=0 f n
    for f in "${prefix}${stem}".*"${ext}"; do
        n="${f#"${prefix}${stem}."}"; n="${n%"${ext}"}"
        if [[ "${n}" =~ ^[0-9]+$ ]] && (( 10#${n} > max )); then
            max=$(( 10#${n} ))
        fi
    done
    local backup="${prefix}${stem}.$(( max + 1 ))${ext}"
    mv "${file}" "${backup}"
    echo "Previous commands file kept as: ${backup}"
}

if [[ -n "${save_script}" ]]; then
    backup_existing_file "${save_script}"
    {
        echo "#!/bin/bash -e"
        echo "# Commands executed by gmx_analysis.sh"
        echo "# Invocation: ${original_invocation}"
        echo "# Generated: $(date)"
    } > "${save_script}"
    chmod +x "${save_script}"
fi

echo "=== GROMACS Trajectory Analysis ==="
if [[ "${single_run}" == true ]]; then
    echo "Run:             ${basename_arg}  (single trajectory)"
    echo "Trajectory:      ${single_xtc}"
    echo "Energy file:     ${single_edr:-(none: energy analyses skipped)}"
else
    echo "Basename:        ${basename_arg}  (replica ${replica})"
    echo "Segments dir:    ${segments_dir}"
    echo "Segments found:  ${#seg_xtc[@]}"
    for i in "${!seg_xtc[@]}"; do
        echo "  - $(basename "${seg_xtc[$i]}")"
    done
fi
echo "Structure:       ${structure_file}"
echo "Group:           ${group} (fit group for RMSD: ${fit_group})"
[[ -n "${index_file}" ]] && echo "Index:           ${index_file}"
[[ -n "${begin_ps}" ]] && echo "Begin:           ${begin_ps} ps"
[[ -n "${end_ps}" ]]   && echo "End:             ${end_ps} ps"
[[ -n "${dt_ps}" ]]    && echo "dt:              ${dt_ps} ps"
echo "Plot time unit:  ${plot_time_unit}"
if [[ "${concat_only}" == true ]]; then
    echo "Analyses:        none (--concat-only)"
else
    echo "Energy terms:    ${energy_terms}"
    [[ -n "${measure_file}" ]] && echo "Measurements:    ${#measure_types[@]} from ${measure_file}"
fi
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
# With --capture <file>, the tool's output is also saved to <file> (overwritten), e.g. to
# parse the summary table gmx energy prints; it is still appended to log_file as well.
run_piped_cmd() {
    local capture_file=""
    if [[ "$1" == "--capture" ]]; then
        capture_file="$2"; shift 2
    fi
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
            [[ -n "${capture_file}" ]] && printf '> %q 2>&1' "${capture_file}"
            echo ""
        } >> "${save_script}"
    fi
    if [[ "${dry_run}" == true ]]; then
        echo "(dry run: not executed)"
        return 0
    fi
    local status=0
    if [[ -n "${capture_file}" ]]; then
        printf '%s\n' "${stdin_data}" | "$@" > "${capture_file}" 2>&1 || status=$?
        cat "${capture_file}" >> "${log_file}"
    else
        printf '%s\n' "${stdin_data}" | "$@" >> "${log_file}" 2>&1 || status=$?
    fi
    if [[ ${status} -eq 0 ]]; then
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

# Returns success if gmx <tool> accepts option <opt>, parsed from its help text
# (looked up once per tool). Option sets differ between tools and GROMACS versions,
# e.g. gmx energy has no -tu/-dt. If the help can't be read (e.g. --dry-run without
# gmx available), every option is assumed to be supported.
declare -A tool_options_cache=()
tool_supports() {
    local tool="$1" opt="$2"
    if [[ -z "${tool_options_cache[${tool}]+set}" ]]; then
        tool_options_cache[${tool}]="$( { "${gmx_bin}" "${tool}" -h 2>/dev/null || true; } \
            | grep -oE '^ +-[A-Za-z0-9_]+' | tr -d ' ' | tr '\n' ' ' || true)"
    fi
    [[ -z "${tool_options_cache[${tool}]}" ]] && return 0
    [[ " ${tool_options_cache[${tool}]} " == *" ${opt} "* ]]
}

# Sets the time_flags array for gmx <tool> from --begin/--end/--dt/--plot-tu,
# passing only options the tool supports. Tools with -tu get the --plot-tu unit,
# and -b/-e/-dt in that unit (gmx reads them in the -tu unit); tools without -tu
# (e.g. gmx energy, gmx rmsf) read -b/-e/-dt in ps.
# Pass "range-only" for tools whose time axis is not set with -tu: gmx rmsf (no
# time axis) and gmx energy (rescaled afterwards, see --energy-tu). These always
# get -b/-e/-dt in ps and never -tu.
set_time_flags() {
    local tool="$1" mode="$2" unit="ps"
    time_flags=()
    if [[ "${mode}" != "range-only" ]] && tool_supports "${tool}" -tu; then
        unit="${plot_time_unit}"
        time_flags+=(-tu "${unit}")
    elif [[ "${mode}" != "range-only" && "${plot_time_unit}" != "ps" ]]; then
        echo "Note: gmx ${tool} has no -tu option; its time axis stays in ps"
    fi
    local opt value_ps
    for opt in -b -e -dt; do
        case "${opt}" in
            -b)  value_ps="${begin_ps}" ;;
            -e)  value_ps="${end_ps}" ;;
            -dt) value_ps="${dt_ps}" ;;
        esac
        [[ -z "${value_ps}" ]] && continue
        if tool_supports "${tool}" "${opt}"; then
            time_flags+=("${opt}" "$(convert_from_ps "${value_ps}" "${unit}")")
        else
            echo "Note: gmx ${tool} has no ${opt} option; it is not applied to this analysis"
        fi
    done
}

# Post-processes an .xvg file written by gmx energy, which always titles its plots
# "GROMACS Energies" and has no -tu option:
#   - sets the title to the energy term name, taken from the file's legend
#     ('@ s0 legend "Total Energy"' -> '@ title "Total Energy"')
#   - if <unit> is not ps, rescales the time column (first column) to <unit> and
#     updates the x-axis label to match; otherwise data lines are left untouched
#   - adds the statistics gmx energy printed for the term (read from <stats_file>, its
#     captured output) as a header comment line, e.g.
#       # gmx_analysis stats: Average=299.99 Err.Est.=0.098 RMSD=3.53537 Tot-Drift=0.6368 unit=K
#     which plot_xvg.py --stats shows in a text box; <stats_file> is removed afterwards
# Runs through run_cmd, so it is logged to --save-script and honors --dry-run.
postprocess_energy_xvg() {
    local xvg="$1" unit="$2" stats_file="$3" factor desc
    case "${unit}" in
        ps) factor="1" ;;
        ns) factor="0.001" ;;
        us) factor="0.000001" ;;
        fs) factor="1000" ;;
    esac
    # Reads three inputs, told apart by the role=... assignments on the command line:
    # the gmx energy output (stats table), then the .xvg twice (the first pass only reads
    # the legend, which follows the title), writing the result on the second pass.
    local prog='function is_num(v) { return v ~ /^[-+]?([0-9]+[.]?[0-9]*|[.][0-9]+)([eE][-+]?[0-9]+)?$/ }
    role == "stats" {
        # Summary table: a header line, a line of dashes, then the row
        #   "<name>  Average  Err.Est.  RMSD  Tot-Drift  (unit)"
        # where <name> may contain spaces, so fields are counted from the end. The output
        # also holds other text (e.g. the GROMACS quote, other dashed lines), so only the
        # first line after that exact header + dashes is used, and only if its values are
        # numbers (Err.Est. is "--" when there is too little data to estimate it).
        if (stats != "" || done_table) next
        if (!in_header && $0 ~ /^Energy[ \t]+Average[ \t]+Err[.]Est[.][ \t]+RMSD[ \t]+Tot-Drift/) {
            in_header = 1; next
        }
        if (in_header && !in_table) {
            if ($0 ~ /^-+[ \t]*$/) in_table = 1; else in_header = 0
            next
        }
        if (in_table) {
            done_table = 1
            if (NF >= 6 && $NF ~ /^\(.*\)$/ && is_num($(NF-4)) && is_num($(NF-2)) \
                    && is_num($(NF-1)) && (is_num($(NF-3)) || $(NF-3) == "--")) {
                u = $NF; gsub(/[()]/, "", u)
                stats = "# gmx_analysis stats: Average=" $(NF-4) " Err.Est.=" $(NF-3) \
                        " RMSD=" $(NF-2) " Tot-Drift=" $(NF-1) " unit=" u
            }
        }
        next
    }
    role == "scan" {
        if (name == "" && $0 ~ /^@ +s0 +legend/) {
            name = $0; sub(/^@ +s0 +legend +/, "", name); gsub(/"/, "", name)
        }
        next
    }
    # role == "out": stats line goes at the end of the leading # comment block
    !stats_done && !/^#/ {
        if (stats != "") print stats > out
        stats_done = 1
    }
    /^[@#]/ || NF == 0 {
        if (name != "" && $0 ~ /^@ +title /) $0 = "@    title \"" name "\""
        if (factor != 1 && $0 ~ /^@ *xaxis +label/) sub(/\(ps\)/, "(" unit ")")
        print > out; next
    }
    { if (factor != 1) $1 = sprintf("%.6f", $1 * factor); print > out }'
    desc="Set title of $(basename "${xvg}") to the energy term name, add its statistics"
    [[ "${unit}" != "ps" ]] && desc+=", convert its time axis to ${unit}"
    run_cmd "${desc}" \
        awk -v factor="${factor}" -v unit="${unit}" -v out="${xvg}.tmp" "${prog}" \
            role=stats "${stats_file}" role=scan "${xvg}" role=out "${xvg}" &&
    run_cmd "Replace $(basename "${xvg}") with the processed file" mv "${xvg}.tmp" "${xvg}" &&
    run_cmd "Remove temporary gmx energy output" rm -f "${stats_file}"
}

# Sets the '@ title' of an .xvg file (e.g. a measurement's name instead of gmx's generic
# "Distance"). Runs through run_cmd, so it is logged to --save-script and honors --dry-run.
set_xvg_title() {
    local xvg="$1" title="$2"
    run_cmd "Set title of $(basename "${xvg}") to '${title}'" \
        awk -v title="${title}" -v out="${xvg}.tmp" \
            '/^@ +title / { print "@    title \"" title "\"" > out; next } { print > out }' "${xvg}" &&
    run_cmd "Replace $(basename "${xvg}") with the processed file" mv "${xvg}.tmp" "${xvg}"
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

# --- Check the --measure selections before doing any work ------------------------
# Each selection must give exactly one position (one atom, or "com of"/"cog of" a group).
# One gmx select call per measurement, on the structure file (+ index): its -os output has
# one column per selection with the number of positions it matches.
if [[ ${#measure_types[@]} -gt 0 ]] && [[ "${concat_only}" == false ]] && { should_run distance || should_run angle; }; then
    if ! command -v "${gmx_bin%% *}" &> /dev/null; then
        echo ""
        echo "Note: gmx not available, measurement selections not checked (dry run)"
    else
        echo ""
        echo "--- Checking ${#measure_types[@]} measurement(s) in ${measure_file} ---"
        size_xvg="$(mktemp --suffix=.xvg)"
        select_err="$(mktemp)"
        failed=0
        for i in "${!measure_types[@]}"; do
            IFS=$'\x1f' read -ra m_parts <<< "${measure_sels[i]}"
            rm -f "${size_xvg}"
            "${gmx_bin}" select -s "${structure_file}" "${index_flag[@]}" -select "${m_parts[@]}" \
                -os "${size_xvg}" -xvg none -quiet > /dev/null 2> "${select_err}" || true
            # First data line: time, then one count per selection
            counts=()
            [[ -s "${size_xvg}" ]] && read -ra counts <<< "$(grep -v '^[@#]' "${size_xvg}" | head -n 1)"
            if [[ ${#counts[@]} -ne $(( ${#m_parts[@]} + 1 )) ]]; then
                echo "✗ ${measure_types[i]} ${measure_names[i]}: gmx select rejected a selection:"
                # GROMACS' explanation sits between a heading ("Error in user input",
                # "Inconsistency in user input", "Fatal error") and "For more information"
                reason="$(sed -n '/in user input:\|Fatal error/,/For more information/p' "${select_err}" \
                    | grep -vE '^[[:space:]]*$|in user input:|Fatal error|For more information|Invalid command-line options|In command-line option')"
                if [[ -z "${reason}" ]]; then
                    reason="$(grep -vE '^[[:space:]]*$|^-+$|^(Program|Source file|Function|website):' "${select_err}" | tail -n 5)"
                fi
                echo "${reason}" | sed 's/^[[:space:]]*/      /'
                failed=1
                continue
            fi
            for j in "${!m_parts[@]}"; do
                count="${counts[j + 1]}"
                if ! [[ "${count}" =~ ^[0-9]+(\.0+)?$ ]]; then
                    echo "✗ ${measure_types[i]} ${measure_names[i]}: unexpected gmx select output '${count}'"
                    failed=1
                elif [[ "${count%%.*}" -ne 1 ]]; then
                    echo "✗ ${measure_types[i]} ${measure_names[i]}: '${m_parts[j]}' matches ${count%%.*} atom(s); it needs exactly 1 (one atom, or \"com of\"/\"cog of\" a group)"
                    failed=1
                fi
            done
        done
        rm -f "${size_xvg}" "${select_err}"
        if [[ "${failed}" -ne 0 ]]; then
            echo "Error: Fix the measurements in ${measure_file} (nothing has been run)"
            exit 1
        fi
        echo "✓ All selections valid"
    fi
fi

# --- 1. Concatenate segments (segmented runs) ----------------------------------
# traj_xtc/traj_edr: the trajectory and energy file the analyses read
if [[ "${single_run}" == true ]]; then
        traj_xtc="${single_xtc}"
        traj_edr="${single_edr}"
else
    concat_xtc="${output_dir}/${base_name}_concat.xtc"
    concat_edr="${output_dir}/${base_name}_concat.edr"
    traj_xtc="${concat_xtc}"
    traj_edr="${concat_edr}"

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
fi

if [[ "${concat_only}" == true ]]; then
    echo ""
    echo "=== Concatenation complete (--concat-only) ==="
    if [[ "${dry_run}" == true ]]; then
        echo "Dry run: no files were written."
    else
        echo "Concatenated files:"
        echo "  - ${concat_xtc}"
        echo "  - ${concat_edr}"
        echo "Analyze them later with --skip-concat. Full command log: ${log_file}"
    fi
    exit 0
fi

# --- 2. Energy terms (temperature, pressure, potential, total energy, ...) --
if should_run energy && [[ -z "${traj_edr}" ]]; then
    echo ""
    echo "--- Energy terms ---"
    echo "No energy file (${single_prefix}.edr): energy analyses skipped."
elif should_run energy; then
    set_time_flags energy range-only
    IFS=',' read -ra terms <<< "${energy_terms}"
    for term in "${terms[@]}"; do
        term="$(echo "${term}" | xargs)"
        [[ -z "${term}" ]] && continue
        slug="$(echo "${term}" | tr '[:upper:]' '[:lower:]' | tr '-' '_')"
        out="${output_dir}/${slug}.xvg"
        cmd=("${gmx_bin}" energy -f "${traj_edr}" -s "${structure_file}" -o "${out}" -xvg "${xvg_format}" "${time_flags[@]}")
        # gmx energy's output is captured to a temporary file, to read its statistics table
        stats_file="${out}.energy.txt"
        if run_piped_cmd --capture "${stats_file}" "Energy term: ${term}" "${term}" "${cmd[@]}"; then
            postprocess_energy_xvg "${out}" "${energy_time_unit}" "${stats_file}" || true
        else
            rm -f "${stats_file}"
        fi
    done
fi

# --- 3. RMSD to reference structure -----------------------------------------
if should_run rmsd; then
    out="${output_dir}/rmsd.xvg"
    set_time_flags rms
    cmd=("${gmx_bin}" rms -s "${structure_file}" -f "${traj_xtc}" -o "${out}" -xvg "${xvg_format}" "${index_flag[@]}" "${time_flags[@]}")
    stdin_data="${fit_group}"$'\n'"${group}"
    run_piped_cmd "RMSD (fit=${fit_group}, calc=${group})" "${stdin_data}" "${cmd[@]}" || true
fi

# --- 4. Radius of gyration ---------------------------------------------------
if should_run gyration; then
    out="${output_dir}/gyration_radius.xvg"
    set_time_flags gyrate
    cmd=("${gmx_bin}" gyrate -s "${structure_file}" -f "${traj_xtc}" -o "${out}" -sel "${group}" -xvg "${xvg_format}" "${index_flag[@]}" "${time_flags[@]}")
    run_cmd "Radius of gyration (${group})" "${cmd[@]}" || true
fi

# --- 5. Solvent accessible surface area --------------------------------------
if should_run sasa; then
    out="${output_dir}/sasa.xvg"
    set_time_flags sasa
    cmd=("${gmx_bin}" sasa -s "${structure_file}" -f "${traj_xtc}" -o "${out}" -surface "${group}" -probe "${probe_radius}" -xvg "${xvg_format}" "${index_flag[@]}" "${time_flags[@]}")
    run_cmd "Solvent accessible surface area (${group})" "${cmd[@]}" || true
fi

# --- 6. Secondary structure over time (gmx dssp, GROMACS >= 2023) ----------
if should_run secondary-structure; then
    if [[ "${dry_run}" == true ]] || has_gmx_dssp; then
        out_dat="${output_dir}/secondary_structure.dat"
        out_num="${output_dir}/secondary_structure_counts.xvg"
        set_time_flags dssp
        cmd=("${gmx_bin}" dssp -s "${structure_file}" -f "${traj_xtc}" -sel "${group}" \
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
    # rmsf's x-axis is residue/atom index, not time; -b/-e/-dt only select the frames to average over
    set_time_flags rmsf range-only
    cmd=("${gmx_bin}" rmsf -s "${structure_file}" -f "${traj_xtc}" -o "${out}" -xvg "${xvg_format}" "${index_flag[@]}" "${time_flags[@]}")
    [[ "${per_residue}" == true ]] && cmd+=(-res)
    run_piped_cmd "RMS fluctuations (${group})" "${group}" "${cmd[@]}" || true
fi

# --- 8. Distances, angles and dihedrals (--measure) ------------------------------
for i in "${!measure_types[@]}"; do
    m_type="${measure_types[i]}"; m_name="${measure_names[i]}"
    if [[ "${m_type}" == distance ]]; then
        should_run distance || continue
    else
        should_run angle || continue
    fi
    IFS=$'\x1f' read -ra m_parts <<< "${measure_sels[i]}"
    m_selection="${m_parts[0]}"
    for m_part in "${m_parts[@]:1}"; do m_selection+=" plus ${m_part}"; done
    out="${output_dir}/${m_type}_${m_name}.xvg"
    if [[ "${m_type}" == distance ]]; then
        set_time_flags distance
        cmd=("${gmx_bin}" distance -s "${structure_file}" -f "${traj_xtc}" -select "${m_selection}" \
             -oall "${out}" -xvg "${xvg_format}" "${index_flag[@]}" "${time_flags[@]}")
    else
        set_time_flags gangle
        cmd=("${gmx_bin}" gangle -s "${structure_file}" -f "${traj_xtc}" -g1 "${m_type}" -group1 "${m_selection}" \
             -oall "${out}" -xvg "${xvg_format}" "${index_flag[@]}" "${time_flags[@]}")
    fi
    if run_cmd "${m_type^} ${m_name}" "${cmd[@]}"; then
        set_xvg_title "${out}" "${m_name} (${m_type})" || true
    fi
done

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

#!/usr/bin/env python3
"""
bfactors_vs_rmsf.py - compare MD fluctuations (RMSF) with crystallographic B-factors.

Reads a per-residue C-alpha RMSF file written by gmx rmsf (-res, C-alpha group) and the
C-alpha B-factors of a crystal structure (PDB file, or a plain list "residue B"), converts the
B-factors to RMSF units with  RMSF = sqrt(3 B / (8 pi^2))  (B in A^2, RMSF in A, then nm), pairs
residues by number, and writes an .xvg file with both curves and their Pearson correlation.

Used by gmx_analysis.sh --bfactors, and usable on its own:
  bfactors_vs_rmsf.py rmsf_calpha.xvg crystal.pdb -o rmsf_vs_bfactors.xvg [--scale mean]
  bfactors_vs_rmsf.py --check crystal.pdb          (only read and check the B-factor file)
"""
import argparse
import math
import sys
from datetime import datetime

SCALES = ("none", "mean", "zscore")
MIN_MATCHED = 3


class InputError(Exception):
    """A problem with an input file, reported without a traceback."""


def b_to_rmsf_nm(b):
    """Isotropic B-factor (A^2) -> RMSF (nm):  B = 8 pi^2/3 <dr^2>."""
    return math.sqrt(3.0 * b / (8.0 * math.pi ** 2)) / 10.0


def read_rmsf_xvg(path):
    """Per-residue RMSF from gmx rmsf: {residue number: RMSF in nm}."""
    rmsf = {}
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            if line.startswith(("#", "@")) or not line.strip():
                continue
            fields = line.split()
            try:
                if len(fields) != 2:
                    raise ValueError
                res, value = float(fields[0]), float(fields[1])
            except ValueError:
                raise InputError(f"{path}, line {lineno}: expected 'residue RMSF', got: {line.strip()}")
            if res != int(res):
                raise InputError(f"{path}, line {lineno}: residue number {fields[0]} is not an integer "
                                 f"(per-residue RMSF from gmx rmsf -res expected)")
            rmsf[int(res)] = value
    if not rmsf:
        raise InputError(f"{path}: no RMSF values found")
    return rmsf


def read_pdb_bfactors(path):
    """C-alpha B-factors from a PDB file: {residue number: B in A^2}.

    One chain only; for alternative conformations the highest occupancy is used. Several
    MODELs (e.g. NMR ensembles, whose B-factors are meaningless) and insertion codes are errors.
    """
    best = {}       # residue -> (occupancy, B)
    chains = set()
    models = 0
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            if line.startswith("MODEL"):
                models += 1
                continue
            if not line.startswith(("ATOM  ", "HETATM")) or line[12:16].strip() != "CA":
                continue
            where = f"{path}, line {lineno}"
            if len(line.rstrip("\n")) < 66:
                raise InputError(f"{where}: line too short to hold a B-factor (columns 61-66)")
            if line[26] != " ":
                raise InputError(f"{where}: insertion code '{line[26]}' is not supported "
                                 f"(renumber the residues)")
            try:
                res = int(line[22:26])
                occ = float(line[54:60]) if line[54:60].strip() else 1.0
                b = float(line[60:66])
            except ValueError:
                raise InputError(f"{where}: cannot read residue number, occupancy or B-factor")
            if b < 0:
                raise InputError(f"{where}: negative B-factor {b}")
            chains.add(line[21])
            if res not in best or occ > best[res][0]:
                best[res] = (occ, b)
    if models > 1:
        raise InputError(f"{path}: {models} MODELs (an NMR ensemble?); use a single crystal structure")
    if len(chains) > 1:
        raise InputError(f"{path}: C-alpha atoms in several chains ({', '.join(sorted(chains))}); "
                         f"keep only one chain")
    if not best:
        raise InputError(f"{path}: no C-alpha atoms (ATOM/HETATM records named CA) found")
    return {res: b for res, (occ, b) in best.items()}


def read_list_bfactors(path):
    """B-factors from a plain list: 'residue B' per line, '#' starts a comment line."""
    bfac = {}
    with open(path) as f:
        for lineno, line in enumerate(f, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            try:
                if len(fields) != 2:
                    raise ValueError
                res, b = int(fields[0]), float(fields[1])
            except ValueError:
                raise InputError(f"{path}, line {lineno}: expected 'residue B-factor', got: {line.strip()}")
            if b < 0:
                raise InputError(f"{path}, line {lineno}: negative B-factor {b}")
            if res in bfac:
                raise InputError(f"{path}, line {lineno}: residue {res} listed twice")
            bfac[res] = b
    if not bfac:
        raise InputError(f"{path}: no B-factors found")
    return bfac


def read_bfactors(path):
    """B-factors from a PDB file (recognized by ATOM/HETATM records) or a plain list."""
    with open(path) as f:
        is_pdb = any(line.startswith(("ATOM  ", "HETATM")) for line in f)
    return read_pdb_bfactors(path) if is_pdb else read_list_bfactors(path)


def floor_warning(bfac):
    """Warning text if several residues share exactly the lowest B-factor, which usually means
    a lower limit applied during refinement (not real motion), else None."""
    low = min(bfac.values())
    at_floor = [r for r, b in bfac.items() if b == low]
    if len(at_floor) < 3:
        return None
    return (f"Warning: {len(at_floor)} residues have exactly the lowest B-factor ({low:.2f} A^2): "
            f"{ranges(at_floor)}. This is usually a lower limit set during refinement, not real "
            f"motion; those residues flatten the crystal curve.")


def ranges(numbers):
    """Compact text for a set of residue numbers: '1-5, 9, 12-14'."""
    numbers = sorted(numbers)
    parts, start = [], None
    for i, n in enumerate(numbers):
        if start is None:
            start = n
        if i == len(numbers) - 1 or numbers[i + 1] != n + 1:
            parts.append(str(start) if start == n else f"{start}-{n}")
            start = None
    return ", ".join(parts)


def mean(values):
    return sum(values) / len(values)


def stdev(values):
    m = mean(values)
    return math.sqrt(sum((v - m) ** 2 for v in values) / len(values))


def pearson(x, y):
    mx, my = mean(x), mean(y)
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    if sxx == 0 or syy == 0:
        return float("nan")
    return sxy / math.sqrt(sxx * syy)


def compare(rmsf_path, bfac_path, out_path, scale):
    md = read_rmsf_xvg(rmsf_path)
    bfac = read_bfactors(bfac_path)

    matched = sorted(set(md) & set(bfac))
    only_md, only_crystal = set(md) - set(bfac), set(bfac) - set(md)
    if len(matched) < MIN_MATCHED:
        raise InputError(
            f"only {len(matched)} residue numbers in common (MD: {ranges(md)}; crystal: "
            f"{ranges(bfac)}); the residue numbering of the two must match")
    if only_md:
        print(f"Note: {len(only_md)} residue(s) without a crystal B-factor, left out: {ranges(only_md)}")
    if only_crystal:
        print(f"Note: {len(only_crystal)} crystal residue(s) not in the MD, left out: {ranges(only_crystal)}")
    if len(matched) < 0.5 * min(len(md), len(bfac)):
        print(f"Warning: only {len(matched)} of {len(md)} MD residues matched; is the residue "
              f"numbering of the crystal structure the same as in the topology?")

    warning = floor_warning({r: bfac[r] for r in matched})
    if warning:
        print(warning)

    y_md = [md[r] for r in matched]
    y_cr = [b_to_rmsf_nm(bfac[r]) for r in matched]
    r = pearson(y_md, y_cr)
    if scale == "mean":
        factor = mean(y_md) / mean(y_cr)
        y_cr = [v * factor for v in y_cr]
        ylabel, cr_label = "RMSF (nm)", f"Crystal B-factors (scaled x{factor:.2f} to the MD mean)"
    elif scale == "zscore":
        y_md = [(v - mean(y_md)) / stdev(y_md) for v in y_md]
        y_cr = [(v - mean(y_cr)) / stdev(y_cr) for v in y_cr]
        ylabel, cr_label = "RMSF (z-score)", "Crystal B-factors (z-score)"
    else:
        ylabel, cr_label = "RMSF (nm)", "Crystal B-factors"

    with open(out_path, "w") as f:
        f.write(f"# Created by bfactors_vs_rmsf.py on {datetime.now():%Y-%m-%d %H:%M:%S}\n")
        f.write(f"# MD RMSF (C-alpha): {rmsf_path}\n")
        f.write(f"# Crystal B-factors: {bfac_path}, converted with RMSF = sqrt(3 B / (8 pi^2))\n")
        f.write(f"# Scale: {scale}; Pearson r = {r:.4f} over {len(matched)} residues\n")
        f.write('@    title "RMSF: MD vs crystal B-factors"\n')
        f.write(f'@    subtitle "Pearson r = {r:.2f} ({len(matched)} residues, scale: {scale})"\n')
        f.write('@    xaxis  label "Residue"\n')
        f.write(f'@    yaxis  label "{ylabel}"\n')
        f.write("@TYPE xy\n")
        f.write('@ s0 legend "MD (C-alpha)"\n')
        f.write(f'@ s1 legend "{cr_label}"\n')
        for res, a, b in zip(matched, y_md, y_cr):
            f.write(f"{res:6d} {a:10.5f} {b:10.5f}\n")
    print(f"Pearson r = {r:.3f} over {len(matched)} residues (scale: {scale}) -> {out_path}")


def main():
    ap = argparse.ArgumentParser(
        description="Compare MD C-alpha RMSF (gmx rmsf -res) with crystal B-factors.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s rmsf_calpha.xvg 2ckd.pdb -o rmsf_vs_bfactors.xvg
  %(prog)s rmsf_calpha.xvg bfactors.txt -o cmp.xvg --scale zscore
  %(prog)s --check 2ckd.pdb

The B-factor file is a PDB file (C-alpha atoms of one chain; for alternative conformations the
highest occupancy is used) or a plain list with one "residue B-factor" pair per line.
Residues are paired by number, so the crystal numbering must match the topology.
""")
    ap.add_argument("rmsf", nargs="?", help="per-residue C-alpha RMSF .xvg from gmx rmsf -res")
    ap.add_argument("bfactors", nargs="?", help="PDB file, or list of 'residue B-factor' lines")
    ap.add_argument("-o", "--output", help="output .xvg file")
    ap.add_argument("--scale", choices=SCALES, default="mean",
                    help="none: crystal B converted to RMSF (nm); mean: crystal curve scaled to "
                         "the MD mean (default); zscore: both as z-scores")
    ap.add_argument("--check", metavar="FILE",
                    help="only read and check a B-factor file, then exit")
    args = ap.parse_args()

    try:
        if args.check:
            bfac = read_bfactors(args.check)
            print(f"{args.check}: {len(bfac)} C-alpha B-factors, residues {ranges(bfac)}")
            warning = floor_warning(bfac)
            if warning:
                print(warning)
            return 0
        if not (args.rmsf and args.bfactors and args.output):
            ap.error("give the RMSF file, the B-factor file and -o OUTPUT (or --check FILE)")
        compare(args.rmsf, args.bfactors, args.output, args.scale)
    except InputError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except OSError as e:
        print(f"Error: {e.filename}: {e.strerror}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())

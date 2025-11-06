#!/usr/bin/env python3
"""
GROMACS PCA automation script

This script automates a typical Principal Component Analysis (PCA) workflow
for MD trajectories using GROMACS tools and produces Matplotlib plots.

Pipeline (default):
1) Build an index using a GROMACS selection (default: Backbone)
2) Fit and center the trajectory on the selected atoms (gmx trjconv)
3) Compute covariance matrix and eigenvectors (gmx covar)
4) Analyze eigenvectors and project trajectory (gmx anaeig)
5) Parse .xvg outputs and generate plots

Requirements:
- GROMACS in PATH (gmx) or provide --gmx-bin
- Python packages: numpy, matplotlib

Tested with GROMACS 2021+ selection syntax.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np

# Matplotlib is optional at import to allow --no-plots on headless systems
try:
    import matplotlib
    matplotlib.use("Agg")  # headless-safe backend
    import matplotlib.pyplot as plt
except Exception as e:
    plt = None  # defer error until we actually try to plot


class CommandError(RuntimeError):
    pass


def run_cmd(cmd: List[str], stdin: Optional[str] = None, cwd: Optional[Path] = None) -> subprocess.CompletedProcess:
    """Run a command and raise on non-zero return code. Capture stdout/stderr.

    Args:
        cmd: List of command tokens
        stdin: Optional string to feed to stdin (e.g., group selections "0\n0\n")
        cwd: Working directory
    Returns:
        CompletedProcess
    Raises:
        CommandError on failure
    """
    proc = subprocess.run(
        cmd,
        input=stdin.encode() if stdin is not None else None,
        cwd=str(cwd) if cwd else None,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if proc.returncode != 0:
        raise CommandError(f"Command failed ({proc.returncode}): {' '.join(cmd)}\nOutput:\n{proc.stdout.decode(errors='ignore')}")
    return proc


def which_gmx(user_bin: Optional[str]) -> str:
    """Find the GROMACS binary to use (gmx or gmx_mpi).

    Priority:
    1) --gmx-bin argument
    2) $GMX_BIN env var
    3) 'gmx' in PATH
    4) 'gmx_mpi' in PATH
    """
    candidates: List[Optional[str]] = [
        user_bin,
        os.environ.get("GMX_BIN"),
        shutil.which("gmx"),
        shutil.which("gmx_mpi"),
    ]
    for c in candidates:
        if c and shutil.which(c) or (c and Path(c).exists() and os.access(c, os.X_OK)):
            return c
    raise FileNotFoundError(
        "No GROMACS binary found. Install GROMACS or provide --gmx-bin or set $GMX_BIN."
    )


def ensure_outdir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def parse_xvg(xvg_path: Path) -> np.ndarray:
    """Parse a GROMACS .xvg file into a numpy array.

    Ignores lines starting with '@' and '#'.
    Returns an (N, M) float array.
    """
    data: List[List[float]] = []
    with xvg_path.open("r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line[0] in ("@", "#"):
                continue
            # Split on whitespace; tolerate multiple spaces/tabs
            parts = line.split()
            try:
                row = [float(x) for x in parts]
            except ValueError:
                # skip any non-numeric line silently
                continue
            data.append(row)
    if not data:
        raise ValueError(f"No numeric data parsed from {xvg_path}")
    return np.array(data, dtype=float)


def plot_eigenvalues(eigvals: np.ndarray, out: Path) -> None:
    if plt is None:
        raise RuntimeError("Matplotlib not available to create plots.")
    # eigvals expected shape (N, 2?) or 1D; robustly handle columns
    vals = eigvals[:, -1] if eigvals.ndim == 2 and eigvals.shape[1] > 1 else eigvals.squeeze()
    idx = np.arange(1, len(vals) + 1)
    plt.figure(figsize=(6, 4))
    plt.plot(idx, vals, marker='o')
    plt.xlabel('Principal Component')
    plt.ylabel('Eigenvalue (variance)')
    plt.title('Scree plot (Eigenvalues)')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out.with_suffix('.png'), dpi=200)
    plt.savefig(out.with_suffix('.pdf'))
    plt.close()


def plot_cumulative_variance(eigvals: np.ndarray, out: Path) -> None:
    if plt is None:
        raise RuntimeError("Matplotlib not available to create plots.")
    vals = eigvals[:, -1] if eigvals.ndim == 2 and eigvals.shape[1] > 1 else eigvals.squeeze()
    total = np.sum(vals)
    if total <= 0:
        total = 1.0
    frac = vals / total
    cum = np.cumsum(frac)
    idx = np.arange(1, len(vals) + 1)
    plt.figure(figsize=(6, 4))
    plt.plot(idx, cum, marker='o')
    plt.xlabel('Number of PCs')
    plt.ylabel('Cumulative variance explained')
    plt.title('Cumulative variance')
    plt.ylim(0, 1.05)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out.with_suffix('.png'), dpi=200)
    plt.savefig(out.with_suffix('.pdf'))
    plt.close()


def plot_pc_scatter(proj: np.ndarray, out: Path, color_by_time: bool = True) -> None:
    if plt is None:
        raise RuntimeError("Matplotlib not available to create plots.")
    # Expect columns: time, PC1, PC2 when -first 1 -last 2 used
    if proj.shape[1] < 3:
        raise ValueError("Projection file must have at least 3 columns: time, PC1, PC2")
    t = proj[:, 0]
    pc1 = proj[:, 1]
    pc2 = proj[:, 2]
    plt.figure(figsize=(5, 5))
    if color_by_time:
        sc = plt.scatter(pc1, pc2, c=t, cmap='viridis', s=10, linewidths=0, alpha=0.8)
        cbar = plt.colorbar(sc)
        cbar.set_label('Time (ps)')
    else:
        plt.scatter(pc1, pc2, s=10, linewidths=0, alpha=0.8)
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PC1 vs PC2 projection')
    plt.tight_layout()
    plt.savefig(out.with_suffix('.png'), dpi=200)
    plt.savefig(out.with_suffix('.pdf'))
    plt.close()


def plot_pc_timeseries(proj: np.ndarray, out: Path) -> None:
    if plt is None:
        raise RuntimeError("Matplotlib not available to create plots.")
    t = proj[:, 0]
    pc1 = proj[:, 1]
    plt.figure(figsize=(6, 3.5))
    plt.plot(t, pc1, lw=1)
    plt.xlabel('Time (ps)')
    plt.ylabel('PC1')
    plt.title('PC1 over time')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out.with_suffix('.png'), dpi=200)
    plt.savefig(out.with_suffix('.pdf'))
    plt.close()


def build_index_with_select(gmx: str, s: Path, selection: str, out_ndx: Path, log_dir: Path) -> None:
    cmd = [gmx, "select", "-s", str(s), "-select", selection, "-on", str(out_ndx)]
    proc = run_cmd(cmd)
    (log_dir / "select.log").write_bytes(proc.stdout)


def trj_fit_center(gmx: str, s: Path, f: Path, ndx: Path, out_xtc: Path, pbc: str, center: bool, fit: str, log_dir: Path) -> None:
    cmd = [gmx, "trjconv", "-s", str(s), "-f", str(f), "-o", str(out_xtc), "-n", str(ndx), "-pbc", pbc, "-ur", "compact"]
    if center:
        cmd.append("-center")
    if fit:
        cmd.extend(["-fit", fit])
    # For trjconv, we need to answer: group for centering (if -center) and output group.
    # Our ndx has exactly one selection (index 0), so answer "0" twice if centering, else once.
    answers = ["0"]
    if center:
        answers = ["0", "0"]
    stdin = "\n".join(answers) + "\n"
    proc = run_cmd(cmd, stdin=stdin)
    (log_dir / "trjconv.log").write_bytes(proc.stdout)


def run_covar(gmx: str, s: Path, f: Path, ndx: Path, out_dir: Path, log_dir: Path) -> Tuple[Path, Path, Path]:
    ev_xvg = out_dir / "eigenvalues.xvg"
    ev_trr = out_dir / "eigenvectors.trr"
    avg_pdb = out_dir / "average.pdb"
    cmd = [
        gmx,
        "covar",
        "-s",
        str(s),
        "-f",
        str(f),
        "-n",
        str(ndx),
        "-o",
        str(ev_xvg),
        "-v",
        str(ev_trr),
        "-av",
        str(avg_pdb),
    ]
    # covar asks for analysis and fitting groups; answer 0 twice
    proc = run_cmd(cmd, stdin="0\n0\n")
    (log_dir / "covar.log").write_bytes(proc.stdout)
    return ev_xvg, ev_trr, avg_pdb


def run_anaeig(gmx: str, s: Path, f: Path, ndx: Path, eig_trr: Path, out_dir: Path, first: int, last: int, log_dir: Path) -> Tuple[Path, Path, Path, Path]:
    proj_xvg = out_dir / "projection.xvg"
    eig_xvg = out_dir / "eig.xvg"
    rmsf_xvg = out_dir / "rmsf.xvg"
    comp_xvg = out_dir / "comp.xvg"
    cmd = [
        gmx,
        "anaeig",
        "-s",
        str(s),
        "-f",
        str(f),
        "-n",
        str(ndx),
        "-v",
        str(eig_trr),
        "-first",
        str(first),
        "-last",
        str(last),
        "-proj",
        str(proj_xvg),
        "-eig",
        str(eig_xvg),
        "-rmsf",
        str(rmsf_xvg),
        "-comp",
        str(comp_xvg),
    ]
    # One group selection for analysis
    proc = run_cmd(cmd, stdin="0\n")
    (log_dir / "anaeig.log").write_bytes(proc.stdout)
    return proj_xvg, eig_xvg, rmsf_xvg, comp_xvg


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(description="Automate GROMACS PCA and produce plots.")
    p.add_argument("-s", "--structure", required=True, help="Structure/topology file (TPR recommended; can be TPR/GRO/PDB)")
    p.add_argument("-f", "--trajectory", required=True, help="Trajectory file (XTC/TRR)")
    p.add_argument("-o", "--outdir", default="pca_out", help="Output directory (default: pca_out)")
    p.add_argument("--selection", default="Backbone", help="GROMACS selection text to define atoms (default: Backbone)")
    p.add_argument("--first", type=int, default=1, help="First eigenvector index to analyze (default: 1)")
    p.add_argument("--last", type=int, default=2, help="Last eigenvector index to analyze (default: 2)")
    p.add_argument("--pbc", default="mol", choices=["no", "mol", "res", "atom"], help="PBC handling for trjconv (default: mol)")
    p.add_argument("--no-center", action="store_true", help="Do not center trajectory")
    p.add_argument("--fit", default="rot+trans", choices=["", "none", "rot+trans", "progressive"], help="Fitting mode for trjconv (default: rot+trans). Use 'none' or '' to disable.")
    p.add_argument("--gmx-bin", default=None, help="GROMACS binary to use (e.g., gmx or full path)")
    p.add_argument("--no-plots", action="store_true", help="Skip plotting (still runs GROMACS steps)")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing outputs if present")

    args = p.parse_args(argv)

    # Resolve paths
    s = Path(args.structure).expanduser().resolve()
    f = Path(args.trajectory).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    ensure_outdir(outdir)
    log_dir = outdir / "logs"
    ensure_outdir(log_dir)

    # Basic checks
    if not s.exists():
        print(f"Error: structure file not found: {s}", file=sys.stderr)
        return 2
    if not f.exists():
        print(f"Error: trajectory file not found: {f}", file=sys.stderr)
        return 2

    try:
        gmx = which_gmx(args.gmx_bin)
    except FileNotFoundError as e:
        print(str(e), file=sys.stderr)
        return 2

    # Prepare selection index
    ndx = outdir / "selection.ndx"
    if ndx.exists() and not args.overwrite:
        print(f"[skip] Using existing index: {ndx}")
    else:
        print(f"[1/5] Building selection index with: {args.selection}")
        build_index_with_select(gmx, s, args.selection, ndx, log_dir)

    # Fit/center trajectory
    fit_xtc = outdir / "fit.xtc"
    do_center = not args.no_center
    fit_mode = None if args.fit in ("", "none") else args.fit

    if fit_xtc.exists() and not args.overwrite:
        print(f"[skip] Using existing fitted trajectory: {fit_xtc}")
    else:
        print("[2/5] Fitting and centering trajectory (trjconv)...")
        trj_fit_center(gmx, s, f, ndx, fit_xtc, args.pbc, do_center, fit_mode, log_dir)

    # PCA: covar
    print("[3/5] Running covariance analysis (covar)...")
    ev_xvg, ev_trr, avg_pdb = run_covar(gmx, s, fit_xtc, ndx, outdir, log_dir)

    # PCA: anaeig
    print("[4/5] Analyzing eigenvectors and projecting (anaeig)...")
    proj_xvg, eig_xvg, rmsf_xvg, comp_xvg = run_anaeig(
        gmx, s, fit_xtc, ndx, ev_trr, outdir, args.first, args.last, log_dir
    )

    # Plotting
    if args.no_plots:
        print("[5/5] Skipping plotting as requested (--no-plots)")
        return 0

    if plt is None:
        print("Matplotlib not available; cannot create plots. Install matplotlib or use --no-plots.", file=sys.stderr)
        return 3

    print("[5/5] Parsing XVG files and generating plots...")
    try:
        eig_data = parse_xvg(ev_xvg)
        plot_eigenvalues(eig_data, outdir / "plot_eigenvalues")
        plot_cumulative_variance(eig_data, outdir / "plot_cumulative_variance")
    except Exception as e:
        print(f"Warning: failed eigenvalue plots: {e}")

    try:
        proj_data = parse_xvg(proj_xvg)
        plot_pc_scatter(proj_data, outdir / "plot_pc1_pc2_scatter", color_by_time=True)
        plot_pc_timeseries(proj_data, outdir / "plot_pc1_timeseries")
    except Exception as e:
        print(f"Warning: failed projection plots: {e}")

    print("Done. Outputs in:", outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
gmx_panel - Web dashboard for GROMACS analysis XVG files.

Scans a directory for XVG files (or uses a user-provided list) and serves
an HTML panel displaying all plots via a local web server.
Plots are rendered by the plot_xvg.py tool located in the same directory.

Usage:
  gmx_panel.py [--dir DIR] [--files FILE [FILE ...]] [--port PORT] [OPTIONS]
"""
import argparse
import base64
import io
import os
import re
import socket
import sys
from datetime import datetime
from functools import partial
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

# Set Agg backend before any pyplot import (required for off-screen rendering)
import matplotlib
matplotlib.use('Agg')

# Allow importing plot_xvg from the same bin/ directory
_SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPT_DIR))

from plot_xvg import parse_xvg, plot_xvg  # noqa: E402 (must follow matplotlib.use)
import matplotlib.pyplot as plt            # noqa: E402


# ---------------------------------------------------------------------------
# XVG file discovery
# ---------------------------------------------------------------------------

def find_xvg_files(directory: Path, names: list[str] | None = None) -> list[Path]:
    """Return XVG files from *directory*.

    If *names* is provided, resolve each name relative to *directory* and warn
    about missing ones.  Otherwise, return all *.xvg files sorted by name.
    """
    if names:
        found = []
        for name in names:
            p = (directory / name) if not Path(name).is_absolute() else Path(name)
            if p.exists():
                found.append(p)
            else:
                print(f"Warning: {p} not found, skipping.", file=sys.stderr)
        return found
    return sorted(directory.glob("*.xvg"))


# ---------------------------------------------------------------------------
# Plot rendering
# ---------------------------------------------------------------------------

def render_xvg_to_png_b64(xvg_path: Path, style: str, figsize: tuple[float, float],
                            dpi: int) -> tuple[str | None, str | None]:
    """Render *xvg_path* to a base64-encoded PNG.

    Returns ``(b64_string, None)`` on success or ``(None, error_message)`` on failure.
    """
    try:
        fig, _ax = plot_xvg(str(xvg_path), style=style)
        fig.set_size_inches(*figsize)
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        return base64.b64encode(buf.read()).decode("ascii"), None
    except Exception as exc:
        plt.close("all")
        return None, str(exc)


def get_xvg_display_title(xvg_path: Path) -> str:
    """Return a human-readable title extracted from the XVG header, or the stem."""
    try:
        _, _, labels = parse_xvg(str(xvg_path))
        title = labels.get("title", "")
        subtitle = labels.get("subtitle", "")
        parts = [t for t in (title, subtitle) if t]
        return " — ".join(parts) if parts else xvg_path.stem
    except Exception:
        return xvg_path.stem


# ---------------------------------------------------------------------------
# HTML generation
# ---------------------------------------------------------------------------

_CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: system-ui, -apple-system, sans-serif;
    background: #12121f;
    color: #d8d8e8;
    min-height: 100vh;
}
header {
    background: #1a1a30;
    border-bottom: 2px solid #2a2a50;
    padding: 0.9rem 1.5rem;
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 1rem;
    flex-wrap: wrap;
}
header h1 { font-size: 1.3rem; color: #e06c75; font-weight: 700; }
.header-meta { font-size: 0.8rem; color: #888; font-family: monospace; }
.toolbar {
    display: flex;
    align-items: center;
    gap: 0.75rem;
}
.btn {
    background: #2a2a50;
    color: #d8d8e8;
    border: 1px solid #3a3a70;
    border-radius: 4px;
    padding: 0.3rem 0.8rem;
    cursor: pointer;
    font-size: 0.8rem;
    text-decoration: none;
    transition: background 0.15s;
}
.btn:hover { background: #3a3a70; }
.grid {
    display: grid;
    grid-template-columns: {COLUMNS};
    gap: 1.25rem;
    padding: 1.25rem;
}
.card {
    background: #1a1a30;
    border: 1px solid #2a2a50;
    border-radius: 8px;
    overflow: hidden;
    transition: transform 0.15s, box-shadow 0.15s;
}
.card:hover {
    transform: translateY(-2px);
    box-shadow: 0 4px 18px rgba(224, 108, 117, 0.18);
}
.card-header {
    padding: 0.55rem 0.9rem;
    background: #222240;
    border-bottom: 1px solid #2a2a50;
}
.card-header h3 { font-size: 0.88rem; font-weight: 600; color: #d8d8e8; }
.card-header .fname { font-size: 0.72rem; color: #666; font-family: monospace; }
.card img { width: 100%; height: auto; display: block; }
.card-error { border-color: #e06c75; }
.error-body {
    padding: 0.9rem;
    font-size: 0.8rem;
    color: #e06c75;
    font-family: monospace;
    white-space: pre-wrap;
    word-break: break-word;
}
@media (max-width: 800px) { .grid { grid-template-columns: 1fr; } }
"""

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  {REFRESH_META}
  <title>{TITLE}</title>
  <style>{CSS}</style>
</head>
<body>
  <header>
    <div>
      <h1>{TITLE}</h1>
      <div class="header-meta">{META}</div>
    </div>
    <div class="toolbar">
      <a class="btn" href="/">&#8635; Reload</a>
    </div>
  </header>
  <div class="grid">
    {CARDS}
  </div>
</body>
</html>
"""

_CARD_OK = """\
<div class="card">
  <div class="card-header">
    <h3>{TITLE}</h3>
    <div class="fname">{FNAME}</div>
  </div>
  <img src="data:image/png;base64,{IMG}" alt="{TITLE}" loading="lazy" />
</div>"""

_CARD_ERR = """\
<div class="card card-error">
  <div class="card-header">
    <h3>{TITLE}</h3>
    <div class="fname">{FNAME}</div>
  </div>
  <div class="error-body">Error rendering plot:\n{ERR}</div>
</div>"""


def build_html_page(xvg_files: list[Path], title: str, style: str,
                    figsize: tuple[float, float], dpi: int,
                    columns: int, refresh: int) -> str:
    """Render all XVG files and assemble the dashboard HTML."""
    cards = []
    for xvg in xvg_files:
        display_title = get_xvg_display_title(xvg)
        img_b64, err = render_xvg_to_png_b64(xvg, style=style, figsize=figsize, dpi=dpi)
        if img_b64:
            cards.append(_CARD_OK.format(TITLE=display_title, FNAME=xvg.name, IMG=img_b64))
        else:
            cards.append(_CARD_ERR.format(TITLE=display_title, FNAME=xvg.name,
                                          ERR=err or "unknown error"))

    refresh_meta = (f'<meta http-equiv="refresh" content="{refresh}">'
                    if refresh > 0 else "")
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    meta = (f"{len(xvg_files)} plot(s) | rendered {now}"
            + (f" | auto-refresh {refresh}s" if refresh else ""))
    grid_cols = f"repeat({columns}, 1fr)"

    return _HTML_TEMPLATE.format(
        TITLE=title,
        REFRESH_META=refresh_meta,
        CSS=_CSS.replace("{COLUMNS}", grid_cols),
        META=meta,
        CARDS="\n    ".join(cards),
    )


# ---------------------------------------------------------------------------
# HTTP server
# ---------------------------------------------------------------------------

class PanelHandler(BaseHTTPRequestHandler):
    """Serve the XVG dashboard; re-renders on every GET /."""

    def __init__(self, xvg_files, panel_args, *args, **kwargs):
        self._xvg_files = xvg_files
        self._args = panel_args
        super().__init__(*args, **kwargs)

    def do_GET(self):
        # Strip query string for path matching
        path = self.path.split("?")[0]
        if path not in ("/", "/index.html"):
            self.send_error(404)
            return

        # On each request, re-scan the directory so new files appear automatically
        args = self._args
        if args.files:
            xvg_files = self._xvg_files
        else:
            xvg_files = find_xvg_files(Path(args.dir))

        html = build_html_page(
            xvg_files,
            title=args.title,
            style=args.style,
            figsize=(args.figsize[0], args.figsize[1]),
            dpi=args.dpi,
            columns=args.columns,
            refresh=args.refresh,
        )
        body = html.encode("utf-8")
        self.send_response(200)
        self.send_header("Content-Type", "text/html; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, fmt, *args):
        # Only log non-2xx responses to avoid noise
        if args and len(args) >= 2 and not str(args[1]).startswith("2"):
            super().log_message(fmt, *args)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        prog="gmx_panel",
        description="Serve a web dashboard of GROMACS XVG analysis plots.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Scan the current directory and serve all .xvg files
  %(prog)s

  # Serve from a specific analysis output directory
  %(prog)s --dir ./MnMT4_apo_0_analysis

  # Show only selected files (searched inside --dir)
  %(prog)s --dir ./run_analysis --files rmsd.xvg temperature.xvg rmsf.xvg

  # Serve on a different port with auto-refresh every 60 s
  %(prog)s --dir ./analysis --port 9090 --refresh 60

  # Three-column grid with line+dot style
  %(prog)s --columns 3 --style lines+dots

  # Wider, higher-resolution plot cards
  %(prog)s --figsize 10 5 --dpi 120
        """,
    )
    parser.add_argument(
        "--dir", "-d", type=str, default=".",
        help="Directory to scan for XVG files (default: current directory)",
    )
    parser.add_argument(
        "--files", "-f", nargs="+", type=str, default=None,
        help="Explicit list of XVG filenames/paths to include (names resolved "
             "relative to --dir; overrides auto-scan)",
    )
    parser.add_argument(
        "--port", "-p", type=int, default=8080,
        help="Port to serve the dashboard on (default: 8080)",
    )
    parser.add_argument(
        "--host", type=str, default="0.0.0.0",
        help="Host/interface to bind to (default: 0.0.0.0 — all interfaces)",
    )
    parser.add_argument(
        "--title", type=str, default="GROMACS Analysis Dashboard",
        help='Dashboard title shown in the browser (default: "GROMACS Analysis Dashboard")',
    )
    parser.add_argument(
        "--style", "-s", type=str, default="lines",
        choices=["dots", "lines", "lines+dots"],
        help="Plot style passed to plot_xvg (default: lines)",
    )
    parser.add_argument(
        "--columns", "-c", type=int, default=2,
        help="Number of columns in the plot grid (default: 2)",
    )
    parser.add_argument(
        "--figsize", nargs=2, type=float, default=[9.0, 4.5],
        metavar=("W", "H"),
        help="Plot figure size in inches (default: 9 4.5)",
    )
    parser.add_argument(
        "--dpi", type=int, default=100,
        help="Image resolution in DPI (default: 100)",
    )
    parser.add_argument(
        "--refresh", "-r", type=int, default=0,
        help="Auto-refresh interval in seconds (default: 0 = disabled)",
    )

    args = parser.parse_args()

    directory = Path(args.dir).resolve()
    if not directory.is_dir():
        sys.exit(f"Error: '{args.dir}' is not a directory.")

    # Initial file list (for --files mode the list is fixed; for scan mode it
    # is re-built on every request inside the handler)
    xvg_files = find_xvg_files(directory, args.files)
    if not xvg_files:
        msg = f"No XVG files found in '{directory}'."
        if not args.files:
            msg += " Use --files to specify files or point --dir at an analysis directory."
        sys.exit(f"Error: {msg}")

    print(f"Found {len(xvg_files)} XVG file(s) in {directory}:")
    for f in xvg_files:
        print(f"  {f.name}")

    try:
        hostname = socket.gethostname()
    except Exception:
        hostname = "localhost"

    handler_factory = partial(PanelHandler, xvg_files, args)
    server = HTTPServer((args.host, args.port), handler_factory)

    print(f"\nDashboard available at:")
    print(f"  http://localhost:{args.port}/")
    if hostname != "localhost":
        print(f"  http://{hostname}:{args.port}/")
    if args.refresh:
        print(f"  (auto-refresh every {args.refresh}s)")
    print("\nPress Ctrl+C to stop.\n")

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        server.shutdown()


if __name__ == "__main__":
    main()

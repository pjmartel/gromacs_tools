# gmx_panel.py

Web dashboard for GROMACS analysis XVG files.

## Overview

Scans a directory for `.xvg` files (or uses a list you provide) and serves an HTML page
showing every plot as a card in a grid, through a small local web server. Plots are drawn by
[plot_xvg.py](plot_xvg.md) from the same `bin/` directory, so they look the same as
command-line plots. It pairs naturally with the output directory of
[gmx_analysis.sh](gmx_analysis.md).

It is especially handy on remote machines/clusters: start it next to the analysis
output and open the dashboard in a local browser instead of copying files around.

## Features

- ✅ **Zero-config**: `gmx_panel.py` in an analysis directory serves all its `.xvg` files
- ✅ **Live**: files are re-scanned and plots re-rendered on every page load, so new or
  updated `.xvg` files appear on reload (or automatically with `--refresh`)
- ✅ **Titles from the XVG header**: each card shows the file's `@ title` (and `@ subtitle`), plus the file name
- ✅ **Fault tolerant**: a file that can't be plotted shows an error card instead of breaking the page
- ✅ **Analysis order**: plots from a directory scan are shown in a fixed, meaningful order
  (RMSD, secondary structure, RMSF, gyration, SASA, energies; see [Plot Order](#plot-order))
- ✅ **Scalable plots**: a scale factor (`--scale`, or the −/+ buttons in the page) sets the plot size,
  and the grid fits as many plots per row as the window allows
- ✅ **Energy statistics**: optionally show the `gmx energy` average, RMSD and drift on each energy
  plot (`--stats`, or the **Stats** button; see [Energy Statistics](#energy-statistics))
- ✅ **Configurable layout**: fixed grid columns, figure size, resolution, plot style
- ✅ **No extra dependencies**: uses Python's built-in `http.server`, plus numpy/matplotlib

## Usage

### Basic Syntax

```bash
gmx_panel.py [--dir DIR] [--files FILE [FILE ...]] [--port PORT] [OPTIONS]
```

Then open `http://localhost:8080/` in a browser. Stop the server with `Ctrl+C`.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `-d, --dir <dir>` | `.` | Directory to scan for `*.xvg` files |
| `-f, --files <file> [...]` | - | Explicit list of XVG files to show, in this order. Names are resolved relative to `--dir`; absolute paths are used as is. Disables auto-scan |
| `-p, --port <port>` | `8080` | Port to serve the dashboard on |
| `--host <addr>` | `0.0.0.0` | Interface to bind to (`0.0.0.0` = all interfaces) |
| `--title <text>` | `GROMACS Analysis Dashboard` | Page title |
| `-s, --style <style>` | `lines` | Plot style passed to `plot_xvg`: `dots`, `lines` or `lines+dots` |
| `--stats` | off | Show the `gmx energy` statistics box on energy plots (see [Energy Statistics](#energy-statistics)). Can also be toggled from the page |
| `--scale <factor>` | `1.0` | Plot size factor, from `0.25` to `4`: each plot card is 600 px wide times this factor, and as many fit in a row as the window allows. Can also be changed from the page (see [Plot Size](#plot-size)) |
| `-c, --columns <n>` | automatic | Fixed number of columns in the plot grid, overriding the automatic layout from `--scale` (collapses to one column on narrow screens) |
| `--figsize <W> <H>` | `9 4.5` | Figure size of each plot, in inches |
| `--dpi <n>` | `100` | Image resolution |
| `-r, --refresh <s>` | `0` | Auto-refresh interval in seconds (`0` = disabled) |
| `-h, --help` | - | Show help and exit |

## How It Works

1. On startup, the file list is built once to check there is something to show; the server
   exits with an error if no `.xvg` files are found.
2. On every `GET /`, the directory is re-scanned (unless `--files` fixed the list), the files
   are put in the [default order](#plot-order), and each one is rendered with
   `plot_xvg.plot_xvg()` using the matplotlib `Agg` backend.
3. Each plot is embedded in the page as a base64 PNG, so the page is a single
   self-contained HTML document (you can also save it from the browser).

Plots use the default single-file mode of `plot_xvg`: all data columns, with legends and axis
labels from the XVG header. Directory scans only pick up `*.xvg` files; other outputs such
as `secondary_structure.dat` are ignored.

> [!NOTE]
> Rendering happens on every request, so a directory with many large files (e.g. long
> trajectories at full resolution) can take a few seconds to load. Use `--files` to limit the
> set, or a lower `--dpi`.

## Plot Order

When scanning a directory, plots are shown in this order, based on the file names written by
[gmx_analysis.sh](gmx_analysis.md):

| # | Group | Files |
|---|-------|-------|
| 1 | RMSD | `rmsd*.xvg` |
| 2 | Secondary structure (DSSP) | `secondary_structure*.xvg`, `dssp*.xvg` |
| 3 | RMSF | `rmsf*.xvg` |
| 4 | Radius of gyration | `gyrat*.xvg` (e.g. `gyration_radius.xvg`) |
| 5 | SASA | `sasa*.xvg` |
| 6 | Energy terms | `temperature`, `pressure`, `potential`, `total_energy` (in that order), then any other `gmx energy` output (e.g. extra `--energy-terms`), recognized from the command line recorded in its header, or from the default `GROMACS Energies` title |
| 7 | Everything else | alphabetically |

File names are matched case-insensitively, and files within a group are sorted alphabetically.
With `--files`, the plots are shown in the order given instead.

## Plot Size

Each plot card is 600 px wide at scale 1.0. The grid fills each row with as many cards as fit
in the browser window, stretching them to use the full width. Increasing the scale gives bigger
plots and fewer per row; decreasing it fits more plots per row. Resizing the window reflows
the grid.

Set the initial scale with `--scale`, or change it while viewing with the toolbar buttons:

| Button | Action |
|--------|--------|
| **−** / **+** | Smaller / larger plots (steps of ×1.25) |
| **Reset** | Back to the `--scale` value |
| **Stats** | Show/hide the [energy statistics](#energy-statistics) (highlighted when shown) |
| **Reload** | Re-render the plots, keeping the current settings |

The current scale and statistics setting are part of the page URL (e.g.
`http://localhost:8080/?scale=1.5&stats=1`), so they survive auto-refresh and can be bookmarked. When enlarged (scale above 1), plots are rendered
at `--dpi` × scale so they stay sharp.

With `--columns`, the number of plots per row is fixed and the zoom buttons are hidden (the
**Stats** button remains).

## Energy Statistics

For energy terms, [gmx_analysis.sh](gmx_analysis.md) stores the statistics `gmx energy`
prints (average, error estimate, RMSD and total drift) in each `.xvg` file. With `--stats`,
or the **Stats** toolbar button, each energy plot shows the average, RMSD and total drift in a
small box in its top-left corner, and its legend moves to the top-right corner. The box is
drawn by [plot_xvg.py --stats](plot_xvg.md#statistics-box---stats), which also sets which
values are shown.

The statistics are off by default. Plots without statistics (RMSD, RMSF, ..., or energy files
produced before this feature) are unaffected; re-run the energy step of `gmx_analysis.sh` to
add statistics to older files.

## Remote Machines

By default the server listens on all interfaces, and on startup it prints both a
`localhost` URL and one using the machine's hostname.

The dashboard has no authentication. On a shared or public machine, bind to localhost only
and reach it through an SSH tunnel:

```bash
# On the remote machine
python bin/gmx_panel.py --dir MnMT4_apo_0_analysis --host 127.0.0.1 --port 8080

# On your local machine, then open http://localhost:8080/
ssh -L 8080:localhost:8080 user@remote-host
```

## Examples

### Current Directory

```bash
# Serve every .xvg file in the current directory
python bin/gmx_panel.py
```

### gmx_analysis.sh Output

```bash
python bin/gmx_panel.py --dir ./MnMT4_apo_0_analysis --title "MnMT4 apo, replica 0"
```

### Selected Files Only

```bash
# Shown in the order given
python bin/gmx_panel.py --dir ./run_analysis --files rmsd.xvg temperature.xvg rmsf.xvg
```

### Monitoring a Running Analysis

```bash
# Different port, reload every 60 s
python bin/gmx_panel.py --dir ./analysis --port 9090 --refresh 60
```

### Layout

```bash
# Larger plots, fewer per row (also adjustable with the -/+ buttons in the page)
python bin/gmx_panel.py --scale 1.5

# Small plots, many per row: a quick overview of many files
python bin/gmx_panel.py --scale 0.6

# Fixed three-column grid with line+dot style
python bin/gmx_panel.py --columns 3 --style lines+dots

# Wider, higher-resolution plot cards
python bin/gmx_panel.py --figsize 10 5 --dpi 120
```

## Troubleshooting

### "No XVG files found"
Point `--dir` at the directory that contains the `.xvg` files (e.g. the `--output-dir` of
`gmx_analysis.sh`), or list them with `--files`.

### "Address already in use"
Another process (or a previous dashboard) is using the port. Pick another one with `--port`.

### A card shows "Error rendering plot"
The file could not be parsed or plotted by `plot_xvg`. The error message is shown on the
card. Check the file with `python bin/plot_xvg.py <file>` for more detail.

#!/usr/bin/env python3
"""
Command-line tool to plot GROMACS XVG files with Matplotlib.
Supports single or multiple XVG files, moving averages, and customization options.
Includes an option to select the Matplotlib backend (e.g., Qt5Agg) from CLI.
"""
import argparse
import sys
from pathlib import Path

# Configure Matplotlib backend BEFORE importing pyplot
import matplotlib

def _get_backend_from_argv(argv):
    """Extract --backend/--mpl-backend from argv without importing pyplot."""
    for i, arg in enumerate(argv):
        if arg in ("--backend", "--mpl-backend"):
            if i + 1 < len(argv):
                return argv[i + 1]
        elif arg.startswith("--backend="):
            return arg.split("=", 1)[1]
        elif arg.startswith("--mpl-backend="):
            return arg.split("=", 1)[1]
    return None

_cli_backend = _get_backend_from_argv(sys.argv[1:])
if _cli_backend:
    try:
        matplotlib.use(_cli_backend, force=True)
    except Exception as _be:
        # Defer error to later user-visible message; keep default backend
        print(f"Warning: Failed to set Matplotlib backend to '{_cli_backend}': {_be}", file=sys.stderr)

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def moving_average(data, window_size):
  if window_size < 1:
    raise ValueError("Window size must be at least 1.")
  return np.convolve(data, np.ones(window_size) / window_size, mode='valid')


def parse_xvg(filename):
    legends = {}
    labels = {}
    data_lines = []

    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            if line.startswith('@'):
                parts = line.split()
                if parts[1] == 'title':
                    labels['title'] = ' '.join(parts[2:]).strip('"')
                elif parts[1] == 'xaxis' and parts[2] == 'label':
                    labels['xlabel'] = ' '.join(parts[3:]).strip('"')
                elif parts[1] == 'yaxis' and parts[2] == 'label':
                    labels['ylabel'] = ' '.join(parts[3:]).strip('"')
                elif parts[1].startswith('s') and parts[2] == 'legend':
                    index = int(parts[1][1:])
                    legends[index] = ' '.join(parts[3:]).strip('"')
            else:
                data_lines.append(line)

    # Parse data into columns
    data = [list(map(float, line.split())) for line in data_lines]
    data_by_columns = list(zip(*data))  # transpose rows to columns

    return data_by_columns, legends, labels

def plot_xvg(filename, show_moving_avg=False, window_size=10, style='dots', 
             scatter_colormap='viridis', use_scatter=False, use_histogram=False, 
             hist_bins=50, markersize=3):
    """Plot a single XVG file and return (fig, ax) without displaying.

    The previous implementation called plt.show() internally which prevented
    external customization (e.g. --title) before display and led to an extra
    blank window. This refactored version leaves display to the caller.
    
    Parameters
    ----------
    use_histogram : bool
        If True, plot histogram of the second column (y-values) instead of xy plot.
    hist_bins : int
        Number of bins for histogram (default: 50).
    markersize : float
        Size of markers/dots in plots (default: 3).
    """
    data_columns, legends, labels = parse_xvg(filename)
    x = data_columns[0]
    num_datasets = len(data_columns) - 1

    fig, ax = plt.subplots(figsize=(10, 6))

    if use_histogram:
        # Histogram mode: plot distribution of second column (y values)
        for i in range(num_datasets):
            y = data_columns[i + 1]
            label = legends.get(i, f'Dataset {i}')
            ax.hist(y, bins=hist_bins, alpha=0.7, label=label, edgecolor='black')
        
        ax.set_xlabel(labels.get('ylabel', 'Value'))  # histogram x-axis is the y-values
        ax.set_ylabel('Frequency')
        ax.set_title(labels.get('title', 'Distribution'))
        if num_datasets > 1 or legends:
            ax.legend()
    else:
        # Standard xy plot mode
        for i in range(num_datasets):
            y = data_columns[i + 1]
            label = legends.get(i, f'Dataset {i}')

            if use_scatter:
                colors = np.arange(len(x))
                scatter = ax.scatter(x, y, c=colors, cmap=scatter_colormap,
                                     s=markersize*7, alpha=0.7, label=label)
                if num_datasets == 1:  # Only add colorbar for single dataset
                    cbar = fig.colorbar(scatter, ax=ax)
                    cbar.set_label('Frame / Time order', rotation=270, labelpad=20)
            else:
                if style == 'dots':
                    ax.plot(x, y, 'o', markersize=markersize, label=label)
                elif style == 'lines':
                    ax.plot(x, y, '-', label=label)
                elif style == 'lines+dots':
                    ax.plot(x, y, 'o-', markersize=markersize, label=label)
                else:
                    ax.plot(x, y, style, label=label)

            if show_moving_avg and num_datasets == 1 and not use_scatter:
                y_avg = moving_average(y, window_size)
                x_avg = x[:len(y_avg)]
                ax.plot(x_avg, y_avg, label=f'{label} (Moving Avg, window={window_size})',
                        linestyle='--', linewidth=2)

        ax.set_xlabel(labels.get('xlabel', 'X-axis'))
        ax.set_ylabel(labels.get('ylabel', 'Y-axis'))
        ax.set_title(labels.get('title', ''))
        if num_datasets > 1 or legends:
            ax.legend()
    
    ax.grid(True)
    fig.tight_layout()
    return fig, ax

# This version supports chaining multiple plots on the same axes, based on a
# shared ax variable (should replace the plot_xvg, but it is not completely tested)
def plot_xvg_multi(filename, show_moving_avg=False, window_size=10, ax=None, 
                   custom_legend=None, style='dots', scatter_colormap='viridis', 
                   use_scatter=False, use_histogram=False, hist_bins=50, markersize=3):

    data_columns, legends, labels = parse_xvg(filename)
    x = data_columns[0]
    num_datasets = len(data_columns) - 1

    # Create new Axes if none provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    if use_histogram:
        # Histogram mode: plot distribution of second column
        for i in range(num_datasets):
            y = data_columns[i + 1]
            label = custom_legend if custom_legend and num_datasets == 1 else legends.get(i, f'Dataset {i}')
            ax.hist(y, bins=hist_bins, alpha=0.7, label=label, edgecolor='black')
        
        ax.set_xlabel(labels.get('ylabel', 'Value'))
        ax.set_ylabel('Frequency')
        ax.set_title(labels.get('title', 'Distribution'))
    else:
        # Standard xy plot mode
        for i in range(num_datasets):
            y = data_columns[i + 1]
            label = custom_legend if custom_legend and num_datasets == 1 else legends.get(i, f'Dataset {i}')
            
            if use_scatter:
                # Scatter mode: color by order (time-like)
                colors = np.arange(len(x))
                scatter = ax.scatter(x, y, c=colors, cmap=scatter_colormap, 
                                   s=markersize*7, alpha=0.7, label=label)
                # Note: colorbar is tricky with shared axes, user should add manually if needed
            else:
                # Line/marker mode
                if style == 'dots':
                    ax.plot(x, y, 'o', markersize=markersize, label=label)
                elif style == 'lines':
                    ax.plot(x, y, '-', label=label)
                elif style == 'lines+dots':
                    ax.plot(x, y, 'o-', markersize=markersize, label=label)
                else:  # allow custom matplotlib format strings
                    ax.plot(x, y, style, label=label)

            if show_moving_avg and num_datasets == 1 and not use_scatter:
                y_avg = moving_average(y, window_size)
                x_avg = x[:len(y_avg)]
                ax.plot(x_avg, y_avg, linestyle='--', linewidth=2,
                        label=f'{label} (Moving Avg, window={window_size})')

        ax.set_xlabel(labels.get('xlabel', 'X-axis'))
        ax.set_ylabel(labels.get('ylabel', 'Y-axis'))
        ax.set_title(labels.get('title', ''))

    if num_datasets > 1 or legends or custom_legend:
        ax.legend()
    ax.grid(True)

    return ax


def main():
    """Main function for command-line interface."""
    parser = argparse.ArgumentParser(
        description='Plot GROMACS XVG files with optional moving averages.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Plot a single XVG file (default: dots)
  %(prog)s energy.xvg
  
  # Plot with line style
  %(prog)s energy.xvg --style lines
  
  # Plot with moving average
  %(prog)s energy.xvg --moving-avg --window 50 --style lines
  
  # Scatter plot colored by frame/time order (like PCA plots)
  %(prog)s pc_projection.xvg --scatter --colormap viridis
  
  # Histogram of values (e.g., RMSD distribution)
  %(prog)s rmsd.xvg --histogram --bins 100
  
  # XY correlation: plot y-values from two files against each other
  %(prog)s rmsd1.xvg rmsd2.xvg --xy-correlation --title "RMSD1 vs RMSD2"
  
  # 3D correlation: plot y-values from three files as x, y, z
  %(prog)s pc1.xvg pc2.xvg pc3.xvg --xy-correlation --scatter \\
    --title "PC1 vs PC2 vs PC3" --colormap viridis
  
  # Plot multiple XVG files on the same axes
  %(prog)s file1.xvg file2.xvg file3.xvg --multi
  
  # Plot multiple files with custom legends and dots+lines
  %(prog)s file1.xvg file2.xvg --multi --legends "Run 1" "Run 2" --style lines+dots
  
  # Save plot to file instead of displaying
  %(prog)s energy.xvg --output energy_plot.png
  
  # Combine multiple options
  %(prog)s rmsd.xvg --moving-avg --window 100 --output rmsd.pdf --style lines
        """)
    
    parser.add_argument('files', nargs='+', type=str,
                        help='XVG file(s) to plot')
    
    parser.add_argument('--moving-avg', '--ma', action='store_true',
                        help='Apply moving average filter (only for single dataset)')
    
    parser.add_argument('--window', '-w', type=int, default=10,
                        help='Window size for moving average (default: 10)')
    
    parser.add_argument('--style', '-s', type=str, default='dots',
                        choices=['dots', 'lines', 'lines+dots'],
                        help='Plot style: dots (scatter), lines, or lines+dots (default: dots). '
                             'You can also use custom matplotlib format strings (e.g., "o-", ".-", "s")')
    
    parser.add_argument('--scatter', action='store_true',
                        help='Use scatter plot mode with points colored by frame/time order. '
                             'Creates plots similar to PCA projections with color gradient.')
    
    parser.add_argument('--histogram', '--hist', action='store_true',
                        help='Plot histogram of the second column (y-values) instead of xy plot. '
                             'Useful for distribution analysis (e.g., RMSD, energy, distances).')
    
    parser.add_argument('--xy-correlation', '--xycorr', action='store_true',
                        help='Plot second column of first file vs second column of second file. '
                             'Requires exactly 2 files with the same number of data rows. '
                             'Useful for correlation analysis (e.g., RMSD1 vs RMSD2, PC1 vs PC2). '
                             'With 3 files, creates a 3D correlation plot (x, y, z from files 1, 2, 3).')
    
    parser.add_argument('--bins', type=int, default=50,
                        help='Number of bins for histogram mode (default: 50)')
    
    parser.add_argument('--markersize', '--ms', type=float, default=3.0,
                        help='Size of markers/dots in plots (default: 3.0). '
                             'For scatter plots, the actual size is markersize*7 for better visibility.')
    
    parser.add_argument('--colormap', '--cmap', type=str, default='viridis',
                        help='Colormap for scatter mode (default: viridis). '
                             'Common options: viridis, plasma, inferno, magma, coolwarm, rainbow')
    
    parser.add_argument('--multi', '-m', action='store_true',
                        help='Plot multiple XVG files on the same axes')
    
    parser.add_argument('--legends', '-l', nargs='+', type=str,
                        help='Custom legend labels for multiple files (must match number of files)')
    
    parser.add_argument('--output', '-o', type=str,
                        help='Output file for saving the plot (e.g., plot.png, plot.pdf). '
                             'If not specified, plot will be displayed interactively.')
    
    parser.add_argument('--title', '-t', type=str,
                        help='Custom title for the plot (overrides XVG title)')
    
    parser.add_argument('--xlabel', type=str,
                        help='Custom x-axis label (overrides XVG label)')
    
    parser.add_argument('--ylabel', type=str,
                        help='Custom y-axis label (overrides XVG label)')
    
    parser.add_argument('--figsize', nargs=2, type=float, default=[10, 6],
                        metavar=('WIDTH', 'HEIGHT'),
                        help='Figure size in inches (default: 10 6)')
    
    parser.add_argument('--aspect', type=str, default=None,
                        help='Aspect ratio of the plot. Options: "equal" (1:1), "auto" (default), '
                             'or a number (e.g., 1.5 for width/height ratio). '
                             'Useful for correlation plots or when you need square axes.')
    
    parser.add_argument('--dpi', type=int, default=300,
                        help='DPI for saved figures (default: 300)')

    parser.add_argument('--backend', '--mpl-backend', type=str,
                        help='Matplotlib backend to use (e.g., Qt5Agg, TkAgg, Agg). '
                             'Note: must be provided before other options to take effect; '
                             'this script reads it early to configure Matplotlib before importing pyplot.')
    
    args = parser.parse_args()
    
    # Validate inputs
    if args.legends and len(args.legends) != len(args.files):
        parser.error(f"Number of legends ({len(args.legends)}) must match number of files ({len(args.files)})")
    
    if args.window < 1:
        parser.error("Window size must be at least 1")
    
    if args.xy_correlation:
        if len(args.files) not in (2, 3):
            parser.error("--xy-correlation requires exactly 2 files (2D) or 3 files (3D)")
        if args.multi:
            parser.error("--xy-correlation cannot be used with --multi")
        if args.histogram:
            parser.error("--xy-correlation cannot be used with --histogram")
    
    # Check that all files exist
    for file in args.files:
        if not Path(file).exists():
            print(f"Error: File not found: {file}", file=sys.stderr)
            return 1
    
    try:
        if args.xy_correlation:
            if len(args.files) == 2:
                # 2D correlation mode: plot y-values from file1 vs file2
                data1_columns, legends1, labels1 = parse_xvg(args.files[0])
                data2_columns, legends2, labels2 = parse_xvg(args.files[1])
                
                # Validate that both files have the same number of data points
                y1 = data1_columns[1]  # second column (y-values) from file 1
                y2 = data2_columns[1]  # second column (y-values) from file 2
                
                if len(y1) != len(y2):
                    print(f"Error: Files have different number of data rows!", file=sys.stderr)
                    print(f"  {args.files[0]}: {len(y1)} rows", file=sys.stderr)
                    print(f"  {args.files[1]}: {len(y2)} rows", file=sys.stderr)
                    return 1
                
                # Create the 2D correlation plot
                fig, ax = plt.subplots(figsize=tuple(args.figsize))
                
                if args.scatter:
                    # Scatter mode with color by order
                    colors = np.arange(len(y1))
                    scatter = ax.scatter(y1, y2, c=colors, cmap=args.colormap, 
                                       s=args.markersize*7, alpha=0.7)
                    cbar = fig.colorbar(scatter, ax=ax)
                    cbar.set_label('Frame / Time order', rotation=270, labelpad=20)
                else:
                    # Simple dots or specified style
                    if args.style == 'dots':
                        ax.plot(y1, y2, 'o', markersize=args.markersize)
                    elif args.style == 'lines':
                        ax.plot(y1, y2, '-')
                    elif args.style == 'lines+dots':
                        ax.plot(y1, y2, 'o-', markersize=args.markersize)
                    else:
                        ax.plot(y1, y2, args.style)
                
                # Set labels (use ylabel from each file as axis labels)
                xlabel_default = labels1.get('ylabel', f'Values from {Path(args.files[0]).name}')
                ylabel_default = labels2.get('ylabel', f'Values from {Path(args.files[1]).name}')
                
                ax.set_xlabel(args.xlabel if args.xlabel else xlabel_default)
                ax.set_ylabel(args.ylabel if args.ylabel else ylabel_default)
                ax.set_title(args.title if args.title else 'XY Correlation')
                ax.grid(True)
                
                # Set aspect ratio if specified
                if args.aspect:
                    if args.aspect.lower() == 'equal':
                        ax.set_aspect('equal', adjustable='box')
                    elif args.aspect.lower() == 'auto':
                        ax.set_aspect('auto')
                    else:
                        try:
                            ax.set_aspect(float(args.aspect), adjustable='box')
                        except ValueError:
                            print(f"Warning: Invalid aspect ratio '{args.aspect}', using auto", file=sys.stderr)
                
                fig.tight_layout()
            
            else:  # len(args.files) == 3
                # 3D correlation mode: plot y-values from file1, file2, file3 as x, y, z
                data1_columns, legends1, labels1 = parse_xvg(args.files[0])
                data2_columns, legends2, labels2 = parse_xvg(args.files[1])
                data3_columns, legends3, labels3 = parse_xvg(args.files[2])
                
                # Extract second columns
                y1 = data1_columns[1]
                y2 = data2_columns[1]
                y3 = data3_columns[1]
                
                # Validate that all files have the same number of data points
                if not (len(y1) == len(y2) == len(y3)):
                    print(f"Error: Files have different number of data rows!", file=sys.stderr)
                    print(f"  {args.files[0]}: {len(y1)} rows", file=sys.stderr)
                    print(f"  {args.files[1]}: {len(y2)} rows", file=sys.stderr)
                    print(f"  {args.files[2]}: {len(y3)} rows", file=sys.stderr)
                    return 1
                
                # Create 3D plot
                fig = plt.figure(figsize=tuple(args.figsize))
                ax = fig.add_subplot(111, projection='3d')
                
                if args.scatter:
                    # Scatter mode with color by order
                    colors = np.arange(len(y1))
                    scatter = ax.scatter(y1, y2, y3, c=colors, cmap=args.colormap,
                                       s=args.markersize*7, alpha=0.7)
                    cbar = fig.colorbar(scatter, ax=ax, pad=0.1)
                    cbar.set_label('Frame / Time order', rotation=270, labelpad=20)
                else:
                    # Line or marker plot in 3D
                    if args.style == 'dots':
                        ax.scatter(y1, y2, y3, s=args.markersize*7, alpha=0.7)
                    elif args.style == 'lines':
                        ax.plot(y1, y2, y3, '-', alpha=0.7)
                    elif args.style == 'lines+dots':
                        ax.plot(y1, y2, y3, 'o-', markersize=args.markersize, alpha=0.7)
                    else:
                        ax.plot(y1, y2, y3, args.style, alpha=0.7)
                
                # Set labels (use ylabel from each file)
                xlabel_default = labels1.get('ylabel', f'Values from {Path(args.files[0]).name}')
                ylabel_default = labels2.get('ylabel', f'Values from {Path(args.files[1]).name}')
                zlabel_default = labels3.get('ylabel', f'Values from {Path(args.files[2]).name}')
                
                ax.set_xlabel(args.xlabel if args.xlabel else xlabel_default)
                ax.set_ylabel(args.ylabel if args.ylabel else ylabel_default)
                ax.set_zlabel(zlabel_default)  # z-label doesn't have CLI override (could add if needed)
                ax.set_title(args.title if args.title else '3D Correlation')
                
                # Set equal aspect ratio for 3D if requested
                if args.aspect and args.aspect.lower() == 'equal':
                    # Get data ranges - convert to numpy arrays first
                    y1_arr = np.array(y1)
                    y2_arr = np.array(y2)
                    y3_arr = np.array(y3)
                    max_range = np.array([y1_arr.max()-y1_arr.min(), 
                                         y2_arr.max()-y2_arr.min(), 
                                         y3_arr.max()-y3_arr.min()]).max() / 2.0
                    mid_x = (y1_arr.max()+y1_arr.min()) * 0.5
                    mid_y = (y2_arr.max()+y2_arr.min()) * 0.5
                    mid_z = (y3_arr.max()+y3_arr.min()) * 0.5
                    ax.set_xlim(mid_x - max_range, mid_x + max_range)
                    ax.set_ylim(mid_y - max_range, mid_y + max_range)
                    ax.set_zlim(mid_z - max_range, mid_z + max_range)
                
                fig.tight_layout()
            
        elif len(args.files) == 1 and not args.multi:
            # Single file mode (deferred display)
            fig, ax = plot_xvg(args.files[0], show_moving_avg=args.moving_avg, window_size=args.window,
                               style=args.style, scatter_colormap=args.colormap, use_scatter=args.scatter,
                               use_histogram=args.histogram, hist_bins=args.bins, markersize=args.markersize)

            # Apply custom labels if provided (override XVG metadata)
            if args.title:
                ax.set_title(args.title)
            if args.xlabel:
                ax.set_xlabel(args.xlabel)
            if args.ylabel:
                ax.set_ylabel(args.ylabel)
            
            # Set aspect ratio if specified
            if args.aspect:
                if args.aspect.lower() == 'equal':
                    ax.set_aspect('equal', adjustable='box')
                elif args.aspect.lower() == 'auto':
                    ax.set_aspect('auto')
                else:
                    try:
                        ax.set_aspect(float(args.aspect), adjustable='box')
                    except ValueError:
                        print(f"Warning: Invalid aspect ratio '{args.aspect}', using auto", file=sys.stderr)

        else:
            # Multiple files mode
            fig, ax = plt.subplots(figsize=tuple(args.figsize))
            
            for i, file in enumerate(args.files):
                custom_legend = args.legends[i] if args.legends else None
                plot_xvg_multi(file, show_moving_avg=args.moving_avg, 
                             window_size=args.window, ax=ax, custom_legend=custom_legend,
                             style=args.style, scatter_colormap=args.colormap, 
                             use_scatter=args.scatter, use_histogram=args.histogram,
                             hist_bins=args.bins, markersize=args.markersize)
            
            # Apply custom labels if provided
            if args.title:
                ax.set_title(args.title)
            if args.xlabel:
                ax.set_xlabel(args.xlabel)
            if args.ylabel:
                ax.set_ylabel(args.ylabel)
            
            # Set aspect ratio if specified
            if args.aspect:
                if args.aspect.lower() == 'equal':
                    ax.set_aspect('equal', adjustable='box')
                elif args.aspect.lower() == 'auto':
                    ax.set_aspect('auto')
                else:
                    try:
                        ax.set_aspect(float(args.aspect), adjustable='box')
                    except ValueError:
                        print(f"Warning: Invalid aspect ratio '{args.aspect}', using auto", file=sys.stderr)
            
            plt.tight_layout()
        
        # Save or show the plot
        if args.output:
            fig.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
            print(f"Plot saved to: {args.output}")
        else:
            # Use the pyplot-level show which blocks and opens a window
            # consistently across backends (Qt5Agg, TkAgg, etc.).
            plt.show()
        
        return 0
    
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())

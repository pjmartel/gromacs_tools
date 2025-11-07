#!/usr/bin/env python3
"""
Command-line tool to plot GROMACS XVG files with Matplotlib.
Supports single or multiple XVG files, moving averages, and customization options.
"""
import argparse
import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

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
             scatter_colormap='viridis', use_scatter=False):

    data_columns, legends, labels = parse_xvg(filename)
    x = data_columns[0]
    num_datasets = len(data_columns) - 1

    plt.figure(figsize=(10, 6))

    for i in range(num_datasets):
        y = data_columns[i + 1]
        label = legends.get(i, f'Dataset {i}')
        
        if use_scatter:
            # Scatter mode: color by order (time-like)
            colors = np.arange(len(x))
            scatter = plt.scatter(x, y, c=colors, cmap=scatter_colormap, 
                                s=20, alpha=0.7, label=label)
            if num_datasets == 1:  # Only add colorbar for single dataset
                cbar = plt.colorbar(scatter)
                cbar.set_label('Frame / Time order', rotation=270, labelpad=20)
        else:
            # Line/marker mode
            if style == 'dots':
                plt.plot(x, y, 'o', markersize=3, label=label)
            elif style == 'lines':
                plt.plot(x, y, '-', label=label)
            elif style == 'lines+dots':
                plt.plot(x, y, 'o-', markersize=3, label=label)
            else:  # allow custom matplotlib format strings
                plt.plot(x, y, style, label=label)

        # Apply moving average only if there's a single dataset and the user requested it
        if show_moving_avg and num_datasets == 1 and not use_scatter:
            y_avg = moving_average(y, window_size)
            x_avg = x[:len(y_avg)]  # match lengths
            plt.plot(x_avg, y_avg, label=f'{label} (Moving Avg, window={window_size})', 
                    linestyle='--', linewidth=2)

    plt.xlabel(labels.get('xlabel', 'X-axis'))
    plt.ylabel(labels.get('ylabel', 'Y-axis'))
    plt.title(labels.get('title', ''))
    if num_datasets > 1 or legends:
        plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# This version supports chaining multiple plots on the same axes, based on a
# shared ax variable (should replace the plot_xvg, but it is not completely tested)
def plot_xvg_multi(filename, show_moving_avg=False, window_size=10, ax=None, 
                   custom_legend=None, style='dots', scatter_colormap='viridis', 
                   use_scatter=False):

    data_columns, legends, labels = parse_xvg(filename)
    x = data_columns[0]
    num_datasets = len(data_columns) - 1

    # Create new Axes if none provided
    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 6))

    for i in range(num_datasets):
        y = data_columns[i + 1]
        #label = legends.get(i, f'Dataset {i}')
        label = custom_legend if custom_legend and num_datasets == 1 else legends.get(i, f'Dataset {i}')
        
        if use_scatter:
            # Scatter mode: color by order (time-like)
            colors = np.arange(len(x))
            scatter = ax.scatter(x, y, c=colors, cmap=scatter_colormap, 
                               s=20, alpha=0.7, label=label)
            # Note: colorbar is tricky with shared axes, user should add manually if needed
        else:
            # Line/marker mode
            if style == 'dots':
                ax.plot(x, y, 'o', markersize=3, label=label)
            elif style == 'lines':
                ax.plot(x, y, '-', label=label)
            elif style == 'lines+dots':
                ax.plot(x, y, 'o-', markersize=3, label=label)
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
    
    parser.add_argument('--dpi', type=int, default=300,
                        help='DPI for saved figures (default: 300)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if args.legends and len(args.legends) != len(args.files):
        parser.error(f"Number of legends ({len(args.legends)}) must match number of files ({len(args.files)})")
    
    if args.window < 1:
        parser.error("Window size must be at least 1")
    
    # Check that all files exist
    for file in args.files:
        if not Path(file).exists():
            print(f"Error: File not found: {file}", file=sys.stderr)
            return 1
    
    try:
        if len(args.files) == 1 and not args.multi:
            # Single file mode
            plot_xvg(args.files[0], show_moving_avg=args.moving_avg, window_size=args.window,
                    style=args.style, scatter_colormap=args.colormap, use_scatter=args.scatter)
            
            # Apply custom labels if provided
            if args.title:
                plt.title(args.title)
            if args.xlabel:
                plt.xlabel(args.xlabel)
            if args.ylabel:
                plt.ylabel(args.ylabel)
            
        else:
            # Multiple files mode
            fig, ax = plt.subplots(figsize=tuple(args.figsize))
            
            for i, file in enumerate(args.files):
                custom_legend = args.legends[i] if args.legends else None
                plot_xvg_multi(file, show_moving_avg=args.moving_avg, 
                             window_size=args.window, ax=ax, custom_legend=custom_legend,
                             style=args.style, scatter_colormap=args.colormap, 
                             use_scatter=args.scatter)
            
            # Apply custom labels if provided
            if args.title:
                ax.set_title(args.title)
            if args.xlabel:
                ax.set_xlabel(args.xlabel)
            if args.ylabel:
                ax.set_ylabel(args.ylabel)
            
            plt.tight_layout()
        
        # Save or show the plot
        if args.output:
            plt.savefig(args.output, dpi=args.dpi, bbox_inches='tight')
            print(f"Plot saved to: {args.output}")
        else:
            plt.show()
        
        return 0
    
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())

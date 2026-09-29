# gmxl

GROMACS command wrapper with automatic logging.

## Overview

`gmxl` is a drop-in replacement for `gmx` for commands you run by hand. It runs
`gmx <subcommand> [args...]` exactly as given, shows the output in the terminal as usual, and
for commands that succeed:

- saves the full output (stdout + stderr) to `<subcommand>_<timestamp>.log`
- appends the command line, with a timestamp, to a running `GMX_COMMANDS` file

This keeps a record of the manual analyses done in a directory, the same way the other tools
in this repository record their commands in `<tool>_commands.sh`.

## Features

- ✅ **Transparent**: same arguments as `gmx`, output still shown live, exit status preserved
- ✅ **Interactive prompts still work** (e.g. group selection menus), and piped input too
- ✅ **Per-command output logs**, named by subcommand and time
- ✅ **Running command history** in `GMX_COMMANDS`
- ✅ **Only successful commands are recorded**; a failed command leaves no log file and no history entry

## Usage

### Basic Syntax

```bash
gmxl <subcommand> [args...]
```

Run without arguments to print a short help message. `gmx` must be on `PATH` (e.g. after
sourcing `GMXRC`).

### Files Written

Both files are written to the current directory:

| File | Content |
|------|---------|
| `<subcommand>_<YYYYmmdd_HHMMSS>.log` | Everything the command printed (stdout and stderr) |
| `GMX_COMMANDS` | Two lines per successful command: the command line, and where its output was saved |

Example `GMX_COMMANDS`:

```
[2026-09-29 14:02:11] gmx rms -s md.tpr -f md.xtc -o rmsd.xvg
[2026-09-29 14:02:11] ✓ Success - Output saved to: rms_20260929_140211.log
[2026-09-29 14:03:40] gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg
[2026-09-29 14:03:40] ✓ Success - Output saved to: gyrate_20260929_140340.log
```

## Examples

```bash
# RMSD (answer the group prompts interactively as usual)
gmxl rms -s md.tpr -f md.xtc -o rmsd.xvg

# Radius of gyration
gmxl gyrate -s md.tpr -f md.xtc -o gyrate.xvg

# Energy term, with the selection piped in
echo "Temperature" | gmxl energy -f md.edr -o temperature.xvg

# Review what was run in this directory
cat GMX_COMMANDS
```

## Notes

- `GMX_COMMANDS` is a history log, not a runnable script: lines start with a timestamp, and
  interactive group choices aren't recorded. For a reproducible script, pipe the selections
  in (as in the `energy` example above) and copy the logged command line.
- The exit status is that of `gmx`, so `gmxl` works in scripts and `&&` chains.
- To use it from anywhere, add `bin/` to `PATH` (see the main [README](../README.md)).

# Tab completion

Bash Tab completion for gromacs_tools commands: type `plot_xvg.py --st` and press Tab to get
`--start --stats --style`, `plot_xvg.py --style <Tab>` for `dots lines lines+dots`, and so on.

Currently available for: `plot_xvg.py`.

## How it works

1. `generate.py` reads a tool's command-line options from its argparse parser and writes a
   bash completion file, `<tool>.bash`, using [shtab](https://github.com/iterative/shtab).
   It can add value suggestions argparse doesn't define itself (e.g. Matplotlib colormaps for
   `--colormap`).
2. Each `<tool>.bash` defines a completion function and registers it with bash's `complete`
   builtin for the command name (e.g. `complete -F _shtab_plot_xvg_py plot_xvg.py`).
3. The shell loads every `completions/*.bash` file at startup (in the dotfiles setup:
   `bash/common.sh`). Without dotfiles, add this to `~/.bashrc`:

   ```bash
   for f in ~/github/gromacs_tools/completions/*.bash; do . "$f"; done
   ```

Pressing Tab only runs that small bash function, never the tool itself, so completion is
instant. Completion is tied to the command name, so it works for `plot_xvg.py ...` (the tool on
`PATH`), not for `python bin/plot_xvg.py ...`.

## After changing a tool's options

Regenerate, and commit the updated `.bash` file together with the change (requires
[uv](https://docs.astral.sh/uv/)):

```bash
completions/generate.py           # rewrite the .bash files that changed
completions/generate.py --check   # exit status 1 if any .bash file is out of date
```

Open a new terminal (or `source ~/.bashrc`) to use the new completion.

## Adding another tool

In `generate.py`, write a function returning the tool's configured `ArgumentParser` (with
`parser.prog` set to the command name, file arguments marked with `action.complete = shtab.FILE`)
and add it to `TOOLS`. Tools that build their parser inside `main()` first need it moved into a
separate function, as in `plot_xvg.py`'s `setup_argument_parser()`. Bash tools (which parse their
options by hand) need a different approach.

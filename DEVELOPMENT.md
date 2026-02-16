# Development Setup Guide

This guide explains how to set up a development environment for gromacs_tools.

## Quick Setup (Recommended)

Run the automated setup script:

```bash
./setup_dev_env.sh
```

This will:
1. Create a virtual environment in `.venv/`
2. Install all dependencies from `requirements.txt`
3. Activate the environment

## Manual Setup

If you prefer to set up manually:

### 1. Create Virtual Environment

```bash
python3 -m venv .venv
```

### 2. Activate Virtual Environment

**On Linux/Mac:**
```bash
source .venv/bin/activate
```

**On Windows:**
```bash
.venv\Scripts\activate
```

### 3. Install Dependencies

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

## Working with the Virtual Environment

### Activating

Each time you open a new terminal session:
```bash
source .venv/bin/activate
```

You'll see `(.venv)` prefix in your prompt.

### Deactivating

When you're done:
```bash
deactivate
```

### Using with VSCode

VSCode should automatically detect the `.venv` folder. You can:

1. **Select the interpreter**: 
   - Press `Ctrl+Shift+P` (or `Cmd+Shift+P` on Mac)
   - Type "Python: Select Interpreter"
   - Choose `.venv/bin/python`

2. **Auto-activation**: VSCode will automatically activate the virtual environment in integrated terminals

## Using with AI Coding Assistants (GitHub Copilot, etc.)

When working with AI assistants in VSCode:

1. **Create local venv** (already gitignored):
   ```bash
   ./setup_dev_env.sh
   ```

2. **Run Python code with venv**:
   The assistant should use:
   ```bash
   source .venv/bin/activate && python script.py
   ```
   Or directly:
   ```bash
   .venv/bin/python script.py
   ```

3. **Conda users**: If you prefer conda:
   ```bash
   conda create -n gromacs_tools python=3.9
   conda activate gromacs_tools
   pip install -r requirements.txt
   ```

## Why Virtual Environments?

- **Isolation**: Dependencies don't conflict with system Python or other projects
- **Reproducibility**: Everyone uses the same package versions
- **Safety**: Can be deleted/recreated without affecting system
- **Git-friendly**: `.venv/` is gitignored, won't bloat repository

## Installing New Packages

When adding new dependencies:

1. Install in your venv:
   ```bash
   pip install package-name
   ```

2. Update requirements.txt:
   ```bash
   pip freeze > requirements.txt
   ```
   
   Or manually add to `requirements.txt` with version constraints:
   ```
   package-name>=1.0.0
   ```

3. Commit the updated `requirements.txt` (but NOT the `.venv/` folder)

## Troubleshooting

### "python3: command not found"
Install Python 3.8+ from python.org or your package manager.

### "pip: command not found"
```bash
python3 -m ensurepip --upgrade
```

### VSCode doesn't detect .venv
1. Reload window: `Ctrl+Shift+P` → "Developer: Reload Window"
2. Manually select interpreter: `Ctrl+Shift+P` → "Python: Select Interpreter"

### Packages not found when running scripts
Make sure the venv is activated:
```bash
source .venv/bin/activate
python script.py
```

## Repository Structure

```
gromacs_tools/
├── .venv/              # Local virtual environment (gitignored)
├── .gitignore          # Includes .venv/, venv/, etc.
├── requirements.txt    # Python dependencies (committed)
├── setup_dev_env.sh    # Setup script (committed)
├── bin/                # Tools
├── docs/               # Documentation
└── tests/              # Test files
```

## Best Practices

✅ **DO**:
- Create `.venv/` in project root
- Activate venv before working
- Update `requirements.txt` when adding packages
- Commit `requirements.txt` changes

❌ **DON'T**:
- Commit `.venv/` folder to git
- Install packages globally for project-specific work
- Forget to activate venv before running scripts
- Share compiled Python files (`.pyc`, `__pycache__/`)

## Alternative: Conda Environment

If you use conda instead:

```bash
# Create environment
conda create -n gromacs_tools python=3.9

# Activate
conda activate gromacs_tools

# Install dependencies
pip install -r requirements.txt

# Deactivate
conda deactivate
```

Note: Conda environments are stored globally (not in project folder), which is also fine! They won't be committed to git.

---
name: python-virtual-env
description: Select and activate the correct existing Python environment before running Python tools, tests, or build scripts. Prefer Conda environments under ~/miniconda3; ask the user when the application does not identify a reliable environment.
---

# Select a Python environment

Choose and activate the correct existing Python environment before running project Python commands.
Do not assume that `python`, `pip`, or `pytest` on the default `PATH` belongs to the application being
worked on.

Use this skill when running Python tools, tests, documentation builders, or build-maintenance scripts,
and when diagnosing imports or dependency versions.
Skip environment discovery only when the user already supplied an absolute interpreter or the active
environment has been verified for the current command session.

## Selection policy

Use evidence in this order:

1. An environment or interpreter explicitly named by the user.
2. Repository instructions such as `AGENTS.md`, `README.md`, environment files, or earlier context in
   the current task.
3. An already active environment whose interpreter and dependencies have been verified.
4. The names and package contents of existing Conda environments.

Never choose an environment merely because its Python version looks plausible.
Names such as `env3.14`, `buildbot_3.14`, and `mkdocs` may target different applications despite
sharing packages or Python versions.

If more than one environment remains plausible and using the wrong one could change results, ask the
user which environment to use.
Do not create a new environment or install packages just to avoid this question.

## Common pitfall

Many systems have multiple Python installations:

* system Python
* Conda base environment
* one or more project environments
* virtual environments (venv)

Running

python script.py

or

pip install package

may use the wrong interpreter.

Determine the intended environment before executing environment-dependent commands.

## Discover existing Conda environments

The user normally installs Miniconda under `$HOME/miniconda3`.
In a non-interactive shell, initialize Conda with the known installation before querying it:

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda env list
```

If `$HOME/miniconda3` does not exist, locate Conda with `command -v conda` or ask the user where it is
installed.
Do not run `conda init`; it modifies shell configuration and is unnecessary for an individual command.

The active environment is marked with `*`.
Environment names are clues, not proof of application compatibility.

## Activate and run

Activation and the dependent command must execute in the same shell invocation because activation does
not persist across independent agent command sessions:

```bash
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate ENVIRONMENT_NAME
python --version
python -c 'import sys; print(sys.executable)'
python path/to/tool.py
```

For ABINIT, use `env3.14` only when the user or repository context identifies it as appropriate:

```bash
source /Users/giantomassi/miniconda3/etc/profile.d/conda.sh
conda activate env3.14
```

`conda run` is an alternative when shell activation is impractical:

```bash
conda run -n ENVIRONMENT_NAME python path/to/tool.py
```

Do not silently switch from an explicitly requested environment because activation fails.
Report the failure and investigate the environment name or Conda installation.

## Verify suitability

Before a long or mutating workflow, verify both the interpreter and the important application import:

```bash
which python
python --version
python -c 'import sys; print(sys.executable)'
python -c 'import REQUIRED_PACKAGE; print(REQUIRED_PACKAGE.__file__)'
```

Use `python -m pytest` when it is important that pytest belongs to the selected interpreter.
For project CLIs, check whether the executable resolves inside the active environment.

## Installing packages

Install or update packages only when the user requested it or when installation is an authorized,
necessary part of the task.
Prefer Conda when the required package is available in the project's established channels:

```bash
conda install package_name
```

when the package is available from Conda.

Otherwise use

```bash
python -m pip install package_name
```

instead of

```bash
pip install package_name
```

Using

```bash
python -m pip
```

guarantees that pip belongs to the currently selected Python interpreter.

Never install into Conda `base` by default.
Avoid mixing Conda and pip installations without checking how the environment is already managed.

## Exporting or changing an environment

Environment creation, export, and bulk updates are separate mutations.
Perform them only when explicitly requested.

Export an environment with:

To reproduce an environment:

```bash
conda env export > environment.yml
```

To recreate it:

```bash
conda env create -f environment.yml
```

Update installed packages with:

```bash
conda update --all
```

or update a single package with

```bash
conda update package_name
```

Bulk `conda update --all` can materially change a working application environment.
Do not use it as routine troubleshooting.

## Completion information

When environment selection affected the task, report the environment name and interpreter path.
If the required environment could not be identified, state the candidates and ask one concise question
instead of repeatedly trying environments.

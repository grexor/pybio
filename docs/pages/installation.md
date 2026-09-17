# Installation

The easiest way to install **pybio** is by running:

```bash
pip install pybio
```

!!! note
    On some systems, **pip** installs the executable scripts under `~/.local/bin`. However this folder is not in the `PATH`, which will result in `command not found` if you try to run `pybio` on the command line. To fix this, run `export PATH="$PATH:~/.local/bin"` (and add this to your `.profile`). Another option is to install inside a virtual environment (using `virtualenv`).

If you would instead like to **install the latest development version** from this repository:

```bash
# clone pybio GitHub repository
git clone https://github.com/grexor/pybio.git

# install
pip install .
```

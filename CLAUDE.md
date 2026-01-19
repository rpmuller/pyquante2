# CLAUDE.md

This file provides guidance for Claude Code when working with the pyquante2 codebase.

## Project Overview

PyQuante2 is a Python quantum chemistry suite for developing and testing quantum chemistry methods. It includes implementations of Hartree-Fock (RHF, UHF), DFT, and MP2, with optional Cython extensions for performance-critical routines.

**Status:** Published on PyPI as `pyquante2` (current version: 0.2.0)

## Project Structure

```
pyquante2/
├── src/pyquante2/       # Main source code
│   ├── basis/           # Gaussian basis sets (PGBF, CGBF)
│   ├── cints/           # Cython/C integral routines
│   ├── dft/             # Density functional theory
│   ├── geo/             # Molecular geometry (atoms, molecules)
│   ├── graphics/        # Visualization tools
│   ├── grid/            # Numerical grids for DFT
│   ├── hylleraas/       # Hylleraas variational methods
│   ├── ints/            # Integral evaluation (one/two-electron)
│   ├── pt/              # Perturbation theory (MP2)
│   ├── scf/             # SCF solvers and Hamiltonians
│   └── viewer/          # Molecular viewer
├── tests/               # Test suite
├── attic/               # Archived/deprecated code
├── docs/                # Documentation
├── CHANGELOG.md         # Release history
├── VERSION_MANAGEMENT.md # Version bump checklist
├── FUTURE_TASKS.md      # Planned improvements
├── CODE_REVIEW.md       # Code quality findings
└── TODO.md              # Original task list
```

## Development Commands

```bash
# Install with dev dependencies
uv sync --extra dev

# Run tests (preferred)
make test

# Or directly with pytest
uv run pytest

# Run tests with doctests
uv run pytest --doctest-modules --ignore=attic

# Build package for PyPI
make clean          # Clean old builds
make build          # Build wheel + source distribution
uv run twine check dist/*  # Verify package

# Publish to PyPI (see VERSION_MANAGEMENT.md first!)
make install_test   # Upload to TestPyPI
make install        # Upload to PyPI
```

## Testing

- Tests are in `tests/` directory using pytest
- Many modules include doctests in docstrings
- Run `make test` to run all tests including doctests
- 3 DFT tests are skipped (conditional functionality)

## Code Style

- Python 3.11+ required
- Uses numpy extensively for numerical operations
- Classes use lowercase names (e.g., `molecule`, `atom`, `cgbf`)
- Doctests serve as both documentation and tests
- No type annotations currently used

## Key Concepts

- **PGBF**: Primitive Gaussian Basis Function
- **CGBF**: Contracted Gaussian Basis Function
- **HGP**: Head-Gordon-Pople integral algorithm
- **RHF/UHF**: Restricted/Unrestricted Hartree-Fock
- **DFT**: Density Functional Theory (LDA, GGA functionals)

## Dependencies

- Core: numpy, ipykernel
- Dev/Test: pytest, matplotlib, build, twine
- Build: setuptools, Cython (for C extensions)

## Notes

- Cython extensions in `src/pyquante2/cints/` are optional but improve performance
- The `attic/` directory contains archived code - ignore when searching
- Uses `uv` for package management (uv.lock present)
- **Before version bumps:** Review `VERSION_MANAGEMENT.md` checklist
- **For future tasks:** See `FUTURE_TASKS.md`
- **PyPI package:** Available at https://pypi.org/project/pyquante2/

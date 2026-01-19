# PyQuante2 Future Tasks

This file tracks upcoming tasks and improvements to work on with Claude Code.

## PyPI & Distribution

- [ ] **Cross-platform wheel building with GitHub Actions**
  - Set up cibuildwheel in GitHub Actions
  - Build wheels for Linux (x86_64, aarch64)
  - Build wheels for macOS (Intel + ARM64)
  - Build wheels for Windows (x86_64)
  - Support Python 3.11, 3.12, 3.13
  - See: https://cibuildwheel.readthedocs.io/

- [ ] **Create project-specific PyPI tokens**
  - After first release, switch from account-wide to project-specific tokens
  - Update .pypirc with new tokens

## Code Quality (from CODE_REVIEW.md)

### High Priority

- [ ] **Fix skipped DFT tests** (Issue #7)
  - `tests/test_dft.py:30` - CVWN test consistently fails
  - `tests/test_dft.py:77` - CLYP doesn't work
  - `tests/test_dft.py:94` - CPBE doesn't work
  - Either fix functionals or document as unsupported

- [ ] **Add missing test coverage** (Issue #8)
  - `src/pyquante2/graphics/` - visualization modules
  - `src/pyquante2/grid/` - grid generation
  - `src/pyquante2/io/` - CML parsing
  - `src/pyquante2/viewer/` - molecular viewer
  - `src/pyquante2/hylleraas/` - Hylleraas methods

### Medium Priority

- [ ] **Add type hints** (Issue #11)
  - Start with core modules:
    - `src/pyquante2/utils.py`
    - `src/pyquante2/geo/molecule.py`
    - `src/pyquante2/geo/atom.py`
    - `src/pyquante2/basis/pgbf.py`
    - `src/pyquante2/basis/cgbf.py`
    - `src/pyquante2/scf/iterators.py`
  - Use gradual typing (don't require 100% coverage)
  - Add mypy to dev dependencies

- [ ] **Documentation improvements** (Issue #12)
  - Add module-level docstrings to all modules
  - Add NumPy-style docstrings to public functions
  - Add `__all__` to `__init__.py` files
  - Create API documentation with Sphinx

### Low Priority

- [ ] **Remove bare returns in __init__** (Issue #9)
  - `src/pyquante2/geo/molecule.py:52`
  - `src/pyquante2/basis/basisset.py:47`
  - `src/pyquante2/geo/atom.py:26`
  - `src/pyquante2/basis/cgbf.py:42`

- [ ] **Clean up dead/commented code** (Issue #14)
  - `src/pyquante2/utils.py:52-65` - commented scipy code
  - `src/pyquante2/ints/integrals.py:23-81` - old `twoe_integrals_compressed`

- [ ] **Dependency cleanup** (Issue #13)
  - Consider moving `ipykernel` to optional dependencies only
  - Add version constraints to numpy if needed

## Feature Enhancements (from TODO.md)

### DFT Improvements

- [ ] Get cvwn working
- [ ] Get clyp working
- [ ] Test cpbe
- [ ] Wrap libxc for more functionals
- [ ] Fix Slater exchange high error (2e-6 in test suite)
- [ ] Fix cvwn bug with dfb value at nb=0

### Integral Routines

- [ ] Test routine for all Python and C integral routines for contracted bfs
- [ ] Add rys/crys support
- [ ] Add libint support

### UI/UX Improvements

- [ ] Units for `molecule._repr_html_`
- [ ] `__repr__` and `_repr_html_` for basisset
- [ ] Use numpy in shapes/viewer stuff

### Code Modernization

- [ ] More use of itertools
- [ ] More use of NamedTuples
- [ ] More Cython optimization
- [ ] More einsum usage
- [ ] Settings to ConfigParser

### Performance

- [ ] Profile and optimize hot paths
- [ ] Is there redundancy in Becke reweighting? (loops over atoms in Ps loop AND becke_atomic_...)

### Architecture Questions

- [ ] Move xyz readers into IO module?
- [ ] Should we require scipy for incomplete gamma functions?
- [ ] Should we require scipy for Legendre/Lebedev stuff?

## Breaking API Changes (Consider for v1.0.0)

- [ ] **Rename classes to PEP 8 conventions** (Issue #10)
  - `atom` → `Atom`
  - `molecule` → `Molecule`
  - `pgbf` → `PGBF`
  - `cgbf` → `CGBF`
  - `basisset` → `BasisSet`
  - Note: This is a breaking change - plan migration strategy
  - Consider deprecation warnings in 0.x versions first

## Release Planning

### v0.3.0 (Next minor release)
- Focus on fixing DFT tests
- Add type hints to core modules
- Improve documentation

### v0.4.0
- Cross-platform wheel building
- Expanded test coverage

### v1.0.0 (Stable API)
- Complete API documentation
- PEP 8 class naming (breaking changes)
- Full test coverage
- Performance benchmarks

---

## How to Use This File

1. **Pick a task** - Choose from high priority items or areas of interest
2. **Start a Claude Code session** - Reference this file: "Let's work on <task from FUTURE_TASKS.md>"
3. **Check off completed items** - Use `[x]` when done
4. **Add new tasks** - As you discover improvements, add them here
5. **Review regularly** - Update priorities as the project evolves

---

Last updated: 2026-01-19

# PyQuante2 Code Review

## Summary

| Category | Count | Severity |
|----------|-------|----------|
| ~~Bare exception handlers~~ | ~~9~~ | ~~High~~ DONE |
| Missing type hints | All modules | Medium |
| Test coverage gaps | Multiple modules | High |
| Print statements (should use logging) | 20+ | Medium |
| ~~Mutable default arguments~~ | ~~3~~ | ~~Medium~~ DONE |
| Documentation gaps | Multiple | Medium |
| API design inconsistencies | 5 | Low |

---

## Issues to Address

### 1. ~~Bare Exception Handlers~~ (DONE)

- [x] `src/pyquante2/__init__.py:19` - `except ImportError:`
- [x] `src/pyquante2/ints/integrals.py:6,12` - `except ImportError:`
- [x] `src/pyquante2/grid/grid.py:12` - `except ImportError:`
- [x] `src/pyquante2/graphics/lineplot.py:4` - `except ImportError:`
- [x] `src/pyquante2/graphics/contourplot.py:6` - `except ImportError:`
- [x] `src/pyquante2/viewer/trackball_camera.py:89` - `except ImportError:`
- [x] `src/pyquante2/viewer/viewer.py:12` - `except ImportError:`
- [x] `src/pyquante2/scf/mcscf.py:70` - `except TypeError:`
- [x] `src/pyquante2/dft/reference.py:3621` - `except StopIteration:`

---

### 2. ~~Mutable Default Arguments~~ (DONE)

- [x] `src/pyquante2/geo/molecule.py:39` - `atomlist=None`
- [x] `src/pyquante2/basis/cgbf.py:25` - `exps=None, coefs=None`

---

### 3. ~~Resource Leak in File I/O~~ (DONE)

- [x] `src/pyquante2/geo/molecule.py:214-220` - now uses `with open(fname) as f:`

---

### 4. ~~Deprecated Type Checking~~ (DONE)

- [x] `src/pyquante2/io/cml.py:271` - now uses `isinstance(inp, str)`

---

### 5. ~~Wildcard Imports~~ (DONE)

- [x] `src/pyquante2/__init__.py:6` - `from pyquante2.geo.samples import *` - already fixed with explicit imports
- [x] `src/pyquante2/viewer/viewer.py:10` - `from pyglet.gl import *` - replaced with explicit GL imports
- [x] `src/pyquante2/viewer/trackball_camera.py:88` - `from pyglet.gl import *` - replaced with specific imports

**Fix:** Import specific names or add `__all__` to samples module ✅

---

### 6. ~~Replace Print Statements with Logging~~ (DONE)

- [x] `src/pyquante2/scf/hamiltonians.py:122` - replaced with `logger.info()`
- [x] `src/pyquante2/grid/grid.py:13` - already fixed with `logger.debug()`
- [x] `src/pyquante2/ints/integrals.py:7,13` - already fixed with `logger.debug()`
- [x] `src/pyquante2/utils.py:103,130` - no print statements found (already removed or not present)
- [x] Multiple instances in `src/pyquante2/scf/mcscf.py` - replaced 16+ print statements with logging

**Fix:** ✅ All print statements replaced with logging
```python
import logging
logger = logging.getLogger(__name__)
logger.warning("Couldn't find cython becke routine")
```

---

### 7. Skipped/Broken Tests (High Priority)

- [ ] `tests/test_dft.py:30` - CVWN test consistently fails
- [ ] `tests/test_dft.py:77` - CLYP doesn't work
- [ ] `tests/test_dft.py:94` - CPBE doesn't work

**Fix:** Either fix the functionals or document as unsupported

---

### 8. Missing Test Coverage (Medium Priority)

Modules without dedicated tests:
- [ ] `src/pyquante2/graphics/` - visualization modules
- [ ] `src/pyquante2/grid/` - grid generation (partial)
- [ ] `src/pyquante2/io/` - CML parsing
- [ ] `src/pyquante2/viewer/` - molecular viewer
- [ ] `src/pyquante2/hylleraas/` - Hylleraas methods

---

### 9. Bare Return in __init__ Methods (Low Priority)

- [ ] `src/pyquante2/geo/molecule.py:52`
- [ ] `src/pyquante2/basis/basisset.py:47`
- [ ] `src/pyquante2/geo/atom.py:26`
- [ ] `src/pyquante2/basis/cgbf.py:42`

**Fix:** Remove unnecessary `return` statements from `__init__`

---

### 10. Inconsistent Class Naming (Low Priority)

Classes use lowercase (non-PEP 8):
- [ ] `atom` -> `Atom`
- [ ] `molecule` -> `Molecule`
- [ ] `pgbf` -> `PGBF`
- [ ] `cgbf` -> `CGBF`
- [ ] `basisset` -> `BasisSet`

**Note:** This is a breaking API change - consider carefully

---

### 11. Missing Type Hints (Medium Priority)

No type hints in any module. Priority files:
- [ ] `src/pyquante2/utils.py` - utility functions
- [ ] `src/pyquante2/geo/molecule.py` - core data structures
- [ ] `src/pyquante2/geo/atom.py` - core data structures
- [ ] `src/pyquante2/basis/pgbf.py` - basis functions
- [ ] `src/pyquante2/basis/cgbf.py` - basis functions
- [ ] `src/pyquante2/scf/iterators.py` - SCF algorithms

---

### 12. Documentation Improvements (Low Priority)

- [ ] Update README.md (references Python 2.7, outdated install instructions)
- [ ] Add module-level docstrings to all modules

- [ ] Add NumPy-style docstrings to public functions
- [ ] Add `__all__` to `__init__.py`



---

### 13. Dependency Cleanup (Low Priority)

- [ ] `pyproject.toml` - Remove `pip` from dependencies
- [ ] `pyproject.toml` - Consider moving `ipykernel` to optional
- [ ] Add version constraints to numpy

---

### 14. Dead/Commented Code (Low Priority)

- [ ] `src/pyquante2/utils.py:52-65` - commented scipy code
- [ ] `src/pyquante2/ints/integrals.py:23-81` - `twoe_integrals_compressed` marked as "old"

---

## Positive Observations

- Core quantum chemistry algorithms are correct and well-tested
- Good use of doctests for documentation and testing
- Clean separation between geo, basis, scf, dft modules
- Proper use of numpy vectorized operations and einsum
- Cython extensions provide good performance where needed

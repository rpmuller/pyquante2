# Version Management Reminder

## IMPORTANT: Always Update Version Numbers Before Release!

This document serves as a reminder to update version numbers across all relevant files before building and releasing to PyPI.

## Critical Files to Update

### 1. `pyproject.toml` (PRIMARY SOURCE OF TRUTH)
**Location**: Line 7
```toml
version = "0.2.0"  # <-- UPDATE THIS!
```

### 2. `CHANGELOG.md`
- Update the version number in the latest release section
- Replace "2026-01-XX" with the actual release date
- Move items from `[Unreleased]` to the new version section

Example:
```markdown
## [0.2.0] - 2026-01-19  # <-- UPDATE DATE!
```

## Semantic Versioning Guide

Follow [Semantic Versioning](https://semver.org/): **MAJOR.MINOR.PATCH**

- **MAJOR** (1.0.0): Breaking changes, incompatible API changes
- **MINOR** (0.2.0): New features, backward-compatible functionality
- **PATCH** (0.2.1): Bug fixes, backward-compatible fixes

### Current Status: 0.2.0 (Alpha/Development)
- Use 0.x.y for development/pre-release versions
- Move to 1.0.0 when API is stable and production-ready

## Pre-Release Checklist

Before running `make build`:

- [ ] Update version in `pyproject.toml`
- [ ] Update version in `CHANGELOG.md`
- [ ] Update release date in `CHANGELOG.md` (replace "2026-01-XX")
- [ ] Add release notes describing changes
- [ ] Commit version changes: `git commit -m "Bump version to X.Y.Z"`
- [ ] Run `make clean` to remove old builds
- [ ] Run `make build` to create new distribution
- [ ] Run `twine check dist/*` to verify metadata
- [ ] Tag the release: `git tag vX.Y.Z`
- [ ] Push tags: `git push --tags`

## Common Version Bump Scenarios

### Bug Fix (Patch)
0.2.0 → 0.2.1
- Fixed calculation errors
- Fixed typos in documentation
- Performance improvements without API changes

### New Feature (Minor)
0.2.0 → 0.3.0
- Added new basis sets
- New calculation methods
- New convenience functions

### Breaking Change (Major)
0.2.0 → 1.0.0
- Changed function signatures
- Removed deprecated features
- Restructured package layout

## Post-Release Checklist

After successful PyPI upload:

- [ ] Verify package appears on PyPI: https://pypi.org/project/pyquante2/
- [ ] Test installation: `pip install pyquante2==X.Y.Z`
- [ ] Update CHANGELOG.md with new `[Unreleased]` section
- [ ] Commit and push changes

## Quick Command Reference

```bash
# View current version
grep '^version' pyproject.toml

# Clean build artifacts
make clean

# Build package
make build

# Check distribution
twine check dist/*

# Upload to TestPyPI first (recommended)
make install_test

# Upload to PyPI (after testing)
make install

# Tag release
git tag v0.2.0
git push --tags
```

## Automation Option (Future)

Consider using tools to automate version management:
- `bump2version` / `bumpversion`: Automatically update version numbers
- `setuptools-scm`: Derive version from git tags
- GitHub Actions: Automate releases on tag push

## WARNING: Version Number Mistakes

**Common mistakes to avoid:**
1. Forgetting to update version before building
2. Reusing the same version number (PyPI won't allow re-upload)
3. Uploading to PyPI with wrong version
4. Not updating CHANGELOG.md
5. Not tagging releases in git

**If you make a mistake:**
- You CANNOT delete or replace a PyPI release
- You MUST bump to a new version (e.g., 0.2.1) and re-upload
- This is why testing on TestPyPI first is critical!

## Current Version: 0.2.0

Last updated: 2026-01-19

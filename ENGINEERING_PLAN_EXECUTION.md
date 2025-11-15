# Sheriff Professionalization: Engineering Execution Plan
## Complete Implementation Roadmap

**Date:** 2025-11-15
**Status:** IN PROGRESS
**Goal:** Execute professionalization to scvi-tools standards

---

## Phase 1: CI/CD Infrastructure (HIGH PRIORITY)

### 1.1 GitHub Actions Workflows

**Files to Create:**

#### `.github/workflows/tests.yml` - Main Test Suite
```yaml
Purpose: Run Python + Rust tests on every push/PR
Platforms: Ubuntu, macOS (Windows future)
Python versions: 3.10, 3.11, 3.12
Includes: Code coverage upload to Codecov
```

#### `.github/workflows/rust-tests.yml` - Rust-Specific Tests
```yaml
Purpose: Dedicated Rust testing (cargo test, clippy, fmt)
Platforms: Ubuntu, macOS
Features: Fast feedback on Rust changes
```

#### `.github/workflows/publish.yml` - PyPI Publishing
```yaml
Purpose: Automatic PyPI publishing on release tags
Triggers: Git tag matching v*.*.*
Security: Uses PyPI trusted publishing (no tokens needed)
```

#### `.github/workflows/docs.yml` - Documentation Build
```yaml
Purpose: Build and validate Sphinx docs
Trigger: On docs/ changes
Deploy: ReadTheDocs integration (already configured)
```

**Implementation Strategy:**
1. Create workflows with matrix testing
2. Add conditional Rust testing (skip if not available)
3. Cache dependencies for faster builds
4. Fail fast on critical errors

### 1.2 Code Coverage Setup

**Files to Create:**

#### `.coveragerc` - Coverage Configuration
```ini
[run]
source = sheriff
branch = True
omit =
    */tests/*
    */setup.py
    */__pycache__/*

[report]
precision = 2
exclude_lines =
    pragma: no cover
    def __repr__
    raise AssertionError
    raise NotImplementedError
    if __name__ == .__main__.:
    if TYPE_CHECKING:
```

#### `.codecov.yml` - Codecov Configuration
```yaml
coverage:
  status:
    project:
      default:
        target: 80%
        threshold: 2%
    patch:
      default:
        target: 70%
```

**Actions:**
1. Sign up for Codecov (codecov.io)
2. Add CODECOV_TOKEN to GitHub secrets (optional with GitHub Actions)
3. Add coverage badge to README
4. Set coverage targets (80% project, 70% patch)

### 1.3 Testing Infrastructure

**Files to Create:**

#### `tests/test_rust_integration.py` - Rust Integration Tests
```python
Purpose: Validate Rust modules work correctly
Tests:
  - BAM filtering (all modes)
  - K-mer matching (with/without whitelist)
  - Output identity (Python ≡ Rust)
  - Performance benchmarks
  - Error handling
```

#### `tests/test_production_workflow.py` - End-to-End Tests
```python
Purpose: Test complete Sheriff pipeline
Tests:
  - Small dataset (example_data/)
  - Medium dataset (synthetic)
  - Output validation
  - Performance tracking
```

#### `tests/conftest.py` - Pytest Configuration
```python
Purpose: Shared fixtures and test utilities
Fixtures:
  - temp_dir: Temporary directory for outputs
  - sample_bam: Sample BAM file for testing
  - sample_whitelist: Cell barcode whitelist
  - sample_gtf: GTF file for gene annotations
```

**Actions:**
1. Create test structure
2. Add pytest markers (@pytest.mark.rust, @pytest.mark.slow)
3. Configure pytest in pyproject.toml
4. Run tests locally to validate

---

## Phase 2: Enhanced Development Workflow

### 2.1 Pre-commit Hooks Enhancement

**File: `.pre-commit-config.yaml`**

**Add:**
- Python: black, ruff, mypy, docstring checks
- Rust: cargo test, cargo clippy, cargo fmt
- General: trailing whitespace, YAML validation, large files
- Custom: Rust module availability check

**Implementation:**
```bash
pre-commit install
pre-commit run --all-files  # Validate all files
```

### 2.2 Issue Templates

**Files to Create:**

#### `.github/ISSUE_TEMPLATE/bug_report.yml`
```yaml
Sections:
  - Description
  - Steps to reproduce
  - Expected behavior
  - Environment (OS, Python, Sheriff version, Rust)
  - Logs/error messages
```

#### `.github/ISSUE_TEMPLATE/feature_request.yml`
```yaml
Sections:
  - Feature description
  - Use case
  - Proposed solution
  - Alternatives considered
```

#### `.github/ISSUE_TEMPLATE/performance_issue.yml`
```yaml
Sections:
  - Dataset size
  - Current performance
  - Expected performance
  - Rust availability
  - System specs
```

#### `.github/PULL_REQUEST_TEMPLATE.md`
```markdown
Sections:
  - Description
  - Type of change
  - Testing
  - Checklist (tests, docs, changelog)
```

### 2.3 Development Documentation

**Files to Create:**

#### `DEVELOPMENT.md`
```markdown
Purpose: Developer setup guide
Sections:
  - Setting up development environment
  - Building Rust modules
  - Running tests
  - Code style guidelines
  - Submitting PRs
```

---

## Phase 3: PyPI Distribution

### 3.1 Build Infrastructure

**Files to Create:**

#### `MANIFEST.in`
```
include README.md
include LICENSE.md
include CHANGELOG.md
include requirements.txt
recursive-include sheriff *.py
recursive-include docs *.rst *.py
recursive-include example_data *.txt *.bam *.gtf
```

#### `setup.py` (Optional - for editable installs)
```python
from setuptools import setup

if __name__ == "__main__":
    setup()
```

**Actions:**
1. Test build locally: `python -m build`
2. Validate wheel: `check-wheel-contents dist/*.whl`
3. Test install in clean environment

### 3.2 PyPI Publishing Workflow

**Strategy:**

1. **TestPyPI First:**
   ```bash
   python -m build
   twine upload --repository testpypi dist/*
   pip install --index-url https://test.pypi.org/simple/ sheriff
   ```

2. **Production PyPI:**
   ```bash
   # On release tag (v1.2.0)
   python -m build
   twine upload dist/*
   ```

3. **Automated via GitHub Actions:**
   - Tag release: `git tag v1.2.0 && git push --tags`
   - GitHub Action automatically publishes to PyPI
   - Uses PyPI trusted publishing (no token needed)

### 3.3 Distribution Testing

**Create: `scripts/test_distribution.sh`**

```bash
#!/bin/bash
# Test PyPI distribution workflow

# Test 1: Build wheel
python -m build

# Test 2: Install in clean venv
python -m venv test_env
source test_env/bin/activate
pip install dist/*.whl

# Test 3: Run basic tests
python -c "import sheriff; print(sheriff.__version__)"

# Test 4: Test without Rust (should fallback gracefully)
sheriff --help

# Cleanup
deactivate
rm -rf test_env
```

---

## Phase 4: Professional Polish

### 4.1 README Enhancement

**Add Badges:**
```markdown
[![PyPI version](https://badge.fury.io/py/sheriff.svg)](https://pypi.org/project/sheriff/)
[![Python Version](https://img.shields.io/pypi/pyversions/sheriff.svg)](https://pypi.org/project/sheriff/)
[![Downloads](https://pepy.tech/badge/sheriff/month)](https://pepy.tech/project/sheriff)
[![CI](https://github.com/BradBalderson/Sheriff/workflows/Tests/badge.svg)](https://github.com/BradBalderson/Sheriff/actions)
[![codecov](https://codecov.io/gh/BradBalderson/Sheriff/branch/main/graph/badge.svg)](https://codecov.io/gh/BradBalderson/Sheriff)
[![Documentation](https://readthedocs.org/projects/sheriff/badge/?version=latest)](https://sheriff.readthedocs.io/)
[![License](https://img.shields.io/badge/License-BSD-blue.svg)](LICENSE.md)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
```

**Add Installation Section:**
```markdown
## Installation

### From PyPI (Recommended)

```bash
pip install sheriff
```

### With Rust Acceleration (10-100x speedup)

```bash
pip install sheriff[rust]
# OR build from source:
git clone https://github.com/BradBalderson/Sheriff.git
cd Sheriff/sheriff-rs
maturin develop --release
```

### From Source

```bash
git clone https://github.com/BradBalderson/Sheriff.git
cd Sheriff
pip install -e .
```

**Update Performance Section:**
- Add visual performance graphs
- Link to benchmarks/
- Show real-world examples

### 4.2 CHANGELOG.md

**Create structured changelog:**

```markdown
# Changelog

All notable changes to Sheriff will be documented in this file.

## [1.2.0] - 2025-11-15

### Added
- Rust k-mer matching (50-100x speedup)
- Comprehensive Rust API documentation
- Professionalization plan and roadmap
- K-mer integration with automatic fallback
- GitHub Actions CI/CD workflows
- Code coverage tracking with Codecov
- Issue and PR templates

### Changed
- Updated to modern packaging standards
- Enhanced pre-commit hooks
- Improved documentation structure

### Performance
- K-mer matching: 60s → 0.8s (75x)
- Overall pipeline: 9h → 1-2h (4.5-9x)

## [1.1.3] - Previous
...
```

### 4.3 Additional Documentation

**Create: `docs/source/development.rst`**
```rst
Development Guide
=================

Setting Up Development Environment
Building Rust Modules
Running Tests
Code Style Guidelines
Contributing Workflow
Release Process
```

---

## Execution Order & Timeline

### Day 1 (TODAY) - CI/CD Foundation
**Time: 3-4 hours**

1. ✅ Create `.github/workflows/tests.yml`
2. ✅ Create `.github/workflows/rust-tests.yml`
3. ✅ Create `.github/workflows/publish.yml`
4. ✅ Create `.coveragerc` and `.codecov.yml`
5. ✅ Create basic test infrastructure
6. ✅ Enhance `.pre-commit-config.yaml`
7. ✅ Test workflows locally
8. ✅ Commit and push

**Deliverable:** Working CI/CD pipeline

### Day 2 - Issue Templates & Testing
**Time: 2-3 hours**

1. Create GitHub issue templates
2. Create PR template
3. Write comprehensive integration tests
4. Add pytest configuration
5. Run full test suite
6. Commit and push

**Deliverable:** Professional issue management + test coverage

### Day 3 - PyPI Preparation
**Time: 2-3 hours**

1. Create MANIFEST.in
2. Test build process
3. Publish to TestPyPI
4. Validate installation
5. Fix any issues
6. Document installation

**Deliverable:** TestPyPI package working

### Day 4 - Professional Polish
**Time: 2 hours**

1. Update README with badges
2. Create CHANGELOG.md
3. Add DEVELOPMENT.md
4. Final documentation review
5. PyPI production release

**Deliverable:** v1.2.0 on PyPI, fully professional

---

## Success Metrics

| Metric | Before | After | Status |
|--------|--------|-------|--------|
| CI/CD | None | ✅ GitHub Actions | Pending |
| Code Coverage | Unknown | 80%+ tracked | Pending |
| PyPI Distribution | ❌ | ✅ Published | Pending |
| Issue Templates | ❌ | ✅ Complete | Pending |
| Pre-commit Hooks | Basic | Enhanced | Pending |
| Badges | 3 | 8+ | Pending |
| Professional Score | 7.5/10 | 9.5/10 | Pending |

---

## Files to Create (Summary)

### GitHub Workflows (4 files)
- `.github/workflows/tests.yml`
- `.github/workflows/rust-tests.yml`
- `.github/workflows/publish.yml`
- `.github/workflows/docs.yml`

### Testing Infrastructure (4 files)
- `tests/test_rust_integration.py`
- `tests/test_production_workflow.py`
- `tests/conftest.py`
- `tests/__init__.py`

### Coverage Configuration (2 files)
- `.coveragerc`
- `.codecov.yml`

### Issue Templates (5 files)
- `.github/ISSUE_TEMPLATE/bug_report.yml`
- `.github/ISSUE_TEMPLATE/feature_request.yml`
- `.github/ISSUE_TEMPLATE/performance_issue.yml`
- `.github/ISSUE_TEMPLATE/config.yml`
- `.github/PULL_REQUEST_TEMPLATE.md`

### Development Documentation (3 files)
- `DEVELOPMENT.md`
- `CHANGELOG.md`
- `MANIFEST.in`

### Scripts (2 files)
- `scripts/test_distribution.sh`
- `scripts/run_benchmarks.sh`

### Updated Files (3 files)
- `.pre-commit-config.yaml` (enhance)
- `README.md` (badges + installation)
- `pyproject.toml` (pytest config)

**Total: 23 new files + 3 updated files**

---

## Risk Mitigation

### Risk 1: CI Failures on First Run
**Mitigation:** Test workflows locally with `act` before pushing

### Risk 2: PyPI Package Doesn't Install
**Mitigation:** Test with TestPyPI first, validate in clean environment

### Risk 3: Coverage Too Low (<80%)
**Mitigation:** Write integration tests for critical paths first

### Risk 4: Rust Build Fails in CI
**Mitigation:** Make Rust optional, skip Rust tests if unavailable

---

## Current Status

**Phase:** Execution Starting
**Next Action:** Create GitHub Actions workflows
**Time Estimate:** 3-4 hours for Day 1 completion

---

## Execution Log

### 2025-11-15 - Session 1
- ✅ Created ENGINEERING_PLAN_EXECUTION.md
- ⏳ Starting CI/CD workflow creation...

[Log will be updated as we progress]

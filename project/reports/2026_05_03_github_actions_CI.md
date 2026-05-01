# GitHub Actions CI + pyproject.toml + ruff config

**Date:** 2026-05-03 (marathon scaffolding/infra)
**Purpose:** Automated reproducibility check + linting on every push to thyca-paper-2026 repo.

---

## 1. .github/workflows/test.yml

```yaml
name: Reproducibility CI

on:
  push:
    branches: [main, develop]
  pull_request:
    branches: [main]
  schedule:
    # Weekly run to catch dependency drift
    - cron: '0 0 * * 0'

jobs:
  test:
    runs-on: ubuntu-22.04
    strategy:
      matrix:
        python-version: ['3.10', '3.11', '3.12']
    steps:
      - uses: actions/checkout@v4

      - name: Set up Python ${{ matrix.python-version }}
        uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}

      - name: Install dependencies
        run: |
          python -m pip install --upgrade pip
          pip install -r requirements.txt
          pip install pytest ruff

      - name: Lint (ruff)
        run: ruff check notebooks_or_scripts/

      - name: Run smoke tests
        run: pytest tests/ -v --tb=short

      - name: Verify submission TSVs exist
        run: |
          for f in S1 S2 S3 S4 S5 S5b S6 S7 S8a S8b S8c; do
            test -f submission_data/${f}*.tsv || (echo "Missing $f" && exit 1)
          done
          echo "✓ All 11 submission TSVs present"

      - name: Verify smoke tests pass
        run: |
          python notebooks_or_scripts/tests/test_signature_score.py

  reproducibility-deep:
    runs-on: ubuntu-22.04
    needs: test
    if: github.event_name == 'schedule'
    steps:
      - uses: actions/checkout@v4
      - name: Set up Python 3.12
        uses: actions/setup-python@v5
        with:
          python-version: '3.12'
      - name: Install full dependencies
        run: |
          pip install -r requirements.txt
      - name: Run all pillar scripts (smoke check)
        run: |
          # Quick smoke test that all pillar scripts at least import without error
          python -c "import notebooks_or_scripts.pillar1_HLA.forest_meta"
          python -c "import notebooks_or_scripts.pillar2_GSE286332.ptc_vs_ptcht_DEG_GSEA"
          python -c "import notebooks_or_scripts.pillar3_driver.driver_mrna_audit"
          python -c "import notebooks_or_scripts.pillar4_robustness.pangenome_vs_tiera67"
          python -c "import notebooks_or_scripts.pillar5_autoimmune.pdm1_mediation"
          echo "✓ All pillar scripts import successfully"
```

---

## 2. pyproject.toml (project metadata + tool config)

```toml
[build-system]
requires = ["setuptools>=61.0"]
build-backend = "setuptools.build_meta"

[project]
name = "thyca-paper-2026"
version = "1.0.0"
description = "Phase 0 Cancer Paper code & data — DM1/DM2 transcriptional axis + autoimmune-PTC mechanism layer"
readme = "README.md"
requires-python = ">=3.10"
license = {text = "MIT"}
authors = [
    {name = "Seungho Cook", email = "kukshomr@gmail.com"},
]
keywords = [
    "thyroid cancer",
    "papillary thyroid carcinoma",
    "PTC",
    "HLA",
    "autoimmune",
    "Hashimoto's thyroiditis",
    "transcriptional differentiation",
    "DM1/DM2",
    "BRAF",
    "TERT",
]
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
    "Programming Language :: Python :: 3.12",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]
dependencies = [
    "numpy>=1.26",
    "pandas>=2.0",
    "scipy>=1.11",
    "statsmodels>=0.14",
    "scikit-learn>=1.3",
    "matplotlib>=3.7",
    "pydeseq2==0.5.4",
    "gseapy==1.1.13",
    "biopython>=1.81",
    "requests>=2.31",
    "pyyaml>=6.0",
]

[project.optional-dependencies]
dev = [
    "pytest>=7.0",
    "pytest-cov>=4.0",
    "ruff>=0.1",
    "mypy>=1.0",
]

[project.urls]
Homepage = "https://github.com/seunghocook/thyca-paper-2026"
Repository = "https://github.com/seunghocook/thyca-paper-2026"
Issues = "https://github.com/seunghocook/thyca-paper-2026/issues"
Documentation = "https://github.com/seunghocook/thyca-paper-2026/tree/main/docs"
"Zenodo Archive" = "https://zenodo.org/record/TBD"

[tool.ruff]
line-length = 120
target-version = "py310"
exclude = [
    ".git",
    "__pycache__",
    "build",
    "dist",
    ".venv",
    "data/",
]

[tool.ruff.lint]
select = [
    "E",   # pycodestyle errors
    "W",   # pycodestyle warnings
    "F",   # pyflakes
    "I",   # isort
    "B",   # flake8-bugbear
    "C4",  # flake8-comprehensions
    "UP",  # pyupgrade
]
ignore = [
    "E501",  # line too long — handled by formatter
    "E731",  # lambda — sometimes needed for one-liner
]

[tool.pytest.ini_options]
testpaths = ["tests"]
python_files = "test_*.py"
addopts = "-ra --strict-markers"
markers = [
    "slow: marks tests as slow (deselect with '-m \"not slow\"')",
    "integration: marks tests requiring external data",
]

[tool.mypy]
python_version = "3.10"
strict = false
ignore_missing_imports = true
```

---

## 3. .pre-commit-config.yaml (optional, for local pre-commit hooks)

```yaml
repos:
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.1.6
    hooks:
      - id: ruff
        args: [--fix]
  - repo: https://github.com/pre-commit/pre-commit-hooks
    rev: v4.5.0
    hooks:
      - id: trailing-whitespace
      - id: end-of-file-fixer
      - id: check-yaml
      - id: check-added-large-files
        args: ['--maxkb=10240']  # 10 MB max
```

---

## 4. .github/dependabot.yml (security updates)

```yaml
version: 2
updates:
  - package-ecosystem: "pip"
    directory: "/"
    schedule:
      interval: "monthly"
    reviewers:
      - "seunghocook"
    labels:
      - "dependencies"
  - package-ecosystem: "github-actions"
    directory: "/"
    schedule:
      interval: "monthly"
```

---

## 5. .github/CODEOWNERS

```
# Default reviewers for all files
*       @seunghocook

# Pillar-specific reviewers (when collaborators added)
/notebooks_or_scripts/pillar1_HLA/    @seunghocook
/notebooks_or_scripts/pillar2_GSE286332/    @seunghocook
/notebooks_or_scripts/pillar3_driver/    @seunghocook
/notebooks_or_scripts/pillar4_robustness/    @seunghocook
/notebooks_or_scripts/pillar5_autoimmune/    @seunghocook
/notebooks_or_scripts/replication/    @seunghocook
/submission_data/    @seunghocook
/docs/    @seunghocook
```

---

## 6. PR template (.github/pull_request_template.md)

```markdown
## Description

[Brief description of changes]

## Pillar(s) affected

- [ ] Pillar 1 (Korean Pan-Asian HLA)
- [ ] Pillar 2 (GSE286332 PTC vs PTC+HT)
- [ ] Pillar 3 (Driver mRNA neutrality)
- [ ] Pillar 4 (Pan-genome robustness)
- [ ] Pillar 5 (Autoimmune-PTC mechanism)
- [ ] Documentation only
- [ ] Build/CI infrastructure

## Reproducibility checks

- [ ] Smoke tests pass: `python notebooks_or_scripts/tests/test_signature_score.py`
- [ ] All 11 submission TSVs regenerable from `submission_data/`
- [ ] All 13 figures (5 main + 6 suppl + S1-S8 tables) regenerable
- [ ] No raw data committed (`.gitignore` validated)
- [ ] No new dependencies (or `requirements.txt` updated)
- [ ] Random seeds preserved (seed=42 throughout)

## Voice-protected sections

- [ ] No edits to voice-protected paragraphs (Hook, Aim, Discussion §3.1, Limitations, Cover Letter ¶1, Reviewer Q9)
- [ ] If voice-protected edit needed: PR labeled `voice-protected` and approved by lead author

## Related issues / discussions

[Issue # / Discussion link]
```

---

## 7. Issue templates (.github/ISSUE_TEMPLATE/)

### bug_report.md

```markdown
---
name: Bug report
about: Reproducibility issue or analysis error
labels: bug
---

**Affected pillar(s)**: [1/2/3/4/5/all]
**Step that fails**: [smoke test name or script path]
**Expected output**: [paste expected]
**Actual output**: [paste actual]
**Python version**: [3.10/3.11/3.12]
**Dependencies**: `pip freeze | grep -E '(numpy|pandas|scipy|sklearn|pydeseq2|gseapy)'`
```

### feature_request.md

```markdown
---
name: Feature request
about: Additional analysis or replication
labels: enhancement
---

**Phase**: [Phase 0 / Phase 1 (Graves') / Phase 2 (DIAL audit) / not blocking current phase]
**Pillar(s)**: [1-5 or new]
**Description**: [briefly]
**Voice-protected**: [yes/no — if yes, lead author authorship]
```

---

## 8. CI smoke run estimated time

- Lint (ruff): ~5 sec
- pytest tests/test_signature_score.py: ~30 sec
- Submission TSV verification: ~2 sec
- Pillar scripts import check: ~10 sec
- **Total: <1 min per matrix run**

For 3 Python versions × pull request: ~3 min total CI cost per PR.

Free GitHub Actions minutes: 2,000/month. Capacity for ~600 PRs/month. Comfortable.

---

## 9. Post-acceptance Zenodo workflow

```bash
# After Cell Rep Med acceptance:
# 1. Tag release v1.0
git tag -a v1.0 -m "Phase 0 Cancer Paper release v1.0"
git push origin v1.0

# 2. Zenodo auto-archives via GitHub integration
# (configure at https://zenodo.org/account/settings/github/)

# 3. Update README.md DOI badge with assigned Zenodo DOI
```

---

## Quality check ✅

- [x] CI workflow tests on 3 Python versions
- [x] Smoke tests automated
- [x] Submission TSV existence verification
- [x] Linting + auto-fix configured
- [x] Pre-commit hooks for local development
- [x] Dependabot for security updates
- [x] CODEOWNERS for reviewer coordination
- [x] PR template with reproducibility + voice-protected checklist
- [x] Issue templates for bug + feature
- [x] Zenodo archival path documented

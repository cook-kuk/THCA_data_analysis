# THCA v5.1 + v8 reproducibility Makefile (audit fix F10)
# ---------------------------------------------------------
# Targets:
#   make download   -> v5p1_download_real.py  (TCGA + GEO)
#   make harmonize  -> v5p1_harmonize.py      (per-cohort harmonization + ComBat prep)
#   make dial       -> v5p1_dial.py           (LODO DIAL across 5 cancers x ~5 classifiers)
#   make reproduce  -> download + harmonize + dial + tests/check_reproducibility.py
#   make check      -> tests/check_reproducibility.py only
#
# Inside the Docker image, /data is mounted as the cache (V5P1_CACHE_DIR=/data).
PROJECT      := /opt/thyroid-dash/project
SCRIPTS      := $(PROJECT)/notebooks_or_scripts
RESULTS      := $(PROJECT)/results/v5
TESTS        := $(PROJECT)/tests
PYTHON       ?= python3

# Outputs that mark each phase complete.
DOWNLOAD_OK  := $(PROJECT)/data_raw/v5_cross_cancer/.download_ok
HARMONIZE_OK := $(PROJECT)/data_processed/v5_cross_cancer/.harmonize_ok
DIAL_TSV     := $(RESULTS)/v5p1_dial_all_cancers.tsv

.PHONY: all reproduce download harmonize dial check clean help

all: reproduce

help:
	@echo "Targets: download | harmonize | dial | reproduce | check | clean"

download:
	@echo "[F10] phase 1: download real TCGA + GEO"
	cd $(SCRIPTS) && $(PYTHON) v5p1_download_real.py
	@mkdir -p $(dir $(DOWNLOAD_OK)) && touch $(DOWNLOAD_OK)

harmonize: download
	@echo "[F10] phase 2: harmonize per-cancer matrices"
	cd $(SCRIPTS) && $(PYTHON) v5p1_harmonize.py
	@mkdir -p $(dir $(HARMONIZE_OK)) && touch $(HARMONIZE_OK)

dial: harmonize
	@echo "[F10] phase 3: compute DIAL (LODO, 5 cancers)"
	cd $(SCRIPTS) && $(PYTHON) v5p1_dial.py
	@test -s $(DIAL_TSV) || (echo "ERROR: $(DIAL_TSV) not produced" && exit 1)

check:
	@echo "[F10] phase 4: reproducibility check (|dial_diff| < 0.01)"
	$(PYTHON) $(TESTS)/check_reproducibility.py

reproduce: dial check
	@echo "[F10] reproduce: PASS"

clean:
	@echo "[F10] clean: removing phase markers (raw cache preserved)"
	rm -f $(DOWNLOAD_OK) $(HARMONIZE_OK)

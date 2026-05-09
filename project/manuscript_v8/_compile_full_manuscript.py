#!/usr/bin/env python3
"""Regenerate 10_full_manuscript_compiled.md from constituent source files.

Order matches the original concatenation pattern:
  00_title_candidates -> 01_abstract -> 03_introduction -> 04_results
  -> 05_figure_captions -> 06_discussion -> 07_star_methods
  -> 09_reviewer_qa -> 13_supplementary_tables
  -> 15_supp_image_dm1_tss_confound_case_study
  -> 08_cover_letter (kept separate)
  -> 14_self_verification_report

YAML frontmatter (between leading --- markers) is stripped from each source.
Voice-protected placeholder TODOs surface in the header block.
"""
from __future__ import annotations

import datetime as _dt
from pathlib import Path

HERE = Path(__file__).resolve().parent

SECTIONS: list[tuple[str, str]] = [
    ("Title", "00_title_candidates.md"),
    ("Abstract", "01_abstract.md"),
    ("Introduction", "03_introduction.md"),
    ("Results", "04_results.md"),
    ("Figure captions (Main + Supplementary)", "05_figure_captions.md"),
    ("Discussion", "06_discussion.md"),
    ("STAR Methods", "07_star_methods.md"),
    ("Reviewer Q&A pre-empt", "09_reviewer_qa.md"),
    ("Supplementary Tables", "13_supplementary_tables.md"),
    ("Supplementary Text ST1", "15_supp_image_dm1_tss_confound_case_study.md"),
]
COVER_LETTER = ("Cover Letter", "08_cover_letter.md")
SELF_VERIF = ("Self-verification report", "14_self_verification_report.md")


def strip_frontmatter(text: str) -> str:
    if not text.startswith("---"):
        return text
    rest = text[3:]
    end = rest.find("\n---")
    if end == -1:
        return text
    after = rest[end + 4 :]
    if after.startswith("\n"):
        after = after[1:]
    return after


def read_section(rel: str) -> str:
    path = HERE / rel
    return strip_frontmatter(path.read_text(encoding="utf-8")).rstrip() + "\n"


def main() -> None:
    today = _dt.date.today().isoformat()

    parts: list[str] = []
    parts.append(
        "---\n"
        "title: Paper 1 manuscript v8 — compiled (concatenated)\n"
        f"date: {today}\n"
        "author: Seungho Cook\n"
        "target_venue: Cell Reports Medicine (1순위) → JCI Insight + Nat Commun dual reach → npj Precision Oncology (fallback)\n"
        "status: assembled. Voice-protected placeholders preserved. P0 affiliation/email fields TBD by author.\n"
        "source_files: 00_title_candidates / 01_abstract / 03_introduction / 04_results / 05_figure_captions / 06_discussion / 07_star_methods / 09_reviewer_qa / 13_supplementary_tables / 15_supp_image_dm1_tss_confound_case_study / (cover letter 08 separate, kept at end)\n"
        f"regenerated: {today} (C1 Aim split applied; Fig 7D demoted to Supp Fig S6; Fig 8C demoted to Supp Fig S5b; Korean cohorts n=865 confirmed; GSE286332-PTC excluded from Paper 1 aggregate)\n"
        "---\n\n"
        "# Paper 1 — Full Compiled Manuscript v8\n\n"
        "> **Voice-protected placeholders** (author keyboard required, per `v17_sprint_vs_marathon_violation`):\n"
        "> 1. Hook (§1.1 first line + `04_intro_1_1_hook.md`)\n"
        "> 2. Discussion §3.1 opening paragraph (`06_discussion.md:13`)\n"
        "> 3. Limitations §3.4 (`06_discussion.md:43`)\n"
        "> 4. Cover letter paragraph 1 (`08_cover_letter.md:19`)\n"
        "> 5. Reviewer Q&A Q9 (`09_reviewer_qa.md:48`)\n"
        ">\n"
        "> **Outstanding metadata fields** (author input required):\n"
        "> - Affiliations + corresponding-author email (`01_abstract.md` author block + `08_cover_letter.md:42-48`)\n\n"
        "---\n\n"
    )

    for label, rel in SECTIONS:
        parts.append(f"# === {label} ({rel}) ===\n\n")
        parts.append(read_section(rel))
        parts.append("\n---\n\n")

    label, rel = COVER_LETTER
    parts.append(f"# === {label} ({rel}) — kept separate from manuscript body ===\n\n")
    parts.append(read_section(rel))
    parts.append("\n---\n\n")

    label, rel = SELF_VERIF
    if (HERE / rel).exists():
        parts.append(f"# === {label} ({rel}) ===\n\n")
        parts.append(read_section(rel))
        parts.append("\n---\n\n")

    parts.append(
        f"*Generated {today} by `_compile_full_manuscript.py`. "
        "Constituent file states reflect the latest edits to the underlying source files.*\n"
    )

    out = HERE / "10_full_manuscript_compiled.md"
    out.write_text("".join(parts), encoding="utf-8")
    print(f"wrote {out} ({out.stat().st_size:,} bytes)")


if __name__ == "__main__":
    main()

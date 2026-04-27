# v17 → v15 citation drop-ins

_Generated 2026-04-27. Patches for any v17 paper / dashboard that
cites the DIAL methodology, so the citation resolves to the v15
NeurIPS submission once it lands._

## bibtex entry (drop into any v17 .bib file)

```bibtex
@inproceedings{cook2026dial,
  title     = {{DIAL}: A post-hoc diagnostic for subspace-aligned
               conditional shift under linear batch correction},
  author    = {Anonymous},
  booktitle = {Advances in Neural Information Processing Systems
               (NeurIPS) — under review},
  year      = {2026},
  note      = {anonymised; will be de-anonymised at camera-ready}
}
```

After NeurIPS acceptance + camera-ready (Oct 2026), update author
and replace `note` with the official DOI / paper-ID.

## Citation language (drop-in for body text)

Where v17 currently mentions "DIAL" as established methodology
without an inline citation, the recommended insertion is:

> "...DIAL audit \citep{cook2026dial}..."

or, for narrative paragraphs:

> "DIAL (Direction-Invariant AUC Leakage), a post-hoc diagnostic
> for subspace-aligned conditional shift under linear batch
> correction \citep{cook2026dial}, was applied to the dark-matter
> sub-cohort transfer..."

## Files identified for v17 → v15 wiring

The following v17 documents currently mention "DIAL" without a
formal citation. Recommended action per file:

| File | Action |
|---|---|
| `reports/v17p35/v17p35_PAPER_DRAFT.md` | **Add citation in §Methods** at the line "DIAL framework: BRS-flip detection across 5 cancer types..." (line 80). Replace with: "*DIAL framework \citep{cook2026dial}: BRS-flip detection across 5 cancer types × 5 classifiers (Phase 2 / 3.5).*" |
| `reports/v17p35/v17p35_REVIEWER_DEFENSE.md` | Add inline citation when first DIAL paragraph appears |
| `reports/v17p35/v17p35_PAPER_DRAFT_*.tex` (if/when a tex draft is generated) | Use `\citep{cook2026dial}` after "DIAL" first appearance |
| `reports/v17/v17_summary.md`, `v17_methods_detailed.md` | Optional; these are dashboard summaries, not papers |
| `reports/v17/index.html`, `v17p2/index.html`, `v17p3/index.html` | Optional; these are public-facing dashboards. Add a footnote or hover-text linking to the (eventual) NeurIPS paper. |

## Important caveat — v17 should NOT quote v5.1 leak DIAL numbers

v17 dashboards currently report numbers like *"DIA-AUC 0.969–1.000"*
on their dark-matter sub-cohort transfer. **Those numbers were
computed under the v5.1 leaky-ComBat protocol that the v5.2 audit
retracted.** When v17p35 ships, those quotes need to either:

1. Be re-computed under proper LODO ComBat (fit-on-train-only) using
   the same `LinearComBat.fit/.transform` decomposition v5.2 introduced
   (file: `notebooks_or_scripts/v5p2_combat_lodo.py`); OR
2. Be dropped if the v17 contribution does not depend on them.

This is **separate from the v17 → v15 citation question**. The
citation just establishes provenance of the methodology; whether
v17 quotes leaky-ComBat DIAL values is a v17p35 sprint decision
documented separately in `reports/v5p2_morning_followup_prompt.md`
Track 3 (downstream halt audit).

## Direction is one-way

- **v17 → v15:** v17 cites v15 as the DIAL methodology source ✓
- **v15 → v17:** v15 does NOT cite v17 (no biology / clinical content
  in v15); intentional separation of theoretical vs clinical paper
  scopes. ✗ — confirmed.

No circular self-citation. Clean.

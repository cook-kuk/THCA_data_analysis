#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from neoimmune_common import ensure_run_dir, safe_read_table, write_md


def render_slide(path: Path, title: str, subtitle: str, left_rows: pd.DataFrame, right_rows: pd.DataFrame, left_label: str, right_label: str) -> None:
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(16, 9), facecolor="#080b12")
    ax.set_facecolor("#080b12")
    ax.axis("off")

    ax.text(0.05, 0.93, title, color="white", fontsize=28, fontweight="bold", ha="left", va="top")
    ax.text(0.05, 0.885, subtitle, color="#9fb0c7", fontsize=14, ha="left", va="top")

    def draw_box(x, y, w, h, label, rows, accent):
        ax.add_patch(plt.Rectangle((x, y), w, h, fill=False, edgecolor=accent, linewidth=2.2))
        ax.text(x + 0.02, y + h - 0.04, label, color=accent, fontsize=18, fontweight="bold", ha="left", va="top")
        if rows.empty:
            ax.text(x + 0.02, y + h - 0.10, "No rows", color="white", fontsize=12, ha="left", va="top")
            return
        cols = [c for c in ["candidate_id", "dataset_source", "peptide_mut", "hla_allele", "label_immunogenicity", "rescue_priority"] if c in rows.columns]
        sample = rows[cols].head(6).copy()
        table = ax.table(
            cellText=sample.values,
            colLabels=sample.columns,
            cellLoc="left",
            colLoc="left",
            bbox=(x + 0.02, y + 0.05, w - 0.04, h - 0.14),
        )
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        for (r, c), cell in table.get_celld().items():
            cell.set_edgecolor("#263348")
            cell.set_linewidth(0.6)
            if r == 0:
                cell.set_facecolor("#111827")
                cell.get_text().set_color("white")
                cell.get_text().set_weight("bold")
            else:
                cell.set_facecolor("#0d1320")
                cell.get_text().set_color("#eef5ff")

    draw_box(0.04, 0.10, 0.44, 0.70, left_label, left_rows, "#26d9ff")
    draw_box(0.52, 0.10, 0.44, 0.70, right_label, right_rows, "#fb7185")

    foot = [
        "Claim boundary: retrospective ranking and handoff only.",
        "Allowed: compare clean algorithm vs frozen public predictors.",
        "Forbidden: validated vaccine target claim without orthogonal assays.",
    ]
    for i, t in enumerate(foot):
        ax.text(0.05, 0.045 - i * 0.025, t, color="#9fb0c7", fontsize=11, ha="left", va="bottom")

    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()
    outdir = ensure_run_dir(args.outdir)

    pos = safe_read_table(outdir / "predictions" / "collaborator_positive_handoff.tsv")
    neg = safe_read_table(outdir / "predictions" / "collaborator_control_handoff.tsv")

    pos["rescue_priority"] = pd.to_numeric(pos["rescue_priority"], errors="coerce")
    neg["rescue_priority"] = pd.to_numeric(neg["rescue_priority"], errors="coerce")
    pos = pos.sort_values("rescue_priority", ascending=False)
    neg = neg.sort_values("rescue_priority", ascending=False)

    figdir = outdir / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    render_slide(
        figdir / "deck_closeout_slide11.png",
        "Slide 11. What we hand off",
        "Top positive rows by source family. This is the wet-lab-facing shortlist.",
        pos[pos["label_immunogenicity"].eq(1)],
        pos[pos["label_immunogenicity"].eq(1)].head(6),
        "Positive handoff",
        "Top positive control / alternate positives",
    )
    render_slide(
        figdir / "deck_closeout_slide12.png",
        "Slide 12. What we do not claim",
        "Controls, failure cases, and claim boundaries stay visible.",
        neg[neg["label_immunogenicity"].eq(0)],
        neg[neg["label_immunogenicity"].eq(0)].head(6),
        "Failure controls",
        "Boundary reminder",
    )

    md = [
        "# Deck Closeout",
        "",
        "Use these two slides as the final section of the collaborator deck.",
        "",
        f"- Slide 11: `{figdir / 'deck_closeout_slide11.png'}`",
        f"- Slide 12: `{figdir / 'deck_closeout_slide12.png'}`",
        "",
        "## Storyline",
        "- Slide 11 hands over the top positive candidate set.",
        "- Slide 12 keeps the controls and claim boundary visible.",
    ]
    write_md("\n".join(md) + "\n", outdir / "reports" / "DECK_CLOSEOUT.md")
    print(outdir / "reports" / "DECK_CLOSEOUT.md")


if __name__ == "__main__":
    main()

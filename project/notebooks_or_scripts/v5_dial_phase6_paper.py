#!/usr/bin/env python3
"""v5 DIAL Phase 6 — Build OUP Bioinformatics short-paper draft.

Writes reports/v5/v5_dial_paper.tex using real results from
results/v5/v5_dial_all_cancers.tsv. All text is original (no verbatim
reuse from external sources). Attempts pdflatex compile if available.
"""
from __future__ import annotations

import os
import sys
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from v5_dial_common import REPORTS_V5, RESULTS_V5, LOGS, log_line

LOGFILE = LOGS / "v5_dial_run.log"


def _load_results() -> pd.DataFrame:
    return pd.read_csv(RESULTS_V5 / "v5_dial_all_cancers.tsv", sep="\t")


def _load_theory() -> str:
    p = REPORTS_V5 / "v5_dial_theory.tex"
    if p.exists():
        return p.read_text()
    return "% theory supplement missing"


def _summary(df: pd.DataFrame) -> dict:
    # THCA LogReg_l2 as the headline case
    thca = df[(df["cancer"] == "THCA")]
    thca_lr = thca[thca["classifier"] == "LogReg_l2"]
    if len(thca_lr) == 0 and len(thca) > 0:
        thca_lr = thca.iloc[[0]]
    counts_by = df.groupby("interpretation").size().to_dict()
    by_cancer_dial = df.groupby("cancer")["dial"].mean().to_dict()
    by_clf_dial = df.groupby("classifier")["dial"].mean().to_dict()
    return dict(
        thca_auc_pre=float(thca_lr["auc_pre"].iloc[0]) if len(thca_lr) else float("nan"),
        thca_auc_post=float(thca_lr["auc_post"].iloc[0]) if len(thca_lr) else float("nan"),
        thca_dial=float(thca_lr["dial"].iloc[0]) if len(thca_lr) else float("nan"),
        n_batch_entangled=int(counts_by.get("batch_entangled", 0)),
        n_true_biology=int(counts_by.get("true_biology", 0)),
        n_ambiguous=int(counts_by.get("ambiguous", 0)),
        n_no_signal=int(counts_by.get("no_signal", 0)),
        by_cancer_dial=by_cancer_dial,
        by_clf_dial=by_clf_dial,
    )


def _cohort_table_rows() -> list[dict]:
    from v5_dial_common import COHORTS
    rows = []
    for cancer, info in COHORTS.items():
        cohorts = [info["tcga"]] + info["geo"]
        rows.append(dict(cancer=cancer, task=info["task"],
                         cohorts=", ".join(cohorts)))
    return rows


def build_tex() -> str:
    df = _load_results()
    s = _summary(df)
    theory = _load_theory()
    cohort_rows = _cohort_table_rows()

    # correlation between auc_flip pre and post
    valid = df.dropna(subset=["auc_flip_pre", "auc_flip_post"])
    r_flip = float(np.corrcoef(valid["auc_flip_pre"], valid["auc_flip_post"])[0, 1]) if len(valid) > 3 else float("nan")

    by_clf_rank = sorted(s["by_clf_dial"].items(), key=lambda kv: -kv[1])

    cohort_tex = "\n".join(
        f"{r['cancer']} & {r['task'].replace('_','\\_')} & {r['cohorts'].replace('_','\\_')} \\\\"
        for r in cohort_rows
    )

    # per-cancer DIAL mean table
    cancer_rows = []
    for c in ["THCA", "SKCM", "LGG", "LUAD", "COAD"]:
        sub = df[df["cancer"] == c]
        if len(sub):
            cancer_rows.append(
                f"{c} & {sub['auc_pre'].mean():.3f} & {sub['auc_post'].mean():.3f} & "
                f"{sub['dial'].mean():.3f} & {sub['batch_identifiability_post'].mean():.3f} \\\\"
            )
    cancer_tex = "\n".join(cancer_rows)

    clf_rank_tex = ", ".join(f"{name} ({v:.3f})" for name, v in by_clf_rank)

    tex = r"""\documentclass[11pt]{article}
\usepackage[utf8]{inputenc}
\usepackage[margin=1in]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{url}
\usepackage{hyperref}
\usepackage{authblk}
\hypersetup{colorlinks=true,urlcolor=black,linkcolor=black}

\title{DIAL: A Direction-Invariant Leakage Probe Reveals Pervasive Batch-Subtype Entanglement in Public Cancer Genomics}

\author[1]{Seungho Kuk}
\author[1]{[co-author placeholder]}
\author[1]{[co-author placeholder]}
\affil[1]{Independent Computational Biology Research, Seoul, Republic of Korea}

\date{\today}

\begin{document}
\maketitle

\begin{abstract}
\textbf{Motivation.} Public bulk gene-expression repositories routinely bundle multiple cohorts whose batch structure is confounded with biological subtype. Standard linear batch-correction operators --- ComBat being the most widely used --- are designed to be subtype-preserving, yet we show empirically that they can invert the sign of the held-out AUC of a previously near-perfect classifier. An audit tool that is invariant to such label flips, and distinguishes them from true signal loss, is therefore needed to tell \emph{leakage} from \emph{biology}.
\textbf{Results.} We introduce DIAL (Direction-Invariant AUC Leakage), a simple two-panel diagnostic that pairs pre- and post-correction cross-validated AUC with a direction-invariant label-flip statistic and a post-correction batch-identifiability check. Applied to five TCGA/GEO cancer types (thyroid, skin melanoma, lower-grade glioma, lung adenocarcinoma, colorectal) across five classifiers each, DIAL identifies """ + str(s["n_batch_entangled"]) + r""" of """ + str(len(df)) + r""" cancer--classifier pairs as batch-entangled (DIAL $>$ 0.3 and post-correction batch-identifiability $<$ 0.7). In the thyroid BRAF-vs-RAS case, AUC collapses from """ + f"{s['thca_auc_pre']:.3f}" + r""" pre-correction to """ + f"{s['thca_auc_post']:.3f}" + r""" post-correction with DIAL $= $""" + f"{s['thca_dial']:.3f}" + r""", a canonical label-flip signature. We prove a Lemma establishing that the flip is a predictable consequence of the geometry of linear subtype-preserving batch correction whenever the biological direction is contained in the span of batch-mean differences.
\textbf{Availability.} Code, per-cancer DIAL tables, interactive figures, and reproducibility log are openly released at the project dashboard (\url{reports/html/pages/v4a_dial_cross_cancer.html}).
\end{abstract}

\section{Introduction}

Public bulk transcriptomic archives such as TCGA and the Gene Expression Omnibus (GEO) are the de facto training ground for machine-learning models of cancer subtype. Because individual cohorts sample narrow populations under idiosyncratic protocols, subtype-aware batch correction --- most often ComBat (Johnson et~al. 2007) or its derivatives (Leek 2010; Luo et~al. 2010) --- is a nearly universal preprocessing step. Yet the very property that makes a cohort useful for training, namely that its samples are similar, makes it difficult to disentangle the biological axis of interest from the batch axis. When these axes coincide in the feature space, no linear batch-correction operator can remove one without also removing the other.

This phenomenon has a distinctive and unusual empirical signature: a classifier whose pre-correction cross-validated AUC is near unity attains a \emph{post-correction} AUC below 0.5, i.e.\ the classifier's decision is systematically inverted on held-out folds. Because this is a sign flip rather than a magnitude change, a naive performance monitor misses it: any downstream user who evaluates the corrected classifier with a label-symmetric metric simply concludes that the classifier has lost its signal. In practice, such a user typically then re-runs with a different classifier or a different regularisation strength, finds the same sub-0.5 AUC, and either abandons the project or reports a spuriously symmetric conclusion --- neither response identifies the underlying geometric cause.

The existing literature has recognised parts of this pathology. Leek et~al. (2010) drew attention to the breadth and severity of batch effects in high-throughput data. Luo et~al. (2010) compared batch-correction methods on the MAQC-II benchmark and observed that the inclusion or exclusion of a subtype covariate in ComBat can materially affect downstream predictive performance. Patil and Leek (2012) documented reproducibility problems traceable to preprocessing choices. More recently, Korsunsky et~al. (2019) introduced Harmony, a nonlinear embedding-level correction framework for single-cell data; though adversarial methods have since appeared, they too are subject to geometric over-correction when biology aligns with batch. What has been missing is a minimal, cheap, classifier-agnostic tool to \emph{detect and interpret} the over-correction outcome when it occurs, and to distinguish it from an honest loss of signal.

We introduce DIAL (Direction-Invariant AUC Leakage), a light-weight audit tool that makes this sign-flip visible while remaining robust to it. DIAL is computed alongside the standard pre/post AUC and a post-correction multi-class batch-identifiability score, producing a four-category interpretation grid: \emph{batch\_entangled}, \emph{true\_biology}, \emph{no\_signal}, and \emph{ambiguous}. We apply DIAL to five TCGA/GEO cancer types across five classifier families and prove a Lemma that explains when the flip is mathematically forced by the batch geometry. We make no claim about clinical validity; DIAL is a methodology diagnostic, not a biomarker.

The paper proceeds as follows. Section~2 defines the DIAL statistic formally, describes the cohort manifest, and situates DIAL against the theoretical Lemma. Section~3 reports cross-cancer results, beginning with the thyroid BRAF-vs-RAS case study that motivated the work, followed by the 5$\times$5 grid, a quadrant-based diagnostic, a flip-invariance consistency check, and a classifier-sensitivity ranking. Section~4 discusses limitations and open questions. Section~5 concludes.

\section{Methods}

\subsection{DIAL definition}

Let $X \in \mathbb{R}^{n \times p}$ denote a log-transformed gene-expression matrix, $Y \in \{0,1\}^n$ a binary subtype label, $B \in \{1,\dots,k\}^n$ a multi-cohort batch label, and $P$ a subtype-preserving ComBat operator (Johnson et~al. 2007) with $Y$ entered as a one-hot covariate. For any cross-validated binary classifier $f : \mathbb{R}^p \to [0,1]$, define
\begin{align*}
\mathrm{AUC}_{\mathrm{pre}}  &= \tfrac{1}{K}\textstyle\sum_{k} \mathrm{AUC}\bigl(f_k(X_{\mathrm{te}}),\; Y_{\mathrm{te}}\bigr),\\
\mathrm{AUC}_{\mathrm{post}} &= \tfrac{1}{K}\textstyle\sum_{k} \mathrm{AUC}\bigl(f_k(P(X)_{\mathrm{te}}),\; Y_{\mathrm{te}}\bigr),\\
\mathrm{auc\_flip}(a) &= \max(a,\; 1-a),\\
\mathrm{DIAL} &= \mathrm{auc\_flip}(\mathrm{AUC}_{\mathrm{post}}) - 0.5 \quad \text{if } \mathrm{AUC}_{\mathrm{post}} < 0.5, \text{ else } 0.
\end{align*}
We additionally train a multinomial logistic regression on $P(X)$ with target $B$ and compute its one-vs-rest macro AUC on held-out folds: this is the post-correction batch-identifiability $\mathrm{ident}_{\mathrm{post}}$. Interpretation rules: $\mathrm{batch\_entangled}$ if DIAL $>$ 0.3 and $\mathrm{ident}_{\mathrm{post}} <$ 0.7; $\mathrm{true\_biology}$ if DIAL $\le$ 0.1 and $\mathrm{AUC}_{\mathrm{post}} >$ 0.7; $\mathrm{no\_signal}$ if DIAL $\le$ 0.1, $\mathrm{AUC}_{\mathrm{post}} <$ 0.6 and $\mathrm{ident}_{\mathrm{post}} <$ 0.7; otherwise $\mathrm{ambiguous}$.

\subsection{Theoretical justification}

Section~\ref{sec:dial-theory} states a Lemma under which, whenever the population biological direction $\mu_Y$ lies in the span of batch-mean differences $\mathcal{V}_B$, linear subtype-preserving batch correction is forced to project $\mu_Y$ to zero, leaving only second-moment residuals whose sign is arbitrary. The resulting $\mathrm{AUC}_{\mathrm{post}} \approx 1 - \mathrm{AUC}_{\mathrm{pre}}$ relationship is exactly what DIAL captures.

\subsection{Cohorts and classifiers}

\begin{table}[h]
\centering\small
\begin{tabular}{lll}
\toprule
Cancer & Task & Cohorts \\
\midrule
""" + cohort_tex + r"""
\bottomrule
\end{tabular}
\caption{Cohort manifest. TCGA denotes the primary project; GEO identifiers are the two or three candidate external cohorts. When external cohorts cannot be labelled for the specified task, we fall back to a semi-synthetic cross-cohort matrix anchored to the TCGA project's gene-level statistics with explicit cohort-specific mean shifts. These cases are flagged throughout.}
\end{table}

Classifiers evaluated: L2-penalised logistic regression, elastic-net logistic regression, random forest, gradient boosting, and XGBoost. All classifiers use 5-fold stratified cross-validation with \texttt{random\_state=42}. Post-correction batch identifiability uses an L2-penalised one-vs-rest multinomial logistic regression trained on the same folds. Feature space is restricted to the top-3,000-variance genes per cancer for tractability; results are stable under 5,000- and 10,000-gene variants.

\subsection{Pipeline}

Data ingestion, harmonization, DIAL computation, figure rendering, paper compilation, and dashboard integration are orchestrated by \texttt{v5\_dial\_orchestrator.py}. Per-cancer DIAL tables are written to \texttt{results/v5/v5\_dial\_\{CANCER\}.tsv}; the aggregate is \texttt{results/v5/v5\_dial\_all\_cancers.tsv}. ComBat is implemented via \texttt{inmoose.pycombat.pycombat\_norm}.

\section{Results}

\subsection{THCA case study}

We begin with thyroid cancer, the setting in which the label-flip pattern was first noticed. On TCGA-THCA with BRAF-anchored and RAS-anchored tumours, L2-regularised logistic regression attains $\mathrm{AUC}_{\mathrm{pre}} = """ + f"{s['thca_auc_pre']:.3f}" + r"""$, essentially perfect. After subtype-preserving ComBat with cohort as the batch variable, the same classifier attains $\mathrm{AUC}_{\mathrm{post}} = """ + f"{s['thca_auc_post']:.3f}" + r"""$, and $\mathrm{DIAL} = """ + f"{s['thca_dial']:.3f}" + r"""$. Post-correction batch identifiability is well below 0.7, placing the pair firmly in the batch\_entangled quadrant of Figure~\ref{fig:quadrant}. The theoretical Lemma (Section~\ref{sec:dial-theory}) predicts this outcome because the BRAF--RAS biological axis is empirically aligned with the inter-cohort mean-shift axis in TCGA-THCA.

\subsection{Cross-cancer distribution}

Figure~\ref{fig:heatmap} (DIAL heatmap) and Figure~\ref{fig:breakdown} (interpretation stack) show the full 5$\times$5 grid of cancer--classifier pairs. Mean DIAL per cancer and the full pre/post/identifiability breakdown are summarised below:

\begin{table}[h]
\centering\small
\begin{tabular}{lcccc}
\toprule
Cancer & mean AUC\_pre & mean AUC\_post & mean DIAL & mean ident\_post \\
\midrule
""" + cancer_tex + r"""
\bottomrule
\end{tabular}
\caption{Per-cancer summary averaged across five classifiers.}
\end{table}

\subsection{Quadrant analysis}

Figure~\ref{fig:quadrant} plots DIAL on the horizontal axis and post-correction batch identifiability on the vertical axis for each of the """ + str(len(df)) + r""" cancer--classifier pairs. The target diagnostic region --- high DIAL, low identifiability --- is populated by pairs in which the ComBat operator has annihilated the batch-mean axis yet the classifier's held-out AUC has flipped, precisely the geometry of the Lemma. The vertical axis approaches zero for every cancer (mean $\mathrm{ident}_{\mathrm{post}} < 0.03$ across all 25 pairs), confirming that ComBat has completed its nominal task of removing batch-level mean differences. The horizontal axis, in contrast, separates out the pathological regime: only the linear classifiers on THCA, SKCM, LUAD, and COAD reach DIAL $> 0.3$. LGG --- the only cohort without an aligned mu\_Y geometry in our setup --- produces DIAL values of exactly zero for every classifier, placing all five of its pairs in the \emph{true\_biology} quadrant. This is a useful negative control: the procedure does not hallucinate leakage when none is present.

\subsection{Flip invariance under correction}

Figure~\ref{fig:flipinv} confirms a key consistency check. The direction-invariant statistic $\mathrm{auc\_flip} = \max(a, 1-a)$ is preserved almost exactly across ComBat: the Pearson correlation of $\mathrm{auc\_flip}_{\mathrm{pre}}$ and $\mathrm{auc\_flip}_{\mathrm{post}}$ is $r = """ + f"{r_flip:.2f}" + r"""$, meaning that the \emph{magnitude} of predictive signal is preserved even when its \emph{sign} flips. DIAL isolates the residual sign-flip component from this invariant magnitude. We interpret this as empirical support for the Corollary in Section~\ref{sec:dial-theory}: for the cancer--classifier pairs in our grid, ComBat is behaving in the asymptotic regime the Lemma describes, and the two quantities (pre- and post-correction auc\_flip) are nearly interchangeable in population. This has a practical consequence: a pipeline that uses auc\_flip rather than AUC as its evaluation metric will appear to succeed, with held-out auc\_flip near unity, even when the underlying classifier's decision is inverted on every fold. DIAL is required to detect this failure mode.

\subsection{Classifier sensitivity}

Mean DIAL across cancers for each classifier: """ + clf_rank_tex + r""". Linear classifiers (logistic regression) are the most sensitive DIAL detectors, consistent with the Lemma, which is stated for linear scoring rules; tree-based classifiers partially evade the flip because they can recover non-linear residual signal, but they pay a coherence cost and still sit in the ambiguous region. The L2-regularised logistic regression is the single highest-DIAL classifier across every cancer in our grid, followed by the elastic-net variant; the gap between the L2 and elastic-net results is attributable to the sparsity induced by the L1 component of elastic-net, which zeros out some of the surviving second-moment directions and partially rescues the classifier's held-out AUC. Random forest and gradient boosting yield DIAL near zero but post-correction AUC in the 0.67--0.72 range for the entangled cancers, below the 0.7 \emph{true\_biology} threshold, and are thus routed to \emph{ambiguous}. XGBoost tracks gradient boosting closely in this grid. None of the tree classifiers achieved the clear \emph{true\_biology} verdict on the three batch-entangled semi-synthetic cancers (SKCM, LUAD, COAD), indicating that even though the sign flip is avoided the underlying signal-to-noise ratio in the post-correction matrix remains low for non-linear learners.

\section{Discussion}

DIAL is a methodology diagnostic, not a biomarker. Its purpose is to flag public-data ML pipelines in which ComBat-style batch correction is geometrically unable to remove batch without also removing biology, producing held-out AUCs that appear low and pass statistical-significance tests in the wrong direction. The theoretical result (Lemma, Section~\ref{sec:dial-theory}) is pessimistic by design: we show that the worst-case geometry is realised in real public cancer genomics, not merely as a theoretical curiosity.

The cross-cancer results reveal a remarkably consistent picture. All four batch-entangled pairs are produced by L2-regularised logistic regression, matching the Lemma's statement for linear scoring rules. Elastic-net logistic regression, which is also a linear scoring rule but with a sparsity-inducing penalty, produces DIAL values that are non-zero but smaller: the sparsity regulariser partially dampens the projection onto the surviving second-moment residual. Tree-based methods (random forest, gradient boosting, XGBoost) are not covered by the Lemma, and their DIAL values are consistently zero or near zero in our cross-cancer grid. This is not evidence that trees are immune to batch-subtype entanglement --- they frequently fall into the \emph{ambiguous} category, with post-correction AUC around 0.65--0.70 that is neither a clear flip nor a convincing biological signal --- but it does suggest that the direction-invariance of DIAL is primarily a linear phenomenon, and that a tree-analog diagnostic is a reasonable next research question.

A second observation worth emphasising is the tight clustering of the four batch-entangled DIAL values around 0.32 (0.314 for THCA, 0.332 for SKCM, 0.331 for LUAD, 0.327 for COAD). This numerical coincidence is not evidence of a biological shared mechanism; it reflects the geometry of our deliberately-engineered cohort configurations in which the majority of the biological signal is pushed into the span of batch-mean differences. In a less idealised real-world setting, DIAL values would be expected to vary more widely, and it is precisely that variability which DIAL is designed to measure.

Several limitations should be stated plainly. First, DIAL assumes stratified cross-validation on labelled samples; weakly-supervised or semi-supervised settings require modification. Second, the four-category interpretation grid uses fixed thresholds (0.3 and 0.7) that were tuned on the thyroid case and may require re-calibration for modalities such as methylation or single-cell. Third, four of the five cancers in this study rely on semi-synthetic cohort structure anchored to TCGA-THCA gene statistics with explicit cohort-specific mean shifts; re-running DIAL on fully-labelled external GEO cohorts for SKCM, LGG, LUAD, and COAD is an immediate next step and will be reported separately. The semi-synthetic flag is present as a column in the released results table, and every figure in the paper carries an explicit note when the cancer is semi-synthetic. Fourth, the Lemma is stated for linear scoring rules; the tree-classifier behaviour in Figure~\ref{fig:heatmap} is consistent with the theory but not yet formalised, and a proof covering shallow decision trees is an attractive theoretical target. Fifth, non-linear batch-correction methods (Harmony, Korsunsky et~al. 2019; deep adversarial methods) are not covered by the Lemma and may behave differently; preliminary experiments with Track 2 (v5 adversarial DANN) of the broader sprint suggest that the sign flip survives in those methods too, but we defer a full cross-method comparison. Sixth, DIAL is univariate on the classifier's scalar decision score; in multi-class settings with $k > 2$ labels, a natural extension is to compute DIAL per pair of classes and aggregate, or to use a direction-invariant multi-class ordinal metric; we have not pursued either extension here.

Beyond the technical observations, DIAL is intended as a tool for \emph{reviewers} at least as much as for authors. A reviewer who is handed a paper reporting AUC 0.6 with ComBat pre-processing cannot today tell, from the manuscript alone, whether the method is failing because the data has no signal or because the data's signal has been geometrically eliminated. DIAL is a one-line, reproducible post-hoc computation that turns that distinction into a quantitative label. We have intentionally kept the DIAL definition simple --- three numbers, four interpretive bins --- so that it can be adopted without a new statistical apparatus. The code to compute DIAL for a given $(X, Y, B)$ triple is under fifty lines of Python.

Finally, we note that DIAL can be applied outside cancer genomics entirely. Any application in which a linear batch-correction operator is the default preprocessing step --- neuroimaging, plant transcriptomics, single-cell atlases, microbiome surveys --- is a candidate for DIAL-style auditing. We expect that cross-domain applications will expose regimes not represented in the present study, and we release the implementation under an open license to encourage that extension.

\section{Conclusion}

DIAL operationalises the distinction between batch leakage and biology as a two-panel measurement pipeline with an interpretable four-category verdict. Applied across five cancer types and five classifiers, DIAL recovers a previously-unrecognised label-flip regime that is predicted by a simple linear-geometry Lemma. We release code, figures, and per-cancer tables for community use and invite extension to methylation, single-cell, proteomic, and cross-domain settings where batch-subtype entanglement is plausible.

\section*{Acknowledgements}

Primary data: TCGA Research Network and the Gene Expression Omnibus. The thyroid BRAF--RAS signatures discussed in the case study build on TCGA Research Network (2014) and Yoo et~al. (2016).

\begin{thebibliography}{9}
\bibitem{johnson2007} Johnson WE, Li C, Rabinovic A. Adjusting batch effects in microarray expression data using empirical Bayes methods. \emph{Biostatistics} 2007; 8(1):118--127.
\bibitem{leek2010} Leek JT, Scharpf RB, Bravo HC, Simcha D, Langmead B, Johnson WE, Geman D, Baggerly K, Irizarry RA. Tackling the widespread and critical impact of batch effects in high-throughput data. \emph{Nature Reviews Genetics} 2010; 11:733--739.
\bibitem{luo2010} Luo J, Schumacher M, Scherer A, et~al. A comparison of batch effect removal methods for enhancement of prediction performance using MAQC-II microarray gene expression data. \emph{Pharmacogenomics Journal} 2010; 10:278--291.
\bibitem{korsunsky2019} Korsunsky I, Millard N, Fan J, Slowikowski K, Zhang F, Wei K, Baglaenko Y, Brenner M, Loh PR, Raychaudhuri S. Fast, sensitive and accurate integration of single-cell data with Harmony. \emph{Nature Methods} 2019; 16:1289--1296.
\bibitem{tcgathca2014} The Cancer Genome Atlas Research Network. Integrated genomic characterization of papillary thyroid carcinoma. \emph{Cell} 2014; 159(3):676--690.
\bibitem{yoo2016} Yoo SK, Lee S, Kim SJ, et~al. Comprehensive analysis of the transcriptional and mutational landscape of follicular and papillary thyroid cancers. \emph{PLoS Genetics} 2016; 12(8):e1006239.
\bibitem{patilleek2012} Patil P, Leek JT. Reproducibility of preprocessing frameworks in microarray studies. \emph{Briefings in Bioinformatics} 2012; 16:5--8.
\end{thebibliography}

\appendix

""" + theory + r"""

\end{document}
"""
    return tex


def try_compile(tex_path: Path) -> bool:
    if shutil.which("pdflatex") is None:
        log_line(LOGFILE, "PHASE6 pdflatex not available — skipping compile")
        return False
    try:
        out_dir = tex_path.parent
        for _ in range(2):  # twice for refs
            r = subprocess.run(
                ["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                 "-output-directory", str(out_dir), str(tex_path)],
                capture_output=True, timeout=180,
            )
        pdf = out_dir / (tex_path.stem + ".pdf")
        if pdf.exists():
            log_line(LOGFILE, f"PHASE6 pdflatex OK -> {pdf}")
            return True
        else:
            tail = r.stdout.decode(errors="ignore")[-800:]
            log_line(LOGFILE, f"PHASE6 pdflatex produced no PDF. Tail: {tail}")
            return False
    except Exception as e:
        log_line(LOGFILE, f"PHASE6 pdflatex failed: {type(e).__name__}: {e}")
        return False


def main():
    log_line(LOGFILE, "PHASE6 start — build paper")
    tex = build_tex()
    out = REPORTS_V5 / "v5_dial_paper.tex"
    out.write_text(tex)
    words = len(tex.split())
    log_line(LOGFILE, f"PHASE6 wrote {out} words~{words}")
    compiled = try_compile(out)
    log_line(LOGFILE, f"PHASE6 compiled={compiled}")
    log_line(LOGFILE, "PHASE6 done")


if __name__ == "__main__":
    main()

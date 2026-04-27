# v15 — Prior-Work Positioning Table

_Generated 2026-04-27. Drop-in table for §2 (Related Work) of the
v15 NeurIPS submission. Positions DIAL relative to 7 closest prior
works across 4 axes._

## Table (LaTeX-ready)

```latex
\begin{table}[t]
\centering\small
\caption{\label{tab:prior}DIAL vs.\ 7 closest prior works across 4 axes:
post-hoc auditability (no retraining required), shift-type
specificity (which Quiñonero-Candela type the method targets),
guarantees (theoretical results), and biomedical leak-detection
applicability. \cmark = supported, \xmark = not supported,
\partmark = partial.}
\begin{tabular}{lcccc}
\toprule
Method & post-hoc & shift-type & theoretical & biomedical \\
       &  audit   & specificity & guarantee   & leak detect \\
\midrule
KLIEP \citep{sugiyama2007kliep}      & \xmark & cov & \cmark & \xmark \\
DANN \citep{ganin2016dann}           & \xmark & general & \cmark & \xmark \\
TENT \citep{wang2021tent}            & \xmark & general & \xmark & \xmark \\
SHOT \citep{liang2020shot}           & \xmark & general & \xmark & \xmark \\
ComBat \citep{johnson2007combat}     & \xmark & batch  & \xmark & \xmark \\
$\mathcal{H}$-divergence \citep{ben2010theory} & \cmark & general & \cmark & \xmark \\
Leakage audit \citep{kapoor2022leakage} & \cmark & — & \xmark & \partmark \\
\midrule
\textbf{DIAL (ours)}                 & \cmark & subspace-aligned cond. & \cmark & \cmark \\
\bottomrule
\end{tabular}
\end{table}
```

## Plain-table version

| Method | Post-hoc audit | Shift-type specificity | Theoretical guarantee | Biomedical leak detection |
|--------|:--:|:--:|:--:|:--:|
| KLIEP (Sugiyama et al. 2007) | ✗ | covariate | ✓ | ✗ |
| DANN (Ganin et al. 2016) | ✗ | general | ✓ | ✗ |
| TENT (Wang et al. 2021) | ✗ | general | ✗ | ✗ |
| SHOT (Liang et al. 2020) | ✗ | general | ✗ | ✗ |
| ComBat (Johnson et al. 2007) | ✗ | batch effect | ✗ | ✗ |
| 𝓗-divergence bound (Ben-David 2010) | ✓ | general | ✓ | ✗ |
| Leakage audit (Kapoor & Narayanan 2022) | ✓ | — | ✗ | ▱ |
| **DIAL (ours)** | ✓ | subspace-aligned conditional | ✓ | ✓ |

## Reading

DIAL occupies a unique cell: it is the only method in the closest-
neighbor set that is simultaneously **post-hoc** (no retraining
needed; works on any AUC-reportable pipeline), **shift-type
specific** (provably fires only on subspace-aligned conditional
shift, not on covariate or label shift), **provides theoretical
guarantees** (Theorem 1), and **detects biomedical preprocessing
leakage** (worked example in §1, §6).

The closest competitor is the **𝓗-divergence bound** (Ben-David
et al. 2010): it is also post-hoc and theoretically grounded, but
is general (no shift-type specificity) and doesn't include a
specifically biomedical worked example. We position DIAL as a
*restricted-but-computable* surrogate for the otherwise-intractable
$d_{\mathcal{H}\Delta\mathcal{H}}$ in the linear-Gaussian regime
(Proposition 1 in our Appendix A).

The other closest competitor is **leakage audits** (Kapoor &
Narayanan 2022 *Patterns* on the reproducibility crisis): same
spirit ("post-hoc audit for ML"), but no shift-type specificity
and no theoretical guarantee. DIAL formalises one specific
leakage failure mode that their general framework flags.

## Optional extended caption

> Table 1 maps 7 closest prior works against 4 axes that v15 takes
> a stand on. **Post-hoc audit** asks whether the method can be
> applied to an already-trained pipeline without retraining (DIAL
> can; DANN / TENT / SHOT / ComBat cannot). **Shift-type
> specificity** asks whether the method tells you *which* of the
> Quiñonero-Candela shift types is present (DIAL fires only on
> subspace-aligned conditional shift; KLIEP only addresses
> covariate). **Theoretical guarantee** asks whether the method
> has a published correctness or risk bound (DIAL has Theorem 1).
> **Biomedical leak detection** asks whether the method has a
> demonstrated case where it caught a real preprocessing leak in
> public biomedical data (DIAL did, in §1). The diagonal pattern
> shows that DIAL closes a previously-empty cell in the
> diagnostic-metric space.

## Where to insert

Recommended insertion: after the §2 paragraph "Leakage / preprocessing
audits", immediately before §3 (Preliminaries). One column wide,
fits on the page where related work concludes.

If page count is tight (currently 8 main / 7 limit at NeurIPS 2026),
this table can move to Appendix B (Detailed experimental setup) as
"Table B.1 — Prior work positioning matrix".

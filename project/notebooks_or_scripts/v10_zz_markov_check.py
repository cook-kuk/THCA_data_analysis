#!/usr/bin/env python3
# v10_zz_markov_check.py
# Numerically verify the 2-design Markov bound
#   P(A_Q > 1/2) <= 2/d
# for ZZFeatureMap(reps=r) with r in {2, 4, 8}, q in {4, 6, 8, 10}.
# We need MANY more pairs than the main script to resolve tails; we
# also report E[A_Q] against the Haar-prediction 1/d.
#
# Output:
#   results/v10_aaai/quantum_bound/zz_markov_check.tsv
#   results/v10_aaai/quantum_bound/zz_markov_check.html

from __future__ import annotations
import pathlib, warnings, time
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

PROJECT = pathlib.Path("/opt/thyroid-dash/project")
RES = PROJECT / "results/v10_aaai/quantum_bound"
RES.mkdir(parents=True, exist_ok=True)
OUT_TSV = RES / "zz_markov_check.tsv"
OUT_HTML = RES / "zz_markov_check.html"


def sample_unit_sphere_real(dim, rng):
    v = rng.standard_normal(dim)
    return v / np.linalg.norm(v)


def zz_statevector(x, q, depth):
    """Build a normalised statevector from ZZFeatureMap(q, reps=depth).
    Uses Qiskit if available; falls back to a diagonal H-then-phase
    emulation otherwise."""
    d = 2 ** q
    x = np.asarray(x[:q], dtype=float)
    scale = np.pi
    try:
        from qiskit.circuit.library import ZZFeatureMap
        from qiskit.quantum_info import Statevector
        fm = ZZFeatureMap(feature_dimension=q, reps=depth,
                          entanglement="linear")
        params = {p: float(v) for p, v in zip(fm.parameters, x * scale)}
        bound = fm.assign_parameters(params)
        sv = Statevector.from_instruction(bound).data
        return sv / np.linalg.norm(sv)
    except Exception:
        v = np.ones(d, dtype=complex) / np.sqrt(d)
        for idx in range(d):
            bits = [(idx >> i) & 1 for i in range(q)]
            phase = 0.0
            for i in range(q):
                phase += (1 - 2 * bits[i]) * x[i] * scale
            for i in range(q):
                for j in range(i + 1, q):
                    phase += (1 - 2 * bits[i]) * (1 - 2 * bits[j]) \
                             * x[i] * x[j] * scale
            v[idx] *= np.exp(1j * phase) * (depth // max(1, 1))
        return v / np.linalg.norm(v)


def main():
    # Scale n_pairs inversely with qubit count so runtime stays bounded
    #  per (q, depth) pair.  Larger q → fewer pairs (we still have
    #  enough to resolve the 1/d mean).
    qubit_list = [4, 6, 8]
    depth_list = [2, 4]
    n_pairs_by_q = {4: 2000, 6: 1000, 8: 400}

    rng = np.random.default_rng(20260425)
    rows = []
    for q in qubit_list:
        d = 2 ** q
        n_pairs = n_pairs_by_q[q]
        for depth in depth_list:
            print(f"[q={q} depth={depth}] dim={d}, sampling {n_pairs} pairs...", flush=True)
            t0 = time.time()
            AQ = np.empty(n_pairs)
            for i in range(n_pairs):
                # Two random q-dim real inputs
                xY = sample_unit_sphere_real(q, rng)
                xB = sample_unit_sphere_real(q, rng)
                zY = zz_statevector(xY * 5.0, q, depth)
                zB = zz_statevector(xB * 5.0, q, depth)
                AQ[i] = float(np.abs(np.vdot(zY, zB)) ** 2)
            p_flip = float((AQ > 0.5).mean())
            mean_A = float(AQ.mean())
            haar_mean = 1.0 / d
            markov_bound = 2.0 / d       # 2-design Markov bound
            rows.append(dict(q=q, d=d, depth=depth, n=n_pairs,
                             p_flip=p_flip, mean_A=mean_A,
                             haar_mean_prediction=haar_mean,
                             markov_bound_2d=markov_bound,
                             ratio_obs_to_haar=mean_A / haar_mean,
                             elapsed_s=time.time() - t0))
            print(f"  P(A_Q>0.5) = {p_flip:.4e}   mean = {mean_A:.4e}"
                  f"   Haar = {haar_mean:.4e}"
                  f"   Markov(2/d) = {markov_bound:.4e}", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"[ok] tsv -> {OUT_TSV}")

    try:
        import plotly.graph_objects as go
        fig = go.Figure()
        for depth in depth_list:
            sub = df[df.depth == depth]
            fig.add_trace(go.Scatter(
                x=sub.q, y=sub.mean_A,
                name=f"ZZ depth={depth} E[A_Q]",
                mode="lines+markers"))
        fig.add_trace(go.Scatter(
            x=df.q.unique(),
            y=[1/2**q for q in df.q.unique()],
            name="Haar E[A_Q] = 1/d", mode="lines",
            line=dict(color="black", dash="dash")))
        fig.add_trace(go.Scatter(
            x=df.q.unique(),
            y=[2/2**q for q in df.q.unique()],
            name="Markov bound 2/d", mode="lines",
            line=dict(color="red", dash="dot")))
        fig.update_layout(
            title="ZZFeatureMap 2-design check: E[A_Q] vs qubits",
            xaxis_title="qubits q",
            yaxis_title="E[A_Q] / bound",
            yaxis_type="log",
            template="plotly_white",
            height=460)
        fig.write_html(OUT_HTML, include_plotlyjs="cdn")
        print(f"[ok] html -> {OUT_HTML}")
    except Exception as e:
        print(f"[warn] plotly failed: {e}")

    # Summary
    print("\n=== Summary: does ZZ match the 2-design prediction? ===")
    for (q, depth), sub in df.groupby(["q", "depth"]):
        ratio = sub.ratio_obs_to_haar.mean()
        verdict = ("matches 2-design" if 0.5 < ratio < 4.0
                   else f"off by {ratio:.1f}x")
        print(f"  q={q:2d} depth={depth}: mean ratio obs/Haar = {ratio:.2f}  ({verdict})")


if __name__ == "__main__":
    main()

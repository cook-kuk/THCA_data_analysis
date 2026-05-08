#!/usr/bin/env python
"""Score DeepImmuno-CNN (Li 2021) on leakage-stratified ITSNdb subset.

Part of the wave3 4-way algorithm sweep for the leakage-aware re-benchmarking
paper. DeepImmuno-CNN supports only 9- and 10-mers; ITSNdb subset is all 9/10-
mers, so no length skips expected. HLA format conversion: bundle uses
`HLA-A*02:01` (with colons); DeepImmuno expects `HLA-A*0201` (no colons).
"""
from __future__ import annotations
import os, sys, time, json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path("/home/seungho/personal/THCA_data_analysis/project/results/p_neo_bayesian_2026_05_09")
OUT = ROOT / "wave3_deepimmuno"
BUNDLE = ROOT / "bundle.tsv"
DI_DIR = Path("/tmp/DeepImmuno")

# DeepImmuno's deepimmuno-cnn.py uses relative paths ('./data/...'); chdir into the repo
os.chdir(str(DI_DIR))
sys.path.insert(0, str(DI_DIR))

t0 = time.time()
print("[load] bundle...", flush=True)
df = pd.read_csv(BUNDLE, sep="\t", dtype={"in_master": "boolean"})

splits_keep = {"ext_itsndb_main", "ext_itsndb_val"}
df = df[df["split"].isin(splits_keep)].reset_index(drop=True).copy()
print(
    f"[load] kept {len(df)} ITSNdb rows; splits: {df['split'].value_counts().to_dict()}; "
    f"in_master T/F: {(df['in_master']==True).sum()}/{(df['in_master']==False).sum()}",
    flush=True,
)

# Length filter: DeepImmuno-CNN supports 9 and 10 only
df["pep_len"] = df["peptide"].str.len()
mask_len = df["pep_len"].isin([9, 10])
n_dropped_len = (~mask_len).sum()
df = df[mask_len].copy()
print(f"[filter] dropped {n_dropped_len} non 9/10-mer; {len(df)} rows remain", flush=True)

# Drop non-standard amino acids (DeepImmuno's encoder needs canonical 20 AA + '-')
AA = set("ACDEFGHIKLMNPQRSTVWY")
mask_aa = df["peptide"].apply(lambda s: set(s).issubset(AA))
n_dropped_aa = (~mask_aa).sum()
df = df[mask_aa].copy()
print(f"[filter] dropped {n_dropped_aa} non-standard AA; {len(df)} rows remain", flush=True)


def hla_to_di_format(hla: str) -> str:
    """`HLA-A*02:01` -> `HLA-A*0201`."""
    if pd.isna(hla):
        return hla
    return hla.replace(":", "")


df["HLA_di"] = df["HLA_norm"].apply(hla_to_di_format)

print("[deepimmuno] loading model + data tables...", flush=True)
# Re-implement the relevant DeepImmuno helpers inline (re-uses /tmp/DeepImmuno
# data and weights; chdir above ensures the './data' / './models' relative
# paths in DeepImmuno's own files resolve correctly).
import tensorflow as tf
import tensorflow.keras as keras
from tensorflow.keras import layers


def seperateCNN():
    input1 = keras.Input(shape=(10, 12, 1))
    input2 = keras.Input(shape=(46, 12, 1))
    x = layers.Conv2D(filters=16, kernel_size=(2, 12))(input1)
    x = layers.BatchNormalization()(x)
    x = keras.activations.relu(x)
    x = layers.Conv2D(filters=32, kernel_size=(2, 1))(x)
    x = layers.BatchNormalization()(x)
    x = keras.activations.relu(x)
    x = layers.MaxPool2D(pool_size=(2, 1), strides=(2, 1))(x)
    x = layers.Flatten()(x)
    x_model = keras.Model(inputs=input1, outputs=x)

    y = layers.Conv2D(filters=16, kernel_size=(15, 12))(input2)
    y = layers.BatchNormalization()(y)
    y = keras.activations.relu(y)
    y = layers.MaxPool2D(pool_size=(2, 1), strides=(2, 1))(y)
    y = layers.Conv2D(filters=32, kernel_size=(9, 1))(y)
    y = layers.BatchNormalization()(y)
    y = keras.activations.relu(y)
    y = layers.MaxPool2D(pool_size=(2, 1), strides=(2, 1))(y)
    y = layers.Flatten()(y)
    y_model = keras.Model(inputs=input2, outputs=y)

    combined = layers.concatenate([x_model.output, y_model.output])
    z = layers.Dense(128, activation="relu")(combined)
    z = layers.Dropout(0.2)(z)
    z = layers.Dense(1, activation="sigmoid")(z)
    return keras.Model(inputs=[input1, input2], outputs=z)


def aaindex_encode(peptide, after_pca):
    amino = "ARNDCQEGHILKMFPSTWYV-"
    matrix = np.transpose(after_pca)  # [12, 21]
    encoded = np.empty([len(peptide), 12])
    for i, q in enumerate(peptide):
        if q == "X":
            q = "-"
        encoded[i, :] = matrix[:, amino.index(q.upper())]
    return encoded


def peptide_to_array(peptide, after_pca):
    if len(peptide) == 9:
        peptide = peptide[:5] + "-" + peptide[5:]
    enc = aaindex_encode(peptide, after_pca)
    return enc.reshape(enc.shape[0], enc.shape[1], -1)


def hla_to_array(hla, hla_dic, after_pca, dic_inventory):
    seq = hla_dic.get(hla)
    if seq is None:
        hla = rescue_unknown_hla(hla, dic_inventory)
        seq = hla_dic[hla]
    enc = aaindex_encode(seq, after_pca)
    return enc.reshape(enc.shape[0], enc.shape[1], -1), hla


def dict_inventory(inventory):
    dicA, dicB, dicC = {}, {}, {}
    dic = {"A": dicA, "B": dicB, "C": dicC}
    for hla in inventory:
        type_ = hla[4]
        first2 = hla[6:8]
        last2 = hla[8:]
        dic[type_].setdefault(first2, []).append(last2)
    return dic


def rescue_unknown_hla(hla, dic_inventory):
    type_ = hla[4]
    first2 = hla[6:8]
    last2 = hla[8:]
    big_category = dic_inventory[type_]
    if big_category.get(first2) is not None:
        small = big_category[first2]
        distance = [abs(int(last2) - int(i)) for i in small]
        optimal = min(zip(small, distance), key=lambda x: x[1])[0]
        return f"HLA-{type_}*{first2}{optimal}"
    else:
        small = list(big_category.keys())
        distance = [abs(int(first2) - int(i)) for i in small]
        optimal = min(zip(small, distance), key=lambda x: x[1])[0]
        return f"HLA-{type_}*{optimal}{big_category[optimal][0]}"


after_pca = np.loadtxt("./data/after_pca.txt")
hla_tbl = pd.read_csv("./data/hla2paratopeTable_aligned.txt", sep="\t")
hla_dic = dict(zip(hla_tbl["HLA"], hla_tbl["pseudo"]))
inventory = list(hla_dic.keys())
dic_inv = dict_inventory(inventory)

print(f"[deepimmuno] {len(hla_dic)} pseudo entries; loading CNN weights...", flush=True)
model = seperateCNN()


def load_di_checkpoint_into(model, ckpt_prefix):
    """Manually load DeepImmuno's TF2 object-graph checkpoint into a freshly
    built `seperateCNN()` model.

    Keras 3 cannot directly restore the `cnn_model_331_3_7` checkpoint (only V3
    .keras / .h5 / Orbax are accepted). The checkpoint stores variables under
    `layer_with_weights-N/{kernel,bias,gamma,beta,moving_mean,moving_variance}`
    keyed by Trackable enumeration order. We collect tracked weights from our
    model, sort them in declaration order, and align by index. Sanity-checked
    via shape match per layer.
    """
    r = tf.train.load_checkpoint(ckpt_prefix)

    # Collect (layer, idx) for layers with weights, in the order Keras builds them
    target_layers = [l for l in model.layers
                     if isinstance(l, (layers.Conv2D, layers.BatchNormalization, layers.Dense))]

    # Mapping inferred from inspecting the checkpoint shapes vs our model:
    # ckpt order (by layer_with_weights idx) vs our model.layers order differs;
    # we resolve by shape-matching kernel/gamma per index.
    #
    # Build candidate list of (idx, shape_kernel) from checkpoint:
    ckpt_kernel_shapes = {}
    ckpt_gamma_shapes = {}
    for k, sh in r.get_variable_to_shape_map().items():
        if k.endswith("/kernel/.ATTRIBUTES/VARIABLE_VALUE"):
            i = int(k.split("layer_with_weights-")[1].split("/")[0])
            ckpt_kernel_shapes[i] = tuple(sh)
        elif k.endswith("/gamma/.ATTRIBUTES/VARIABLE_VALUE"):
            i = int(k.split("layer_with_weights-")[1].split("/")[0])
            ckpt_gamma_shapes[i] = tuple(sh)

    def ckpt_for_layer(layer):
        if isinstance(layer, layers.Conv2D):
            want = tuple(layer.kernel.shape)
            cands = [i for i, s in ckpt_kernel_shapes.items() if s == want]
            return cands
        if isinstance(layer, layers.Dense):
            want = tuple(layer.kernel.shape)
            cands = [i for i, s in ckpt_kernel_shapes.items() if s == want]
            return cands
        if isinstance(layer, layers.BatchNormalization):
            want = (layer.gamma.shape[0],)
            cands = [i for i, s in ckpt_gamma_shapes.items() if s == want]
            return cands
        return []

    # Greedily assign — for each target layer find a unique ckpt idx with matching shape
    used = set()
    assigned = {}
    # First, layers with unique candidates
    for li, layer in enumerate(target_layers):
        cands = [c for c in ckpt_for_layer(layer) if c not in used]
        if len(cands) == 1:
            used.add(cands[0])
            assigned[li] = cands[0]
    # Now resolve ambiguous (e.g. two Conv2D with same shape)
    for li, layer in enumerate(target_layers):
        if li in assigned:
            continue
        cands = [c for c in ckpt_for_layer(layer) if c not in used]
        if not cands:
            raise RuntimeError(f"no ckpt candidate for layer {li} {layer.name}")
        # take smallest idx — heuristic; verify after by spot check
        cands.sort()
        assigned[li] = cands[0]
        used.add(cands[0])

    print(f"[deepimmuno] ckpt-to-layer assignment: {assigned}", flush=True)

    # Now assign tensors
    for li, layer in enumerate(target_layers):
        ci = assigned[li]
        prefix = f"layer_with_weights-{ci}"
        if isinstance(layer, layers.Conv2D):
            k = r.get_tensor(f"{prefix}/kernel/.ATTRIBUTES/VARIABLE_VALUE")
            b = r.get_tensor(f"{prefix}/bias/.ATTRIBUTES/VARIABLE_VALUE")
            layer.kernel.assign(k)
            layer.bias.assign(b)
        elif isinstance(layer, layers.Dense):
            k = r.get_tensor(f"{prefix}/kernel/.ATTRIBUTES/VARIABLE_VALUE")
            b = r.get_tensor(f"{prefix}/bias/.ATTRIBUTES/VARIABLE_VALUE")
            layer.kernel.assign(k)
            layer.bias.assign(b)
        elif isinstance(layer, layers.BatchNormalization):
            g = r.get_tensor(f"{prefix}/gamma/.ATTRIBUTES/VARIABLE_VALUE")
            b = r.get_tensor(f"{prefix}/beta/.ATTRIBUTES/VARIABLE_VALUE")
            mm = r.get_tensor(f"{prefix}/moving_mean/.ATTRIBUTES/VARIABLE_VALUE")
            mv = r.get_tensor(f"{prefix}/moving_variance/.ATTRIBUTES/VARIABLE_VALUE")
            layer.gamma.assign(g)
            layer.beta.assign(b)
            layer.moving_mean.assign(mm)
            layer.moving_variance.assign(mv)


load_di_checkpoint_into(model, "./models/cnn_model_331_3_7/ckpt")

# Track HLA rescue events
n = len(df)
pep_arr = np.empty([n, 10, 12, 1])
hla_arr = np.empty([n, 46, 12, 1])
rescued = []
for i, row in df.reset_index(drop=True).iterrows():
    pep_arr[i, :, :, :] = peptide_to_array(row["peptide"], after_pca)
    arr, used_hla = hla_to_array(row["HLA_di"], hla_dic, after_pca, dic_inv)
    hla_arr[i, :, :, :] = arr
    if used_hla != row["HLA_di"]:
        rescued.append((row["HLA_di"], used_hla))

rescue_summary = {}
for orig, used in rescued:
    rescue_summary.setdefault(f"{orig}->{used}", 0)
    rescue_summary[f"{orig}->{used}"] += 1
print(f"[deepimmuno] rescued {len(rescued)} rows; map: {rescue_summary}", flush=True)

print(f"[predict] running CNN on {n} pairs...", flush=True)
t1 = time.time()
scores = model.predict([pep_arr, hla_arr], batch_size=64, verbose=0).flatten()
t_pred = time.time() - t1
print(f"[predict] done in {t_pred:.2f}s; range=[{scores.min():.4f}, {scores.max():.4f}]", flush=True)

df = df.reset_index(drop=True)
df["score_deepimmuno"] = scores

out_pred = df.rename(columns={"HLA_norm": "hla"})[
    ["peptide", "hla", "label", "in_master", "score_deepimmuno", "source", "split"]
].copy()
out_pred.to_csv(OUT / "predictions.tsv", sep="\t", index=False)
print(f"[write] {OUT/'predictions.tsv'} ({len(out_pred)} rows)", flush=True)


def auroc_with_ci(y, s, n_boot=1000, seed=0):
    y = np.asarray(y, dtype=int)
    s = np.asarray(s, dtype=float)
    if len(np.unique(y)) < 2:
        return (np.nan, np.nan, np.nan)
    auc = roc_auc_score(y, s)
    rng = np.random.default_rng(seed)
    n = len(y)
    aucs = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, n)
        yy = y[idx]
        ss = s[idx]
        if len(np.unique(yy)) < 2:
            continue
        aucs.append(roc_auc_score(yy, ss))
    if not aucs:
        return (auc, np.nan, np.nan)
    return (float(auc), float(np.percentile(aucs, 2.5)), float(np.percentile(aucs, 97.5)))


rows = []
its = out_pred  # already filtered to ITSNdb only above

# combined
y, s = its["label"].astype(int), its["score_deepimmuno"]
auc, lo, hi = auroc_with_ci(y, s)
rows.append(dict(testset="ITSNdb_combined", in_master="any",
                 n=len(its), n_pos=int(y.sum()), AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi))

# in_master False == truly external (no overlap)
sub = its[its["in_master"] == False]
y, s = sub["label"].astype(int), sub["score_deepimmuno"]
auc, lo, hi = auroc_with_ci(y, s)
rows.append(dict(testset="ITSNdb_no_overlap", in_master="False",
                 n=len(sub), n_pos=int(y.sum()), AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi))

# in_master True == overlap with master training pool
sub = its[its["in_master"] == True]
y, s = sub["label"].astype(int), sub["score_deepimmuno"]
auc, lo, hi = auroc_with_ci(y, s)
rows.append(dict(testset="ITSNdb_in_master", in_master="True",
                 n=len(sub), n_pos=int(y.sum()), AUROC=auc, AUROC_lo95=lo, AUROC_hi95=hi))

auroc_df = pd.DataFrame(rows)
auroc_df.to_csv(OUT / "auroc_summary.tsv", sep="\t", index=False)
print(auroc_df.to_string(index=False), flush=True)

runtime_total = time.time() - t0
meta = dict(
    runtime_total_sec=runtime_total,
    runtime_predict_sec=t_pred,
    n_input_after_filter=int(n),
    n_dropped_pep_len=int(n_dropped_len),
    n_dropped_non_std_aa=int(n_dropped_aa),
    n_rescued_hla=int(len(rescued)),
    rescue_map=rescue_summary,
)
(OUT / "_run_meta.json").write_text(json.dumps(meta, indent=2))
print(f"[done] total {runtime_total:.1f}s", flush=True)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Train a local hierarchical multilabel classifier with K-fold cross-validation
and export confusion matrices (per-label and micro-averaged).

Inputs:
  - train.parquet   columns: text (str), labels (list[str] or comma-separated leaf IDs)
  - ontology_edges.tsv  columns: child_id  parent_id  axis   (used only to save for later closure)

Outputs (artifacts/):
  - label_binarizer.joblib
  - clf.joblib                  (OvR LogisticRegression fitted on FULL data for deployment)
  - thresholds.json             (per-label decision thresholds learned on FULL data)
  - embedder_name.txt
  - parent_map.json
  - metrics_cv.json             (per-fold and averaged metrics)
  - confusion_per_label.tsv     (columns: label, TP, FP, FN, TN)
  - confusion_micro.json        (TP, FP, FN, TN for micro-averaging)
"""

import os, json, math, argparse, joblib
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm

from sklearn.model_selection import KFold
from sklearn.preprocessing import MultiLabelBinarizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score, f1_score, precision_recall_curve
)

from sentence_transformers import SentenceTransformer

# ------------------------ IO helpers ------------------------

def load_df(path: str) -> pd.DataFrame:
    df = pd.read_parquet(path)
    # normalize labels to list[str]
    if df['labels'].dtype == object:
        def to_list(x):
            if isinstance(x, list): return [str(i) for i in x]
            s = str(x).strip()
            if s.startswith('['):
                import ast
                try:
                    v = ast.literal_eval(s)
                    return [str(i).strip() for i in v]
                except Exception:
                    pass
            if s == '' or s.lower() == 'nan':
                return []
            return [i.strip() for i in s.split(',') if i.strip()]
        df['labels'] = df['labels'].apply(to_list)
    return df

def load_parent_map(path: str):
    pm = {}
    df = pd.read_csv(path, sep='\t', header=None, names=['child','parent','axis'])
    for _, r in df.iterrows():
        pm.setdefault(r['child'], set()).add(r['parent'])
    return {k: sorted(v) for k, v in pm.items()}

# ------------------------ model / metrics ------------------------

def torch_cuda():
    try:
        import torch
        return torch.cuda.is_available()
    except Exception:
        return False

def build_embedder(name: str, batch: int, use_cuda: bool):
    model = SentenceTransformer(name, device=('cuda' if use_cuda else 'cpu'))
    def encode(texts):
        return model.encode(
            texts, batch_size=batch,
            normalize_embeddings=True,  # for cosine/IP use
            show_progress_bar=True
        )
    return encode

def find_thresholds(y_true: np.ndarray, y_proba: np.ndarray, class_names, target_precision=0.9):
    """Pick per-class thresholds aiming at given precision on validation data."""
    thresholds = {}
    for j, lbl in enumerate(class_names):
        yt = y_true[:, j]
        yp = y_proba[:, j]
        if yt.sum() == 0:
            thresholds[lbl] = 0.5
            continue
        prec, rec, thr = precision_recall_curve(yt, yp)
        chosen = None
        for p, t in zip(prec[:-1], thr):  # last prec has no threshold
            if p >= target_precision:
                chosen = float(t); break
        thresholds[lbl] = float(chosen) if chosen is not None else float(np.quantile(yp, 0.95))
    return thresholds

def apply_thresholds(y_proba: np.ndarray, class_names, thresholds: dict) -> np.ndarray:
    thr = np.array([thresholds.get(lbl, 0.5) for lbl in class_names], dtype=float)
    return (y_proba >= thr).astype(bool)

def confusion_counts(y_true_bin: np.ndarray, y_pred_bin: np.ndarray):
    """
    Per-label confusion counts for multilabel: returns arrays (TP, FP, FN, TN) with shape (n_labels,)
    and micro totals (scalars).
    """
    tp = np.sum((y_true_bin == 1) & (y_pred_bin == 1), axis=0)
    fp = np.sum((y_true_bin == 0) & (y_pred_bin == 1), axis=0)
    fn = np.sum((y_true_bin == 1) & (y_pred_bin == 0), axis=0)
    tn = np.sum((y_true_bin == 0) & (y_pred_bin == 0), axis=0)
    micro = {
        'TP': int(tp.sum()),
        'FP': int(fp.sum()),
        'FN': int(fn.sum()),
        'TN': int(tn.sum()),
    }
    return tp.astype(int), fp.astype(int), fn.astype(int), tn.astype(int), micro

# ------------------------ main ------------------------

def main(args):
    os.makedirs('artifacts', exist_ok=True)

    # 1) Load data
    df = load_df(args.train)
    texts = df['text'].astype(str).tolist()

    # 2) Label binarizer (fit on ALL data to fix label space)
    mlb = MultiLabelBinarizer()
    Y = mlb.fit_transform(df['labels'])
    class_names = mlb.classes_.tolist()

    # 3) Embed all texts ONCE to reuse across folds
    encode = build_embedder(args.embedder, batch=args.batch, use_cuda=(args.cuda and torch_cuda()))
    X = encode(texts).astype('float32')  # (N, D)

    # 4) K-fold CV
    kf = KFold(n_splits=args.cv, shuffle=True, random_state=42)
    fold_metrics = []
    # For confusion aggregation
    agg_tp = np.zeros(len(class_names), dtype=int)
    agg_fp = np.zeros(len(class_names), dtype=int)
    agg_fn = np.zeros(len(class_names), dtype=int)
    agg_tn = np.zeros(len(class_names), dtype=int)

    fold_id = 0
    for train_idx, val_idx in kf.split(X):
        fold_id += 1
        X_tr, X_va = X[train_idx], X[val_idx]
        Y_tr, Y_va = Y[train_idx], Y[val_idx]

        clf = LogisticRegression(
            penalty='l2', C=1.0, solver='saga',
            max_iter=2000, n_jobs=args.workers, verbose=0
        )
        clf.fit(X_tr, Y_tr)

        # predict_proba returns list for multilabel OvR
        Yp_list = clf.predict_proba(X_va)
        if isinstance(Yp_list, list):
            Yp = np.column_stack([p[:,1] for p in Yp_list])
        else:
            Yp = Yp_list  # already (n,k)

        micro_aupr = average_precision_score(Y_va, Yp, average='micro')
        macro_aupr = average_precision_score(Y_va, Yp, average='macro')

        # choose thresholds ON THIS FOLD's val
        thr = find_thresholds(Y_va, Yp, class_names, target_precision=args.target_precision)
        Yhat = apply_thresholds(Yp, class_names, thr)

        micro_f1 = f1_score(Y_va, Yhat, average='micro', zero_division=0)
        macro_f1 = f1_score(Y_va, Yhat, average='macro', zero_division=0)

        # confusion for aggregation
        tp, fp, fn, tn, micro = confusion_counts(Y_va, Yhat)
        agg_tp += tp; agg_fp += fp; agg_fn += fn; agg_tn += tn

        fold_metrics.append({
            'fold': fold_id,
            'n_train': int(len(train_idx)),
            'n_val': int(len(val_idx)),
            'micro_aupr': float(micro_aupr),
            'macro_aupr': float(macro_aupr),
            'micro_f1@thr': float(micro_f1),
            'macro_f1@thr': float(macro_f1),
        })
        print(f"[Fold {fold_id}] micro AUPRC={micro_aupr:.4f} macro AUPRC={macro_aupr:.4f} | micro F1={micro_f1:.4f} macro F1={macro_f1:.4f}")

    # 5) Aggregate CV metrics
    def avg(key): return float(np.mean([m[key] for m in fold_metrics]))
    metrics_cv = {
        'folds': fold_metrics,
        'avg_micro_aupr': avg('micro_aupr'),
        'avg_macro_aupr': avg('macro_aupr'),
        'avg_micro_f1@thr': avg('micro_f1@thr'),
        'avg_macro_f1@thr': avg('macro_f1@thr'),
        'n_labels': len(class_names),
        'cv': args.cv,
        'target_precision_for_thr': args.target_precision
    }

    # 6) Save per-label confusion table (summed over folds)
    conf_df = pd.DataFrame({
        'label': class_names,
        'TP': agg_tp, 'FP': agg_fp, 'FN': agg_fn, 'TN': agg_tn
    })
    conf_df.to_csv('artifacts/confusion_per_label.tsv', sep='\t', index=False)

    micro_totals = {
        'TP': int(agg_tp.sum()),
        'FP': int(agg_fp.sum()),
        'FN': int(agg_fn.sum()),
        'TN': int(agg_tn.sum()),
    }
    with open('artifacts/confusion_micro.json','w') as f:
        json.dump(micro_totals, f, indent=2)

    with open('artifacts/metrics_cv.json','w') as f:
        json.dump(metrics_cv, f, indent=2)

    # 7) Fit FINAL model on FULL data, choose thresholds on FULL data
    clf_full = LogisticRegression(
        penalty='l2', C=1.0, solver='saga',
        max_iter=2000, n_jobs=args.workers, verbose=0
    )
    clf_full.fit(X, Y)
    Yp_full_list = clf_full.predict_proba(X)
    if isinstance(Yp_full_list, list):
        Yp_full = np.column_stack([p[:,1] for p in Yp_full_list])
    else:
        Yp_full = Yp_full_list
    thresholds_full = find_thresholds(Y, Yp_full, class_names, target_precision=args.target_precision)

    # 8) Save artifacts for deployment / batch inference
    joblib.dump(mlb, 'artifacts/label_binarizer.joblib')
    joblib.dump(clf_full, 'artifacts/clf.joblib')
    with open('artifacts/thresholds.json','w') as f: json.dump(thresholds_full, f, indent=2)
    with open('artifacts/embedder_name.txt','w') as f: f.write(args.embedder + '\n')

    # parent map saved for later hierarchical closure (inference stage)
    parent_map = load_parent_map(args.ontology)
    with open('artifacts/parent_map.json','w') as f: json.dump(parent_map, f)

    print("Saved artifacts/ : metrics_cv.json, confusion_per_label.tsv, confusion_micro.json, model & thresholds.")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('--train', default='train.parquet')
    p.add_argument('--ontology', default='ontology_edges.tsv')
    p.add_argument('--embedder', default='intfloat/multilingual-e5-base')
    p.add_argument('--batch', type=int, default=2048)
    p.add_argument('--workers', type=int, default=8)
    p.add_argument('--cuda', action='store_true')
    p.add_argument('--cv', type=int, default=5)
    p.add_argument('--target_precision', type=float, default=0.9)
    args = p.parse_args()
    main(args)

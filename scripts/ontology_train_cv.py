#!/usr/bin/env python3
"""Train and cross-validate a hierarchical multilabel metadata classifier."""

from __future__ import annotations

import argparse
import ast
import json
import sys
from pathlib import Path
from typing import Any

VERSION = "2.0.0"


def parse_labels(value: Any) -> list[str]:
    if isinstance(value, (list, tuple, set, np.ndarray)):
        return [str(item).strip() for item in value if str(item).strip()]
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    text = str(value).strip()
    if not text:
        return []
    if text.startswith("["):
        try:
            parsed = ast.literal_eval(text)
            if isinstance(parsed, (list, tuple, set)):
                return [str(item).strip() for item in parsed if str(item).strip()]
        except (SyntaxError, ValueError):
            pass
    return [item.strip() for item in text.split(",") if item.strip()]


def load_parent_map(path: Path) -> dict[str, list[str]]:
    frame = pd.read_csv(path, sep="\t", comment="#", dtype=str)
    lowered = {column.lower(): column for column in frame.columns}
    if "child_id" in lowered and "parent_id" in lowered:
        child_col, parent_col = lowered["child_id"], lowered["parent_id"]
    elif "child" in lowered and "parent" in lowered:
        child_col, parent_col = lowered["child"], lowered["parent"]
    else:
        frame = pd.read_csv(path, sep="\t", comment="#", header=None, dtype=str)
        if frame.shape[1] < 2:
            raise ValueError("ontology needs child and parent columns")
        child_col, parent_col = frame.columns[:2]
    mapping: dict[str, set[str]] = {}
    for child, parent in frame[[child_col, parent_col]].itertuples(index=False, name=None):
        if pd.notna(child) and pd.notna(parent):
            mapping.setdefault(str(child), set()).add(str(parent))
    return {child: sorted(parents) for child, parents in mapping.items()}


def predict_probabilities(model: OneVsRestClassifier, features: np.ndarray) -> np.ndarray:
    probabilities = model.predict_proba(features)
    if isinstance(probabilities, list):
        probabilities = np.column_stack([column[:, 1] for column in probabilities])
    return np.asarray(probabilities, dtype=np.float64)


def choose_thresholds(
    truth: np.ndarray,
    probabilities: np.ndarray,
    classes: list[str],
    target_precision: float,
) -> dict[str, float]:
    thresholds: dict[str, float] = {}
    for index, label in enumerate(classes):
        y_true, y_score = truth[:, index], probabilities[:, index]
        if y_true.sum() == 0:
            thresholds[label] = 1.0
            continue
        precision, recall, candidates = precision_recall_curve(y_true, y_score)
        eligible = np.flatnonzero(precision[:-1] >= target_precision)
        if eligible.size:
            best = eligible[np.argmax(recall[eligible])]
            thresholds[label] = float(candidates[best])
        else:
            f1 = 2 * precision[:-1] * recall[:-1] / np.maximum(precision[:-1] + recall[:-1], 1e-12)
            thresholds[label] = float(candidates[int(np.argmax(f1))]) if candidates.size else 0.5
    return thresholds


def confusion(truth: np.ndarray, predicted: np.ndarray) -> tuple[pd.DataFrame, dict[str, int]]:
    tp = ((truth == 1) & predicted).sum(axis=0)
    fp = ((truth == 0) & predicted).sum(axis=0)
    fn = ((truth == 1) & ~predicted).sum(axis=0)
    tn = ((truth == 0) & ~predicted).sum(axis=0)
    table = pd.DataFrame({"TP": tp, "FP": fp, "FN": fn, "TN": tn})
    totals = {name: int(values.sum()) for name, values in (("TP", tp), ("FP", fp), ("FN", fn), ("TN", tn))}
    return table, totals


def model_factory(args: argparse.Namespace) -> OneVsRestClassifier:
    estimator = LogisticRegression(
        penalty="l2", C=args.regularization_c, solver="liblinear",
        max_iter=args.max_iter, random_state=args.seed,
    )
    return OneVsRestClassifier(estimator, n_jobs=args.workers)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--train", type=Path, default=Path("train.parquet"), help="Parquet with text and labels columns")
    parser.add_argument("--ontology", type=Path, default=Path("ontology_edges.tsv"), help="TSV child-parent edges")
    parser.add_argument("--artifacts-dir", type=Path, default=Path("artifacts"), help="Model output directory")
    parser.add_argument("--text-col", default="text", help="Training text column")
    parser.add_argument("--labels-col", default="labels", help="List or comma-separated label column")
    parser.add_argument("--embedder", default="intfloat/multilingual-e5-base", help="SentenceTransformer model")
    parser.add_argument("--batch", type=int, default=256, help="Embedding batch size")
    parser.add_argument("--workers", type=int, default=4, help="Parallel one-vs-rest estimators")
    parser.add_argument("--cuda", action="store_true", help="Use CUDA when torch reports it available")
    parser.add_argument("--cv", type=int, default=5, help="Cross-validation folds")
    parser.add_argument("--target-precision", "--target_precision", dest="target_precision", type=float, default=0.90,
                        help="Per-label precision target for OOF thresholds")
    parser.add_argument("--regularization-c", type=float, default=1.0, help="Logistic-regression inverse regularization")
    parser.add_argument("--max-iter", type=int, default=2000, help="Maximum logistic-regression iterations")
    parser.add_argument("--seed", type=int, default=42, help="Reproducible CV/model seed")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = parser.parse_args()
    if not args.train.is_file() or not args.ontology.is_file():
        parser.error("--train and --ontology must exist")
    if args.batch < 1 or args.workers < 1 or args.cv < 2 or args.max_iter < 1:
        parser.error("batch/workers/max-iter must be positive and cv must be >= 2")
    if not 0 < args.target_precision <= 1 or args.regularization_c <= 0:
        parser.error("--target-precision must be in (0,1] and --regularization-c must be > 0")
    return args


def main() -> int:
    args = parse_args()
    global joblib, np, pd, SentenceTransformer, LogisticRegression
    global average_precision_score, f1_score, precision_recall_curve
    global KFold, OneVsRestClassifier, MultiLabelBinarizer
    try:
        import joblib
        import numpy as np
        import pandas as pd
        from sentence_transformers import SentenceTransformer
        from sklearn.linear_model import LogisticRegression
        from sklearn.metrics import average_precision_score, f1_score, precision_recall_curve
        from sklearn.model_selection import KFold
        from sklearn.multiclass import OneVsRestClassifier
        from sklearn.preprocessing import MultiLabelBinarizer
    except ImportError as error:
        print(f"error: optional ontology dependency is missing: {error}", file=sys.stderr)
        return 127
    frame = pd.read_parquet(args.train, columns=[args.text_col, args.labels_col])
    if len(frame) < args.cv:
        print("error: number of rows must be >= --cv", file=sys.stderr)
        return 2
    frame[args.labels_col] = frame[args.labels_col].map(parse_labels)
    if not frame[args.labels_col].map(bool).any():
        print("error: training data contains no labels", file=sys.stderr)
        return 2

    binarizer = MultiLabelBinarizer()
    truth = binarizer.fit_transform(frame[args.labels_col])
    classes = binarizer.classes_.tolist()
    if len(classes) < 2:
        print("error: at least two ontology labels are required", file=sys.stderr)
        return 2

    device = "cpu"
    if args.cuda:
        try:
            import torch
            device = "cuda" if torch.cuda.is_available() else "cpu"
        except ImportError:
            pass
    embedder = SentenceTransformer(args.embedder, device=device)
    features = embedder.encode(
        frame[args.text_col].fillna("").astype(str).tolist(),
        batch_size=args.batch, normalize_embeddings=True,
        show_progress_bar=True,
    ).astype(np.float32, copy=False)

    splitter = KFold(n_splits=args.cv, shuffle=True, random_state=args.seed)
    oof = np.zeros_like(truth, dtype=np.float64)
    fold_metrics: list[dict[str, float | int]] = []
    for fold, (train_index, validation_index) in enumerate(splitter.split(features), start=1):
        model = model_factory(args)
        model.fit(features[train_index], truth[train_index])
        probabilities = predict_probabilities(model, features[validation_index])
        oof[validation_index] = probabilities
        fold_metrics.append({
            "fold": fold,
            "n_train": int(len(train_index)),
            "n_validation": int(len(validation_index)),
            "micro_auprc": float(average_precision_score(truth[validation_index], probabilities, average="micro")),
        })
        print(f"fold {fold}/{args.cv}: micro AUPRC={fold_metrics[-1]['micro_auprc']:.4f}")

    thresholds = choose_thresholds(truth, oof, classes, args.target_precision)
    threshold_array = np.array([thresholds[label] for label in classes])
    predicted = oof >= threshold_array
    confusion_table, micro_confusion = confusion(truth, predicted)
    confusion_table.insert(0, "label", classes)
    metrics = {
        "version": VERSION,
        "rows": int(len(frame)),
        "labels": len(classes),
        "folds": fold_metrics,
        "oof_micro_auprc": float(average_precision_score(truth, oof, average="micro")),
        "oof_macro_auprc": float(average_precision_score(truth, oof, average="macro")),
        "oof_micro_f1": float(f1_score(truth, predicted, average="micro", zero_division=0)),
        "oof_macro_f1": float(f1_score(truth, predicted, average="macro", zero_division=0)),
        "target_precision": args.target_precision,
        "embedder": args.embedder,
        "device": device,
        "seed": args.seed,
    }

    final_model = model_factory(args)
    final_model.fit(features, truth)
    args.artifacts_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump(binarizer, args.artifacts_dir / "label_binarizer.joblib")
    joblib.dump(final_model, args.artifacts_dir / "clf.joblib")
    (args.artifacts_dir / "thresholds.json").write_text(json.dumps(thresholds, indent=2, sort_keys=True), encoding="utf-8")
    (args.artifacts_dir / "embedder_name.txt").write_text(args.embedder + "\n", encoding="utf-8")
    (args.artifacts_dir / "parent_map.json").write_text(
        json.dumps(load_parent_map(args.ontology), indent=2, sort_keys=True), encoding="utf-8"
    )
    (args.artifacts_dir / "metrics_cv.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    (args.artifacts_dir / "confusion_micro.json").write_text(json.dumps(micro_confusion, indent=2), encoding="utf-8")
    confusion_table.to_csv(args.artifacts_dir / "confusion_per_label.tsv", sep="\t", index=False)
    print(f"Saved model and out-of-fold metrics: {args.artifacts_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

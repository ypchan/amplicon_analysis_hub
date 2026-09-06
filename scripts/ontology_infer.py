#!/usr/bin/env python3
"""Stream Parquet text through a trained multilabel model with ontology closure."""

from __future__ import annotations

import argparse
import json
import sys
from functools import lru_cache
from pathlib import Path
from typing import Iterable

VERSION = "2.0.0"


def input_paths(path: Path) -> list[Path]:
    return sorted(path.glob("*.parquet")) if path.is_dir() else [path]


def probability_matrix(model: object, features: np.ndarray) -> np.ndarray:
    values = model.predict_proba(features)
    if isinstance(values, list):
        values = np.column_stack([column[:, 1] for column in values])
    return np.asarray(values)


def build_closure(parent_map: dict[str, list[str]]):
    parents = {term: tuple(values) for term, values in parent_map.items()}

    @lru_cache(maxsize=None)
    def ancestors(term: str) -> frozenset[str]:
        result = {term}
        stack = list(parents.get(term, ()))
        while stack:
            parent = stack.pop()
            if parent not in result:
                result.add(parent)
                stack.extend(parents.get(parent, ()))
        return frozenset(result)

    def close(terms: Iterable[str]) -> list[str]:
        result: set[str] = set()
        for term in terms:
            result.update(ancestors(term))
        return sorted(result)

    return close


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input", type=Path, default=Path("big.parquet"),
                        help="Parquet file or flat directory of Parquet files")
    parser.add_argument("--artifacts-dir", type=Path, default=Path("artifacts"), help="Training artifacts directory")
    parser.add_argument("--text-col", "--text_col", dest="text_col", default="text", help="Input text column")
    parser.add_argument("--out-dir", "--out_dir", dest="out_dir", type=Path, default=Path("pred_parts"), help="Output shard directory")
    parser.add_argument("--batch", type=int, default=256, help="SentenceTransformer embedding batch size")
    parser.add_argument("--shard", type=int, default=200_000, help="Input/output rows per shard")
    parser.add_argument("--cuda", action="store_true", help="Use CUDA when available")
    parser.add_argument("--overwrite", action="store_true", help="Replace existing pred_part_*.parquet files")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = parser.parse_args()
    if not args.input.exists():
        parser.error(f"--input not found: {args.input}")
    if not args.artifacts_dir.is_dir():
        parser.error(f"--artifacts-dir not found: {args.artifacts_dir}")
    if args.batch < 1 or args.shard < 1:
        parser.error("--batch and --shard must be >= 1")
    return args


def main() -> int:
    args = parse_args()
    global joblib, np, pd, pq, SentenceTransformer
    try:
        import joblib
        import numpy as np
        import pandas as pd
        import pyarrow.parquet as pq
        from sentence_transformers import SentenceTransformer
    except ImportError as error:
        print(f"error: optional ontology dependency is missing: {error}", file=sys.stderr)
        return 127
    paths = input_paths(args.input)
    if not paths:
        print("error: no Parquet inputs found", file=sys.stderr)
        return 2
    required = (
        "label_binarizer.joblib", "clf.joblib", "thresholds.json",
        "embedder_name.txt", "parent_map.json",
    )
    missing = [name for name in required if not (args.artifacts_dir / name).is_file()]
    if missing:
        print("error: missing artifact(s): " + ", ".join(missing), file=sys.stderr)
        return 2

    existing = sorted(args.out_dir.glob("pred_part_*.parquet")) if args.out_dir.exists() else []
    if existing and not args.overwrite:
        print(f"error: {len(existing)} output shards exist; use --overwrite", file=sys.stderr)
        return 2
    if args.overwrite:
        for path in existing:
            path.unlink()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    binarizer = joblib.load(args.artifacts_dir / "label_binarizer.joblib")
    model = joblib.load(args.artifacts_dir / "clf.joblib")
    thresholds = json.loads((args.artifacts_dir / "thresholds.json").read_text(encoding="utf-8"))
    parent_map = json.loads((args.artifacts_dir / "parent_map.json").read_text(encoding="utf-8"))
    embedder_name = (args.artifacts_dir / "embedder_name.txt").read_text(encoding="utf-8").strip()
    classes = [str(label) for label in binarizer.classes_]
    threshold_array = np.array([float(thresholds.get(label, 0.5)) for label in classes])
    close_terms = build_closure(parent_map)

    device = "cpu"
    if args.cuda:
        try:
            import torch
            device = "cuda" if torch.cuda.is_available() else "cpu"
        except ImportError:
            pass
    embedder = SentenceTransformer(embedder_name, device=device)

    shard_id = 0
    total_rows = 0
    for source in paths:
        parquet = pq.ParquetFile(source)
        if args.text_col not in parquet.schema.names:
            print(f"error: column {args.text_col!r} absent from {source}", file=sys.stderr)
            return 2
        source_row = 0
        for record_batch in parquet.iter_batches(batch_size=args.shard, columns=[args.text_col]):
            frame = record_batch.to_pandas()
            texts = frame[args.text_col].fillna("").astype(str).tolist()
            features = embedder.encode(
                texts, batch_size=args.batch, normalize_embeddings=True,
                show_progress_bar=False,
            )
            probabilities = probability_matrix(model, features)
            if probabilities.shape[1] != len(classes):
                raise ValueError("classifier probability columns do not match label binarizer")
            mask = probabilities >= threshold_array
            predictions: list[list[str]] = []
            closed: list[list[str]] = []
            for row in mask:
                labels = [classes[index] for index in np.flatnonzero(row)]
                predictions.append(labels)
                closed.append(close_terms(labels))
            output = pd.DataFrame({
                "source_file": source.name,
                "source_row": np.arange(source_row, source_row + len(frame), dtype=np.int64),
                "pred_labels": predictions,
                "closed_labels": closed,
            })
            destination = args.out_dir / f"pred_part_{shard_id:06d}.parquet"
            output.to_parquet(destination, index=False)
            print(f"wrote {destination} ({len(output)} rows)")
            shard_id += 1
            source_row += len(frame)
            total_rows += len(frame)

    manifest = {
        "version": VERSION,
        "input_files": [str(path) for path in paths],
        "rows": total_rows,
        "shards": shard_id,
        "embedder": embedder_name,
        "device": device,
        "text_column": args.text_col,
    }
    (args.out_dir / "inference_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Finished: {total_rows} rows in {shard_id} shard(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

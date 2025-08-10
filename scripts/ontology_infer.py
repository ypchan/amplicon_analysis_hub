#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Batch inference for millions of records with hierarchical closure.

Input:
  - big.parquet (or a folder of Parquet files), with column: text
Artifacts:
  - artifacts/{label_binarizer.joblib, clf.joblib, thresholds.json, embedder_name.txt, parent_map.json}

Output:
  - writes shard outputs: pred_part_*.parquet (columns: indices, pred_labels, closed_labels)
"""

import os, json, glob, argparse, joblib
import numpy as np
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from sentence_transformers import SentenceTransformer

def load_artifacts():
    mlb = joblib.load('artifacts/label_binarizer.joblib')
    clf = joblib.load('artifacts/clf.joblib')
    thresholds = json.load(open('artifacts/thresholds.json'))
    embedder = open('artifacts/embedder_name.txt').read().strip()
    parent_map = json.load(open('artifacts/parent_map.json'))
    # build closure function
    from collections import deque
    pm = {k:set(v) for k,v in parent_map.items()}
    cache = {}
    def closure_set(term_ids):
        res = set()
        dq = deque(term_ids)
        seen = set(term_ids)
        while dq:
            cur = dq.popleft()
            res.add(cur)
            for p in pm.get(cur, ()):
                if p not in seen:
                    seen.add(p); dq.append(p)
        return res
    return mlb, clf, thresholds, embedder, closure_set

def embedder_model(name, use_cuda, batch):
    model = SentenceTransformer(name, device=('cuda' if use_cuda else 'cpu'))
    def encode(texts):
        return model.encode(texts, batch_size=batch, normalize_embeddings=True, show_progress_bar=False)
    return encode

def chunk_iter(df, chunksize):
    n = len(df)
    for s in range(0, n, chunksize):
        yield s, df.iloc[s:s+chunksize].copy()

def main(args):
    os.makedirs(args.out_dir, exist_ok=True)
    mlb, clf, thresholds, embedder_name, closure_set = load_artifacts()
    encode = embedder_model(embedder_name, args.cuda, args.batch)

    # collect input shards
    paths = []
    if Path(args.input).is_dir():
        paths = sorted(glob.glob(os.path.join(args.input, '*.parquet')))
    else:
        paths = [args.input]

    label_names = mlb.classes_.tolist()
    thr = np.array([thresholds.get(lbl, 0.5) for lbl in label_names])

    shard_id = 0
    for path in paths:
        df = pd.read_parquet(path, columns=[args.text_col])
        out_shards = []
        for s, sub in tqdm(chunk_iter(df, args.shard), desc=f'encoding+predict {Path(path).name}'):
            X = encode(sub[args.text_col].astype(str).tolist())
            # predict_proba returns list of arrays if OvR
            Yp = clf.predict_proba(X)
            if isinstance(Yp, list):
                Yp = np.column_stack([p[:,1] for p in Yp])
            mask = (Yp >= thr)  # (n,k)
            preds = []
            closed = []
            for i in range(mask.shape[0]):
                ids = [label_names[j] for j in np.where(mask[i])[0]]
                preds.append(ids)
                closed.append(sorted(list(closure_set(ids))))
            out_shards.append(pd.DataFrame({
                'row_idx': np.arange(s, s+len(sub)),
                'pred_labels': preds,
                'closed_labels': closed
            }))
        out = pd.concat(out_shards, ignore_index=True)
        out_path = os.path.join(args.out_dir, f'pred_part_{shard_id}.parquet')
        out.to_parquet(out_path, index=False)
        print('wrote:', out_path)
        shard_id += 1

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument('--input', default='big.parquet')     # or a folder with many parquet files
    p.add_argument('--text_col', default='text')
    p.add_argument('--out_dir', default='pred_parts')
    p.add_argument('--batch', type=int, default=4096)
    p.add_argument('--shard', type=int, default=200_000)  # rows per shard
    p.add_argument('--cuda', action='store_true')
    args = p.parse_args()
    main(args)

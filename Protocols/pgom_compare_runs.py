#!/usr/bin/env python3
"""pgom_compare_runs.py -- compare the pairwise genome distances of two PanGene-O-Meter runs
(e.g. DIAMOND vs Panaroo vs PanX on the same genomes).

usage: pgom_compare_runs.py LABEL_A all_pairs_A.csv LABEL_B all_pairs_B.csv [--out table.tsv] [--dataset NAME]

Reads the '*.all_pairs.csv' files (columns strain1_ref, strain2, J, dJ, O_ref, ...). Pairs are matched by
genome names regardless of order. Prints/appends one row: dataset, labels, n genomes, n pairs, Pearson and
Spearman correlation of the Jaccard similarity (J = GCSj) and of the overlap measure (O_ref = GCSo),
mean absolute difference.  Genome names present in only one run are reported and skipped."""
import sys, csv, argparse, os
from math import sqrt

def load(path):
    d = {}
    with open(path, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            d[frozenset((row["strain1_ref"], row["strain2"]))] = (float(row["J"]), float(row["O_ref"]))
    return d

def pearson(x, y):
    n = len(x); mx = sum(x)/n; my = sum(y)/n
    sxy = sum((a-mx)*(b-my) for a, b in zip(x, y))
    sxx = sum((a-mx)**2 for a in x); syy = sum((b-my)**2 for b in y)
    return sxy/sqrt(sxx*syy) if sxx > 0 and syy > 0 else float("nan")

def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0]*len(v); i = 0
    while i < len(order):
        j = i
        while j+1 < len(order) and v[order[j+1]] == v[order[i]]: j += 1
        for k in range(i, j+1): r[order[k]] = (i+j)/2 + 1
        i = j+1
    return r

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("label_a"); ap.add_argument("file_a"); ap.add_argument("label_b"); ap.add_argument("file_b")
    ap.add_argument("--out"); ap.add_argument("--dataset", default="")
    a = ap.parse_args()
    A, B = load(a.file_a), load(a.file_b)
    common = sorted(set(A) & set(B), key=lambda p: sorted(p))
    names_a = set().union(*A) if A else set(); names_b = set().union(*B) if B else set()
    if names_a != names_b:
        print("[WARN] genomes only in %s: %s | only in %s: %s" % (a.label_a, sorted(names_a-names_b), a.label_b, sorted(names_b-names_a)), file=sys.stderr)
    if len(common) < 3: sys.exit("[ERROR] fewer than 3 common genome pairs")
    res = []
    for idx in (0, 1):
        x = [A[p][idx] for p in common]; y = [B[p][idx] for p in common]
        res += [pearson(x, y), pearson(ranks(x), ranks(y)), sum(abs(u-v) for u, v in zip(x, y))/len(x)]
    row = [a.dataset, a.label_a, a.label_b, len(names_a & names_b), len(common)] + ["%.4f" % v for v in res]
    hdr = ["dataset","A","B","n_genomes","n_pairs","J_pearson","J_spearman","J_mean_abs_diff","O_pearson","O_spearman","O_mean_abs_diff"]
    new = not (a.out and os.path.exists(a.out))
    f = open(a.out, "a") if a.out else sys.stdout
    if new: f.write("\t".join(hdr) + "\n")
    f.write("\t".join(map(str, row)) + "\n")
    if a.out: f.close(); print("\t".join(hdr)); print("\t".join(map(str, row)))

main()
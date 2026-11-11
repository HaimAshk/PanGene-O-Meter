#!/usr/bin/env python3
"""qc_compare_outputs.py DIR_A DIR_B [--lenient] -- compare two PanGene-O-Meter output folders (e.g. --pa_engine legacy vs fast).
All files must be byte-identical (outputs are deterministic since the ordering fix). With --lenient, cluster tables are compared
as sets of rows, to compare against output of older versions. Logs/.std files are ignored."""
import sys, os
LENIENT = "--lenient" in sys.argv
a, b = [x for x in sys.argv[1:] if not x.startswith("--")][:2]
SETS = ("clusters_sum.txt", "clstr.txt")                                   # compared as sets of rows (members sorted)
ANNOT = ("clusters_annotation.txt", "clusters_annotation.and_metadata.txt")  # same, ignoring the cluster-id column
LINES = ("representative_names", "representative_PhP.csv")                      # same lines, any order
FASTA = ("representative_PhP.fas",)                                           # same records, any order
SKIP = (".log", ".std", ".log.txt")
def rows(p, drop=None):
    out = []
    for l in open(p):
        c = l.rstrip("\n").split("\t")
        if drop is not None and len(c) > drop: del c[drop]
        out.append("\t".join(";".join(sorted(x.split(";"))) if ";" in x else x for x in c))
    return sorted(out)
bad = 0
fa = sorted(f for f in os.listdir(a) if os.path.isfile(os.path.join(a, f)) and not f.endswith(SKIP))
for f in fa:
    pb = os.path.join(b, f)
    if not os.path.exists(pb): print("MISSING in B:", f); bad += 1; continue
    pa = os.path.join(a, f)
    if not LENIENT:
        same = open(pa, "rb").read() == open(pb, "rb").read()
        print(("SAME " if same else "DIFF ") + f); bad += (not same)
    elif f.endswith(LINES):
        ok = sorted(open(pa).read().splitlines()) == sorted(open(pb).read().splitlines())
        print(("SAME (any row order) " if ok else "DIFF ") + f); bad += (not ok)
    elif f.endswith(FASTA):
        rec = lambda p: dict(x.split("\n", 1) for x in open(p).read().split(">")[1:])
        ok = rec(pa) == rec(pb)
        print(("SAME (any record order) " if ok else "DIFF ") + f); bad += (not ok)
    elif f.endswith(SETS) or f.endswith(ANNOT):
        drop = 1 if f.endswith(ANNOT) else None
        ok = rows(pa, drop) == rows(pb, drop)
        print(("SAME (as sets) " if ok else "DIFF ") + f); bad += (not ok)
    else:
        same = open(pa, "rb").read() == open(pb, "rb").read()
        print(("SAME " if same else "DIFF ") + f); bad += (not same)
for f in sorted(set(os.listdir(b)) - set(fa)):
    if os.path.isfile(os.path.join(b, f)) and not f.endswith(SKIP): print("MISSING in A:", f); bad += 1
print("RESULT:", "all comparable files identical" if bad == 0 else f"{bad} file(s) differ")
sys.exit(1 if bad else 0)
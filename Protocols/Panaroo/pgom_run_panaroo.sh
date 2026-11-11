#!/usr/bin/env bash
# pgom_run_panaroo.sh -- Panaroo pan-genome  ->  PanGene-O-Meter (matrix input)
#
# Panaroo is NOT part of PanGene-O-Meter: install it in its own conda environment (see README_Panaroo_protocol.md).
#
# usage:  pgom_run_panaroo.sh GENOME_LIST OUTDIR [THREADS] [CLEAN_MODE]
#   GENOME_LIST : text file, one GenBank file (.gbff/.gbk/.gb) per line. Genome name = file name without extension.
#   OUTDIR      : created; OUTDIR/panaroo (Panaroo output), OUTDIR/pgom (PanGene-O-Meter output), OUTDIR/timing.tsv
#   THREADS     : default 8
#   CLEAN_MODE  : strict (default) | moderate | sensitive
# environment variables:
#   PGOM_PANAROO_ENV : path of the conda environment that contains panaroo (default: panaroo must be in PATH)
#   PGOM_SCRIPT      : path of PanGeneOmeter.pl (default: ../../PanGeneOmeter/PanGeneOmeter.pl relative to this script)
#   PGOM_EXTRA_ARGS  : extra arguments passed to PanGeneOmeter.pl (e.g. "--GCS_clustering_cutoff 0.95")
set -euo pipefail
LIST="${1:?usage: $0 GENOME_LIST OUTDIR [THREADS] [CLEAN_MODE]}"
OUT="${2:?usage: $0 GENOME_LIST OUTDIR [THREADS] [CLEAN_MODE]}"
THREADS="${3:-8}"
MODE="${4:-strict}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PGOM_SCRIPT="${PGOM_SCRIPT:-$HERE/../../PanGeneOmeter/PanGeneOmeter.pl}"

[ -s "$LIST" ] || { echo "[ERROR] genome list '$LIST' is missing or empty" >&2; exit 1; }
[ -f "$PGOM_SCRIPT" ] || { echo "[ERROR] PanGeneOmeter.pl not found at '$PGOM_SCRIPT' (set PGOM_SCRIPT)" >&2; exit 1; }
case "$MODE" in strict|moderate|sensitive) ;; *) echo "[ERROR] CLEAN_MODE must be strict, moderate or sensitive" >&2; exit 1;; esac
# genome names are the file names without extension: they must be unique
dups=$(sed 's#.*/##; s#\.[^.]*$##' "$LIST" | sort | uniq -d)
[ -z "$dups" ] || { echo "[ERROR] duplicate genome names (file names without extension): $dups" >&2; exit 1; }
while read -r g; do [ -z "$g" ] || [ -s "$g" ] || { echo "[ERROR] listed file not found: $g" >&2; exit 1; }; done < "$LIST"

LIST="$(cd "$(dirname "$LIST")" && pwd)/$(basename "$LIST")"
mkdir -p "$OUT"; OUT="$(cd "$OUT" && pwd)"
if [ -e "$OUT/panaroo" ]; then echo "[ERROR] '$OUT/panaroo' already exists; Panaroo needs a new directory (remove it or choose another OUTDIR)" >&2; exit 1; fi

if [ -n "${PGOM_PANAROO_ENV:-}" ]; then PANAROO=(conda run --no-capture-output -p "$PGOM_PANAROO_ENV" panaroo)
else command -v panaroo >/dev/null || { echo "[ERROR] panaroo not found: activate its environment or set PGOM_PANAROO_ENV" >&2; exit 1; }; PANAROO=(panaroo); fi

echo -e "step\tseconds\tversion" > "$OUT/timing.tsv"
PV=$("${PANAROO[@]}" --version 2>&1 | tail -1)
# NOTE: GenBank files are only converted by Panaroo when -i gets ONE list file (several files given directly are read as GFF3)
t0=$(date +%s)
"${PANAROO[@]}" -i "$LIST" -o "$OUT/panaroo" --clean-mode "$MODE" --remove-invalid-genes -t "$THREADS" 2>&1 | tee "$OUT/panaroo.log" | tail -n 5
t1=$(date +%s)
RTAB="$OUT/panaroo/gene_presence_absence.Rtab"
[ -s "$RTAB" ] || { echo "[ERROR] Panaroo did not produce $RTAB; see $OUT/panaroo.log" >&2; exit 1; }
echo -e "panaroo_${MODE}\t$((t1-t0))\t$PV" >> "$OUT/timing.tsv"

# shellcheck disable=SC2086
perl "$PGOM_SCRIPT" --phyletic_pattern_file "$RTAB" --outDir "$OUT/pgom" --threads "$THREADS" ${PGOM_EXTRA_ARGS:-} 2>&1 | tee "$OUT/pgom.log" | tail -n 3
t2=$(date +%s)
echo -e "pangeneometer_matrix\t$((t2-t1))\t-" >> "$OUT/timing.tsv"
echo "Done. Panaroo: $OUT/panaroo   PanGene-O-Meter: $OUT/pgom   timing: $OUT/timing.tsv"
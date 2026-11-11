# Panaroo -> PanGene-O-Meter protocol

Panaroo is an external pan-genome tool. PanGene-O-Meter reads its presence/absence matrix (`gene_presence_absence.Rtab`) through the matrix input; nothing of Panaroo is installed with PanGene-O-Meter.

## 1. Install Panaroo in its own environment (tested with Panaroo 1.8.1)
```bash
conda create -y -p /path/to/panaroo_env -c conda-forge -c bioconda panaroo
conda activate /path/to/panaroo_env
panaroo --version
```

## 2. Prepare a genome list
One GenBank file (`.gbff`, `.gbk` or `.gb`) per line. The genome name is the file name without the extension and must be unique:
```
/data/gbff/Ecoli_K12_MG1655.gbff
/data/gbff/Ecoli_K12_DH10B.gbff
```
**Important:** give Panaroo the list file itself (`-i list.txt`). GenBank files are converted to GFF3 only in this case; if several files are given directly Panaroo reads them as GFF3 and fails with "Error reading prokka input!".

## 3. Run Panaroo and PanGene-O-Meter
```bash
conda activate /path/to/panaroo_env
export PGOM_SCRIPT=/path/to/PanGene-O-Meter/PanGeneOmeter/PanGeneOmeter.pl   # optional if run from the repository
bash Protocols/Panaroo/pgom_run_panaroo.sh genomes.txt out_panaroo 8 strict
```
Arguments: genome list, output directory, threads (default 8), Panaroo clean mode (`strict` default, `moderate`, `sensitive`). The script runs
`panaroo -i genomes.txt -o OUT/panaroo --clean-mode strict --remove-invalid-genes -t 8` and then
`PanGeneOmeter.pl --phyletic_pattern_file OUT/panaroo/gene_presence_absence.Rtab --outDir OUT/pgom`.
`OUT/timing.tsv` holds the running time of each step. The Panaroo messages "invalid gene ... Has stop: True" are expected: `--remove-invalid-genes` drops annotations with internal stop codons (pseudogenes).

Run PanGene-O-Meter alone on an existing Panaroo result:
```bash
perl PanGeneOmeter.pl --phyletic_pattern_file out_panaroo/panaroo/gene_presence_absence.Rtab --outDir out_pgom
```
The matrix may also come from Roary (`gene_presence_absence.Rtab`) or Microbializer (`phyletic_patterns.fas`).

## 4. Comparing two runs
```bash
python3 Protocols/pgom_compare_runs.py DIAMOND out_diamond/*.all_pairs.csv Panaroo out_panaroo/pgom/*.all_pairs.csv --dataset Ecoli8 --out comparison.tsv
```
Reports Pearson/Spearman correlation and mean absolute difference of the Jaccard (GCSj) and overlap (GCSo) values over all genome pairs.
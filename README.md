# PanGene-O-Meter
Quantifying bacterial genomes' gene-content similarity and gene-content-based dereplication

Introduction
============
Bacterial genome evolution is shaped to a great extent by horizontal gene transfer, detectable as genes with a presence-absence pattern of variation that does not follow phylogenetic relationships across the majority of the genome. While average nucleotide identity (ANI) is the most common measure to measure bacterial genome similarity, its ability to capture gene-content differences is limited, and its use therefore often misses a major factor for functional variability. **PanGene-O-Meter provides a complementary gene-content-based method to quantify similarity between bacterial genomes, leveraging knowledge of pan-genome structure and relationships between orthology groups.** PanGene-O-Meter, provides a higher level of granularity than ANI when classifying bacterial isolates of the same species. Furthermore, PanGene-O-Meter is useful for efficient bacterial genomes clustering, allowing for facile selection of one or several 'representative' genomes, effectively dereplicating a set of genomes by their gene content. 

## Table of contents
  * [Pipeline overview](#pipeline-overview)
  * [Quick start](#quick-start)
  * [How to run](#how-to-run)
  * [Command line arguments](#command-line-arguments)
  * [Directory structure and analysis output](#output-files)
  * [Versions](#versions)
  * [Citing and Credit](#Citing-and-Credit)

## Pipeline overview
1.	Orthology group assignment: Assign predicted genes to orthology groups (pangenes) to construct a pan-genome. This step is done by PanGene-O-Meter (proteins from GenBank files, or genes predicted from genome FASTA files, clustered with DIAMOND DeepClust), or it is skipped when you provide the presence/absence matrix of another pan-genome tool.
2.	Phyletic pattern representation: Represent each genome as a binary presence/absence (P/A) vector, where each element corresponds to a pangene and the value is either 1 (present) or 0 (absent).
3.	Gene-content similarity (GCS) calculation: Calculate GCS based on Jaccard similarity \($GCSj$\) between the P/A vectors of two bacterial genomes \($pS1$ and $pS2$\):

$$GCSj(pS1,pS2)=\frac{M11}{(M01+M10+M11)}$$

Where $M11$ represents the total number of elements (orthology-groups) for which both $pS1$ and $pS2$ have a value of 1; $M01$ represents the total number of elements for which the value of $pS1$ is $0$ and $pS2$ is $1$; and $M10$ represents the total number of elements for which the value of $pS1$ is $1$ and the value of $pS2$ is $0$.

Alternatively, the GCS can be defined by the shared pangenes, based on their overlap coefficient \($GCSo$\).  This approach reduces the impact of partial genomes on the GCS metric.

$$GCSo(pS1,pS2)=\frac{M11}{(min(M01+M11,M10+M11))}$$

4. **GeneContRep**, a greedy incremental clustering algorithm to dereplicate a list of genomes. GeneContRep, allows the selection of representative genomes according to specific gene-content similarity cutoffs.

## Quick start
### Installing
PanGene-O-Meter needs Perl, [DIAMOND](https://github.com/bbuchfink/diamond) (tested with **version 2.1.10**) and, for genome FASTA input, [pyrodigal](https://github.com/althonos/pyrodigal). All of them are installed in one conda environment. PanX is **not** needed (see [Optional: PanX](#optional-panx)).

0. If conda is not already installed, we recommend installing Miniforge (community-maintained, uses the free conda-forge channel): https://github.com/conda-forge/miniforge
Any recent conda distribution works; we have tested installation with Miniforge and Miniconda3. PanGene-O-Meter does not depend on a specific conda version.
1. Download PanGene-O-Meter
```
git clone https://github.com/HaimAshk/PanGene-O-Meter.git
cd PanGene-O-Meter
```
2. Create the environment (named `pangeneometer`; `mamba` can be used instead of `conda` for a faster installation)
```
conda env create -f environment.yml
```
**OR** [replace /some/where/ with your desired conda environment path]
```
conda env create -p /some/where/pangeneometer -f environment.yml
```
3. Activate the environment and test the installation
```
conda activate pangeneometer
## OR ##
conda activate /some/where/pangeneometer

perl PanGeneOmeter/PanGeneOmeter.pl --help
diamond version
```

Notes:
* The `PGOM_*.pm` files must stay in the same folder as `PanGeneOmeter.pl` (they are found automatically). Do not move `PanGeneOmeter.pl` alone; to run it from any folder, either call it with its full path, or make it executable and link it: `chmod +x /full/path/PanGeneOmeter/PanGeneOmeter.pl; ln -s /full/path/PanGeneOmeter/PanGeneOmeter.pl $CONDA_PREFIX/bin/PanGeneOmeter.pl`.
* Please use the DIAMOND version of the environment. DIAMOND v2.1.24 crashes in the `recluster` step and v2.1.25 has no `recluster` command, so these versions can not be used at the moment; PanGene-O-Meter stops with an error if DIAMOND returns an empty result.

## How to run
PanGene-O-Meter needs the presence/absence of orthology groups in every genome. There are three ways to provide it. In all of them, `--outDir` is the output directory and `--threads` sets the number of CPU cores \[default=8\].

**1. Annotated genomes in [*GenBank Flat File*](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/) format** (a text file with one genome file per line, `.gz` files are accepted):

``perl PanGeneOmeter/PanGeneOmeter.pl --in_gbk_list_file list_of_gbk --outDir output_directory --threads 16``

**2. Genome sequences in FASTA format** (one file per genome, genes are predicted with [pyrodigal](https://github.com/althonos/pyrodigal)):

``perl PanGeneOmeter/PanGeneOmeter.pl --in_genome_list_file list_of_fna --input_type fna --outDir output_directory --threads 16``

**3. A presence/absence matrix produced by another tool** (e.g. Microbializer, Panaroo, Roary). No proteins are extracted and no clustering of genes is done:

``perl PanGeneOmeter/PanGeneOmeter.pl --phyletic_pattern_file matrix.fas --outDir output_directory``

* To see all supported parameters run: ` perl PanGeneOmeter/PanGeneOmeter.pl --help `

### Presence/absence matrix formats (`--phyletic_pattern_file`)
The format is detected automatically; use `--pattern_format auto|rtab|roary|fasta|csv` to set it.

| Source | File | Format |
|---|---|---|
| Microbializer | `phyletic_patterns.fas` | FASTA, one record per genome; the sequence is a string of `0`/`1` (one character per orthology group). The first word of the header is the genome name. Wrapped lines are fine |
| Panaroo | `gene_presence_absence.Rtab` | tab-delimited, one row per gene family, one column per genome (0/1) |
| Roary | `gene_presence_absence.Rtab` or `gene_presence_absence.csv` | as above; the `.csv` is converted to presence/absence |
| PanGene-O-Meter or other | `*.01.csv` | comma-delimited; the first column is the genome name, the header has the family names, values 0/1 |

Requirements: at least 2 genomes, the same number of orthology groups for every genome, and only `0`/`1` values (paralog counts must be converted to presence/absence first). The orthology-group names of a FASTA matrix are not stored in the file; they are called `OG1 ... OGn`.

With a matrix there are no genome lengths. When representatives are chosen, genomes are ordered by their number of present genes (instead of their genome length), so within a cluster of near-identical genomes the representative can differ from a run on the same data with genome input. The GCS distances and the cluster members are the same.

### Genome names
Whatever the input, genomes are renamed to `S1 ... SN` (in list order) for the analysis and renamed back in the result tables (phyletic patterns, pairwise distances, clusters, representatives and metadata). The mapping is saved in `genome_names.map.txt` (original file, original name, `S<n>`). Intermediate files (`All_proteins.faa`, `protein_faa/`, the DIAMOND `.tsv` files) keep the `S<n>` names. Genome names must be unique: two genomes with the same base name stop the run with an error.

## Command line arguments
### Mandatory parameters:
* `--outDir <Directory>`: Provide the full path to the directory where output files will be saved.
* **One** of the following inputs:
  * `--in_gbk_list_file <File>` (alias `--in_genome_list_file`): Specify the path to the file containing a list of genomes to be analyzed, one file per line. Each genome is in [*GenBank Flat File*](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/) format (default), or in FASTA format with `--input_type fna`.
  * `--phyletic_pattern_file <File>`: A presence/absence matrix (see [formats](#presenceabsence-matrix-formats---phyletic_pattern_file)).

### Additional (optional) parameters:
#### General Parameters:
* `--input_type <gbff|fna>`: The format of the genomes in the list file: GenBank (`gbff`, default) or genome sequences in FASTA (`fna`; genes are predicted with pyrodigal). (The values `faa` and `single_faa`, protein FASTA input, are experimental and not described here.)
* `--prefix_name <Prefix>`: Specify a prefix to be added to the names of the output files. Default is an empty string (""), meaning no prefix will be used.
* `--metadata_file <File>` (alias `--in_metadata_file`): Provide the path to a tab-delimited file containing metadata and descriptions for the genomes being analyzed; the first column is the genome name. Default is an empty string (""). If not provided, a table with the genome source, length and number of proteins is used (with `--pangenome_alg PanX` the metadata extracted by PanX from the GenBank files is used).
* `--ignore_private_genes`: Ignore private genes (present in a single genome) in the GCS scores. Default: private genes are used. The excluded groups are listed in `*.filtered_pos_below_1.txt`.
* `--skip_all_vs_all`: Skip the calculation of all pairwise GCS distances (the distances needed for clustering are calculated anyway).
* `--pa_engine <fast|legacy>`: Implementation of the phyletic-pattern and all-vs-all steps. `fast` (default) uses packed bit vectors and `--threads`; `legacy` is the previous implementation. Both give the same results (see the QC in the repository);
* `--threads <Number>`: Set the number of threads to use for parallel processing. The default is `8`.
* `--help`: Display help information and exit the program.

#### PanGenome Construction Parameters (not used with `--phyletic_pattern_file`):
* `--pangenome_alg <Algorithm>`: Choose the algorithm for PanGenome construction. Options are `DIAMONDClust` or `PanX`. The default option is `DIAMONDClust`. `PanX` needs a separate installation (see [Optional: PanX](#optional-panx)).
* `--diamond_pid <Percentage>`: Set the percentage identity for clustering sequences when using `DIAMONDClust`. The default is `80`.
* `--diamond_cover <Coverage Cutoff>`: Define the mutual coverage cutoff for use in `DIAMONDClust`. The default value is `70`.
* `--protein_extractor <builtin|panx>`: How proteins are taken from the GenBank files for `DIAMONDClust`. The default `builtin` needs no PanX.
* `--reuse_PanGenome`: If this option is specified, the PanGenome construction step will be skipped, and existing files from previous runs in the output directory will be reused. By default, the program will override and execute all steps anew.
* `--pattern_format <auto|rtab|roary|fasta|csv>`: The format of `--phyletic_pattern_file`. The default is `auto`.

#### *GeneContRep* Clustering Parameters:
* `--GCS_method <Method>`: Select the gene-content similarity metric to be used in the *GeneContRep* step. Options are `GCSj` or `GCSo`, with the default being `GCSj`.
* `--GCS_clustering_cutoff <Cutoff>` (alias `--CGS_clustering_cutoff`): Set the cutoff value for gene-content similarity clustering during the *GeneContRep* step. The value should be between `0` and `1`, with the default being `0.9`.

## Output files
Output files are written to `--outDir`. The names start with `DIAMOND_DeepClust.` (genome input), `PanX.` (`--pangenome_alg PanX`) or `InputMatrix.` (matrix input), after the optional `--prefix_name`.
#### Phyletic Patterns
* **.PhyleticPattern.01.fas**: A FASTA format file representing the Phyletic pattern (0/1). Each site (or position in sequence) corresponds to an orthology group, as specified in `PhyleticPattern.clusters_names_and_pos.txt` (matrix input: `PhyleticPattern.og_names.txt`).
* **.PhyleticPattern.clusters_names_and_pos.txt**: Contains the names and positions of orthology groups as detailed in the FASTA Phyletic Pattern files.
* **.PhyleticPattern.01.csv**: Presents the Phyletic pattern as a CSV matrix, indicating the presence (1) or absence (0) of orthology groups across genomes.
* **.PhyleticPattern.csv**: Displays the Phyletic pattern as a CSV table, where multi-copy orthology groups are represented by the number of copies observed in each genome.
* **.PhyleticPattern.clusters_names_and_NumOfSpecies.txt**: Lists orthology groups alongside the number of genomes in which each group is found, and the average number of copies per genome. This file provides insights into the distribution and conservation of gene clusters.

(With matrix input only the `.01.fas`, `.01.csv` and `og_names.txt` files are written.)
#### Pairwise Distances
* **.distance.all_pairs.and_metadata.csv**: This CSV file includes pairwise distance measurements for all genome pairs, calculated using the gene-content similarity metrics `GCSj` and `GCSo`. Additionally, it incorporates metadata for each genome pair, offering quantitative insights into the gene-content similarity along with contextual information for each pair of genomes.
#### Clusters
* **.clstr.txt**: Describes cluster assignments based on a specified similarity cutoff. Each cluster comprises genomes grouped according to the gene-content similarity metric, sharing similarity above the selected cutoff in relation to the cluster representative.
* **.clusters_annotation_and_metadata.txt**: Provides details on the cluster membership of each genome, including associated metadata to present additional contextual information regarding each genome's biological and ecological characteristics.
* **.representative_PhP.csv**: Contains phyletic pattern data formatted as a CSV matrix for the representative genomes selected for each cluster. This file provides insights into orthology group presence or absence across cluster representatives.
* **.representative_PhP.fas**: Provides phyletic pattern data in FASTA format for the representative genomes selected for each cluster. Each position in the FASTA sequence corresponds to an orthology group, with detailed information available in the PhyleticPattern.clusters_names_and_pos.txt file.
* **.clusters_sum.txt**: Indicates the representative sequence ID for each cluster, along with the number of genome members in each cluster and their respective genome IDs. This file provides a concise overview of the composition of each cluster.
#### Other files
* **genome_names.map.txt**: The mapping between the original genome files and names and the internal names `S1 ... SN`.
* **genome_stats.tsv**: For each genome (original name; `internal_id` gives the S1..SN name used inside intermediate files), the source of the proteins (GenBank or pyrodigal), the genome length and the number of proteins (matrix input: the number of present genes). The clustering log (`*.log.txt`) also lists the original names.
* **All_proteins.faa**, **protein_faa/**, **All_proteins.diamond_cluster...tsv**: Intermediate files of the DIAMOND step (proteins are named `S<n>|locus`).

## Optional: PanX
PanX is not needed for the default analysis. `--pangenome_alg PanX` (and `--protein_extractor panx`) use the [PanX](https://github.com/neherlab/pan-genome-analysis/) pan-genome pipeline, which needs Python 2.7 and its old dependencies. Install it in **its own** conda environment, not in the `pangeneometer` one:

<details>
<summary>PanX installation (legacy)</summary>

```
conda create -n PanX python==2.7.*
conda activate PanX
git clone https://github.com/HaimAshk/PanGene-O-Meter.git $CONDA_PREFIX/PanGeneOmeter_github
git clone https://github.com/neherlab/pan-genome-analysis.git $CONDA_PREFIX/panx_github
conda env update --file $CONDA_PREFIX/PanGeneOmeter_github/PanGeneOmeter/PanGeneOmeter-environment.yml --prune
cp $CONDA_PREFIX/PanGeneOmeter_github/PanX_Patch/sf_extract_sequences.py $CONDA_PREFIX/panx_github/scripts/
```
Run PanGene-O-Meter from this environment with the full path of the script (do not move it), for example:
`perl $CONDA_PREFIX/PanGeneOmeter_github/PanGeneOmeter/PanGeneOmeter.pl --in_gbk_list_file list_of_gbk --outDir out --pangenome_alg PanX`.
If PanX is installed elsewhere, set the environment variable `PGOM_PANX` to the full path of `panX.py`.

If pip installation fails \[previous steps ends with error like: `CondaEnvException: Pip failed` or typing `pip` results in an error like: `ImportError: No module named pip._internal.cli.main`\]
```
python2.7 -m ensurepip --default-pip
conda env update --file $CONDA_PREFIX/PanGeneOmeter_github/PanGeneOmeter/PanGeneOmeter-environment.yml --prune
```
PanX dependencies: [PanX](https://github.com/neherlab/pan-genome-analysis/), [DIAMOND](https://github.com/bbuchfink/diamond), [MCL](http://micans.org/mcl/), [mafft](http://mafft.cbrc.jp/alignment/software/), [fasttree](http://www.microbesonline.org/fasttree/), [raxml](https://github.com/stamatak/standard-RAxML), [treetime](http://github.com/neherlab/treetime).
</details>

## Versions
| Tag / branch | Description |
|---|---|
| `v1.2.0` | DIAMOND DeepClust is the default (no PanX needed); GenBank proteins are extracted by PanGene-O-Meter; genome FASTA input (pyrodigal) and presence/absence matrices as input |
| `v1.1.0` | merge of two working copies with start-up fixes; PanX was the default |
| `v1.0.0` | the original release |

Older versions remain available, e.g. `git checkout v1.0.0`.

## Citing and Credit
**When using the tool in published research, please cite:**
-   Ashkenazy H and Weigel D,<br>
    \"PanGene-O-Meter: Intra-Species Diversity Based on Gene-Content\",<br>
    *bioRxiv* 2025.09.01.673544;
    [doi: 10.1101/2025.09.01.673544](https://doi.org/10.1101/2025.09.01.673544)

    **When using `--pangenome_alg DIAMONDClust` (the default) please also cite:**
    -   Buchfink, B.J., Barbé, É., Ashkenazy, H., Reuter K., Kennedy J.A., Drost H.G.,<br>
        \"Clustering the protein universe of life using DIAMOND DeepClust\",<br>
        *Nat Methods* **23, 724–727 (2026)**;
        [[doi: 10.1038/s41592-026-03030-z](https://doi.org/10.1038/s41592-026-03030-z)]

    **When using `--pangenome_alg PanX` please also cite:**
    -   Ding W, Baumdicker F, Neher RA,<br>
        \"panX: pan-genome analysis and exploration\",<br>
        *Nucleic Acids Research*,**46(1):e5**
        [doi: 10.1093/nar/gkx977](https://doi.org/10.1093/nar/gkx977).

    **When using `--input_type fna` please also cite:**
    -   Larralde M,<br>
        \"Pyrodigal: Python bindings and interface to Prodigal, an efficient method for gene prediction in prokaryotes\",<br>
        *Journal of Open Source Software* **7(72), 4296 (2022)**;
        [doi: 10.21105/joss.04296](https://doi.org/10.21105/joss.04296)
    -   Hyatt D, Chen G-L, LoCascio PF, Land ML, Larimer FW, Hauser LJ,<br>
        \"Prodigal: prokaryotic gene recognition and translation initiation site identification\",<br>
        *BMC Bioinformatics* **11, 119 (2010)**;
        [doi: 10.1186/1471-2105-11-119](https://doi.org/10.1186/1471-2105-11-119)

    **When using `--phyletic_pattern_file`, please also cite the tool that produced the matrix.**

Support
=======
PanGene-O-Meter is actively supported and developed software. Please use the [issue tracker](https://github.com/HaimAshk/PanGene-O-Meter/issues) for malfunctions and the [GitHub discussions](https://github.com/HaimAshk/PanGene-O-Meter/discussions) for questions, comments, feature requests, etc.
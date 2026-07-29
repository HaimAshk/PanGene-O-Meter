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
  * [Citing and Credit](#Citing-and-Credit)

## Pipeline overview
1.	Orthology group assignment: Assign predicted genes to orthology groups (pangenes) to construct a pan-genome (different methods can be used for this step).
2.	Phyletic pattern representation: Represent each genome as a binary presence/absence (P/A) vector, where each element corresponds to a pangene and the value is either 1 (present) or 0 (absent).
3.	Gene-content similarity (GCS) calculation: Calculate GCS based on Jaccard similarity \($GCSj$\) between the P/A vectors of two bacterial genomes \($pS1$ and $pS2$\):

$$GCSj(pS1,pS2)=\frac{M11}{(M01+M10+M11)}$$

Where $M11$ represents the total number of elements (orthology-groups) for which both $pS1$ and $pS2$ have a value of 1; $M01$ represents the total number of elements for which the value of $pS1$ is $0$ and $pS2$ is $1$; and $M10$ represents the total number of elements for which the value of $pS1$ is $1$ and the value of $pS2$ is $0$.

Alternatively, the GCS can be defined by the shared pangenes, based on their overlap coefficient \($GCSo$\).  This approach reduces the impact of partial genomes on the GCS metric.

$$GCSo(pS1,pS2)=\frac{M11}{(min(M01+M11,M10+M11))}$$

4. **GeneContRep**, a greedy incremental clustering algorithm to dereplicate a list of genomes. GeneContRep, allows the selection of representative genomes according to specific gene-content similarity cutoffs.

## Quick start
### Installing
0. If conda is not already installed, we recommend installing Miniforge (community-maintained, uses the free conda-forge channel): https://github.com/conda-forge/miniforge
Any recent conda distribution works; we have tested installation with Miniforge and Miniconda3. PanGene-O-Meter does not depend on a specific conda version.
1. Create a new environment for PanGeneOmeter
```
conda create -n PanGeneOmeter python==2.7.*
```
**OR** [replace /some/where/ with your desired conda environment path]
```
conda create -p /some/where/PanGeneOmeter python==2.7.*
```
2. Activate the conda environment
```
conda activate PanGeneOmeter
## OR ##
conda activate /some/where/PanGeneOmeter
```
3. Download and install PanGeneOmeter and dependencies
```
git clone https://github.com/HaimAshk/PanGene-O-Meter.git $CONDA_PREFIX/PanGeneOmeter_github
git clone https://github.com/neherlab/pan-genome-analysis.git $CONDA_PREFIX/panx_github
conda env update --file $CONDA_PREFIX/PanGeneOmeter_github/PanGeneOmeter/PanGeneOmeter-environment.yml --prune 
mv $CONDA_PREFIX/PanGeneOmeter_github/PanGeneOmeter/PanGeneOmeter.pl $CONDA_PREFIX/bin
cp $CONDA_PREFIX/PanGeneOmeter_github/PanX_Patch/sf_extract_sequences.py $CONDA_PREFIX/panx_github/scripts/
```

Note: PanGene-O-Meter itself is Perl-based and does not require Python. Python 2.7 is needed  for PanX, which installs and runs without issue in a dedicated environment (e.g. via conda/mamba). Please follow the instructions above to create a dedicated environment for PanGene-O-Meter.

### Installation troubleshooting
* If pip installation fails \[previous steps ends with error like: `CondaEnvException: Pip failed` or typing `pip` results in an error like: `ImportError: No module named pip._internal.cli.main`\]


  ```
  python2.7 -m ensurepip --default-pip
  conda env update --file $CONDA_PREFIX/PanGeneOmeter_github/PanGeneOmeter/PanGeneOmeter-environment.yml --prune
  ```
  
#### Overview of dependencies (installed in the environment):
  * [PanX](https://github.com/neherlab/pan-genome-analysis/) 
  * [DIAMOND](https://github.com/bbuchfink/diamond)
  * [MCL](http://micans.org/mcl/)
  * [mafft](http://mafft.cbrc.jp/alignment/software/)
  * [fasttree](http://www.microbesonline.org/fasttree/)
  * [raxml](https://github.com/stamatak/standard-RAxML)
  * [treetime](http://github.com/neherlab/treetime)

## How to run
PanGene-O-Meter requires a file containing a list of annotated genomes in [*GenBank Flat File*](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/) format for analysis, 
<br>specified with the `--in_gbk_list_file` argument, and an output directory, designated by the `--outDir` argument. 
<br>`--threads ` sets the number of CPU cores \[default=1\].

**So for example:**

``PanGeneOmeter.pl --in_gbk_list_file list_of_gbk  --outDir  output_directory --threads 64``

* To see all supported parameters run: ` PanGeneOmeter.pl -h `

## Command line arguments
### Mandatory parameters:
* `--in_gbk_list_file <File>`: Specify the path to the file containing a list of genomes to be analyzed, each in [*GenBank Flat File*](https://www.ncbi.nlm.nih.gov/genbank/samplerecord/) format.
* `--outDir <Directory>`: Provide the full path to the directory where output files will be saved.
### Additional (optional) parameters:
#### General Parameters:
* `--prefix_name <Prefix>`: Specify a prefix to be added to the names of the output files. Default is an empty string (""), meaning no prefix will be used.
* `--metadata_file <File>`: Provide the path to a tab-delimited file containing metadata and descriptions for the genomes being analyzed. Default is an empty string (""). If not provided, metadata will be automatically extracted from the GenBank file.
* `--threads <Number>`: Set the number of threads to use for parallel processing. The default is `1`, meaning the program will use a single thread.
* `--help`: Display help information and exit the program.

#### PanGenome Construction Parameters:
* `--pangenome_alg <Algorithm>`: Choose the algorithm for PanGenome construction. Options are `PanX` or `DIAMONDClust`. The default option is `PanX`.
* `--diamond_pid <Percentage>:` Set the percentage identity for clustering sequences when using `DIAMONDClust`. The default is `80`.
* `--diamond_cover <Coverage Cutoff>`: Define the mutual coverage cutoff for use in `DIAMONDClust`. The default value is `70`.
* `--reuse_PanGenome`: If this option is specified, the PanGenome construction step will be skipped, and existing files from previous runs in the output directory will be reused. By default, the program will override and execute all steps anew.

#### *GeneContRep* Clustering Parameters:
* `--GCS_method <Method>`: Select the gene-content similarity metric to be used in the *GeneContRep* step. Options are `GCSj` or `GCSo`, with the default being `GCSj`.
* `--CGS_clustering_cutoff <Cutoff>`: Set the cutoff value for gene-content similarity clustering during the *GeneContRep* step. The value should be between `0` and `1`, with the default being `0.9`.

## Citing and Credit
**When using the tool in published research, please cite:**
-   Ashkenazy H and Weigel D,<br>
    \"PanGene-O-Meter: Intra-Species Diversity Based on Gene-Content\",<br>
    *bioRxiv* 2025.09.01.673544;
    [doi: 10.1101/2025.09.01.673544](https://doi.org/10.1101/2025.09.01.673544)

    **When using `--pangenome_alg PanX` (or for default run) please also cite:**
    -   Ding W, Baumdicker F, Neher RA,<br>
        \"panX: pan-genome analysis and exploration\",<br>
        *Nucleic Acids Research*,**46(1):e5**
        [doi: 10.1093/nar/gkx977](https://doi.org/10.1093/nar/gkx977).
 
    **When using `--pangenome_alg DIAMONDClust` please also cite:**
    -   Buchfink, B.J., Barbé, É., Ashkenazy, H., Reuter K., Kennedy J.A., Drost H.G.,<br>
        \"Clustering the protein universe of life using DIAMOND DeepClust\",<br>
        *Nat Methods* **23, 724–727 (2026)**;
        [[doi: 10.1038/s41592-026-03030-z](https://doi.org/10.1038/s41592-026-03030-z)]

## Output files
#### Phyletic Patterns
* **.PhyleticPattern.fas**: A FASTA format file representing the Phyletic pattern. Each site (or position in sequence) corresponds to an orthology group, as specified in `PhyleticPattern.clusters_names_and_pos.txt`.
* **.PhyleticPattern.clusters_names_and_pos.txt**: Contains the names and positions of orthology groups as detailed in the FASTA Phyletic Pattern files.
* **.PhyleticPattern.01.csv**: Presents the Phyletic pattern as a CSV matrix, indicating the presence (1) or absence (0) of orthology groups across genomes.
* **.PhyleticPattern.csv**: Displays the Phyletic pattern as a CSV table, where multi-copy orthology groups are represented by the number of copies observed in each genome.
* **.PhyleticPattern.clusters_names_and_NumOfSpecies.txt**: Lists orthology groups alongside the number of genomes in which each group is found, and the average number of copies per genome. This file provides insights into the distribution and conservation of gene clusters.
#### Pairwise Distances
* **.distance.all_pairs.and_metadata.csv**: This CSV file includes pairwise distance measurements for all genome pairs, calculated using the gene-content similarity metrics `GCSj` and `GCSo`. Additionally, it incorporates metadata for each genome pair, offering quantitative insights into the gene-content similarity along with contextual information for each pair of genomes.
#### Clusters
* **.clstr.txt**: Describes cluster assignments based on a specified similarity cutoff. Each cluster comprises genomes grouped according to the gene-content similarity metric, sharing similarity above the selected cutoff in relation to the cluster representative.
* **.clusters_annotation_and_metadata.txt**: Provides details on the cluster membership of each genome, including associated metadata to present additional contextual information regarding each genome's biological and ecological characteristics.
* **.representative_PhP.csv**: Contains phyletic pattern data formatted as a CSV matrix for the representative genomes selected for each cluster. This file provides insights into orthology group presence or absence across cluster representatives.
* **.representative_PhP.fas**: Provides phyletic pattern data in FASTA format for the representative genomes selected for each cluster. Each position in the FASTA sequence corresponds to an orthology group, with detailed information available in the PhylleticPattern.clusters_names_and_pos.txt file.
* **.clusters_sum.txt**: Indicates the representative sequence ID for each cluster, along with the number of genome members in each cluster and their respective genome IDs. This file provides a concise overview of the composition of each cluster.


Support
=======
PanGene-O-Meter is actively supported and developed software. Please use the [issue tracker](https://github.com/HaimAshk/PanGene-O-Meter/issues) for malfunctions and the [GitHub discussions](https://github.com/HaimAshk/PanGene-O-Meter/discussions) for questions, comments, feature requests, etc.

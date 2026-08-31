#!/usr/bin/env perl

use strict;
use warnings;
use File::Copy;
use File::Path qw(make_path remove_tree);
use File::Basename;
use List::Util qw(min);
use Getopt::Long;
use FindBin qw($Bin);
use Archive::Tar;
use Cwd;
use POSIX qw(strftime);
# use IO::Zlib;
# use IO::Compress::Gzip qw(gzip $GzipError);

my $PANX_PATH=$ENV{CONDA_PREFIX}."/panx_github/panX.py";

# my $PANX_PATH="/ebio/abt6_projects/small_projects/hashkenazy/Programs/PanX/pan-genome-analysis";
# my $DIAMOND_CLUST_PATH="/ebio/abt6_projects/small_projects/hashkenazy/Programs/DIAMOND_v2.1.9/diamond";
my $DIAMOND_CLUST_PATH="diamond";

# PanGene-O-Meter


my $help_usage="USAGE: perl $0 --in_gbk_list_file <file with list of genomes, each in genebank format> --outDir <full path output direcory>
 \tOptional parameters:
 \t--threads <number of threads> default=1
 \t--prefix_name <Prefix to add to the output files> default=\"\"
 \t--metadata_file  <tab delimited file with metadata and descriptions for the analyzed genomes> default=\"\"
 \t--help --> print help and exit
  \t[PanGenome construction]
 \t\t--pangenome_alg <PanX | DIAMONDClust> default=PanX
 \t\t--diamond_pid <percent identity to cluster sequence with DIAMONDClust> default=80
 \t\t--diamnod_cover <mutual-coverage cutoff to be used in DIAMONDClust>, default=70
 \t\t--reuse_PanGenome --> if provided the PanGenome construction step is skipped and the available files in the out dir (from previous runs) are used. default: override and start all steps.
 \t[GeneContRep clustering]
 \t\t--GCS_method <GCSj |GCSJo>, default=GCSJ --> gene-content similarity metric to use in the GeneContRep step
 \t\t--CGS_clustering_cutoff <{0-1}>, default=0.9 --> cutoff for the GCS clustering for GeneContRep step
 \n";

### TO ADD: use calculated PanGenome flag --> otherwise should delete the PanX directory before starting a new one

if (@ARGV<4) {die $help_usage;}
my @orig_ARGV = @ARGV;
# print "[QA INFO] perl $0 @ARGV\n";
my ($in_gbk_list,$out_dir_base,$ignore_private_genes,$threads,$PanGenomeAlg,$skip_all_vs_all,$Prefix_Name,$GCS_method,$clustering_cutoff,$diamond_pid,$diamnod_cover,$reuse_PanGenome,$print_help,$in_metadata_file);


# set default value for optional arguments
$ignore_private_genes=0;
$threads=1;
$print_help=0;
$PanGenomeAlg="PanX";
$skip_all_vs_all=0;
$Prefix_Name="";
$GCS_method="GCSj";
$clustering_cutoff="0.9";
$diamond_pid=80;
$diamnod_cover=70;
$reuse_PanGenome=0;
$in_metadata_file="NA";
my $getoptResult = GetOptions ( "in_gbk_list_file=s"=>\$in_gbk_list,  # list of the genomes and annotations to analyze (each file is expected to be in a GeneBank format) [= means that this parameter is required, s means string]
                                "outDir=s"=>\$out_dir_base,      # path for the output directory

                                # optional parameters
                                "ignore_private_genes!"=>\$ignore_private_genes,     # Flag to indicate private genes should be ignored in the GCS scores, default: noignore_private_genes [Option does not take an argument and may be negated, i.e. prefixed by "no". E.g. "foo!" will allow --foo (with value 1) and -nofoo (with value 0). The option variable will be set to 1, or 0 if negated.
                                "threads:i"=>\$threads,                              # Number of threads to use, default: 8 [: means that this parameter is optional]
                                "pangenome_alg:s"=>\$PanGenomeAlg,                   # Algorithm to use in the PanGenome construction step, {PanX|DIAMONDClust}, default: PanX
                                "skip_all_vs_all!"=>\$skip_all_vs_all,               # flag to skip all-vs-all calculations [Option does not take an argument and may be negated, i.e. prefixed by "no". E.g. "foo!" will allow --foo (with value 1) and -nofoo (with value 0). The option variable will be set to 1, or 0 if negated.]
                                "prefix_name:s"=>\$Prefix_Name,                      # Prefix to add to the output files, default=""
                                "GCS_method:s"=>\$GCS_method,                        # metric to use in the GeneContRep step {GCSj|GCSJo}, default: GCSJ
                                "CGS_clustering_cutoff:s"=>\$clustering_cutoff,      # cutoff for the GCS clustering for GeneContRep step [0-1], default: 0.9
                                "diamond_pid:i"=>\$diamond_pid,                      # percent identity to cluster sequence in DIAMONDClust, default: 80
                                "diamnod_cover:i"=>\$diamnod_cover,                  # mutual-coverage cutoff to be used in DIAMONDClust, default: 70
                                "reuse_PanGenome!"=>\$reuse_PanGenome,               # if provided will skip the PanGenome construction step, otherwise the relevant PanGenome directories and files are DELETED
                                "in_metadata_file:s"=>\$in_metadata_file,            # tab delimited file with metadata and descriptions for the analyzed genomes> default=\"\". Assumes the first column equals to the genome id
                                "help!"=>\$print_help        
);
if ($print_help==1) {print $help_usage; exit(0);}

# arguemnts validation
if (!defined $in_gbk_list) {die "[FATAL ERROR] --in_gbk_list is obligatory argument. The minimal usage is: perl $0 --in_gbk_list GENOME_LIST --outDir OUTPUT_PATH\n";}
if (!defined $out_dir_base) {die "[FATAL ERROR] --outDir is obligatory argument. The minimal usage is: perl $0 --in_gbk_list GENOME_LIST --outDir OUTPUT_PATH\n";}

# optional 
# my $in_fna_list=shift; #STILL NOT SUPPORTED
# my $ignore_private_genes=""; # {Y|[N]}
my $print_OG_diff="N"; # for now don't print the diff
# my $threads=96;
my $raxml_exec="raxml";
my $PanXExtra_params=""; ### ToDO: support from the cmdline
# my $PanGenomeAlg="DIAMONDClust"; #"PanX";
# my $skip_all_vs_all="NO";
# my $Prefix_Name="";
# my $GCS_method="GCSj";
# my $clustering_cutoff;
# my $diamond_pid;
# my $diamnod_cover;
if (!defined $PanGenomeAlg) {$PanGenomeAlg="PANX"} else {$PanGenomeAlg=uc($PanGenomeAlg);}

if (!defined $threads) {$threads=8;}
if (defined $skip_all_vs_all) {if ($skip_all_vs_all==0) {$skip_all_vs_all="NO";} else {$skip_all_vs_all="YES";}} # to correctly handle params
if (!defined $skip_all_vs_all) {$skip_all_vs_all="NO";} else {$skip_all_vs_all=uc($skip_all_vs_all);}
if (!defined $Prefix_Name) {$Prefix_Name="";}
if (!defined $ignore_private_genes) {$ignore_private_genes="N";} else {$ignore_private_genes=uc($ignore_private_genes);}
if (!defined $GCS_method) {$GCS_method="GCSj"} elsif (uc ($GCS_method) eq "GCSJ") {$GCS_method="GCSj";} elsif (uc($GCS_method) eq "GCSO") {$GCS_method="GCSo";} else {print "[UNKNOWN OPTION] --GCS_method can be only GCSj or GCSo, unknown option provided: '$GCS_method' --> exit\n";exit (1);}
if (!defined $clustering_cutoff) {$clustering_cutoff=0.9;}
if (!defined $diamond_pid) {$diamond_pid=80;}
if (!defined $diamnod_cover) {$diamnod_cover=70;}
if (!defined $reuse_PanGenome) {$reuse_PanGenome="NO"} elsif ($reuse_PanGenome==0) {$reuse_PanGenome="NO";} elsif ($reuse_PanGenome==1) {$reuse_PanGenome="YES"} else {print "[UNKNOWN OPTION] --reuse_PanGenome option was not interpreted correctly, [value='$reuse_PanGenome']' --> exit\n";exit (1);}
# add option to provide tsv metadata file --> otherwise will take the one PanX extract
# validate arguments


# output structure
if ($out_dir_base!~/\/$/){$out_dir_base=$out_dir_base."/";}

my $PanX_out_dir=$out_dir_base."PanX/out/";
my $PanX_input_gbk_dir=$PanX_out_dir."input_GenBank/";


my $PanGenomeAlg_for_file_names="";
if ($PanGenomeAlg eq "PANX") {
    $PanGenomeAlg_for_file_names="PanX";
}
elsif ($PanGenomeAlg eq "DIAMONDCLUST") {
    $PanGenomeAlg_for_file_names="DIAMOND_DeepClust";
} else {
    $PanGenomeAlg_for_file_names=$PanGenomeAlg;
}
my $clustering_cutoff_for_print=$clustering_cutoff*100;
my $log_file="";
if (defined $Prefix_Name and $Prefix_Name ne "") {$log_file=$out_dir_base.$Prefix_Name.".PanGeneOmeter.".$PanGenomeAlg_for_file_names.".".$GCS_method.".".$clustering_cutoff_for_print.".log";} 
else {$log_file=$out_dir_base."PanGeneOmeter.".$PanGenomeAlg_for_file_names.".".$GCS_method.".".$clustering_cutoff_for_print.".log";}
# constatnt
my @gbk_suffix=(".gbk",".gbff","gb");

my @all_output_files=(); # will be populated during the run

if (! -d $out_dir_base) {make_path($out_dir_base) || die "Can't create the output directory '$out_dir_base' $!";}
open (my $LOG,">",$log_file) || die "Can't open LOG file '$log_file' $!";
print $LOG "============ PanGene-O-Meter ============\n";
my $datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
print $LOG "$datestring [INFO] user_cmd: $0 @orig_ARGV\n";
print $LOG "$datestring [INFO] interpreted_cmd: $0 --in_gbk_list_file $in_gbk_list --outDir $out_dir_base --ignore_private_genes $ignore_private_genes --threads $threads --pangenome_alg $PanGenomeAlg --skip_all_vs_all $skip_all_vs_all --prefix_name $Prefix_Name --GCS_method $GCS_method --CGS_clustering_cutoff $clustering_cutoff --diamond_pid $diamond_pid --diamnod_cover $diamnod_cover --reuse_PanGenome $reuse_PanGenome --in_metadata_file $in_metadata_file\n";
if (!defined $in_metadata_file or $in_metadata_file eq "" or !-e $in_metadata_file or -s $in_metadata_file==0) {print  $LOG "$datestring [INFO] Metadata file was not provided or is empty. PanGene-O-Meter Will use the data extracted from the gbk files\n";$in_metadata_file="NA";}



# 1. create an output dir for PanX and run it
# Copy the provided gbk files and change their names such they will be S1...Sn 
# store a map between user provided names and the new one
my %GenomesNames_New2Old=();
my %GenomesNames_Old2New=();

my %distances=();
my %genomes_length=();

my $Out_PhyleticPattern_prefix="";
if ($PanGenomeAlg eq "PANX") {
    if (defined $Prefix_Name and $Prefix_Name ne "") {$Out_PhyleticPattern_prefix=$out_dir_base.$Prefix_Name.".".$PanGenomeAlg_for_file_names.".PhyleticPattern";} 
    else {$Out_PhyleticPattern_prefix=$out_dir_base.$PanGenomeAlg_for_file_names.".PhyleticPattern";}
}
elsif ($PanGenomeAlg eq "DIAMONDCLUST") {
    if (defined $Prefix_Name and $Prefix_Name ne "") {$Out_PhyleticPattern_prefix=$out_dir_base.$Prefix_Name.".".$PanGenomeAlg_for_file_names.".mutual_cover_".$diamnod_cover.".id_".$diamond_pid.".PhyleticPattern";} 
    else {$Out_PhyleticPattern_prefix=$out_dir_base.$PanGenomeAlg_for_file_names.".mutual_cover_".$diamnod_cover.".id_".$diamond_pid.".PhyleticPattern";}
}

my $Out_all_vs_all_J_overlap_distance_file="";

# clean old results
if ($reuse_PanGenome eq "NO")
{
    if (-e $PanX_out_dir and -d $PanX_out_dir) {
        print "== IMPORTANT! previously created PanGenome will be DELETED, to avoid this behavior please provide --reuse_PanGenome argument\n";
        print $LOG "== IMPORTANT! previously created PanGenome will be DELETED, to avoid this behavior please provide --reuse_PanGenome argument\n";
    }
    if ($PanGenomeAlg eq "PANX")
    {
        if (-d $PanX_out_dir) {
            $datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
            print $LOG "$datestring [INFO] DELETING directory: '$PanX_out_dir'\n";
            remove_tree($PanX_out_dir);
        }
    } elsif ($PanGenomeAlg eq "DIAMONDCLUST") 
    {
        if (-d $PanX_out_dir) {
            $datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
            print $LOG "$datestring [INFO] DELETING directory: '$PanX_out_dir'\n";
            remove_tree($PanX_out_dir);
        }
        my $diamond_clusters_results_file="";
        if (defined $Prefix_Name and $Prefix_Name ne "") {$diamond_clusters_results_file=$out_dir_base.$Prefix_Name.".All_proteins.diamond_cluster.mutual_cover_".$diamnod_cover.".id_".$diamond_pid.".tsv";}
        else {$diamond_clusters_results_file=$out_dir_base."All_proteins.diamond_cluster.mutual_cover_".$diamnod_cover.".id_".$diamond_pid.".tsv";}
        if (-e $diamond_clusters_results_file)
        {
            $datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
            print "[INFO] DELETING file: '$diamond_clusters_results_file'\n";
            print $LOG "$datestring [INFO] DELETING file: '$diamond_clusters_results_file'\n";
            unlink($diamond_clusters_results_file)
        }
        my $diamond_REclusters_results_file="";
        if (defined $Prefix_Name and $Prefix_Name ne "") {$diamond_REclusters_results_file=$out_dir_base.$Prefix_Name.".All_proteins.diamond_cluster.mutual_cover_".$diamnod_cover.".id_".$diamond_pid.".recluster.tsv";}
        else {$diamond_REclusters_results_file=$out_dir_base."All_proteins.diamond_cluster.mutual_cover_".$diamnod_cover.".id_".$diamond_pid.".recluster.tsv";}
        if (-e $diamond_REclusters_results_file){
            $datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
            print "[INFO] DELETING file: '$diamond_REclusters_results_file'\n";
            print $LOG "$datestring [INFO] DELETING file: '$diamond_REclusters_results_file'\n";

            unlink($diamond_REclusters_results_file)
        }
    }
}
if (!-d $PanX_input_gbk_dir) {make_path($PanX_input_gbk_dir) or exit_on_error ("[ERROR] Can not create the output directory: '$PanX_input_gbk_dir' $!");} 

# Process input GBK list
my $genomes_name_index_file=$out_dir_base."genome_names.map.txt";
my $genome_list_for_PanX=$out_dir_base."List_of_genomes_for_PanX.txt";
open (my $LIST_OF_GENOMES_PANX,">",$genome_list_for_PanX) || exit_on_error ("[ERROR] Can't open LIST_OF_GENOMES_PANX '$genome_list_for_PanX' $!");
open (my $GENOME_NAME_MAP,">",$genomes_name_index_file) || exit_on_error ("[ERROR] Can't open GENOME_NAME_MAP '$genomes_name_index_file' $!");
open (my $GBK_LIST,"<",$in_gbk_list) || exit_on_error ("[ERROR] Can't open GBK_LIST '$in_gbk_list' $!");
my $genome_number=0;
# read the list of GBK files
my @list_of_gbk=<$GBK_LIST>;
chomp (@list_of_gbk);
close ($GBK_LIST);

# filter GBK without CDS
$datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
print $LOG "$datestring [INFO] before filtering we have ",scalar(@list_of_gbk)," gbk files\n";
my ($good_list_ArrRef,$errors_ArrRef)=filter_gbk_without_cds(\@list_of_gbk);
# my ($good_list_ArrRef,$errors_ArrRef)=validate_CDS_remove_db_xref_lines(\@list_of_gbk,$tmp_dir);
$datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
print $LOG "$datestring [INFO] after filtering we have ",scalar(@{$good_list_ArrRef})," gbk files\n";
if ((scalar @{$errors_ArrRef})>0)
{
    print "\n==== ERRORS in GBKs: =======\n";
    print join("\n",@{$errors_ArrRef}),"\n";
    print $LOG "\n==== ERRORS in GBKs: =======\n";
    print $LOG join("\n",@{$errors_ArrRef}),"\n";
}
if (scalar (@{$good_list_ArrRef})==0) {
    exit_on_error ("[FATAL_ERROR] Could not identify any coding sequences on the provided gbk files. Please make sure they include CDS tags and resubmit.\n");
} else {
    $datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
    print "$datestring [INFO] After filtering genomes that do not have enough CDS, the run continues with ".scalar(@{$good_list_ArrRef})." genomes.\n";
    print $LOG "$datestring [INFO] After filtering genomes that do not have enough CDS, the run continues with ".scalar(@{$good_list_ArrRef})." genomes.\n";
}
#if (scalar (@{$good_list_ArrRef})>0) {
#    print "\n==== GOOD: =======\n";
#    print join("\n",@{$good_list_ArrRef}),"\n";
#}
# copy the good ones and change names
# while (my $gbk_file=<$GBK_LIST>)
foreach my $gbk_file (@{$good_list_ArrRef})
{
#    chomp ($gbk_file);
    $genome_number++;
    my $new_Sname="S".$genome_number;
    my $gbk_file_base=fileparse($gbk_file, @gbk_suffix);
    my $cp_gbk_file=$PanX_input_gbk_dir.$new_Sname.".gbk";
    copy ($gbk_file, $cp_gbk_file) || exit_on_error ("[FATAL_ERROR] Can't copy '$gbk_file' to '$cp_gbk_file' $!");

    print $GENOME_NAME_MAP "$gbk_file\t$gbk_file_base\t$new_Sname\n";
    print $LIST_OF_GENOMES_PANX $new_Sname.".gbk\n";
    
    if (!exists  $GenomesNames_Old2New{$gbk_file_base}) {$GenomesNames_Old2New{$gbk_file_base}=$new_Sname;}
    else {exit_on_error ("[FATAL_ERROR] More than one genome with the same name was provided: '$gbk_file_base' is duplicate!. Please make sure all names are unique and run again...\n");}
    $GenomesNames_New2Old{$new_Sname}=$gbk_file_base;

    # get the genome length
    my ($genome_length)=get_genome_length_from_gbk($cp_gbk_file);
    # print "$gbk_file_base\t$genome_length\n";# QA
    $genomes_length{$new_Sname}=$genome_length;
}
close ($GENOME_NAME_MAP);
close ($LIST_OF_GENOMES_PANX);
    
# Finish preparing the gbk files

# if PanGenome_algorithm == PanX
if ($PanGenomeAlg eq "PANX")
{
    # run PanX
    if ($genome_number>50) {$PanXExtra_params=$PanXExtra_params."-dmdc -dcs 50"}
    my $PanX_cmd="cd $PanX_input_gbk_dir; $PANX_PATH -mo -rxm $raxml_exec -fn $PanX_out_dir -sl $genome_list_for_PanX  $PanXExtra_params -sitr -t $threads -cg 0.7 -ct -st 1 2 3 4 5 6 > ".$PanX_out_dir."PanX.std 2>&1"; ##  7 8 9
    $datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
    print $LOG "$datestring [INFO] $PanX_cmd\n";
    if (!-e "$PanX_out_dir/geneCluster") {
        system ($PanX_cmd);
    }
    # list OGs
    my @OG_aa_aln_files = grep {/_aa_aln\.fa$/} <$PanX_out_dir/geneCluster/*>; ### ToDo: test if it work for very large number of files 
    # print join ("\n",@OG_aa_aln_files),"\n";
    if (@OG_aa_aln_files<1) {exit_on_error ("[ERROR] PanX clustering failed -- no clusters were found at: '$PanX_out_dir/geneCluster/' -- look for errors at: '$PanX_out_dir"."PanX.std'\n");}
    open (my $OG_LIST,">",$out_dir_base."list_of_PanX_aa_aln_for_phylettic_pattern.txt") || exit_on_error ("[ERROR] Can't open 'OG_LIST' '$out_dir_base.list_of_PanX_aa_aln_for_phylettic_pattern.txt' $!");
    print $OG_LIST join ("\n",@OG_aa_aln_files),"\n";
    close ($OG_LIST);
    # ToDo: here we can also extract all the functions from the OGs files.
    PanX_cleanup($PanX_out_dir);
    # make PhyleticPattern
    my ($ans)=Builed_PhyleticPattern_from_PanX_aa_aln_list(\@OG_aa_aln_files,$Out_PhyleticPattern_prefix, "NA","PANX_PANGENOMETER"); # "PANX_PROKKA"
    if ($ans ne "OK") {exit_on_error ("[ERROR] 'Builed_PhyleticPattern_from_PanX_aa_aln_list(OG_aa_aln_files,$Out_PhyleticPattern_prefix, \"NA\",\"PANX_PANGENOMETER\")' --> FAILED\n");}
} elsif ($PanGenomeAlg eq "DIAMONDCLUST")
{
    # use PanX just as method to QA the genebank and extract the protein_faa files used in DIAMONDClust --> https://github.com/neherlab/pan-genome-analysis/blob/master/step-tutorials.md
    my $PanX_cmd="cd $PanX_input_gbk_dir; $PANX_PATH -mo -rxm $raxml_exec -fn $PanX_out_dir -sl $genome_list_for_PanX  $PanXExtra_params -sitr -t $threads -cg 0.7 -ct -st 1 3 4 > ".$PanX_out_dir."PanX_extract_proteins.std 2>&1"; ##  7 8 9
    $datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
    print $LOG "$datestring [INFO] $PanX_cmd\n";
    system ($PanX_cmd);

    # concate all proteins faa
    my @faa_files = grep {/\.faa$/} <$PanX_out_dir/protein_faa/*>; ### ToDo: test if it work for very large number of files 
    my $single_faa_file_for_DIAMOND="";

    if (defined $Prefix_Name and $Prefix_Name ne "") {$single_faa_file_for_DIAMOND=$out_dir_base.$Prefix_Name.".All_proteins.faa";} 
    else {$single_faa_file_for_DIAMOND=$out_dir_base."All_proteins.faa";}

    concatenate_FASTA_files($single_faa_file_for_DIAMOND,\@faa_files);
    # cluster by DIAMOND
    my $diamond_clusters_results_file="";
    if (defined $Prefix_Name and $Prefix_Name ne "") {$diamond_clusters_results_file=$out_dir_base.$Prefix_Name.".All_proteins.diamond_cluster.mutual_cover_".$diamnod_cover.".id_".$diamond_pid.".tsv";}
    else {$diamond_clusters_results_file=$out_dir_base."All_proteins.diamond_cluster.mutual_cover_".$diamnod_cover.".id_".$diamond_pid.".tsv";}

    my $diamond_cluster_cmd="$DIAMOND_CLUST_PATH cluster -d $single_faa_file_for_DIAMOND -o $diamond_clusters_results_file --approx-id $diamond_pid --mutual-cover $diamnod_cover --masking 0 --soft-masking 0 --comp-based-stats 0 -M 64G --threads $threads 1> $diamond_clusters_results_file.std 2>&1";
    $datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
    print $LOG "$datestring [INFO] $diamond_cluster_cmd\n";
    system($diamond_cluster_cmd);

    my $diamond_REclusters_results_file="";
    if (defined $Prefix_Name and $Prefix_Name ne "") {$diamond_REclusters_results_file=$out_dir_base.$Prefix_Name.".All_proteins.diamond_cluster.mutual_cover_".$diamnod_cover.".id_".$diamond_pid.".recluster.tsv";}
    else {$diamond_REclusters_results_file=$out_dir_base."All_proteins.diamond_cluster.mutual_cover_".$diamnod_cover.".id_".$diamond_pid.".recluster.tsv";}

    if (-e $diamond_clusters_results_file) # recluster
    {
        my $diamond_recluster_cmd="$DIAMOND_CLUST_PATH recluster -d $single_faa_file_for_DIAMOND --clusters  $diamond_clusters_results_file -o $diamond_REclusters_results_file --approx-id $diamond_pid --mutual-cover $diamnod_cover --masking 0 --soft-masking 0 --comp-based-stats 0 --threads $threads 1> $diamond_REclusters_results_file.std 2>&1";
        $datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
        print $LOG "$datestring [INFO] $diamond_recluster_cmd\n";
        system($diamond_recluster_cmd);
    } else {
        exit_on_error ("[ERROR] Diamond clustering failed, see '$diamond_clusters_results_file.std' for details.");
    }
    if (!-e $diamond_REclusters_results_file) {
        exit_on_error ("[ERROR] Diamond reclustering failed, see '$diamond_REclusters_results_file.std' for details.");
    }
    # recluster
#    /ebio/abt6_projects7/small_projects/hashkenazy/Programs/DIAMOND_v2.1.9/diamond  recluster -d /ebio/scratch/hashkenazy/H_pylori/DIAMOND_Clust/Analyzed200/H_pylori.Analyzed200.protein_faa.fas --clusters /ebio/scratch/hashkenazy/H_pylori/DIAMOND_Clust/Analyzed200/H_pylori.Analyzed200.protein_faa.diamond_cluster.mutual_cover_80.id_90.tsv -o /ebio/scratch/hashkenazy/H_pylori/DIAMOND_Clust/Analyzed200/H_pylori.Analyzed200.protein_faa.diamond_cluster.mutual_cover_80.id_90.reclust.tsv --mutual-cover 80 --masking 0 --soft-masking 0 --comp-based-stats 0 --approx-id 90 1> /ebio/scratch/hashkenazy/H_pylori/DIAMOND_Clust/Analyzed200/H_pylori.Analyzed200.protein_faa.diamond_cluster.mutual_cover_80.id_90.reclust.std 2>&1

# convert to Phylletic pattern
# perl /ebio/abt6_projects7/small_projects/hashkenazy/Pseudomonas/Scripts/build_phyletic_pattern_based_on_DIAMONDclusters.pl /ebio/scratch/hashkenazy/H_pylori/DIAMOND_Clust/Analyzed200/H_pylori.Analyzed200.protein_faa.diamond_cluster.mutual_cover_80.id_90.reclust.tsv /ebio/scratch/hashkenazy/H_pylori/DIAMOND_Clust/Analyzed200/H_pylori.Analyzed200.protein_faa.diamond_cluster.mutual_cover_80.id_90.reclust.PhyleticPattern "\|" 1

    # Do phylettic pattern
    my ($ans)=build_phyletic_pattern_based_on_DIAMONDclusters($diamond_REclusters_results_file,$Out_PhyleticPattern_prefix,"\\|",1);
    if ($ans ne "OK") {exit_on_error ("[ERROR] 'build_phyletic_pattern_based_on_DIAMONDclusters($diamond_REclusters_results_file,$Out_PhyleticPattern_prefix, \"\|\",\"1\")' --> FAILED\n");}

} else {
    exit_on_error ("[ERROR] possible options for PanGenome construction algorithms are: PANX or DIAMONDCLUST [value provided unknown: '$PanGenomeAlg']\n");
}

# After creating the phyletic patterns calculating the GCS scores is the same
my $PhyleticPattern_Matrix_01_file=$Out_PhyleticPattern_prefix.".01.csv";
# calculate all distances
if ($skip_all_vs_all ne "YES") 
{
    if (defined $Prefix_Name and $Prefix_Name ne "") {$Out_all_vs_all_J_overlap_distance_file=$out_dir_base.$Prefix_Name.".".$PanGenomeAlg_for_file_names.".All_vs_All.All_vs_All_GCSj_GCSo_distance";} 
    else {$Out_all_vs_all_J_overlap_distance_file=$out_dir_base.$PanGenomeAlg_for_file_names.".All_vs_All.All_vs_All_GCSj_GCSo_distance";}

    my ($ans1,$GCSj_distances_HashRef,$GCSo_distances_HashRef)=calc_jaccard_and_overlap_all_pairs ($PhyleticPattern_Matrix_01_file,$Out_all_vs_all_J_overlap_distance_file,$ignore_private_genes,",",$print_OG_diff); 
    if ($ans1 ne "OK") {exit_on_error ("[ERROR] 'calc_jaccard_and_overlap_all_pairs ($PhyleticPattern_Matrix_01_file,$Out_all_vs_all_J_overlap_distance_file,$ignore_private_genes,\",\",$print_OG_diff);'  --> FAILED\n");}
    # keep the distances for the clustering step
    if ($GCS_method eq "GCSj" and $skip_all_vs_all ne "YES") {%distances=%{$GCSj_distances_HashRef};}
    elsif ($GCS_method eq "GCSo" and $skip_all_vs_all ne "YES") {%distances=%{$GCSo_distances_HashRef};}
}
# anyway run GenContRep
# if not all distances were calculated this %distances will be empty and calculated by the function from the phylletic_pattern
my $out_clustering_prefix="";
if (defined $Prefix_Name and $Prefix_Name ne "") {$out_clustering_prefix=$out_dir_base.$Prefix_Name.".".$PanGenomeAlg_for_file_names.".".$GCS_method.".".$clustering_cutoff_for_print;}
else {$out_clustering_prefix=$out_dir_base.$PanGenomeAlg_for_file_names.".".$GCS_method.".".$clustering_cutoff_for_print;}

my ($ans)=cluster_cd_hit_like($PhyleticPattern_Matrix_01_file,$clustering_cutoff,$out_clustering_prefix,\%genomes_length,\%distances,$GCS_method); # ($phyletic_pattern_csv_file,$clustering_cutoff,$out_prefix,$genomes_length_HashRef,$GCS_method)
if ($ans ne "OK") {exit_on_error ($ans);}
# perl /ebio/abt6_projects/small_projects/hashkenazy/Pseudomonas/Scripts/eclipse-workspace/PanGenomeStat/cluster_by_Jaccard_cd_hit_like.for_PATHOCOM.pl  /ebio/abt6_projects/pacoseq/bacterial_strains_collection/Barcoded_Isolates_ONT/20250607_PATHOCOM_BARCODED_GENOMES/list_of_genomes_files.Pseudomonas /ebio/abt6_projects/pacoseq/bacterial_strains_collection/Barcoded_Isolates_ONT/Pseudomonas_viridiflava/PanX/PATHOCOM_Barcoded_Pseudomonas_viridiflava.ONT.Phyletic_Pattern_all_geneClusters_PanX.01.csv $x /ebio/abt6_projects/pacoseq/bacterial_strains_collection/Barcoded_Isolates_ONT/Pseudomonas_viridiflava/PanX/PATHOCOM_Barcoded_Pseudomonas_viridiflava.ONT.Phyletic_Pattern_all_geneClusters_PanX.Clusterd_by_Jaccard_cd_hit_like.$x 1> /ebio/abt6_projects/pacoseq/bacterial_strains_collection/Barcoded_Isolates_ONT/Pseudomonas_viridiflava/PanX/PATHOCOM_Barcoded_Pseudomonas_viridiflava.ONT.Phyletic_Pattern_all_geneClusters_PanX.Clusterd_by_Jaccard_cd_hit_like.$x.std 2>&1 
# perl /ebio/abt6_projects/small_projects/hashkenazy/Pseudomonas/Scripts/eclipse-workspace/PanGenomeStat/calc_jaccard_and_overlap_all_pairs.pl /ebio/abt6_projects/pacoseq/bacterial_strains_collection/Barcoded_Isolates_ONT/Pseudomonas_viridiflava/PanX/PATHOCOM_Barcoded_Pseudomonas_viridiflava.ONT.Phyletic_Pattern_all_geneClusters_PanX.01.csv /ebio/abt6_projects/pacoseq/bacterial_strains_collection/Barcoded_Isolates_ONT/Pseudomonas_viridiflava/PanX/PATHOCOM_Barcoded_Pseudomonas_viridiflava.ONT.Phyletic_Pattern_all_geneClusters_PanX.01.All_vs_All_Jaccard_Overlap_distance N "," "Y" 1> /ebio/abt6_projects/pacoseq/bacterial_strains_collection/Barcoded_Isolates_ONT/Pseudomonas_viridiflava/PanX/PATHOCOM_Barcoded_Pseudomonas_viridiflava.ONT.Phyletic_Pattern_all_geneClusters_PanX.01.All_vs_All_Jaccard_Overlap_distance.std 2>&1

# bring back orig names
my $clusters_annotation_file=""; # will be used later
my $distance_all_pairs_file="NA";     # if exists will be used later
foreach my $file (@all_output_files) {
    if ($file=~/.fas$/) {
        update_tax_name($file,\%GenomesNames_New2Old);
    } elsif ($file=~/PhyleticPattern.01.csv$/) {
        update_tax_name($file,\%GenomesNames_New2Old,[(0)],",","T");
    } elsif ($file=~/.PhyleticPattern.csv$/) {
        update_tax_name($file,\%GenomesNames_New2Old,[(0)],",","T");
    } elsif ($file=~/distance.all_pairs.csv$/) {
        update_tax_name($file,\%GenomesNames_New2Old,[(0,1,9)],",","T");
        $distance_all_pairs_file=$file;
    } elsif ($file=~/\d+.csv$/) { # the rep CSv phylettic pattern matrix
        update_tax_name($file,\%GenomesNames_New2Old,[(0)],",","T");
    } elsif ($file=~/.clstr.txt$/) {
        update_tax_name($file,\%GenomesNames_New2Old,[(2)],"\t","F");
    } elsif ($file=~/.clusters_sum.txt$/) {
        update_tax_name($file,\%GenomesNames_New2Old,[(0,2)],"\t","T","Y");
    } elsif ($file=~/.clusters_annotation.txt$/) {
        $clusters_annotation_file=$file;
        update_tax_name($file,\%GenomesNames_New2Old,[(0)],"\t","T");
    } elsif ($file=~/.rep_names$/) {
        update_tax_name($file,\%GenomesNames_New2Old,[(0)],"\t","F");
    }
}

my $out_metadata_file=$out_dir_base."meta_data_file.txt";
if ($in_metadata_file ne "NA") {copy ($in_metadata_file,$out_metadata_file)}
else {
    copy ($PanX_out_dir."/metainfo.tsv",$out_metadata_file);
    update_tax_name($out_metadata_file,\%GenomesNames_New2Old,[(0)],"\t","T");
}
# merge metadata file with the clusters annotation
my ($clusters_annotation_with_metadata_file)=fileparse($clusters_annotation_file,".txt");
$clusters_annotation_with_metadata_file=$out_dir_base.$clusters_annotation_with_metadata_file.".and_metadata.txt";
merge_files($clusters_annotation_file,$out_metadata_file,1,1,$clusters_annotation_with_metadata_file,"\t","\t","\t","YES");
push (@all_output_files,$clusters_annotation_with_metadata_file);

# merge metadata with the all_vs_pairs
if ($distance_all_pairs_file ne "NA")
{
    my ($dinstance_all_pairs_file_with_metadata)=fileparse($distance_all_pairs_file,".csv");
    $dinstance_all_pairs_file_with_metadata=$out_dir_base.$dinstance_all_pairs_file_with_metadata.".and_metadata.csv";
    merge_files($distance_all_pairs_file,$out_metadata_file,1,1,$dinstance_all_pairs_file_with_metadata.".tmp",",",",","\t","YES","","S1");
    merge_files($dinstance_all_pairs_file_with_metadata.".tmp",$out_metadata_file,2,1,$dinstance_all_pairs_file_with_metadata,",",",","\t","YES","","S2");
    unlink ($dinstance_all_pairs_file_with_metadata.".tmp");
    push (@all_output_files,$dinstance_all_pairs_file_with_metadata);
}

# Run finished update about all files created
$datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
print "$datestring [FINISHED] PanGene-O-Meter run has finished succesfully. The following outputs were created:\n";
print $LOG "$datestring [FINISHED] PanGene-O-Meter run has finished succesfully. The following outputs were created:\n";

my @Phyletic_pattern_files=grep (/PhyleticPattern/,@all_output_files);
print "== Phyletic Patterns\n";
print $LOG "== Phyletic Patterns\n";
foreach my $file (@Phyletic_pattern_files) {
    print "- $file\n";
    print $LOG "- $file\n";
}

my @Pairwise_files=grep (/all_pairs/,@all_output_files);
@Pairwise_files=grep(!/_distance.all_pairs.csv/,@Pairwise_files);
print "== Pairwise distances\n";
print $LOG "== Pairwise distances\n";
foreach my $file (@Pairwise_files) {
    print "- $file\n";
    print $LOG "- $file\n";
}

my @Clustering_files=sort(grep (/(representative)|(clstr)|(clusters)/,@all_output_files));
@Clustering_files=grep (!/(names_and)|(clusters_annotation\.txt)/,@Clustering_files);
print "== Clusters\n";
print $LOG "== Clusters\n";
foreach my $file (@Clustering_files) {
    print "- $file\n";
    print $LOG "- $file\n";
}
close ($LOG);


######## SUBRUTINES #########
#############################
sub filter_gbk_without_cds
{
    my ($gbk_list_ArrRef)=@_;
    my @good_gbk_list=();
    my @errors=();
    my $MINIMAL_CDS=2;
    foreach my $gbk_file (@{$gbk_list_ArrRef})
    {
        my ($open_result)=open(my $GBK,"<",$gbk_file);
        if (!$open_result) { # open_result is 1 if file opened OK and 0 if not
            push (@errors,"[ERROR] cannot open gbk file '$gbk_file' $!"); 
            next;
        } else 
        {
            my @GBK_content=<$GBK>;
            close ($GBK);
            my @CDS_lines= grep { /CDS/i } @GBK_content;
            my $number_of_CDS=scalar(@CDS_lines);
            if ($number_of_CDS>=$MINIMAL_CDS) {
                push (@good_gbk_list,$gbk_file);
            } else {
                push (@errors,"[ERROR] the gbk file '$gbk_file' has $number_of_CDS CDS which is less than the minmal number of CDS [$MINIMAL_CDS] and therefore will be further ignored!")
            }
        }
        
    }
    return (\@good_gbk_list,\@errors);
}

# sub filter_gbk_without_cds
# {
#     my $gbk_list_ArrRef=@_;
#     my @good_gbk_list=();
#     my @errors=();
#     my $MINIMAL_CDS=2;
#     foreach my $gbk_file (@{$gbk_list_ArrRef})
#     {
#         open(my $GBK,"<",$gbk_file) || {push (@errors,"[ERROR] cannot open gbk file '$gbk_file'"); next};
#         my @GBK_content=<$GBK>;
#         close ($GBK);
#         @CDS_lines= grep { /CDS/i } @GBK_content;
#         if (@CDS_lines>=$MINIMAL_CDS) {
#             push (@good_gbk_list,$gbk_file);
#         } else {
#             push (@errors,"[ERROR] the gbk file '$gbk_file' has ",scalar @CDS_lines," CDS which is less than the minmal number of CDS [$MINIMAL_CDS] and therefore will be further ignored!")
#         }
#     }
#     return (\@good_gbk_list,\@errors);
# }

sub merge_files {
    # /ebio/abt6_projects/small_projects/hashkenazy/Scripts/join_files.pl
    # if (@ARGV<5) {die "USAGE: perl $0 <File1 - base file> <File2 - To add> <File1Field> <File2Field> <Out> <new_delimiter [\\t]> <File1Delim [\\s+]> <File2Delim [\\s+]> <ALL {YES|[NO]}?>\n";}
    my ($File1,$File2,$File1Field,$File2Field,$out_file,$new_delimiter,$File1Delim,$File2Delim,$All,$File1_tag,$File2_tag)=@_;
    $File1Field--;
    $File2Field--;
    if (!defined $File1_tag) {$File1_tag="";} # name tag to add to all headers of file1
    if (!defined $File2_tag) {$File2_tag="";} # name tag to add to all headers of file2


    if ((!defined $new_delimiter) or ($new_delimiter eq ""))
    {
        $new_delimiter="\t";
    }
    if ((!defined $File1Delim) or ($File1Delim eq ""))
    {
        $File1Delim='\s+';
    }
    if ((!defined $File2Delim) or ($File2Delim eq ""))
    {
        $File2Delim='\s+';
    }
    if (defined $All) 
    {
        $All=uc($All);
    }
    else
    {
        $All="NO";
    }
    # print "$File1Field,$File2Field\n";<STDIN>;
    my %ToAdd=();


    open (my $FILE_TO_ADD,"<",$File2) || exit_on_error ("[ERROR] Can't open FILE_TO_ADD '$File2' $!");
    my $header_to_add=<$FILE_TO_ADD>;
    chomp ($header_to_add);
    my @header_to_add=split(/$File2Delim/,$header_to_add);
    if ($File2_tag ne "") {foreach my $header_elem (@header_to_add) {$header_elem=$File2_tag."_".$header_elem}}
    my $NumOfFileds2=scalar(@header_to_add);
    $header_to_add=join($new_delimiter,@header_to_add);
    while (my $line=<$FILE_TO_ADD>)
    {
        chomp ($line);
        $line=~s/^\s+|\s+$//g;
        my @line=split(/$File2Delim/,$line);
        if (defined $line[$File2Field])
        {
            $ToAdd{$line[$File2Field]}=join($new_delimiter,@line);
        }
        else
        {
            exit_on_error ("[ERROR] merge_files: can't find field '$File2Field' in line '$line' of file '$File2'\n");
        }
    }
    close ($FILE_TO_ADD);

    open (my $OUT,">",$out_file) || exit_on_error ("[ERROR] merge_files: Can't open OUT '$out_file' $!");
    open (my $FILE,"<",$File1) || exit_on_error ("[ERROR] merge_files: Can't open FILE '$File1' $!");
    my $header=<$FILE>;
    chomp ($header);
    my @header=split(/$File1Delim/,$header);
    if ($File1_tag ne "") {foreach my $header_elem (@header) {$header_elem=$File1_tag."_".$header_elem}}
    $header=join($new_delimiter,@header);
    print $OUT $header,$new_delimiter,$header_to_add,"\n";
    my @missing_array=();
    for (my $i=0;$i<$NumOfFileds2;$i++){push(@missing_array,"NA");}
    while (my $line=<$FILE>)
    {
        chomp ($line);
        $line=~s/^\s+|\s+$//g;
        my @line	=split(/$File1Delim/,$line);
        # print "@line";<STDIN>;
        # print $line[$File1Field];<STDIN>;
        if (exists $ToAdd{$line[$File1Field]})
        {
            print $OUT join ($new_delimiter,@line),$new_delimiter,$ToAdd{$line[$File1Field]},"\n";
        }
        else
        {
            if ($All eq "YES")
            {
                my $missing="$new_delimiter"."NA";
                print $OUT join ($new_delimiter,(@line,@missing_array)),"\n";
            }else {
                print "[ERROR] Could not find '$line[$File1Field]' values in $File2\n"; 
            }
        }
    }
    close ($FILE);
    close ($OUT);
}



sub build_phyletic_pattern_based_on_DIAMONDclusters {
#    use strict;
#    use warnings;
#    use File::Basename;

#     if (@ARGV<2) {die "USAGE: perl $0 <DIAMOND_Cluster> <out_Phyetic_pattern_prefix> <species_delimiter [\|]> <species_field [0]>";}
# original
# build_phyletic_pattern_based_on_DIAMONDclusters.pl /ebio/scratch/hashkenazy/H_pylori/DIAMOND_Clust/Analyzed200/H_pylori.Analyzed200.protein_faa.diamond_cluster.mutual_cover_80.id_90.reclust.tsv /ebio/scratch/hashkenazy/H_pylori/DIAMOND_Clust/Analyzed200/H_pylori.Analyzed200.protein_faa.diamond_cluster.mutual_cover_80.id_90.reclust.PhyleticPattern "\|" 1
# if (@ARGV<2) {die "USAGE: perl $0 <DIAMOND_Cluster> <out_Phyetic_pattern_prefix> <species_delimiter [NA]> <species_field [ALL]>";} # original script
 
    my ($DIAMOND_Clusters,$Out_PhyleticPattern_prefix,$spices_del,$species_field)=@_;

    if (!defined $spices_del) {$spices_del="\\|";}
    if (!defined $species_field) {$species_field=1;} # it consider it later as $species_field-1

    my %phyletic_pattern=(); # key1: species; value: array of Phyletic_pattern according to the @clusters_names; 0, for absence; 1, for presence

    my @clusters_names_orig=();
    my @clusters_names_new=();

    # my @clusters_MSAs=();
    my %NumOfSpecies_and_Copies_Per_Cluster=(); #key: cluster_name; key2: species; value: NumOfProt
    my %SpeciesNames=();

    open (my $DIAMOND_CLUSTERS,"<",$DIAMOND_Clusters) || exit_on_error ("[ERROR] build_phyletic_pattern_based_on_DIAMONDclusters: Could not open the DIAMOND_CLUSTERS '$DIAMOND_Clusters' $!");
    my $cluster_id=-1;
    my $prev_cluster="NA";
    my $cluster_id_formatted="DC_".sprintf("%06d",$cluster_id);
    while (my $line=<$DIAMOND_CLUSTERS>)
    {
        chomp ($line);
        my ($cluster_name,$cluster_mamber)=split(/\t/,$line);
        
        if ($cluster_name ne $prev_cluster){
            $cluster_id++;
            $cluster_id_formatted="DC_".sprintf("%06d",$cluster_id);
            $prev_cluster=$cluster_name;
        }
        $clusters_names_orig[$cluster_id]=$cluster_name;
        $clusters_names_new[$cluster_id]=$cluster_id_formatted;
        
    #	$clusters_MSAs[$cluster_id]=$cluster_file;
        
        my $species_name="";
        if (uc($species_field) eq "ALL")
        {
            $species_name=$cluster_mamber; # new seq name
        }
        else
        {
            my @tmp=split(/$spices_del/,$cluster_mamber);
            $species_name=$tmp[$species_field-1];
        }
        if (exists $SpeciesNames{$species_name})
        {
            $SpeciesNames{$species_name}++;
        }
        else
        {
            $SpeciesNames{$species_name}=0;
        }
        if (!exists $phyletic_pattern{$species_name}) # init
        {
            my @empty_array=();
            $phyletic_pattern{$species_name}=[@empty_array];
        }
        if (exists $phyletic_pattern{$species_name} and defined $phyletic_pattern{$species_name}->[$cluster_id])
        {
            $phyletic_pattern{$species_name}->[$cluster_id]++;
        }
        else
        {
            $phyletic_pattern{$species_name}->[$cluster_id]=1;
        }
        if (exists $NumOfSpecies_and_Copies_Per_Cluster{$cluster_id}{$species_name})
        {
            $NumOfSpecies_and_Copies_Per_Cluster{$cluster_id}{$species_name}++;
        }
        else
        {
            $NumOfSpecies_and_Copies_Per_Cluster{$cluster_id}{$species_name}=1;
        }
        # print "Finish with cluster $cluster_id\n" if ($cluster_id % 1000==0);			
        
    }
    close ($DIAMOND_CLUSTERS);


    # FROM LIST OF PanX geneClusters
    #open (my $LIST,"<",$list_of_clusters_FASTA) || die "Could not open the LIST '$list_of_clusters_FASTA' $!";
    #my $cluster_id=-1;
    #
    #while (my $cluster_file=<$LIST>)
    #{
    #	$cluster_id++;
    #	chomp ($cluster_file);
    #	my ($name,$path)=fileparse($cluster_file,("_aa_aln.fa",".fa",".fas",".aln"));
    #	$clusters_names[$cluster_id]=$name;
    #	$clusters_MSAs[$cluster_id]=$cluster_file;
    #	
    #	open (my $CLUSTER_FASTA,"<",$cluster_file) || die "Could not open CLUSTER_FASTA '$cluster_file' $!";
    #	while (my $line=<$CLUSTER_FASTA>)
    #	{
    #		chomp ($line);
    #		if ($line=~/^>(.*)/)
    #		{
    #			my $species_name="";
    #			if (uc($species_field) eq "ALL")
    #            {
    #            		$species_name=$1; # new seq name
    #            }
    #			else
    #			{
    #				my @tmp=split(/$spices_del/,$1);
    #				$species_name=$tmp[$species_field-1];
    #			}
    #			if (exists $SpeciesNames{$species_name})
    #			{
    #				$SpeciesNames{$species_name}++;
    #			}
    #			else
    #			{
    #				$SpeciesNames{$species_name}=0;
    #			}
    #			if (!exists $phyletic_pattern{$species_name}) # init
    #			{
    #				my @empty_array=();
    #				$phyletic_pattern{$species_name}=[@empty_array];
    #			}
    #			if (exists $phyletic_pattern{$species_name} and defined $phyletic_pattern{$species_name}->[$cluster_id])
    #			{
    #				$phyletic_pattern{$species_name}->[$cluster_id]++;
    #			}
    #			else
    #			{
    #				$phyletic_pattern{$species_name}->[$cluster_id]=1;
    #			}
    #			if (exists $NumOfSpecies_and_Copies_Per_Cluster{$cluster_id}{$species_name})
    #			{
    #				$NumOfSpecies_and_Copies_Per_Cluster{$cluster_id}{$species_name}++;
    #			}
    #			else
    #			{
    #				$NumOfSpecies_and_Copies_Per_Cluster{$cluster_id}{$species_name}=1;
    #			}
    #		}
    #	}
    #	close ($CLUSTER_FASTA);
    #	print "Finish with cluster $cluster_id\n" if ($cluster_id % 1000==0);
    #}
    #close ($LIST);

    my $TotalNumberOfSpecies=scalar(keys %SpeciesNames);
    my $Out_PhyleticPattern_Matrix_01=$Out_PhyleticPattern_prefix.".01.csv";
    my $Out_PhyleticPattern_Matrix=$Out_PhyleticPattern_prefix.".csv";

    open (my $OUT_MATRIX,">",$Out_PhyleticPattern_Matrix) or exit_on_error ("[ERROR] build_phyletic_pattern_based_on_DIAMONDclusters: Could not open OUT_MATRIX '$Out_PhyleticPattern_Matrix' $!");
    print $OUT_MATRIX "Species,",join(",",@clusters_names_new),"\n";

    open (my $OUT_MATRIX_01,">",$Out_PhyleticPattern_Matrix_01) or exit_on_error ("[ERROR] build_phyletic_pattern_based_on_DIAMONDclusters: Could not open OUT_MATRIX_01 '$Out_PhyleticPattern_Matrix_01' $!");
    print $OUT_MATRIX_01 "Species,",join(",",@clusters_names_new),"\n";

    my $Out_PhyleticPattern_FASTA=$Out_PhyleticPattern_prefix.".fas";
    open (my $OUT_FASTA,">",$Out_PhyleticPattern_FASTA) || exit_on_error ("[ERROR] build_phyletic_pattern_based_on_DIAMONDclusters: Could not open OUT_FASTA '$Out_PhyleticPattern_FASTA' $!");
    foreach my $species (sort keys %phyletic_pattern)
    {
        print $OUT_FASTA ">$species\n";
        print $OUT_MATRIX "$species";
        print $OUT_MATRIX_01 "$species";
        for (my $i=0; $i<=$cluster_id; $i++)
        {
            if (defined $phyletic_pattern{$species}->[$i])
            {
                print $OUT_FASTA "1";
                print $OUT_MATRIX ",".$phyletic_pattern{$species}->[$i];
                print $OUT_MATRIX_01 ",1";
            }
            else
            {
                print $OUT_FASTA "0";
                print $OUT_MATRIX ",0";
                print $OUT_MATRIX_01 ",0";
            }
        }
        print $OUT_FASTA "\n";
        print $OUT_MATRIX "\n";
        print $OUT_MATRIX_01 "\n";
    }
    close ($OUT_FASTA);
    close ($OUT_MATRIX);
    close ($OUT_MATRIX_01);

    my $Out_PhyleticPattern_Clusters_Names_Pos=$Out_PhyleticPattern_prefix.".clusters_names_and_pos.txt";
    open (my $OUT_NAMES,">",$Out_PhyleticPattern_Clusters_Names_Pos) or exit_on_error ("[ERROR] build_phyletic_pattern_based_on_DIAMONDclusters: Could not open OUT_NAMES '$Out_PhyleticPattern_Clusters_Names_Pos' $!\n");
    print $OUT_NAMES "pos\tcluster_name_DIAMOND\tcluster_name_new\n";
    for (my $i=0;$i<scalar(@clusters_names_orig);$i++)
    {
        print $OUT_NAMES $i+1,"\t$clusters_names_orig[$i]\t$clusters_names_new[$i]\n";
    }
    close ($OUT_NAMES);

    my $Out_PhyleticPattern_Clusters_Names_and_NumOfSpecies=$Out_PhyleticPattern_prefix.".clusters_names_and_NumOfSpecies.txt";
    open (my $OUT_NAMES_AND_PREVALENCE,">",$Out_PhyleticPattern_Clusters_Names_and_NumOfSpecies) or exit_on_error ("[ERROR] build_phyletic_pattern_based_on_DIAMONDclusters: Could not open OUT_NAMES_AND_PREVALENCE '$Out_PhyleticPattern_Clusters_Names_and_NumOfSpecies' $!\n");
    print $OUT_NAMES_AND_PREVALENCE "pos\tcluster_name_DIAMOND\tcluster_name_new\tPrevalence\tAvg_Copy_Number\n";
    for (my $i=0;$i<scalar(@clusters_names_orig);$i++)
    {
        my $NumOfSpecies_in_cluster=scalar(keys %{$NumOfSpecies_and_Copies_Per_Cluster{$i}});
        my $Prevalence=sprintf("%.3f",$NumOfSpecies_in_cluster/$TotalNumberOfSpecies);
        my @Copies_Per_Species=();
        my $sum=0;
        foreach my $species (keys %{$NumOfSpecies_and_Copies_Per_Cluster{$i}})
        {
            push (@Copies_Per_Species,$NumOfSpecies_and_Copies_Per_Cluster{$i}{$species});
            $sum=$sum+$NumOfSpecies_and_Copies_Per_Cluster{$i}{$species};
        }
        my $AvgCopyNumber=sprintf("%.3f",$sum/scalar(@Copies_Per_Species));
        print $OUT_NAMES_AND_PREVALENCE $i+1,"\t$clusters_names_orig[$i]\t$clusters_names_new[$i]\t$Prevalence\t$AvgCopyNumber\n";
    }
    close ($OUT_NAMES_AND_PREVALENCE);

    print "[INFO] Out_PhyleticPattern_Clusters_Names_Pos was created: '$Out_PhyleticPattern_Clusters_Names_Pos'\n";
    print "[INFO] Out_PhyleticPattern_FASTA was created: '$Out_PhyleticPattern_FASTA'\n";
    print "[INFO] Out_PhyleticPattern_Matrix_01 was created: '$Out_PhyleticPattern_Matrix_01'\n";
    print "[INFO] Out_PhyleticPattern_Matrix was created: '$Out_PhyleticPattern_Matrix'\n";
    print "[INFO] Out_PhyleticPattern_Clusters_Names_and_NumOfSpecies was created: '$Out_PhyleticPattern_Clusters_Names_and_NumOfSpecies'\n";
    push (@all_output_files,($Out_PhyleticPattern_Clusters_Names_Pos,$Out_PhyleticPattern_FASTA,$Out_PhyleticPattern_Matrix_01,$Out_PhyleticPattern_Matrix,$Out_PhyleticPattern_Clusters_Names_and_NumOfSpecies));

    return ("OK");

}

sub Builed_PhyleticPattern_from_PanX_aa_aln_list {
  #  my $OG_list_ArrRef=shift;
  #  my $Out_PhylleticPattern_file=shift;

    # perl /ebio/abt6_projects/small_projects/hashkenazy/Pseudomonas/Scripts/build_phyletic_pattern_based_on_clusters_allow_dash_in_Sname.pl /ebio/abt6_projects/pacoseq/bacterial_strains_collection/Barcoded_Isolates_ONT/Pseudomonas_viridiflava/PanX/List_of_all_PATHOCOM_Barcoded_Pseudomonas_viridiflava.ONT.PanX_geneClusters_aa_aln /ebio/abt6_projects/pacoseq/bacterial_strains_collection/Barcoded_Isolates_ONT/Pseudomonas_viridiflava/PanX/PATHOCOM_Barcoded_Pseudomonas_viridiflava.ONT.Phyletic_Pattern_all_geneClusters_PanX "NA" PANX_PROKKA
    # if (@ARGV<2) {die "USAGE: perl $0 <list_of_clusters> <out_Phyetic_pattern_prefix> <species_delimiter [NA]> <species_field [ALL|PANX_PROKKA]>";}
    my ($list_of_clusters_FASTA_ArrRef,$Out_PhyleticPattern_prefix,$spices_del,$species_field)=@_; # list_of_clusters_FASTA

    if (!defined $spices_del) {$spices_del="NA";}
    if (!defined $species_field) {$species_field="ALL";} # actually deafult in this case should be: PANX_PROKKA

    my %phyletic_pattern=(); # key1: species; value: array of Phyletic_pattern according to the @clusters_names; 0, for absence; 1, for presence

    my @clusters_names=();
    my @clusters_MSAs=();
    my %NumOfSpecies_and_Copies_Per_Cluster=(); #key: cluster_name; key2: species; value: NumOfProt
    my %SpeciesNames=();

    # open (my $LIST,"<",$list_of_clusters_FASTA) || die "Could not open the LIST '$list_of_clusters_FASTA' $!";
    my @list_of_clusters_FASTA=@{$list_of_clusters_FASTA_ArrRef};
    my $cluster_id=-1;

    # while (my $cluster_file=<$LIST>)
    foreach my $cluster_file (@list_of_clusters_FASTA)
    {
        $cluster_id++;
        chomp ($cluster_file);
        my ($name,$path)=fileparse($cluster_file,("_aa_aln.fa",".fa",".fas",".aln"));
        $clusters_names[$cluster_id]=$name;
        $clusters_MSAs[$cluster_id]=$cluster_file;
        
        open (my $CLUSTER_FASTA,"<",$cluster_file) || exit_on_error ("[ERROR] Builed_PhyleticPattern_from_PanX_aa_aln_list: Could not open CLUSTER_FASTA '$cluster_file' $!");
        while (my $line=<$CLUSTER_FASTA>)
        {
            chomp ($line);
            if ($line=~/^>(.*)/)
            {
                my $species_name="";
                if (uc($species_field) eq "ALL")
                        {
                            $species_name=$1; # new seq name
                        }
                elsif (uc($species_field) eq "PANX_PROKKA") # assumes that the header format is as follows: species_name [including -]-PROKKAgenomeID_geneNumber: DW-Xcc-32_2-JEGAEMIA_02207-1-   DW-Xcc-32_2-JEGAEMIA_02207 PANTOEA-0351-KDMFANOJ_04305 
                {
                    # FR_30B7_2-HEIPDFCH_09157-10-ribonuclease_T(2) <unknown description>
                    my $tmp=$1;
                    my $tmp_prokka_id="NA";
                    ($species_name,$tmp_prokka_id)=$tmp=~/^(.+?)\-([^-]+_+\d+)\-\d+\-/;
    #				print "QA: $line -> '$species_name'\t'$tmp_prokka_id'\n";<STDIN>;
                    # in this cpecific case it can even be simpler since we change genomes names: S\d+-
                    # ($species_name,$tmp_prokka_id)=$tmp=~/^(S\d+)\-([^-]+_+\d+)\-\d+\-/;
                }
                #elsif (uc($species_field) eq "PANX_PROKKA") # assumes that the header format is as follows: species_name [including -]-PROKKAgenomeID_geneNumber: PANTOEA-0351-KDMFANOJ_04305 
                #{
                #	my $tmp=$1;
                #	my $tmp_prokka_id="NA";
                #	($species_name,$tmp_prokka_id)=$tmp=~/^(.+?)\-([^-]+_\d+)\-/;
    #				print "QA: $line -> '$species_name'\t'$tmp_prokka_id'\n";<STDIN>;
                #}
                elsif (uc($species_field) eq "PANX_PANGENOMETER") # assumes that the header format is as follows: species_name S\d+-gene_id- [S44-DXY29_RS01245-10-DUF3143_domain-containing_protein <unknown description>]
                {
                    # FR_30B7_2-HEIPDFCH_09157-10-ribonuclease_T(2) <unknown description>
                    # S44-DXY29_RS01245-10-DUF3143_domain-containing_protein <unknown description>
                    my $tmp=$1;
                    my $tmp_prokka_id="NA";
                    # ($species_name,$tmp_prokka_id)=$tmp=~/^(.+?)\-([^-]+_+\d+)\-\d+\-/;
    #				print "QA: $line -> '$species_name'\t'$tmp_prokka_id'\n";<STDIN>;
                    # in this cpecific case it can even be simpler since we change genomes names: S\d+-
                    # ($species_name,$tmp_prokka_id)=$tmp=~/^(S\d+)\-([^_]+?_[^_-])\-\d+\-/;
                    ($species_name)=$tmp=~/^(S\d+)\-/; # we don't use the other arguments so just extract the species name which is gurrenteed to be S\d+ in PANGENEOMETER
                }
                else # assume the species id is the first token of each sequence id
                {
                    my @tmp=split(/$spices_del/,$1);
                    $species_name=$tmp[$species_field-1];
                }
                if (!defined $species_name or $species_name eq "") {print  "[ERROR] Builed_PhyleticPattern_from_PanX_aa_aln_list: Could not extract species name from header: '$line' in '$cluster_file' \n";exit (1)}
                if (exists $SpeciesNames{$species_name})
                {
                    $SpeciesNames{$species_name}++;
                }
                else
                {
                    $SpeciesNames{$species_name}=0;
                }
                if (!exists $phyletic_pattern{$species_name}) # init
                {
                    my @empty_array=();
                    $phyletic_pattern{$species_name}=[@empty_array];
                }
                if (exists $phyletic_pattern{$species_name} and defined $phyletic_pattern{$species_name}->[$cluster_id])
                {
                    $phyletic_pattern{$species_name}->[$cluster_id]++;
                }
                else
                {
                    $phyletic_pattern{$species_name}->[$cluster_id]=1;
                }
                if (exists $NumOfSpecies_and_Copies_Per_Cluster{$cluster_id}{$species_name})
                {
                    $NumOfSpecies_and_Copies_Per_Cluster{$cluster_id}{$species_name}++;
                }
                else
                {
                    $NumOfSpecies_and_Copies_Per_Cluster{$cluster_id}{$species_name}=1;
                }
            }
        }
        close ($CLUSTER_FASTA);
        # print "Finish with cluster $cluster_id\n" if ($cluster_id % 1000==0);
    }
    print "[INFO] Finished reading $cluster_id OGs out of ".scalar(@list_of_clusters_FASTA)."\n";

    # close ($LIST);

    my $TotalNumberOfSpecies=scalar(keys %SpeciesNames);
    my $Out_PhyleticPattern_Matrix_01=$Out_PhyleticPattern_prefix.".01.csv";
    my $Out_PhyleticPattern_Matrix=$Out_PhyleticPattern_prefix.".csv";

    open (my $OUT_MATRIX,">",$Out_PhyleticPattern_Matrix) or exit_on_error ("[ERROR] Builed_PhyleticPattern_from_PanX_aa_aln_list: Could not open OUT_MATRIX '$Out_PhyleticPattern_Matrix' $!");
    print $OUT_MATRIX "Species,",join(",",@clusters_names),"\n";

    open (my $OUT_MATRIX_01,">",$Out_PhyleticPattern_Matrix_01) or exit_on_error ("[ERROR] Builed_PhyleticPattern_from_PanX_aa_aln_list: Could not open OUT_MATRIX_01 '$Out_PhyleticPattern_Matrix_01' $!");
    print $OUT_MATRIX_01 "Species,",join(",",@clusters_names),"\n";

    my $Out_PhyleticPattern_FASTA=$Out_PhyleticPattern_prefix.".fas";
    open (my $OUT_FASTA,">",$Out_PhyleticPattern_FASTA) || exit_on_error ("[ERROR] Builed_PhyleticPattern_from_PanX_aa_aln_list: Could not open OUT_FASTA '$Out_PhyleticPattern_FASTA' $!");
    foreach my $species (sort keys %phyletic_pattern)
    {
        print $OUT_FASTA ">$species\n";
        print $OUT_MATRIX "$species";
        print $OUT_MATRIX_01 "$species";
        for (my $i=0; $i<=$cluster_id; $i++)
        {
            if (defined $phyletic_pattern{$species}->[$i])
            {
                print $OUT_FASTA "1";
                print $OUT_MATRIX ",".$phyletic_pattern{$species}->[$i];
                print $OUT_MATRIX_01 ",1";
            }
            else
            {
                print $OUT_FASTA "0";
                print $OUT_MATRIX ",0";
                print $OUT_MATRIX_01 ",0";
            }
        }
        print $OUT_FASTA "\n";
        print $OUT_MATRIX "\n";
        print $OUT_MATRIX_01 "\n";
    }
    close ($OUT_FASTA);
    close ($OUT_MATRIX);
    close ($OUT_MATRIX_01);

    my $Out_PhyleticPattern_Clusters_Names_Pos=$Out_PhyleticPattern_prefix.".clusters_names_and_pos.txt";
    open (my $OUT_NAMES,">",$Out_PhyleticPattern_Clusters_Names_Pos) or exit_on_error ("[ERROR] Builed_PhyleticPattern_from_PanX_aa_aln_list: Could not open OUT_NAMES '$Out_PhyleticPattern_Clusters_Names_Pos' $!\n");
    print $OUT_NAMES "pos\tcluster_name\n";
    for (my $i=0;$i<scalar(@clusters_names);$i++)
    {
        print $OUT_NAMES $i+1,"\t$clusters_names[$i]\n";
    }
    close ($OUT_NAMES);

    my $Out_PhyleticPattern_Clusters_Names_and_NumOfSpecies=$Out_PhyleticPattern_prefix.".clusters_names_and_NumOfSpecies.txt";
    open (my $OUT_NAMES_AND_PREVALENCE,">",$Out_PhyleticPattern_Clusters_Names_and_NumOfSpecies) or exit_on_error ("[ERROR] Builed_PhyleticPattern_from_PanX_aa_aln_list: Could not open OUT_NAMES_AND_PREVALENCE '$Out_PhyleticPattern_Clusters_Names_and_NumOfSpecies' $!\n");
    print $OUT_NAMES_AND_PREVALENCE "pos\tcluster_name\tcluster_MSA\tPrevalence\tAvg_Copy_Number\n";
    for (my $i=0;$i<scalar(@clusters_names);$i++)
    {
        my $NumOfSpecies_in_cluster=scalar(keys %{$NumOfSpecies_and_Copies_Per_Cluster{$i}});
        my $Prevalence=sprintf("%.3f",$NumOfSpecies_in_cluster/$TotalNumberOfSpecies);
        my @Copies_Per_Species=();
        my $sum=0;
        foreach my $species (keys %{$NumOfSpecies_and_Copies_Per_Cluster{$i}})
        {
            push (@Copies_Per_Species,$NumOfSpecies_and_Copies_Per_Cluster{$i}{$species});
            $sum=$sum+$NumOfSpecies_and_Copies_Per_Cluster{$i}{$species};
        }
        my $AvgCopyNumber=sprintf("%.3f",$sum/scalar(@Copies_Per_Species));
        print $OUT_NAMES_AND_PREVALENCE $i+1,"\t$clusters_names[$i]\t$clusters_MSAs[$i]\t$Prevalence\t$AvgCopyNumber\n";
    }
    close ($OUT_NAMES_AND_PREVALENCE);
    print "[INFO] Out_PhyleticPattern_Clusters_Names_Pos was created: '$Out_PhyleticPattern_Clusters_Names_Pos'\n";
    print "[INFO] Out_PhyleticPattern_FASTA was created: '$Out_PhyleticPattern_FASTA'\n";
    print "[INFO] Out_PhyleticPattern_Matrix_01 was created: '$Out_PhyleticPattern_Matrix_01'\n";
    print "[INFO] Out_PhyleticPattern_Matrix was created: '$Out_PhyleticPattern_Matrix'\n";
    print "[INFO] Out_PhyleticPattern_Clusters_Names_and_NumOfSpecies was created: '$Out_PhyleticPattern_Clusters_Names_and_NumOfSpecies'\n";
    push (@all_output_files,($Out_PhyleticPattern_Clusters_Names_Pos,$Out_PhyleticPattern_FASTA,$Out_PhyleticPattern_Matrix_01,$Out_PhyleticPattern_Matrix,$Out_PhyleticPattern_Clusters_Names_and_NumOfSpecies));
    return ("OK");
}


### START source: calc_jaccard_and_overlap_all_pairs.pl ###
# use strict;
# use List::Util qw(min);
sub calc_jaccard_and_overlap_all_pairs {
# if (@ARGV<2) {die "Usage: perl $0 <Full_PhyleticPattern_CSV> <OutFilePrefix> <ignore private genes [Y|N]> <delimiter [;] <print_diff [Y\[N]]>\n";}
    my ($Full_PhyleticPattern,$OutPrefix,$ignore_private_genes,$delimiter,$print_diff)=@_; # ARGV;
    # use Set::Jaccard::SimilarityCoefficient;
    if (!defined $delimiter){$delimiter=";";}
    if (!defined $ignore_private_genes) {$ignore_private_genes="N"} 
    else {$ignore_private_genes=uc($ignore_private_genes);}
    if (!defined $print_diff) {$print_diff="N"} else {$print_diff=uc($print_diff)}

    my %GCSj_distances=();
    my %GCSo_distances=();

    my $presence_cutoff=1;
    # my $OutMatrix=$OutPrefix.".csv";
    my $OutPairs=$OutPrefix.".all_pairs.csv";

    open (my $PHYLETIC_PATTERN,"<",$Full_PhyleticPattern) ||exit_on_error ("[ERROR] calc_jaccard_and_overlap_all_pairs: Can't open PHYLETIC_PATTERN for reading '$Full_PhyleticPattern' $!\n");
    my %PhyleticPattern=();
    my $PhyleticPatternHeader=<$PHYLETIC_PATTERN>;
    chomp($PhyleticPatternHeader);
    my @OG_IDs=split($delimiter,$PhyleticPatternHeader);
    shift(@OG_IDs);

    my $NumOfPos=0;
    while (my $line=<$PHYLETIC_PATTERN>)
    {
        chomp ($line);
        my ($S,@Pattern)=split($delimiter,$line);
        $PhyleticPattern{$S}=[@Pattern];
        if ($NumOfPos==0) {
            $NumOfPos=scalar(@Pattern);
        }
    }
    close ($PHYLETIC_PATTERN);

    if ($ignore_private_genes eq "Y")
    {
        my $total_pos_filtered=0;
        my $filterred_pos_file=$OutPrefix.".filtered_pos_below_".$presence_cutoff.".txt";
        open (my $FILTERED,">",$filterred_pos_file) || exit_on_error ("[ERROR] calc_jaccard_and_overlap_all_pairs: Can't open FILTERED_OUT '$filterred_pos_file' $!");
        print $FILTERED "OG_id\tOG_Pos\tSum_Presence\n";
        # loop over all pattern postion and delete from the array position below cutoff
        my %sum_per_pos=(); # key: pos in the pattern (col in the phyletic pattern), value: som of this col
        for (my $pos=0;$pos<$NumOfPos; $pos++) {
            my $sum_presence=0;
            foreach my $s (keys %PhyleticPattern) {
                if ($PhyleticPattern{$s}->[$pos]>0)
                {
                    $sum_presence++;
                }
            }
            if ($sum_presence<=$presence_cutoff) {
                foreach my $s (keys %PhyleticPattern) {
                    delete $PhyleticPattern{$s}->[$pos];
                }
                print $FILTERED join ("\t",$OG_IDs[$pos],$pos,$sum_presence),"\n";
                $total_pos_filtered++;
                # probably need to delete the pos also from the OG_ID list --> TODO: cheack!
                delete($OG_IDs[$pos]);
                # $OG_IDs[$pos]=$OG_IDs[$pos].".FILTERED";
            }
            $sum_per_pos{$pos}=$sum_presence;
        }
        close($FILTERED);	
        if ($total_pos_filtered==0) {
            unlink ($filterred_pos_file);
        } else {
            print "[INFO] calc_jaccard_and_overlap_all_pairs: Filtered $total_pos_filtered positions below or equal to the 'presence' cutoff [$presence_cutoff]\n";
        }
    }
    my @all_species=sort(keys(%PhyleticPattern));
    # my $OutPhyl=$Out.".JD_$JaccardD_cutoff.PhyleticPattern.csv";
    # open (my $OUT_PHY,">",$OutPhyl) || die "Can't open OUT_PHY '$OutPhyl' $!\n";
    # print $OUT_PHY $PhyleticPatternHeader;
    # my @PhyleticPatternHeader=split(",",$PhyleticPatternHeader);
    # my $BLANK="NA," x scalar (@PhyleticPatternHeader);
    # chop ($BLANK);

    # quickly make sure all patterns are of 0/1 and not of higher copy number order
    my $num_of_pos_left="NA";
    my $pos_changed="F";
    foreach my $s (keys %PhyleticPattern) {
        my $s_num_of_pos=0;
        for (my $pos=0;$pos<$NumOfPos;$pos++)
        {
            if (defined $PhyleticPattern{$s}->[$pos] and $PhyleticPattern{$s}->[$pos]>1) {
                print "[INFO] calc_jaccard_and_overlap_all_pairs: position $pos '$OG_IDs[$pos]' in tax $s was converted from $PhyleticPattern{$s}->[$pos] to 1\n";
                $PhyleticPattern{$s}->[$pos]=1;
                $pos_changed="T";
            }
            if (defined $PhyleticPattern{$s}->[$pos]) {
                $s_num_of_pos++;
            }
        }
        if ($num_of_pos_left eq "NA") {$num_of_pos_left=$s_num_of_pos}
        if ($num_of_pos_left != $s_num_of_pos) {print "[ERROR] calc_jaccard_and_overlap_all_pairs: The number of position in species '$s' is differnt from the rest! [$s_num_of_pos!=$num_of_pos_left]\n";}
    }
    if ($pos_changed eq "T") {print "[INFO] calc_jaccard_and_overlap_all_pairs: Notice, positions with copy number higher than 1 were treated as 1\n";}
    print "[INFO] calc_jaccard_and_overlap_all_pairs: Total positions (genes) analyzed: $num_of_pos_left\n";
    #my %PairsOfInterest=();
    my $num_of_species=scalar(@all_species);
    #print "### Num Of species= $num_of_species\n";<STDIN>;
    # open (my $OUT,">",$OutMatrix) || die "Can't open OUT '$OutMatrix' $!\n";
    open (my $OUT_PAIRS,">",$OutPairs) || exit_on_error ("[ERROR] calc_jaccard_and_overlap_all_pairs: Can't open OUT_PAIRS '$OutPairs' $!\n");

    # print $OUT "strain;",join(";",@all_species),"\n"; 
    if ($print_diff eq "Y") {
        print $OUT_PAIRS "strain1_ref,strain2,J,dJ,O_ref,Overlap_coeeficient,S1_size,S2_size,smallest_genome_size,smallest_genome_name,M01,M10,M10_diff,M01_diff\n";
    } else {
        print $OUT_PAIRS "strain1_ref,strain2,J,dJ,O_ref,Overlap_coeeficient,S1_size,S2_size,smallest_genome_size,smallest_genome_name,M01,M10\n";
    } 
    for (my $i=0;$i<$num_of_species-1;$i++)
    {
        my $S1=$all_species[$i];
        # print $OUT "$S1",";" x ($i+1);
        for (my $j=$i+1;$j<$num_of_species;$j++)
        {
            my $S2=$all_species[$j];
            # my ($dJ,$J)=calc_Jaccard_similarity_of_asymmetric_binary_attributes($PhyleticPattern{$S1},$PhyleticPattern{$S2});
            # my ($dO,$O)=calc_Overlap_similarity_of_asymmetric_binary_attributes($PhyleticPattern{$S1},$PhyleticPattern{$S2});
            
            my ($dJ,$J,$O_ref,$O,$S1_size,$S2_size,$smallest_genome_size,$smallest_genome_name,$M10,$M01,$M10_indices_ref,$M01_indices_ref)=calc_Jaccard_and_Overlap_similarity_of_asymmetric_binary_attributes($PhyleticPattern{$S1},$PhyleticPattern{$S2},$S1,$S2);
            
            $GCSj_distances{$S1}{$S2}=$J;
            $GCSj_distances{$S2}{$S1}=$J;
            $GCSo_distances{$S1}{$S2}=$O;
            $GCSo_distances{$S2}{$S1}=$O;            
            if ($print_diff eq "Y")
            {
                my @M10_diff=();
                my @M01_diff=();
                # QA
                if 	($M10 != scalar (@{$M10_indices_ref})) {print "[WARNING] calc_jaccard_and_overlap_all_pairs: $S1 and $S2 -- M10 [$M10] is not equal to the number of elements in the array [".scalar (@{$M10_indices_ref})."]\n";}		# print "QA: M10: $M10 -- ",scalar(@{$M10_indices_ref}),"\tM01: $M01 -- ",scalar(@{$M01_indices_ref}),"\n";
                if 	($M01 != scalar (@{$M01_indices_ref})) {print "[WARNING] calc_jaccard_and_overlap_all_pairs: $S1 and $S2 -- M01 [$M01] is not equal to the number of elements in the array [".scalar (@{$M01_indices_ref})."]\n";}		# print "QA: M10: $M10 -- ",scalar(@{$M10_indices_ref}),"\tM01: $M01 -- ",scalar(@{$M01_indices_ref}),"\n";

                for my $pos (@{$M10_indices_ref}) {
                    push(@M10_diff,$OG_IDs[$pos]);
                }
                # print "QA: M10: $M10 -- ",scalar(@{$M10_indices_ref}),"\tM01: $M01 -- ",scalar(@{$M01_indices_ref}),"\n";
                
                for my $pos (@{$M01_indices_ref}) {
                    push(@M01_diff,$OG_IDs[$pos]);
                }
                my $M10_diff="NA";
                my $M01_diff="NA";
                if (@M10_diff>0) {$M10_diff=join(";",@M10_diff);}
                if (@M01_diff>0) {$M01_diff=join(";",@M01_diff);}
                print $OUT_PAIRS "$S1,$S2,$J,$dJ,$O_ref,$O,$S1_size,$S2_size,$smallest_genome_size,$smallest_genome_name,$M10,$M01,$M10_diff,$M01_diff\n";

            }
            else {
                print $OUT_PAIRS "$S1,$S2,$J,$dJ,$O_ref,$O,$S1_size,$S2_size,$smallest_genome_size,$smallest_genome_name,$M10,$M01\n";
            }
            # print $OUT ";$dJ;";
        }
        # print $OUT "\n";
        # print "$i\n";
    }

    print "[INFO] All GCS distances were written to: '$OutPairs'\n";
    push (@all_output_files,($OutPairs));
    return ("OK",\%GCSj_distances,\%GCSo_distances);
    #print $OUT "S1\tS2\tPhylogeneticDistance\tJaccardDist\tJaccard\n";
    #open (my $PAIRS,"<",$ListOfPairs_File) || die "Can't open PAIRS '$ListOfPairs_File' $!\n";
    #while (my $line=<$PAIRS>)
    #{
    #	chomp ($line);
    #	my ($S1,$S2,$Dist)=split('\s+',$line);
    #	if (exists $PhyleticPattern{$S1} and exists $PhyleticPattern{$S2})
    #	{
    #		my ($dJ,$J)=calc_Jaccard_similarity_of_asymmetric_binary_attributes($PhyleticPattern{$S1},$PhyleticPattern{$S2});
    #		print $OUT "$S1\t$S2\t$Dist\t$dJ\t$J\n";
    #		if ($dJ>=$JaccardD_cutoff)
    #		{
    #			print $OUT_PHY "$S1,",join(",",@{$PhyleticPattern{$S1}}),"\n";
    #			print $OUT_PHY "$S2,",join(",",@{$PhyleticPattern{$S2}}),"\n";
    #			print $OUT_PHY $BLANK,"\n";
    #		}
    #	}
    #	else
    #	{
    #		if (!exists $PhyleticPattern{$S1})
    #		{
    #			print "[WARNNING] Could not find Phyletic Patternf for $S1, skip...\n";
    #		}
    #		if (!exists $PhyleticPattern{$S2})
    #		{
    #			print "[WARNNING] Could not find Phyletic Patternf for $S2, skip...\n";
    #		}
    #		print $OUT "$S1\t$S2\t$Dist\tNA\tNA\n";
    #	}
    ##	$PairsOfInterest{$S1}{$S2}=$Dist;
    #}
    #close ($PAIRS);






    # testing....
    # my @v1=(0,1,1,1);
    # my @v2=(1,0,0,0);
    # print "@v1\n@v2\n";
    # my ($dJ,$J)=calc_Jaccard_similarity_of_asymmetric_binary_attributes(\@v1,\@v2);
    # print "dJ: $dJ\tJ: $J\n================================\n";

    # @v1=(1,1,1,1);
    # @v2=(1,0,0,0);
    # print "@v1\n@v2\n";
    # ($dJ,$J)=calc_Jaccard_similarity_of_asymmetric_binary_attributes(\@v1,\@v2);
    # print "dJ: $dJ\tJ: $J\n================================\n";
}
sub calc_Jaccard_and_Overlap_similarity_of_asymmetric_binary_attributes
{
	# calculate Jaccard between 2 binary vectors
	# M_{11} represents the total number of attributes where A and B both have a value of 1.
	# M_{01} represents the total number of attributes where the attribute of A is 0 and the attribute of B is 1.
	# M_{10} represents the total number of attributes where the attribute of A is 1 and the attribute of B is 0.
	# M_{00} represents the total number of attributes where A and B both have a value of 0.
	
	# The Jaccard similarity coefficient, J, is given as
	#{\displaystyle J={M_{11} \over M_{01}+M_{10}+M_{11}}.} J = {M_{11} \over M_{01} + M_{10} + M_{11}}.
	
	# The Jaccard distance, dJ, is given as
	# {\displaystyle d_{J}={M_{01}+M_{10} \over M_{01}+M_{10}+M_{11}}=1-J.} d_J = {M_{01} + M_{10} \over M_{01} + M_{10} + M_{11}} = 1 - J.

	my $v1_ref=shift;
	my $v2_ref=shift;
	my $S1_name=shift;
	my $S2_name=shift;
	
	my $v1_size=scalar(@{$v1_ref});
	my $v2_size=scalar(@{$v2_ref});
	if ($v1_size!=$v2_size) { print "ERROR: calc_Jaccard_similarity_of_asymmetric_binary_attributes: v1 size is not equal to the size ov v2!!! ($v1_size!=$v2_size)\n";return "NA";}
	# calculate
	my $M11=0;
	my $M00=0;
	my $M10=0;
	my $M01=0;
	my @M10_indices=();
	my @M01_indices=();

	for (my $i=0;$i<$v1_size;$i++)
	{
		if (defined $v1_ref->[$i] and defined $v2_ref->[$i]) # not filtered out
		{
			if ($v1_ref->[$i]==1 and $v2_ref->[$i]==1) {$M11++;}
			if ($v1_ref->[$i]==0 and $v2_ref->[$i]==0) {$M00++;}
			if ($v1_ref->[$i]==1 and $v2_ref->[$i]==0) {$M10++;push(@M10_indices,$i);}
			if ($v1_ref->[$i]==0 and $v2_ref->[$i]==1) {$M01++;push(@M01_indices,$i);}
		}
	}
	my ($dJ,$J,$O_ref,$O);

    my $S1_size=($M10+$M11);
    my $S2_size=($M01+$M11);
    my $smallest_genome_size=min($S1_size,$S2_size);
    my $smallest_genome_name=$S2_name;
    if ($S2_size > $S1_size) {$smallest_genome_name=$S1_name}

	if (($M01+$M10+$M11)==0)
	{
		$dJ="NA";
		$J="NA";
		# return ("NA","NA");
	}
	else
	{
		$dJ=sprintf("%.3f",(($M01+$M10)/($M01+$M10+$M11)));
		$J=sprintf("%.3f",(($M11)/($M01+$M10+$M11)));
	}
	
	 if (($M10+$M11)==0)
     {
         $O_ref="NA";
         $O="NA";
     }
     else
     {
		if ($S1_size==0 or $smallest_genome_size==0)
		{
			$O_ref="NA";
			$O="NA";
		} else {
	         $O_ref=sprintf("%.5f",(($M11)/$S1_size)); # overalp of the elements in ref (S1) 
    	     $O=sprintf("%.5f",(($M11)/$smallest_genome_size)); # overlap coefficient: https://en.wikipedia.org/wiki/Overlap_coefficient
		}
     }
	return ($dJ,$J,$O_ref,$O,$S1_size,$S2_size,$smallest_genome_size,$smallest_genome_name,$M10,$M01,\@M10_indices,\@M01_indices);	
}


sub calc_Jaccard_similarity_of_asymmetric_binary_attributes
{
	# calculate Jaccard between 2 binary vectors
	# M_{11} represents the total number of attributes where A and B both have a value of 1.
	# M_{01} represents the total number of attributes where the attribute of A is 0 and the attribute of B is 1.
	# M_{10} represents the total number of attributes where the attribute of A is 1 and the attribute of B is 0.
	# M_{00} represents the total number of attributes where A and B both have a value of 0.
	
	# The Jaccard similarity coefficient, J, is given as
	#{\displaystyle J={M_{11} \over M_{01}+M_{10}+M_{11}}.} J = {M_{11} \over M_{01} + M_{10} + M_{11}}.
	
	# The Jaccard distance, dJ, is given as
	# {\displaystyle d_{J}={M_{01}+M_{10} \over M_{01}+M_{10}+M_{11}}=1-J.} d_J = {M_{01} + M_{10} \over M_{01} + M_{10} + M_{11}} = 1 - J.

	my $v1_ref=shift;
	my $v2_ref=shift;
	
	my $v1_size=scalar(@{$v1_ref});
	my $v2_size=scalar(@{$v2_ref});
	if ($v1_size!=$v2_size) { print "ERROR: calc_Jaccard_similarity_of_asymmetric_binary_attributes: v1 size is not equal to the size ov v2!!! ($v1_size!=$v2_size)\n";return "NA";}
	# calculate
	my $M11=0;
	my $M00=0;
	my $M10=0;
	my $M01=0;
	
	for (my $i=0;$i<$v1_size;$i++)
	{
		if (defined $v1_ref->[$i] and defined $v2_ref->[$i]) # not filtered out
		{
			if ($v1_ref->[$i]==1 and $v2_ref->[$i]==1) {$M11++;}
			if ($v1_ref->[$i]==0 and $v2_ref->[$i]==0) {$M00++;}
			if ($v1_ref->[$i]==1 and $v2_ref->[$i]==0) {$M10++;}
			if ($v1_ref->[$i]==0 and $v2_ref->[$i]==1) {$M01++;}
		}
	}
	if (($M01+$M10+$M11)==0)
	{
		return ("NA","NA");
	}
	else
	{
		my $dJ=sprintf("%.3f",(($M01+$M10)/($M01+$M10+$M11)));
		my $J=sprintf("%.3f",(($M11)/($M01+$M10+$M11)));
		return ($dJ,$J);	
	}
}

sub calc_Overlap_similarity_of_asymmetric_binary_attributes
{
        # calculate Jaccard between 2 binary vectors
        # M_{11} represents the total number of attributes where A and B both have a value of 1.
        # M_{01} represents the total number of attributes where the attribute of A is 0 and the attribute of B is 1.
        # M_{10} represents the total number of attributes where the attribute of A is 1 and the attribute of B is 0.
        # M_{00} represents the total number of attributes where A and B both have a value of 0.
        
        # The Jaccard similarity coefficient, J, is given as
        #{\displaystyle J={M_{11} \over M_{01}+M_{10}+M_{11}}.} J = {M_{11} \over M_{01} + M_{10} + M_{11}}.
        
        # The Jaccard distance, dJ, is given as
        # {\displaystyle d_{J}={M_{01}+M_{10} \over M_{01}+M_{10}+M_{11}}=1-J.} d_J = {M_{01} + M_{10} \over M_{01} + M_{10} + M_{11}} = 1 - J.

        my $v1_ref=shift;
        my $v2_ref=shift;
        
        my $v1_size=scalar(@{$v1_ref});
        my $v2_size=scalar(@{$v2_ref});
        if ($v1_size!=$v2_size) { print "ERROR: calc_Jaccard_similarity_of_asymmetric_binary_attributes: v1 size is not equal to the size ov v2!!! ($v1_size!=$v2_size)\n";return "NA";}
        # calculate
        my $M11=0;
        my $M00=0;
        my $M10=0;
        my $M01=0;
        
        for (my $i=0;$i<$v1_size;$i++)
        {
        	if (defined $v1_ref->[$i] and defined $v2_ref->[$i]) # not filtered out
			{
                if ($v1_ref->[$i]==1 and $v2_ref->[$i]==1) {$M11++;}
                if ($v1_ref->[$i]==0 and $v2_ref->[$i]==0) {$M00++;}
                if ($v1_ref->[$i]==1 and $v2_ref->[$i]==0) {$M10++;}
                if ($v1_ref->[$i]==0 and $v2_ref->[$i]==1) {$M01++;}
			}
        }
        if (($M01+$M10+$M11)==0)
        {
                return ("NA","NA");
        }
        else
        {
                # my $dJ=sprintf("%.5f",(($M01+$M10)/($M01+$M10+$M11)));
                # my $J=sprintf("%.5f",(($M11)/($M01+$M10+$M11)));
                
                my $dJ=sprintf("%.5f",(($M01)/($M10+$M11)));
                my $J=sprintf("%.5f",(($M11)/($M10+$M11)));
                return ($dJ,$J);        
        }
}
### END source: calc_jaccard_and_overlap_all_pairs.pl ###

### START: source: cluster_by_Jaccard_cd_hit_like.for_PATHOCOM.pl  ###
### TO DO -- correctly handle >1 in the PA
# if (@ARGV<3) {die "USAGE: perl $0  <list_of_genomes_fas> <phyletic_pattern_csv_file> <clustering_cutoff [0-1]> <out_prefix>\n";}
# ToDo: cleanup --> get_genomes_length_FASTA is no longer used in PanGenOmeter
sub get_genomes_length_FASTA {
    my ($genomes_files_ArrRef)=@_;
#   open (my $GENOMES_LIST,"<",$list_of_genomes_fas) || die "Can't open GENOMES_LIST '$list_of_genomes_fas' $!";
#    my @genomes_files=<$GENOMES_LIST>;
#    chomp (@genomes_files);
#    close ($GENOMES_LIST);

    my @genomes_files=@{$genomes_files_ArrRef};
    my %genomes_length=();
    # print $LOG "genome_file\tlength_bp\n";
    foreach my $genome (@genomes_files)
    {
        my $length=get_total_genome_length_fasta($genome);
        $genomes_length{$genome}=$length;
#        print "genomes_length{$genome}=$length\n";
#        print $LOG "$genome\t$length\n";
    }
    return (\%genomes_length);
}
sub cluster_cd_hit_like {
    my ($phyletic_pattern_csv_file,$clustering_cutoff,$out_prefix,$genomes_length_HashRef,$distances_HashRef,$GCS_method)=@_;# ARGV;
    my $rep_fasta_file=$out_prefix.".representative_PhP.fas";
    my $rep_CSV_file=$out_prefix.".representative_PhP.csv";
    my $clusters_file=$out_prefix.".clstr.txt";
    my $log_file=$out_prefix.".log.txt";
    my $rep_names_file=$out_prefix.".representative_names";

    my %genomes_length=%{$genomes_length_HashRef};
    my %distances=%{$distances_HashRef};
    open (my $LOG1,">",$log_file) ||  exit_on_error ("[ERROR] cluster_cd_hit_like: Can't open LOG1 '$log_file' $!");
    # print $LOG1 "# running command: perl $0 @ARGV\n\n";
    print $LOG1 "# running cluster_cd_hit_like: $phyletic_pattern_csv_file $clustering_cutoff $out_prefix genomes_length_HashRef distances_HashRef $GCS_method\n\n";
    #my @sorted_list_of_genomes=sort by_total_genome_length_fasta_rev (@genomes_files);
    my @sorted_list_of_genomes=sort {$genomes_length{$b}<=>$genomes_length{$a}} (keys %genomes_length);
    # print join("\n",@sorted_list_of_genomes),"\n"; # QA
    print $LOG1 "Sorted list of genomes by their nt length\n====================================================================================\n",join("\n",@sorted_list_of_genomes),"\n";
    print $LOG1 "\n";
    open (my $PHYLETIC_PATTERN,"<",$phyletic_pattern_csv_file) || exit_on_error ("[ERROR] cluster_cd_hit_like: Can't open PHYLETIC_PATTERN '$phyletic_pattern_csv_file' $!\n");
    my %PhyleticPattern=();
    my $PhyleticPatternHeader=<$PHYLETIC_PATTERN>;
    while (my $line=<$PHYLETIC_PATTERN>)
    {
            chomp ($line);
            my ($S,@Pattern)=split(",",$line);
            $PhyleticPattern{$S}=[@Pattern];
            # print "Reading Phyletic: $S\n"; # QA
            print $LOG1 "Reading Phyletic: $S\n";
    }
    close ($PHYLETIC_PATTERN);

    my %Clusters=();
    my $cluster_number=0;
    open (my $FASTA_REP,">",$rep_fasta_file) || exit_on_error ("[ERROR] cluster_cd_hit_like: Can't open FASTA_REP '$rep_fasta_file' $!");
    open (my $CLUSTERS,">",$clusters_file) || exit_on_error ("[ERROR] cluster_cd_hit_like: Can't open CLUSTERS '$clusters_file' $!");
    open (my $CSV_REP,">",$rep_CSV_file) || exit_on_error ("[ERROR] cluster_cd_hit_like: Can't open CSV '$rep_CSV_file' $!");
    open (my $REP_NAMES,">",$rep_names_file) || exit_on_error ("[ERROR] cluster_cd_hit_like: Can't open REP_NAMES '$rep_names_file' $!");
    print $CSV_REP $PhyleticPatternHeader;
    my $genomes_left_for_clustering=scalar @sorted_list_of_genomes;
    print "[INFO] Start clustering $genomes_left_for_clustering genomes by $GCS_method [cutoff=$clustering_cutoff]\n";
    print $LOG1 "# [INFO] Start clustering $genomes_left_for_clustering genomes by $GCS_method [cutoff=$clustering_cutoff]\n";

    while (@sorted_list_of_genomes>0 and $genomes_left_for_clustering>0)
    {
        my $cluster_size=0;
        my $cluster_head=shift @sorted_list_of_genomes;
        if (!defined $cluster_head) {next;} # print "*'$cluster_head'\n";
        $genomes_left_for_clustering--;
        $cluster_number++;
        my $short_name=fileparse($cluster_head,(".gbk.fna.gz",".fasta.gz",".fna.gz",".genome.fasta",".short.fa",".short.fna",".fasta",".fna",".fa"));
        # print "QA0 - short: $short_name\n";
        if (!exists $PhyleticPattern{$short_name} and $short_name=~/(.*?)_genomic$/)
        {
                $short_name=$1;
        }
        if (!exists $PhyleticPattern{$short_name} and $short_name=~/(.*?).ragtag.scaffolds.Ori.shortH$/)
        {
                $short_name=$1;
        }
        if (!exists $PhyleticPattern{$short_name} and $short_name=~/(.*?)\.SPAdes/i)
        {
            $short_name=$1;
        }
        if (!exists $PhyleticPattern{$short_name} and $short_name=~/(.*?).flye.clean.oriC$/)
        {
            $short_name=$1;
        }
    #        if (!exists $PhyleticPattern{$short_name} and $short_name=~/-/)
    #        {
    #                $short_name=~s/-/_/g;
    #        }
        if (!exists $PhyleticPattern{$short_name} and $short_name=~/(_ASM.+)$/)
        {
            my ($short_name1)=$short_name=~/(.+?)(_ASM.+)$/;
            $short_name=$short_name1;
        }
        if (! exists $PhyleticPattern{$short_name} and $short_name=~/\-/)
        {
            my ($short_name1)=split(/-/,$short_name); # take only the first token to solve some cases caused by PanX phylettic pattern parsing...
            print "[QA] $short_name -> $short_name1\n";
            if (exists $PhyleticPattern{$short_name1}) {
                print "[INFO] $short_name is treated as $short_name1 forwhich we have a phylettic pattern\n";
                $short_name=$short_name1;
            }
        }
        if (! exists $PhyleticPattern{$short_name})
            {
            print "[WARNING] can't find the phyletic pattern of cluster_head '$short_name' ignoring it!\n";
            print $LOG1 "[WARNING] can't find the phyletic pattern of cluster_head '$short_name' ignoring it!\n";
            next;
        }
        # print "== CLUSTER $cluster_number: $short_name\n";
        print $LOG1 "== CLUSTER $cluster_number: $short_name\n";
        
        print $CLUSTERS ">Cluster $cluster_number\n";
        print $CLUSTERS "$cluster_size\t$genomes_length{$cluster_head}\t$short_name\t*\n";
        print $FASTA_REP ">$short_name\n",join("",@{$PhyleticPattern{$short_name}}),"\n";
        print $CSV_REP "$short_name,",join(",",@{$PhyleticPattern{$short_name}}),"\n";
        print $REP_NAMES "$short_name\n";
        # print "TO CLUSTER: ",count_real_elements_in_array (\@sorted_list_of_genomes) ." or $genomes_left_for_clustering\n"; # QA
        for (my $candidate_index=0;$candidate_index<scalar (@sorted_list_of_genomes);$candidate_index++)
        {
            my $candidate_full_name=$sorted_list_of_genomes[$candidate_index];
            if (defined $candidate_full_name)
            {
                my $candidate=fileparse($candidate_full_name,(".gbk.fna.gz",".fasta.gz",".fna.gz",".genome.fasta",".short.fa",".short.fna",".fasta",".fna",".fa"));
                if (!exists $PhyleticPattern{$candidate} and $candidate=~/(.*?)_genomic$/)
                {
                    $candidate=$1;
                }
                if (!exists $PhyleticPattern{$candidate} and $candidate=~/(.*?).ragtag.scaffolds.Ori.shortH$/)
                {
                    $candidate=$1;
                }
    #			if (!exists $PhyleticPattern{$candidate} and $candidate=~/-/)
    #                        {
    #                                $candidate=~s/-/_/g;
    #                        }
                if (!exists $PhyleticPattern{$candidate} and $candidate=~/(_ASM.+)$/)
                {
                    my ($candidate1)=$candidate=~/(.+?)(_ASM.+)$/;
                    $candidate=$candidate1;
                }
                if (!exists $PhyleticPattern{$candidate} and $candidate=~/(.*?)\.SPAdes/i)
                {
                    $candidate=$1;
                }
                if (!exists $PhyleticPattern{$candidate} and $candidate=~/(.*?).flye.clean.oriC$/)
                {
                    $candidate=$1;
                }
                if (! exists $PhyleticPattern{$candidate} and $candidate=~/\-/)
                {
                    my ($candidate1)=split(/-/,$candidate); # take only the first token to solve some cases caused by PanX phylettic pattern parsing...
                    print "[QA] '$candidate'->'$candidate1'\n";
                    $candidate=$candidate1;
                }
                if (! exists $PhyleticPattern{$candidate})
                {
                    print "[WARNING] can't find the phyletic pattern of candidate '$candidate' ignoring it!\n";
                    print $LOG1 "[WARNING] can't find the phyletic pattern of '$candidate' ignoring it!\n";
                    # $sorted_list_of_genomes[$candidate_index]=undef;
                    delete ($sorted_list_of_genomes[$candidate_index]);
                    next;
                }
                my $S2=$PhyleticPattern{$candidate};
                
                # check if distance was already calculated
                my $distance="NA";
                if (exists $distances{$short_name}{$candidate})
                {
                    $distance=$distances{$short_name}{$candidate};
                } elsif (exists $distances{$candidate}{$short_name}) {
                    $distance=$distances{$candidate}{$short_name}
                } else {
                    my ($dJ,$J,$O_ref,$O,$S1_size,$S2_size,$smallest_genome_size,$smallest_genome_name,$M10,$M01,$M10_indices_ref,$M01_indices_ref)=calc_Jaccard_and_Overlap_similarity_of_asymmetric_binary_attributes($PhyleticPattern{$short_name},$PhyleticPattern{$candidate},$short_name,$candidate);
                    # my $jaccard=calc_Jaccard_similarity_of_asymmetric_binary_attributes($PhyleticPattern{$short_name},$PhyleticPattern{$candidate});
                    if ($GCS_method eq "GCSj") {$distance=$J;}
                    else {$distance=$O;}
                }
                if ($distance eq "NA") {return ("[ERROR] Could not calculate the distance between '$short_name' and '$candidate'!");}
                # print "[QA] $short_name\t$candidate\t$jaccard\n";
                # if ($jaccard>=$clustering_cutoff)
                if ($distance>=$clustering_cutoff)
                {
                    # $jaccard=$jaccard*100;
                    $distance=$distance*100;
                    $cluster_size++;
                    # print "- Adding $candidate ($distance)\n";
                    print $LOG1 "- Adding $candidate ($distance)\n";
                    $genomes_left_for_clustering--;
                    print $CLUSTERS "$cluster_size\t$genomes_length{$candidate_full_name}\t$candidate\tat $distance%\n";
                    delete ($sorted_list_of_genomes[$candidate_index]);
                    if (defined $Clusters{$short_name})
                    {
                        push (@{$Clusters{$short_name}},$candidate);
                    }
                    else
                    {
                        $Clusters{$short_name}=[($candidate)];
                    }
                }
            }
        }
        if (!exists $Clusters{$short_name}) # nothing was added, so we have an empty list just with the representative sequence...
        {
            $Clusters{$short_name}=[()]; 
        }
        
        # print "[INFO] $genomes_left_for_clustering genomes left for clustering after finishing with cluster #$cluster_number [".join("*",@sorted_list_of_genomes),"]\n";
        print "[INFO] $genomes_left_for_clustering genomes left for clustering after finishing with cluster #$cluster_number\n";
    }

    my $Final_Clusters_tsv=$out_prefix.".clusters_sum.txt";
    my $Final_Clusters_annotation=$out_prefix.".clusters_annotation.txt";
    open (my $FINAL_CLUSTERS_TSV,">",$Final_Clusters_tsv) || exit_on_error ("[ERROR] cluster_cd_hit_like: Can't open FINAL_CLUSTERS_TSV '$Final_Clusters_tsv' $!");
    open (my $FINAL_CLUSTERS_ANNOTATION,">",$Final_Clusters_annotation) || exit_on_error ("[ERROR] cluster_cd_hit_like: Can't open FINAL_CLUSTERS_ANNOTATION '$Final_Clusters_annotation' $!");
    print $LOG1 "List of clusters (by $GCS_method $clustering_cutoff)\n====================================================================================\ncluster_representative\tcluster_size\n";
    print $FINAL_CLUSTERS_TSV "rep\tcluster_size\tmambers\n";
    print $FINAL_CLUSTERS_ANNOTATION "name\tcluster_id\trep\n";
    my $cluster_id=0;
    foreach my $cluster (sort {scalar @{$Clusters{$b}}<=> scalar @{$Clusters{$a}}} keys %Clusters)
    {
        $cluster_id++;
        # print "$cluster\t",scalar @{$Clusters{$cluster}},"\t",join(";",@{$Clusters{$cluster}}),"\n";
        print $LOG1 "$cluster\t",scalar @{$Clusters{$cluster}},"\t",join(";",@{$Clusters{$cluster}}),"\n";
        print $FINAL_CLUSTERS_TSV "$cluster\t",scalar @{$Clusters{$cluster}}+1,"\t",join(";",sort(($cluster,@{$Clusters{$cluster}}))),"\n";
        print $FINAL_CLUSTERS_ANNOTATION "$cluster\t$cluster_id\t1\n";
        foreach my $member (@{$Clusters{$cluster}})
        {
            print $FINAL_CLUSTERS_ANNOTATION "$member\t$cluster_id\t0\n";
        }
    }
    close ($FINAL_CLUSTERS_TSV);
    close ($CLUSTERS);
    close ($FINAL_CLUSTERS_ANNOTATION);
    close ($FASTA_REP);
    close ($REP_NAMES);

    print $LOG1 "[INFO] Clusters file was created: '$clusters_file'\n";
    print $LOG1 "[INFO] Clusters_tsv file was created: '$Final_Clusters_tsv'\n";
    print $LOG1 "[INFO] Clusters_annotation was created: '$Final_Clusters_annotation'\n";
    print $LOG1 "[INFO] Representatives Phyletic Patterns FASTA was created: '$rep_fasta_file'\n";
    print $LOG1 "[INFO] Representatives Phyletic Patterns CSV was created: '$rep_CSV_file'\n";
    print $LOG1 "[INFO] Representatives names were written to: '$rep_names_file'\n";
    close ($LOG1);

    push (@all_output_files,($Final_Clusters_tsv,$Final_Clusters_annotation,$rep_fasta_file,$clusters_file,$rep_CSV_file,$rep_names_file,));
    return ("OK")
}

# ToDo: cleanup --> by_total_genome_length_fasta_rev is no longer used in PanGenOmeter
sub by_total_genome_length_fasta_rev
{
	my $length_a=get_total_genome_length_fasta($a);
	my $length_b=get_total_genome_length_fasta($b);
	if ($length_a>$length_b){return 1;}
	elsif ($length_a<$length_b){return -1;}
	else {return 0;}
}

# ToDo: cleanup --> get_total_genome_length_fasta is no longer used in PanGenOmeter
#  sub get_total_genome_length_fasta 
#  {
# 	my $fasta_file=shift;
# 	my $IN;
#     if ($fasta_file =~ /.gz$/) 
# 	{
#         $IN = new IO::Zlib;
#         $IN->open($fasta_file, "rb") or exit_on_error ("[ERROR] get_total_genome_length_fasta: can't open IN '$fasta_file' $!\n");     
#         # open(IN, "gunzip -c $fasta_file |") || exit_on_error ("[ERROR] get_total_genome_length_fasta: Can’t open pipe to $fasta_file");
# 	}
# 	else {
# 		open($IN, $fasta_file) or  exit_on_error ("[ERROR] get_total_genome_length_fasta: can't open '$fasta_file': $!\n");
# 	}
	
# 	my $length=0;
# 	while (my $line=<IN>)
# 	{
# 		if ($line!~/^>/)
# 		{
# 			chomp $line;
# 			$length=$length+length($line);
# 		}
# 	}
# 	close (IN);
# 	return $length;
# }

#     OLDER implenetion, does not support filtered positions     #
# sub calc_Jaccard_similarity_of_asymmetric_binary_attributes
# {
#         # calculate Jaccard between 2 binary vectors
#         # M_{11} represents the total number of attributes where A and B both have a value of 1.
#         # M_{01} represents the total number of attributes where the attribute of A is 0 and the attribute of B is 1.
#         # M_{10} represents the total number of attributes where the attribute of A is 1 and the attribute of B is 0.
#         # M_{00} represents the total number of attributes where A and B both have a value of 0.
        
#         # The Jaccard similarity coefficient, J, is given as
#         #{\displaystyle J={M_{11} \over M_{01}+M_{10}+M_{11}}.} J = {M_{11} \over M_{01} + M_{10} + M_{11}}.
        
#         # The Jaccard distance, dJ, is given as
#         # {\displaystyle d_{J}={M_{01}+M_{10} \over M_{01}+M_{10}+M_{11}}=1-J.} d_J = {M_{01} + M_{10} \over M_{01} + M_{10} + M_{11}} = 1 - J.

#         my $v1_ref=shift;
#         my $v2_ref=shift;
        
#         my $v1_size=scalar(@{$v1_ref});
#         my $v2_size=scalar(@{$v2_ref});
#         if ($v1_size!=$v2_size) { print "ERROR: calc_Jaccard_similarity_of_asymmetric_binary_attributes: v1 size is not equal to the size ov v2!!! ($v1_size!=$v2_size)\n";return "NA";}
#         # calculate
#         my $M11=0;
#         my $M00=0;
#         my $M10=0;
#         my $M01=0;
        
#         for (my $i=0;$i<$v1_size;$i++)
#         {
#                 if ($v1_ref->[$i]==1 and $v2_ref->[$i]==1) {$M11++;}
#                 if ($v1_ref->[$i]==0 and $v2_ref->[$i]==0) {$M00++;}
#                 if ($v1_ref->[$i]==1 and $v2_ref->[$i]==0) {$M10++;}
#                 if ($v1_ref->[$i]==0 and $v2_ref->[$i]==1) {$M01++;}
#         }
#         if (($M01+$M10+$M11)==0)
#         {
#                 return ("NA","NA");
#         }
#         else
#         {
#                 my $dJ=sprintf("%.5f",(($M01+$M10)/($M01+$M10+$M11)));
#                 my $J=sprintf("%.5f",(($M11)/($M01+$M10+$M11)));
#                 return ($dJ,$J);        
#         }
# }
sub count_real_elements_in_array
{
	my $arr_ref=shift;
	my @array=@{$arr_ref};
	my $counter=0;
	for (my $i=0;$i<scalar(@array);$i++)
	{
		if (defined $array[$i])
		{
			$counter++;
		}
	}
	return $counter;
}
### END: source: ###

sub get_genome_length_from_gbk
{
    my ($gbk_file)=@_;
    my $genome_length=0;
    open (my $GBK,"<",$gbk_file) || exit_on_error ("[ERROR] get_genome_length_from_gbk: Failed to open GBK file '$gbk_file' $!");
    while (my $line=<$GBK>)
    {
        if ($line=~/LOCUS\s+\S+\s+(\d+) bp/)
        {
                my $contig_length=$1;
                $genome_length=$genome_length+$contig_length;
        }
    }
    close ($GBK);
    return ($genome_length);
}

sub update_tax_name {
    my ($file,$new_to_orig_hashRef,$tax_col_ArrRef,$sep,$header,$is_list_elemnt)=@_;
    if (defined $header) {$header=uc($header);} else {$header="F";}
    open (my $FILE,"<",$file) || exit_on_error ("[ERROR] update_tax_name: Can't open FILE '$file' for reading: $!");
    my @file_content=<$FILE>;
    close ($FILE);

    open (my $OUT_NEW,">",$file) || exit_on_error ("[ERROR] update_tax_name: Can't open OUT_NEW '$file' for writing: $!");
    if ($header eq "T")
    {
        my $line=shift (@file_content);
        print $OUT_NEW $line;
    }
    foreach my $line (@file_content) {
        if ($line=~/>(.*)\n/ and (!defined $tax_col_ArrRef)) # handle true FASTA file
        {
            my $old_name=$1;
            if (exists $new_to_orig_hashRef->{$old_name})
            {
                print $OUT_NEW ">".$new_to_orig_hashRef->{$old_name}."\n";
            }
            else {
                print "[WARN] Could not find the original name of tax '$old_name', can't change file: '$file' --> look at '$genomes_name_index_file' for details...\n";
                print $OUT_NEW $line;
            }
        } elsif (defined $tax_col_ArrRef) { # handle delimited files
            if ($line=~/^>/) {print $OUT_NEW $line;} # to handle the .clstr.txt file correctly
            else {
                chomp ($line);
                my @line=split(/$sep/,$line);
                foreach my $col (@{$tax_col_ArrRef})
                {
                    if (defined $is_list_elemnt and $is_list_elemnt eq "Y")
                    { # the element itself could be a list of ; elemnts
                        my @sub_elements=split(/;/,$line[$col]);
                        my @new_sub_elements=();
                        foreach my $old_name (@sub_elements)
                        {
                            if (exists $new_to_orig_hashRef->{$old_name})
                            {       
                                push (@new_sub_elements,$new_to_orig_hashRef->{$old_name});
                            } else {
                                print "[WARN] Could not find the original name of tax '$old_name', can't change file: '$file' --> look at '$genomes_name_index_file' for details...\n";
                            }
                        }
                        $line[$col]=join(";",@new_sub_elements);
                    } else {
                        my $old_name=$line[$col];
                        if (exists $new_to_orig_hashRef->{$old_name})
                        {
                            $line[$col]=$new_to_orig_hashRef->{$old_name};
                        } else {
                             print "[WARN] Could not find the original name of tax '$old_name', can't change file: '$file' --> look at '$genomes_name_index_file' for details...\n";
                        }
                    }

                }
                print $OUT_NEW join ($sep,@line),"\n";
            }
        } else { # don't change file --> not sure this is needed
            print $OUT_NEW $line;
        }
    }
    close ($OUT_NEW);
    # print "[INFO] Final - tax names update in file: $file\n";
}
sub  concatenate_FASTA_files {
    my ($out_single_fasta_file,$fasta_files_ArrRef)=@_;
    open (my $OUT,">",$out_single_fasta_file) || exit_on_error ("[ERROR] concatenate_FASTA_files can't open OUT '$out_single_fasta_file' for writing $!");
    foreach my $in_file (@{$fasta_files_ArrRef}) {
        open (my $IN,"<",$in_file) || exit_on_error ("[ERROR] concatenate_FASTA_files can't open IN file '$in_file' for reading $!");
        my @in=<$IN>;
        print $OUT @in;
        close ($IN);
        @in=();
    }
    close ($OUT);
    print "[INFO] concatenate_FASTA_files: the output file '$out_single_fasta_file' was created\n";
}

# sub PanX_cleanup
# {
#     my $PanX_out_dir=shift;
#     # list all .nwk
#     my $trees_dir=$PanX_out_dir."/geneCluster/trees/";
#     make_path ($trees_dir);
#     my @PanX_trees=grep {/\.nwk$/} <$PanX_out_dir/geneCluster/*>;
#     # for my $tree (@PanX_trees) {mv ($tree,$trees_dir);}
#     # Create a new tar object:
#     my $tar = Archive::Tar->new();
#     # Add some files:
#     $tar->add_files(@PanX_trees);
#     # Finished:
#     $tar->write( 'file.tar',COMPRESS_GZIP);
# }
sub PanX_cleanup
{
    my $PanX_out_dir=shift;
    
    # all file types as a hash, key=resulted_out_tar_file
    my %PanX_files_to_tar=();
    $PanX_files_to_tar{"PanX_OGs_trees_nwk.tar.gz"}="\.nwk\$";
    $PanX_files_to_tar{"PanX_OGs_trees_json.tar.gz"}="_tree.json\$";
    $PanX_files_to_tar{"PanX_OGs_faa.tar.gz"}="\.faa\$";
    $PanX_files_to_tar{"PanX_OGs_fna.tar.gz"}="\.fna\$";
    $PanX_files_to_tar{"PanX_OGs_na_aln.tar.gz"}="_na_aln.fa\$";
    $PanX_files_to_tar{"PanX_OGs_na_aln_reduced.tar.gz"}="_na_aln_reduced.fa\$";
    $PanX_files_to_tar{"PanX_OGs_aa_aln_reduced.tar.gz"}="_aa_aln_reduced.fa\$";

    my $cwd=getcwd();
    chdir "$PanX_out_dir/geneCluster/";
    
    for my $out_tar_file (keys %PanX_files_to_tar)
    {
        my $files_suffix=$PanX_files_to_tar{$out_tar_file};
        if (defined $files_suffix and $files_suffix ne "") 
        {
            my @files_to_tar=grep {/$files_suffix/} <*>;
 #           print scalar(@files_to_tar)," '$files_suffix' files were found\n";<STDIN>;
            my $num_of_files=@files_to_tar;
            if ($num_of_files>0)
            {
                print "[INFO] Tar: $num_of_files '$files_suffix' files were tar.gz into '$PanX_out_dir/geneCluster/$out_tar_file'\n";
                # Create a new tar object:
                my $tar = Archive::Tar->new();
                # Add all trees files:
                $tar->add_files(@files_to_tar);
                # Finished:
                $tar->write($PanX_out_dir."/geneCluster/$out_tar_file",COMPRESS_GZIP);
                unlink (@files_to_tar); # clean
            }
        }
    }
    chdir $cwd;
#     print "now at: ",getcwd(),"\n";
}
sub exit_on_error {
    my $error_string=shift;
    my $datestring = strftime "%a %b %e %H:%M:%S %Y", localtime;
    print $LOG "$datestring $error_string\n";
    close $LOG;
    print "$datestring $error_string\n";
    exit (1);
}

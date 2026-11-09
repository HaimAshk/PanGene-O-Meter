package PGOM_Names;
# ---------------------------------------------------------------------------
# PGOM_Names.pm -- ONE place that defines how PanGene-O-Meter names genomes
# and proteins, so that every input path (gbk/gbff, fna+pyrodigal, presence/
# absence matrix, PanX) produces identical, trivially parseable names.
#
#   genome   : S1 .. SN   (numbered in input-list order, after dropping bad inputs)
#   protein  : S<n><SEP><locus>          e.g.  S12|WP_011234567   (SEP = '|')
#   map file : genome_names.map.txt      original_input <TAB> original_base_name <TAB> S<n>
#              (same three columns as the current PanGeneOmeter.pl writes)
#
# The separator is the one the existing DIAMOND phyletic-pattern builder already
# splits on:  build_phyletic_pattern_based_on_DIAMONDclusters(..., "\\|", 1).
# To change it, change $SEP below -- nothing else in the modules hard-codes it.
#
# Core Perl modules only.
# ---------------------------------------------------------------------------
use strict;
use warnings;
use File::Basename qw(basename);
use Exporter 'import';

our $VERSION = '0.1';
our $SEP = '|';

our @EXPORT_OK = qw(
    $SEP protein_header split_protein_id clean_locus
    strip_genome_suffix genome_type
    assign_short_names write_name_map read_name_map read_genome_list
);

# suffixes removed to obtain the "original base name" (".gz" is handled separately).
# Superset of the lists in the current script (@gbk_suffix and the one in cluster_cd_hit_like).
my @SUFFIXES = (
    '.gbk.fna', '.genome.fasta', '.short.fasta', '.short.fna', '.short.fa',
    '.gbff', '.gbk', '.gb', '.genbank', '.gbf',
    '.fasta', '.fna', '.fa', '.fas', '.ffn',
);
my @SUFFIXES_LONGEST_FIRST = sort { length($b) <=> length($a) } @SUFFIXES;

# ---- protein ids ----------------------------------------------------------
# Locus text may not contain whitespace (DIAMOND cuts the header there) nor the separator.
sub clean_locus {
    my ($l) = @_;
    $l = '' unless defined $l;
    $l =~ s/^\s+|\s+$//g;
    $l =~ s/\s+/_/g;
    my $sep_re = quotemeta($SEP);
    $l =~ s/$sep_re/_/g;
    return $l;
}

sub protein_header {            # ($S, $locus) -> "S3|locus"   (no leading '>')
    my ($s, $locus) = @_;
    return $s . $SEP . clean_locus($locus);
}

sub split_protein_id {          # "S3|locus" -> ("S3", "locus");  first separator only
    my ($id) = @_;
    my $i = index($id, $SEP);
    return ($id, '') if $i < 0;
    return (substr($id, 0, $i), substr($id, $i + length($SEP)));
}

# ---- genome names ---------------------------------------------------------
sub strip_genome_suffix {       # path -> original base name
    my ($path) = @_;
    my $n = basename($path);
    $n =~ s/\.gz$//i;
    for my $suf (@SUFFIXES_LONGEST_FIRST) {
        last if $n =~ s/\Q$suf\E$//i;
    }
    return $n;
}

# 'gbk' | 'fasta' | undef  (by suffix; unknown suffix -> look at the first characters)
sub genome_type {
    my ($path) = @_;
    my $n = basename($path);
    $n =~ s/\.gz$//i;
    return 'gbk'   if $n =~ /\.(?:gbff|gbk|gb|genbank|gbf)$/i;
    return 'fasta' if $n =~ /\.(?:fasta|fna|fa|fas|ffn)$/i;
    my $fh;
    if ($path =~ /\.gz$/i) { open($fh, '-|', 'gzip', '-dc', $path) or return undef; }
    else                   { open($fh, '<', $path)                 or return undef; }
    my $first = '';
    while (defined(my $l = <$fh>)) { next if $l =~ /^\s*$/; $first = $l; last; }
    close $fh;
    return 'gbk'   if $first =~ /^LOCUS\s/;
    return 'fasta' if $first =~ /^>/;
    return undef;
}

# Assign S1..SN in the order given.
#   items     : arrayref of file paths (base name = suffix stripped) or, with names=>1,
#               of raw genome names (used verbatim, e.g. matrix row names)
#   returns   : { order=>[S1..], old2new=>{base=>S}, new2old=>{S=>base}, rows=>[[orig,base,S],..] }
#   dies      : "[FATAL_ERROR] More than one genome with the same name ..." (wording kept)
sub assign_short_names {
    my ($items, %o) = @_;
    my (%old2new, %new2old, @rows, @order);
    my $n = $o{start} || 0;
    for my $it (@$items) {
        my $base = $o{names} ? $it : strip_genome_suffix($it);
        $n++;
        my $s = "S$n";
        die "[FATAL_ERROR] More than one genome with the same name was provided: '$base' is duplicate!. "
          . "Please make sure all names are unique and run again...\n" if exists $old2new{$base};
        $old2new{$base} = $s;
        $new2old{$s}    = $base;
        push @order, $s;
        push @rows, [$it, $base, $s];
    }
    return { order => \@order, old2new => \%old2new, new2old => \%new2old, rows => \@rows };
}

sub write_name_map {            # rows = [[orig,base,S],...]
    my ($file, $rows) = @_;
    open(my $o, '>', $file) or die "[ERROR] Can't write name map '$file': $!\n";
    print $o join("\t", @$_), "\n" for @$rows;
    close $o;
    return $file;
}

sub read_name_map {             # -> { old2new, new2old, order }
    my ($file) = @_;
    open(my $in, '<', $file) or die "[ERROR] Can't read name map '$file': $!\n";
    my (%o2n, %n2o, @order);
    while (<$in>) {
        chomp; next if /^\s*$/;
        my ($orig, $base, $s) = split /\t/;
        die "[ERROR] Name map '$file' line $.: need 3 tab-separated columns\n" unless defined $s;
        $o2n{$base} = $s; $n2o{$s} = $base; push @order, $s;
    }
    close $in;
    return { old2new => \%o2n, new2old => \%n2o, order => \@order };
}

# list of genome files: one path per line; blank lines and '#' comments ignored
sub read_genome_list {
    my ($file) = @_;
    open(my $in, '<', $file) or die "[ERROR] Can't open genome list '$file': $!\n";
    my @l;
    while (<$in>) {
        s/\r?\n$//; s/^\s+|\s+$//g;
        next if $_ eq '' || /^#/;
        push @l, $_;
    }
    close $in;
    return \@l;
}

1;

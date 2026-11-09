package PGOM_Input;
#
# PGOM_Input.pm -- pluggable input front-ends for PanGene-O-Meter
#
#   1. read_phyletic_pattern()  start from a precomputed presence/absence
#                               matrix (Microbializer .fas, PanGene-O-Meter
#                               .fas, or PanGene-O-Meter .csv), validated.
#   2. run_pyrodigal()          genome .fna -> proteins (.faa) + GFF, with
#                               protein headers formatted S<n>|LOCUS (PGOM_Names).
#      run_pyrodigal_many()     the same for a whole set, plus one combined
#                               .faa ready for DIAMOND.
#
# Core Perl modules only. External programs are called with the list form
# of system() (no shell, so no quoting problems with odd file names).
# The pyrodigal binary can be overridden with $ENV{PGOM_PYRODIGAL}.
#
# Nothing here touches global state, so it can be dropped into the existing
# script with:   use lib '/path/to/dir'; use PGOM_Input;

use strict;
use warnings;
use File::Basename qw(fileparse);
use File::Temp qw(tempfile);
use FindBin qw($Bin);
use lib $Bin;
use PGOM_Names qw($SEP protein_header strip_genome_suffix);
use Exporter 'import';

our @EXPORT_OK = qw(
    read_phyletic_pattern pa_genome_sizes rename_pattern
    sample_name_from_file run_pyrodigal run_pyrodigal_many
);


# ---------------------------------------------------------------------------
# 1. Presence/absence matrix entry point
# ---------------------------------------------------------------------------
# read_phyletic_pattern($file, names_file => $f, min_genomes => 2)
#
# Accepts:
#   FASTA : ">genome" then one line of 0/1 characters (Microbializer and
#           PanGene-O-Meter's own *.representative_PhP.fas)
#   CSV   : header line, then "genome,0,1,0,..." (PanGene-O-Meter's own CSV)
# The format is detected from the first non-blank character ('>' = FASTA).
#
# Returns a hashref:
#   pattern  => { genome => [ '0','1',... ] }   (same structure the existing
#                                                code builds as %PhyleticPattern)
#   order    => [ genome names in file order ]
#   header   => CSV header line ("genome,OG1,...\n") for writing CSV outputs
#   n_og     => number of orthology groups (columns)
#   format   => 'fasta' | 'csv'
#   warnings => [ strings ]  (non-fatal; caller decides whether to print)
# Dies with a message starting "[ERROR] " on any invalid input.
sub read_phyletic_pattern {
    my ($file, %opt) = @_;
    my $min_genomes = defined $opt{min_genomes} ? $opt{min_genomes} : 2;

    open(my $fh, '<', $file) or die "[ERROR] Can't open presence/absence file '$file': $!\n";
    my @lines;
    while (my $l = <$fh>) {
        $l =~ s/[\r\n]+$//;
        next unless $l =~ /\S/;
        push @lines, $l;
    }
    close $fh;
    die "[ERROR] Presence/absence file '$file' is empty.\n" unless @lines;

    my (%pattern, @order, $header, $format);
    my @warnings;

    if ($lines[0] =~ /^>/) {
        $format = 'fasta';
        my ($name, @seq_parts);
        my $flush = sub {
            return unless defined $name;
            my $seq = join('', @seq_parts);
            die "[ERROR] Genome '$name' has an empty pattern in '$file'.\n" unless length $seq;
            die "[ERROR] Duplicate genome ID '$name' in '$file'.\n" if exists $pattern{$name};
            $pattern{$name} = $seq;       # string for now, split after length check
            push @order, $name;
        };
        for my $l (@lines) {
            if ($l =~ /^>\s*(\S+)/) {
                $flush->();
                ($name, @seq_parts) = ($1);
            } else {
                die "[ERROR] Sequence line before the first '>' header in '$file'.\n" unless defined $name;
                (my $s = $l) =~ s/\s+//g;
                push @seq_parts, $s;
            }
        }
        $flush->();
    } else {
        $format = 'csv';
        my $hdr = shift @lines;
        my @cols = split(/,/, $hdr, -1);
        die "[ERROR] CSV header in '$file' has no orthology-group columns.\n" if @cols < 2;
        $header = $hdr . "\n";
        for my $l (@lines) {
            my ($name, @bits) = split(/,/, $l, -1);
            die "[ERROR] CSV row without a genome name in '$file'.\n" unless defined $name && length $name;
            die "[ERROR] Duplicate genome ID '$name' in '$file'.\n" if exists $pattern{$name};
            $pattern{$name} = join('', @bits);
            push @order, $name;
            if (grep { $_ !~ /^[01]$/ } @bits) {
                die "[ERROR] Genome '$name' in '$file' has values other than 0/1 "
                  . "(counts or paralog numbers must be converted to presence/absence first).\n";
            }
        }
        my $expected = scalar(@cols) - 1;
        for my $n (@order) {
            my $len = length $pattern{$n};
            die "[ERROR] Genome '$n' has $len columns but the header declares $expected.\n"
                if $len != $expected;
        }
    }

    die "[ERROR] Only " . scalar(@order) . " genome(s) in '$file'; at least $min_genomes are needed.\n"
        if @order < $min_genomes;

    # --- common validation ---------------------------------------------------
    my $n_og = length $pattern{ $order[0] };
    for my $n (@order) {
        my $s = $pattern{$n};
        die "[ERROR] Genome '$n' has length " . length($s) . " but '$order[0]' has $n_og; "
          . "all genomes must have the same number of orthology groups.\n"
            if length($s) != $n_og;
        if ($s =~ /([^01])/) {
            die "[ERROR] Genome '$n' contains the character '$1'; only 0 and 1 are allowed.\n";
        }
        push @warnings, "Genome '$n' has no genes at all (all zeros)." if $s !~ /1/;
    }

    # names file (optional): one OG name per line, count must match
    my @og_names;
    if (defined $opt{names_file}) {
        open(my $nf, '<', $opt{names_file}) or die "[ERROR] Can't open names file '$opt{names_file}': $!\n";
        while (my $l = <$nf>) { $l =~ s/[\r\n]+$//; push @og_names, $l if length $l; }
        close $nf;
        die "[ERROR] Names file has " . scalar(@og_names) . " entries but the matrix has $n_og columns.\n"
            if @og_names != $n_og;
    }
    if (!defined $header) {
        @og_names = map { "OG$_" } (1 .. $n_og) unless @og_names;
        $header = join(',', 'genome', @og_names) . "\n";
    }

    # columns that carry no information (informational only)
    my ($all0, $all1) = (0, 0);
    my @sums = (0) x $n_og;
    for my $n (@order) {
        my $i = 0;
        $sums[$i++] += $_ for split //, $pattern{$n};
    }
    for my $s (@sums) { $all0++ if $s == 0; $all1++ if $s == @order; }
    push @warnings, "$all0 of $n_og columns are absent from every genome (no effect on GCS)." if $all0;
    push @warnings, "$all1 of $n_og columns are present in every genome (core; they raise all similarities)." if $all1;

    # convert strings to the array-of-chars structure the existing code expects
    my %arr = map { $_ => [ split //, $pattern{$_} ] } @order;

    return {
        pattern  => \%arr,
        order    => \@order,
        header   => $header,
        n_og     => $n_og,
        format   => $format,
        warnings => \@warnings,
    };
}

# Number of present genes per genome. The existing greedy clustering orders
# genomes by nucleotide length to choose representatives; a matrix-only input
# has no genome lengths, so use this as the ordering key instead and say so.
sub pa_genome_sizes {
    my ($pa) = @_;
    my %size;
    for my $n (@{ $pa->{order} }) {
        my $c = 0;
        $c += $_ for @{ $pa->{pattern}{$n} };
        $size{$n} = $c;
    }
    return \%size;
}

# ---------------------------------------------------------------------------
# 2. fna -> pyrodigal -> proteins named S<n>|LOCUS_ID (see PGOM_Names)
# ---------------------------------------------------------------------------

# Original base name from a genome file name (suffix handling lives in PGOM_Names).
sub sample_name_from_file {
    my ($path) = @_;
    return strip_genome_suffix($path);
}

# Re-key a pattern from read_phyletic_pattern() with a map {old_name => 'S<n>'}
# (from PGOM_Names::assign_short_names(\@{$pa->{order}}, names=>1)).
# Everything downstream then works on S1..SN exactly like the gbk path.
sub rename_pattern {
    my ($pa, $old2new) = @_;
    my (%pat, @order);
    for my $g (@{$pa->{order}}) {
        die "[ERROR] rename_pattern: no new name for '$g'\n" unless exists $old2new->{$g};
        $pat{$old2new->{$g}} = $pa->{pattern}{$g};
        push @order, $old2new->{$g};
    }
    return { %$pa, pattern => \%pat, order => \@order };
}

# run_pyrodigal($fna, $sample, $outdir, threads => 1, translation_table => 11)
# Returns { sample, faa, gff, n_genes, mode, total_bp }.
# Writes  $outdir/$sample.faa  (headers: >S<n>|LOCUS_ID, no trailing '*')
#         $outdir/$sample.gff  (as produced by pyrodigal)
sub run_pyrodigal {
    my ($fna, $sample, $outdir, %opt) = @_;
    my $bin     = $ENV{PGOM_PYRODIGAL} || 'pyrodigal';
    my $threads = $opt{threads} || 1;
    my $table   = $opt{translation_table} || 11;

    die "[ERROR] Genome file '$fna' not found.\n" unless -f $fna;
    die "[ERROR] Sample name '$sample' contains the separator '$SEP' or whitespace; "
      . "use the S<n> names from PGOM_Names.\n" if index($sample, $SEP) >= 0 || $sample =~ /\s/;
    die "[ERROR] Sample name is empty for '$fna'.\n" unless length $sample;
    mkdir $outdir unless -d $outdir;
    die "[ERROR] Can't create output directory '$outdir': $!\n" unless -d $outdir;

    # decompress to a temp file if needed (avoids relying on pyrodigal's gz support)
    my $input = $fna;
    my $tmp_in;
    if ($fna =~ /\.gz$/) {
        (my $fh, $tmp_in) = tempfile(SUFFIX => '.fna', UNLINK => 1);
        close $fh;
        system("gzip -dc " . _q($fna) . " > " . _q($tmp_in)) == 0
            or die "[ERROR] Could not decompress '$fna'.\n";
        $input = $tmp_in;
    }

    my $total_bp = _fasta_length($input);
    die "[ERROR] No sequence found in '$fna'.\n" unless $total_bp > 0;

    my $faa_raw = "$outdir/$sample.raw.faa";
    my $gff     = "$outdir/$sample.gff";

    # Single mode trains on the genome itself; pyrodigal itself warns when
    # given less than 100 kb, so below that use meta mode (pre-trained models).
    my $mode = $total_bp >= 100_000 ? 'single' : 'meta';
    my $warn;
    my $run = sub {
        my ($m) = @_;
        my @cmd = ($bin, '-i', $input, '-a', $faa_raw, '-o', $gff, '-f', 'gff',
                   '-p', $m, '-g', $table, '-j', $threads);
        return system(@cmd);
    };
    my $rc = $run->($mode);
    if ($rc != 0 && $mode eq 'single') {
        $warn = "pyrodigal single mode failed for '$sample'; retried with meta mode.";
        $mode = 'meta';
        $rc   = $run->($mode);
    }
    die "[ERROR] pyrodigal failed on '$fna' (is '$bin' installed? set PGOM_PYRODIGAL to its path).\n"
        if $rc != 0;

    # rewrite headers to S<n>|LOCUS_ID, strip trailing '*'
    my $faa = "$outdir/$sample.faa";
    open(my $in,  '<', $faa_raw) or die "[ERROR] Can't read '$faa_raw': $!\n";
    open(my $out, '>', $faa)     or die "[ERROR] Can't write '$faa': $!\n";
    my ($n_genes, $seen) = (0, {});
    while (my $l = <$in>) {
        $l =~ s/[\r\n]+$//;
        if ($l =~ /^>(\S+)/) {
            my $locus = PGOM_Names::clean_locus($1);   # e.g. contig_12 (pyrodigal: contig_geneIndex)
            die "[ERROR] Duplicate locus ID '$locus' in pyrodigal output for '$sample'.\n" if $seen->{$locus}++;
            print $out '>', protein_header($sample, $locus), "\n";
            $n_genes++;
        } else {
            $l =~ s/\*+$//;
            print $out "$l\n";
        }
    }
    close $in;
    close $out;
    unlink $faa_raw;

    return {
        sample => $sample, faa => $faa, gff => $gff,
        n_genes => $n_genes, mode => $mode, total_bp => $total_bp,
        warning => $warn,
    };
}

# run_pyrodigal_many(\%samples, $outdir, combined_faa => $path, threads => N)
#   %samples : sample name => path to .fna
# Returns { results => [ {...}, ... ], combined_faa => $path, warnings => [...] }.
sub run_pyrodigal_many {
    my ($samples, $outdir, %opt) = @_;
    mkdir $outdir unless -d $outdir;
    die "[ERROR] Can't create output directory '$outdir': $!\n" unless -d $outdir;
    my $combined = $opt{combined_faa} || "$outdir/all_proteins.faa";
    my (@results, @warnings);

    open(my $all, '>', $combined) or die "[ERROR] Can't write '$combined': $!\n";
    for my $s (sort keys %$samples) {
        my $r = run_pyrodigal($samples->{$s}, $s, $outdir,
                              threads => $opt{threads}, translation_table => $opt{translation_table});
        push @results, $r;
        push @warnings, $r->{warning} if $r->{warning};
        # very low gene density usually means a bad assembly or wrong translation table
        my $per_mb = $r->{total_bp} ? $r->{n_genes} / ($r->{total_bp} / 1e6) : 0;
        push @warnings, sprintf("'%s' has only %.0f genes/Mb (%d genes in %d bp); check the assembly or translation table.",
                                $s, $per_mb, $r->{n_genes}, $r->{total_bp}) if $r->{total_bp} >= 100_000 && $per_mb < 500;

        open(my $in, '<', $r->{faa}) or die "[ERROR] Can't read '$r->{faa}': $!\n";
        local $/ = \65536;
        while (my $chunk = <$in>) { print $all $chunk; }
        close $in;
    }
    close $all;
    return { results => \@results, combined_faa => $combined, warnings => \@warnings };
}

# ---------------------------------------------------------------------------
sub _fasta_length {
    my ($f) = @_;
    open(my $fh, '<', $f) or die "[ERROR] Can't read '$f': $!\n";
    my $n = 0;
    while (my $l = <$fh>) { next if $l =~ /^>/; $l =~ s/\s+//g; $n += length $l; }
    close $fh;
    return $n;
}
sub _q { my $s = shift; $s =~ s/'/'\\''/g; return "'$s'"; }

1;

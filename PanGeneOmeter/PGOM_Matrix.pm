package PGOM_Matrix;
# ---------------------------------------------------------------------------
# PGOM_Matrix.pm -- read any common presence/absence output into the canonical
# PanGene-O-Meter structure (same as PGOM_Input::read_phyletic_pattern) and write
# it back in PanGene-O-Meter's own files.
#
# Reads (format auto-detected):
#   rtab   Panaroo / Roary  gene_presence_absence.Rtab   (tab, genes = rows, genomes = columns, 0/1)
#   roary  Roary / Panaroo  gene_presence_absence.csv    (comma, quoted; non-empty cell = present)
#   fasta  Microbializer / PanGene-O-Meter .fas          (via PGOM_Input)
#   csv    PanGene-O-Meter .01.csv                       (via PGOM_Input)
#
# Core Perl modules only.
# ---------------------------------------------------------------------------
use strict;
use warnings;
use FindBin qw($Bin);
use lib $Bin;
use PGOM_Input qw(read_phyletic_pattern);
use Exporter 'import';

our $VERSION = '0.1';
our @EXPORT_OK = qw(read_matrix_any detect_matrix_format write_canonical_matrix);

# metadata columns of Roary (14) and Panaroo (3 + optional) gene_presence_absence.csv
my %META = map { lc($_) => 1 } (
    'Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 'No. sequences',
    'Avg sequences per isolate', 'Genome Fragment', 'Order within Fragment',
    'Accessory Fragment', 'Accessory Order with Fragment', 'QC', 'Min group size nuc',
    'Max group size nuc', 'Avg group size nuc',
);

sub _first_line {
    my ($f) = @_;
    open(my $in, '<', $f) or die "[ERROR] Can't open '$f': $!\n";
    my $l;
    while (defined($l = <$in>)) { last if $l =~ /\S/; }
    close $in;
    die "[ERROR] '$f' is empty\n" unless defined $l;
    $l =~ s/\r?\n$//;
    return $l;
}

sub detect_matrix_format {
    my ($f) = @_;
    my $l = _first_line($f);
    return 'fasta' if $l =~ /^>/;
    return 'roary' if $l =~ /^"?Gene"?,/i && $l =~ /Annotation/i;
    return 'rtab'  if $l =~ /\t/ && $l !~ /,.*,.*,/;
    return 'rtab'  if $l =~ /\t/;
    return 'csv';
}

# minimal CSV line parser (quotes, doubled quotes, commas inside quotes)
sub _csv_split {
    my ($l) = @_;
    my @f; my $cur = ''; my $q = 0;
    my @c = split //, $l;
    for (my $i = 0; $i < @c; $i++) {
        my $ch = $c[$i];
        if ($q) {
            if ($ch eq '"') { if ($i + 1 < @c && $c[$i+1] eq '"') { $cur .= '"'; $i++ } else { $q = 0 } }
            else { $cur .= $ch }
        } else {
            if    ($ch eq '"') { $q = 1 }
            elsif ($ch eq ',') { push @f, $cur; $cur = '' }
            else  { $cur .= $ch }
        }
    }
    push @f, $cur;
    return @f;
}

sub _finish {
    my ($genomes, $rows, $ognames, $fmt, $warn) = @_;   # rows: [ [flags per genome], ... ]
    my %pat; my $n_og = scalar @$rows;
    for my $gi (0 .. $#$genomes) { $pat{$genomes->[$gi]} = [ map { $_->[$gi] } @$rows ]; }
    my %seen;
    for my $g (@$genomes) { die "[ERROR] Duplicate genome column '$g'\n" if $seen{$g}++; }
    die "[ERROR] Need at least 2 genomes (found " . scalar(@$genomes) . ")\n" if @$genomes < 2;
    die "[ERROR] No gene rows found\n" unless $n_og;
    for my $g (@$genomes) {
        push @$warn, "genome '$g' has no genes at all" unless grep { $_ eq '1' } @{$pat{$g}};
    }
    return { pattern => \%pat, order => [@$genomes], n_og => $n_og, format => $fmt,
             og_names => $ognames, warnings => $warn, header => "Species," . join(',', @$ognames) . "\n" };
}

sub _read_rtab {
    my ($f) = @_;
    open(my $in, '<', $f) or die "[ERROR] Can't open '$f': $!\n";
    my ($hdr, @genomes, @rows, @names, @warn);
    while (my $l = <$in>) {
        $l =~ s/\r?\n$//; next unless $l =~ /\S/;
        my @c = split /\t/, $l, -1;
        if (!$hdr) { $hdr = 1; shift @c; @genomes = @c; die "[ERROR] Rtab header has no genome columns\n" unless @genomes; next; }
        my $name = shift @c;
        die "[ERROR] Rtab line $.: expected " . scalar(@genomes) . " values, found " . scalar(@c) . "\n" unless @c == @genomes;
        for my $v (@c) { die "[ERROR] Rtab line $.: value '$v' is not 0/1\n" unless $v =~ /^[01]$/; }
        push @names, $name; push @rows, [@c];
    }
    close $in;
    return _finish(\@genomes, \@rows, \@names, 'rtab', \@warn);
}

sub _read_roary {
    my ($f) = @_;
    open(my $in, '<', $f) or die "[ERROR] Can't open '$f': $!\n";
    my ($hdr, @gcols, @genomes, @rows, @names, @warn);
    while (my $l = <$in>) {
        $l =~ s/\r?\n$//; next unless $l =~ /\S/;
        my @c = _csv_split($l);
        if (!$hdr) {
            $hdr = [@c];
            for my $i (0 .. $#c) { next if $META{lc $c[$i]}; push @gcols, $i; push @genomes, $c[$i]; }
            die "[ERROR] No genome columns found in the header (only metadata columns)\n" unless @gcols;
            next;
        }
        push @names, ($c[0] // '');
        push @rows, [ map { (defined $c[$_] && $c[$_] ne '' && $c[$_] ne 'NA') ? 1 : 0 } @gcols ];
    }
    close $in;
    return _finish(\@genomes, \@rows, \@names, 'roary', \@warn);
}

# read_matrix_any($file, format => 'auto', names_file => ...) -> hashref like read_phyletic_pattern(+ og_names)
sub read_matrix_any {
    my ($f, %o) = @_;
    my $fmt = $o{format} && $o{format} ne 'auto' ? $o{format} : detect_matrix_format($f);
    my $pa;
    if    ($fmt eq 'rtab')  { $pa = _read_rtab($f); }
    elsif ($fmt eq 'roary') { $pa = _read_roary($f); }
    elsif ($fmt eq 'fasta' || $fmt eq 'csv') {
        $pa = read_phyletic_pattern($f, names_file => $o{names_file});
        my @h = split /,/, $pa->{header}; chomp $h[-1]; shift @h;
        $pa->{og_names} = \@h;
    } else { die "[ERROR] Unknown matrix format '$fmt' (use auto|rtab|roary|fasta|csv)\n"; }
    return $pa;
}

# writes <prefix>.fas, <prefix>.01.csv and <prefix>.og_names.txt ; returns the three paths
sub write_canonical_matrix {
    my ($pa, $prefix) = @_;
    my ($fas, $csv, $nm) = ("$prefix.fas", "$prefix.01.csv", "$prefix.og_names.txt");
    open(my $F, '>', $fas) or die "[ERROR] Can't write '$fas': $!\n";
    open(my $C, '>', $csv) or die "[ERROR] Can't write '$csv': $!\n";
    print $C "Species,", join(',', @{$pa->{og_names}}), "\n";
    for my $g (@{$pa->{order}}) {
        print $F ">$g\n", join('', @{$pa->{pattern}{$g}}), "\n";
        print $C join(',', $g, @{$pa->{pattern}{$g}}), "\n";
    }
    close $F; close $C;
    open(my $N, '>', $nm) or die "[ERROR] Can't write '$nm': $!\n";
    print $N "pos\tname\n";
    my $i = 0; print $N ++$i, "\t$_\n" for @{$pa->{og_names}};
    close $N;
    return ($fas, $csv, $nm);
}

1;

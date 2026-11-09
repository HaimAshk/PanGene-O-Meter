package PGOM_Convert;
# ---------------------------------------------------------------------------
# PGOM_Convert.pm -- turn user genomes into what the DIAMOND path needs, WITHOUT PanX:
#
#   gbk / gbff / gb (.gz ok) -> protein FASTA (S<n>|locus headers) + genome length
#   fna / fasta     (.gz ok) -> pyrodigal -> protein FASTA (S<n>|contig_idx) + genome length
#
#   prepare_genomes() does the whole front end for a list of files: names S1..SN
#   (PGOM_Names), per-genome protein files, genome lengths, name map, combined .faa.
#
# Replaces, in PanGeneOmeter.pl, the PanX "-st 1 3 4" protein extraction, the
# copy-and-rename of the gbk files, filter_gbk_without_cds(), get_genome_length_from_gbk()
# and concatenate_FASTA_files().
#
# Core Perl modules only (pyrodigal is called through PGOM_Input::run_pyrodigal).
# ---------------------------------------------------------------------------
use strict;
use warnings;
use File::Path qw(make_path);
use FindBin qw($Bin);
use lib $Bin;
use PGOM_Names qw(protein_header assign_short_names write_name_map genome_type clean_locus);
use PGOM_Input qw(run_pyrodigal);
use Exporter 'import';

our $VERSION = '0.1';
our @EXPORT_OK = qw(gbk_to_faa fasta_total_length prepare_genomes concat_files);

my $MIN_PROTEINS = 2;     # same threshold as filter_gbk_without_cds()

sub _open_in {
    my ($f) = @_;
    my $fh;
    if ($f =~ /\.gz$/i) { open($fh, '-|', 'gzip', '-dc', $f) or die "[ERROR] Can't read '$f' (gzip): $!\n"; }
    else                { open($fh, '<', $f)                 or die "[ERROR] Can't open '$f': $!\n"; }
    return $fh;
}

# ---- FASTA length ----------------------------------------------------------
# total number of sequence characters (all records)
sub fasta_total_length {
    my ($f) = @_;
    my $fh = _open_in($f);
    my ($len, $recs) = (0, 0);
    while (<$fh>) {
        if (/^>/) { $recs++; next; }
        s/\s+//g;
        $len += length($_);
    }
    close $fh;
    die "[ERROR] '$f' contains no FASTA records\n" unless $recs;
    return wantarray ? ($len, $recs) : $len;
}

# ---- translation (fallback for CDS without /translation) -------------------
my %CODON;
{
    my @b = qw(T C A G);
    my $aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG';
    my $i = 0;
    for my $x (@b) { for my $y (@b) { for my $z (@b) { $CODON{"$x$y$z"} = substr($aa, $i++, 1); } } }
}
my %TABLE_TWEAK = ( 1 => {}, 11 => {}, 4 => { TGA => 'W' } );
my %ALT_START   = map { $_ => 1 } qw(TTG CTG ATT ATC ATA GTG);

sub _revcomp { my $s = reverse shift; $s =~ tr/ACGTacgt/TGCAtgca/; return $s; }

sub _translate {
    my ($dna, $table, $force_m) = @_;
    my $tw = $TABLE_TWEAK{$table} or return undef;
    $dna = uc $dna;
    my $p = '';
    for (my $i = 0; $i + 2 < length($dna); $i += 3) {
        my $c = substr($dna, $i, 3);
        $p .= exists $tw->{$c} ? $tw->{$c} : ($CODON{$c} // 'X');
    }
    if ($force_m && length($dna) >= 3 && ($ALT_START{substr($dna,0,3)} || substr($dna,0,3) eq 'ATG')) { substr($p,0,1) = 'M'; }
    return $p;
}

# GenBank location -> nucleotide sequence (undef if unsupported: remote refs, '^', unparsable)
sub _location_seq {
    my ($loc, $seqref) = @_;
    $loc =~ s/\s+//g;
    return undef if $loc =~ /\^/ || $loc =~ /[A-Za-z_][A-Za-z0-9_]*(?:\.\d+)?:/;
    $loc =~ s/[<>]//g;
    my $parse;
    $parse = sub {
        my ($s) = @_;
        if ($s =~ /^complement\((.*)\)$/s) { my $x = $parse->($1); return defined $x ? _revcomp($x) : undef; }
        if ($s =~ /^(?:join|order)\((.*)\)$/s) {
            my $inner = $1; my @parts; my ($depth, $cur) = (0, '');
            for my $ch (split //, $inner) {
                if    ($ch eq '(') { $depth++; $cur .= $ch }
                elsif ($ch eq ')') { $depth--; $cur .= $ch }
                elsif ($ch eq ',' && $depth == 0) { push @parts, $cur; $cur = '' }
                else { $cur .= $ch }
            }
            push @parts, $cur;
            my $out = '';
            for my $p (@parts) { my $x = $parse->($p); return undef unless defined $x; $out .= $x; }
            return $out;
        }
        if ($s =~ /^(\d+)\.\.(\d+)$/) {
            return undef if $1 < 1 || $2 < $1 || $2 > length($$seqref);
            return substr($$seqref, $1 - 1, $2 - $1 + 1);
        }
        if ($s =~ /^(\d+)$/) { return undef if $1 < 1 || $1 > length($$seqref); return substr($$seqref, $1 - 1, 1); }
        return undef;
    };
    return $parse->($loc);
}

# ---- GenBank / GBFF -> faa -------------------------------------------------
# gbk_to_faa($gbk, $S, $out_faa)  ->  hashref
#   length          sum of the LOCUS "<n> bp" fields (the value the current code uses)
#   n_records       number of LOCUS records (contigs)
#   n_cds           CDS features seen
#   n_written       proteins written to $out_faa   (headers  >S<n>|locus)
#   n_pseudo        skipped: /pseudo or /pseudogene
#   n_no_protein    skipped: no /translation and it could not be derived from the sequence
#   n_translated    proteins translated by us because /translation was absent
#   id_source       { locus_tag=>n, protein_id=>n, gene=>n, generated=>n }
#   warnings        arrayref
sub gbk_to_faa {
    my ($gbk, $S, $out_faa) = @_;
    my $fh = _open_in($gbk);
    open(my $out, '>', $out_faa) or die "[ERROR] Can't write '$out_faa': $!\n";

    my %st = (length => 0, n_records => 0, n_cds => 0, n_written => 0, n_pseudo => 0,
              n_no_protein => 0, n_translated => 0, id_source => {}, warnings => []);
    my %used;               # locus -> count, to keep ids unique inside a genome
    my ($in_features, $in_origin) = (0, 0);
    my (@cds, @cur, $seq, $locus_bp);

    my $flush_cds = sub {   # turn collected feature text into a hash
        return unless @cur;
        my @lines = @cur; @cur = ();
        my $head = shift @lines;                      # "CDS             location..."
        my ($loc) = $head =~ /^\s{5}\S+\s+(.*)$/;
        my @q; my $loc_done = 0; my %cds; my $last;
        $loc = defined $loc ? $loc : '';
        for my $l (@lines) {
            $l =~ s/^\s{21}//;
            if ($l =~ m{^/([A-Za-z_0-9]+)(?:=(.*))?$}) {
                $loc_done = 1; $last = $1;
                my $v = defined $2 ? $2 : '1';
                if (exists $cds{$last} && $last ne 'translation') { $cds{$last} .= ' ' . $v; } else { $cds{$last} = $v; }
            } elsif (!$loc_done) {
                $loc .= $l;
            } elsif (defined $last) {
                $cds{$last} .= ($last eq 'translation' ? '' : ' ') . $l;
            }
        }
        for my $k (keys %cds) { $cds{$k} =~ s/^"//; $cds{$k} =~ s/"\s*$//; $cds{$k} =~ s/""/"/g; }
        $cds{_loc} = $loc;
        push @cds, \%cds;
    };

    my $finish_record = sub {
        $flush_cds->();
        my $seqlen = defined $seq ? length($seq) : 0;
        if (defined $seq && defined $locus_bp && $seqlen && $seqlen != $locus_bp) {
            push @{$st{warnings}}, "LOCUS says $locus_bp bp but the sequence has $seqlen in a record of '$gbk'";
        }
        for my $c (@cds) {
            $st{n_cds}++;
            if (exists $c->{pseudo} || exists $c->{pseudogene}) { $st{n_pseudo}++; next; }
            my $prot = $c->{translation};
            if (defined $prot) { $prot =~ s/\s+//g; $prot =~ s/\*+$//; }
            if (!defined $prot || $prot eq '') {
                $prot = undef;
                if ($seqlen) {
                    my $dna = _location_seq($c->{_loc}, \$seq);
                    if (defined $dna) {
                        my $cs  = ($c->{codon_start} // 1); $cs = 1 unless $cs =~ /^[123]$/;
                        $dna = substr($dna, $cs - 1);
                        my $t   = $c->{transl_table} // 11; $t =~ s/\D//g; $t ||= 11;
                        my $partial = $c->{_loc} =~ /[<>]/;
                        my $p = _translate($dna, $t, !$partial);
                        if (defined $p) {
                            $p =~ s/\*+$//;
                            if ($p ne '' && $p !~ /\*/) { $prot = $p; $st{n_translated}++; }
                        }
                    }
                }
                if (!defined $prot) { $st{n_no_protein}++; next; }
            }
            # identifier
            my ($id, $src);
            if    (defined $c->{locus_tag} && $c->{locus_tag} ne '') { ($id, $src) = ($c->{locus_tag}, 'locus_tag'); }
            elsif (defined $c->{protein_id} && $c->{protein_id} ne '') { ($id, $src) = ($c->{protein_id}, 'protein_id'); }
            elsif (defined $c->{gene}       && $c->{gene} ne '')       { ($id, $src) = ($c->{gene}, 'gene'); }
            else { ($id, $src) = ('cds' . ($st{n_written} + 1), 'generated'); }
            $id = clean_locus($id);
            $id = 'cds' . ($st{n_written} + 1) if $id eq '';
            $st{id_source}{$src}++;
            if ($used{$id}++) { $id .= '_dup' . $used{$id}; }
            print $out '>', protein_header($S, $id), "\n";
            for (my $i = 0; $i < length($prot); $i += 80) { print $out substr($prot, $i, 80), "\n"; }
            $st{n_written}++;
        }
        @cds = (); $seq = undef; $locus_bp = undef;
    };

    while (my $line = <$fh>) {
        $line =~ s/\r?\n$//;
        if ($line =~ /^LOCUS\s/) {
            $finish_record->() if $st{n_records};
            $st{n_records}++;
            ($in_features, $in_origin) = (0, 0);
            if ($line =~ /^LOCUS\s+\S+\s+(\d+)\s+bp/) { $locus_bp = $1; $st{length} += $1; }
            else { push @{$st{warnings}}, "LOCUS line without 'bp' length: $line"; $locus_bp = undef; }
            next;
        }
        if ($line =~ m{^//}) { $in_features = $in_origin = 0; next; }
        if ($line =~ /^FEATURES\s/) { $in_features = 1; $in_origin = 0; next; }
        if ($line =~ /^ORIGIN/)     { $flush_cds->(); $in_features = 0; $in_origin = 1; $seq = ''; next; }
        if ($line =~ /^CONTIG/)     { $flush_cds->(); $in_features = 0; next; }
        if ($in_origin) { (my $s = $line) =~ s/[^A-Za-z]//g; $seq .= $s if length $s; next; }
        if ($in_features) {
            if ($line =~ /^ {5}(\S+)\s/) {                       # new feature
                $flush_cds->();
                @cur = ($line) if $1 eq 'CDS';
            } elsif (@cur && $line =~ /^ {21}/) {
                push @cur, $line;
            }
        }
    }
    $finish_record->() if $st{n_records};
    close $fh; close $out;

    die "[ERROR] '$gbk' has no LOCUS records -- not a GenBank/GBFF file\n" unless $st{n_records};
    if ($st{n_translated}) {
        push @{$st{warnings}}, "$st{n_translated} CDS had no /translation; translated from the sequence (tables 1/4/11 only)";
    }
    if ($st{n_no_protein}) {
        push @{$st{warnings}}, "$st{n_no_protein} CDS skipped: no /translation and not derivable from the sequence";
    }
    if (($st{id_source}{generated} // 0) > 0) {
        push @{$st{warnings}}, "$st{id_source}{generated} CDS had no locus_tag/protein_id/gene; generated ids used";
    }
    $st{faa} = $out_faa;
    return \%st;
}

# ---- concatenate ------------------------------------------------------------
sub concat_files {              # streaming copy (no slurping)
    my ($out_file, $files) = @_;
    open(my $out, '>', $out_file) or die "[ERROR] concat_files: can't write '$out_file': $!\n";
    for my $f (@$files) {
        open(my $in, '<', $f) or die "[ERROR] concat_files: can't read '$f': $!\n";
        local $/ = \1048576;
        while (my $chunk = <$in>) { print $out $chunk; }
        close $in;
    }
    close $out;
    return $out_file;
}

# ---- whole front end ---------------------------------------------------------
# prepare_genomes(\@files, $outdir, threads=>, map_file=>, combined_faa=>, stats_file=>)
#   * each file is one genome: gbk/gbff/gb(.gz) or fna/fasta/fa(.gz); types may be mixed
#   * genomes that cannot give >= 2 proteins are dropped (reported), and the numbering
#     continues without a gap  -> S1..SN is always dense, exactly like the current script
#   * two inputs with the same base name -> fatal error (same message as the current code)
#   returns { order, old2new, new2old, genomes_length, n_proteins, dropped=>[[file,reason]],
#             faa_files, combined_faa, map_file, stats_file, warnings,
#             inputs=>{S=>{file,type,gff,faa}} }
sub prepare_genomes {
    my ($files, $outdir, %o) = @_;
    my $threads = $o{threads} || 1;
    my $faa_dir = "$outdir/protein_faa";
    make_path($faa_dir) unless -d $faa_dir;
    die "[ERROR] Can't create '$faa_dir'\n" unless -d $faa_dir;

    # names first (so duplicates are caught before any work); numbering is re-done below if files are dropped
    my %seen;
    for my $f (@$files) {
        my $b = PGOM_Names::strip_genome_suffix($f);
        die "[FATAL_ERROR] More than one genome with the same name was provided: '$b' is duplicate!. "
          . "Please make sure all names are unique and run again...\n" if $seen{$b}++;
    }

    my (@kept, @dropped, @warnings, %length, %nprot, @faa_files, @stats, %inputs);
    my $n = 0;
    for my $f (@$files) {
        my $s = 'S' . ($n + 1);
        my $type = genome_type($f);
        my ($ok, $reason, $res) = (0, '');
        if (!-f $f) { $reason = "file not found"; }
        elsif (!defined $type) { $reason = "unrecognised format (neither GenBank nor FASTA)"; }
        else {
            my $faa = "$faa_dir/$s.faa";
            eval {
                if ($type eq 'gbk') {
                    $res = gbk_to_faa($f, $s, $faa);
                    $res->{n_proteins} = $res->{n_written};
                    $res->{source} = 'gbk';
                } else {
                    my $r = run_pyrodigal($f, $s, $faa_dir, threads => $threads);
                    $res = { length => $r->{total_bp}, n_proteins => $r->{n_genes}, faa => $r->{faa}, gff => $r->{gff},
                             source => "pyrodigal-$r->{mode}", warnings => [ $r->{warning} ? ($r->{warning}) : () ] };
                }
                1;
            } or do { $reason = $@; $reason =~ s/^\[ERROR\]\s*//; $reason =~ s/\s+$//; };
            if ($res) {
                if ($res->{n_proteins} < $MIN_PROTEINS) {
                    $reason = "only $res->{n_proteins} protein(s) (< $MIN_PROTEINS)";
                    unlink $res->{faa} if $res->{faa} && -e $res->{faa};
                    $res = undef;
                } else { $ok = 1; }
            }
        }
        if (!$ok) {
            push @dropped, [$f, $reason];
            next;
        }
        $n++;
        push @kept, $f;
        $length{$s} = $res->{length};
        $nprot{$s}  = $res->{n_proteins};
        push @faa_files, $res->{faa};
        $inputs{$s} = { file => $f, type => $type, gff => $res->{gff}, faa => $res->{faa} };
        push @warnings, map { "$s ($f): $_" } @{$res->{warnings} || []};
        push @stats, join("\t", $s, PGOM_Names::strip_genome_suffix($f), $res->{source}, $res->{length}, $res->{n_proteins});
    }
    die "[FATAL_ERROR] None of the provided genomes gave usable proteins (need GenBank with CDS, or FASTA).\n" unless @kept;

    my $names = assign_short_names(\@kept);
    my $map_file = $o{map_file} || "$outdir/genome_names.map.txt";
    write_name_map($map_file, $names->{rows});

    my $stats_file = $o{stats_file} || "$outdir/genome_stats.tsv";
    open(my $so, '>', $stats_file) or die "[ERROR] Can't write '$stats_file': $!\n";
    print $so "genome\tname\tsource\tlength_bp\tn_proteins\n", map { "$_\n" } @stats;
    close $so;

    my $combined = $o{combined_faa} || "$outdir/All_proteins.faa";
    concat_files($combined, \@faa_files);

    return { order => $names->{order}, old2new => $names->{old2new}, new2old => $names->{new2old},
             genomes_length => \%length, n_proteins => \%nprot, dropped => \@dropped,
             faa_files => \@faa_files, combined_faa => $combined, map_file => $map_file,
             stats_file => $stats_file, warnings => \@warnings, inputs => \%inputs };
}

1;

package PGOM_Fast;
# Fast replacements for two steps of PanGene-O-Meter. They write exactly the same files as the legacy subroutines
# (same names, columns, order, rounding); only the internal representation differs:
#   * build_pa_fast     : phyletic pattern from DIAMOND clusters, packed vec() strings instead of hash-of-arrays/hash-of-hashes
#   * all_pairs_fast    : GCSj/GCSo for all genome pairs with packed bit vectors (popcount), optionally in parallel (fork, core Perl only)
use strict; use warnings;
use Exporter 'import';
our @EXPORT_OK = qw(build_pa_fast all_pairs_fast);

# ---------------------------------------------------------------- phyletic pattern from DIAMOND clusters
sub build_pa_fast {
    my ($clusters_file,$prefix,$del,$field)=@_;
    $del   = "\\|" if (!defined $del);
    $field = 1     if (!defined $field);
    my (%pp,%species_seen,@names_orig,@names_new,@n_species,@n_copies);
    open (my $IN,"<",$clusters_file) || return ("[ERROR] build_pa_fast: could not open '$clusters_file' $!");
    my ($cid,$prev,%cur_sp,$cur_copies)=(-1,"NA",(),0);
    my $close_cluster=sub { if ($cid>=0) {$n_species[$cid]=scalar(keys %cur_sp); $n_copies[$cid]=$cur_copies;} };
    while (my $line=<$IN>) {
        chomp ($line);
        my ($cname,$member)=split(/\t/,$line);
        if ($cname ne $prev) {
            $close_cluster->();
            $cid++; $prev=$cname;
            $names_orig[$cid]=$cname; $names_new[$cid]="DC_".sprintf("%06d",$cid);
            %cur_sp=(); $cur_copies=0;
        }
        my $sp;
        if (uc($field) eq "ALL") {$sp=$member;} else {my @t=split(/$del/,$member); $sp=$t[$field-1];}
        $species_seen{$sp}=1; $cur_sp{$sp}++; $cur_copies++;
        # 16 bit counters: up to 65535 copies of one cluster in one genome
        $pp{$sp}='' if (!defined $pp{$sp});
        my $c=vec($pp{$sp},$cid,16); vec($pp{$sp},$cid,16)=($c<65535 ? $c+1 : 65535);
    }
    close ($IN);
    $close_cluster->();
    my $nsp=scalar(keys %species_seen);
    my $f01=$prefix.".01.csv"; my $fcnt=$prefix.".csv"; my $ffas=$prefix.".01.fas";
    my $fpos=$prefix.".clusters_names_and_pos.txt"; my $fnum=$prefix.".clusters_names_and_NumOfSpecies.txt";
    open (my $OM,">",$fcnt) or return ("[ERROR] build_pa_fast: could not open '$fcnt' $!");
    open (my $O01,">",$f01) or return ("[ERROR] build_pa_fast: could not open '$f01' $!");
    open (my $OF,">",$ffas) or return ("[ERROR] build_pa_fast: could not open '$ffas' $!");
    print $OM  "Species,",join(",",@names_new),"\n";
    print $O01 "Species,",join(",",@names_new),"\n";
    foreach my $sp (sort keys %pp) {
        my (@cnt,@b01,@bfa);
        for (my $i=0;$i<=$cid;$i++) {my $c=vec($pp{$sp},$i,16); push @cnt,$c; push @b01,($c>0?1:0);}
        print $OF ">$sp\n", join("",@b01), "\n";
        print $OM  join(",",$sp,@cnt),"\n";
        print $O01 join(",",$sp,@b01),"\n";
    }
    close ($OF); close ($OM); close ($O01);
    open (my $ON,">",$fpos) or return ("[ERROR] build_pa_fast: could not open '$fpos' $!");
    print $ON "pos\tcluster_name_DIAMOND\tcluster_name_new\n";
    for (my $i=0;$i<=$cid;$i++) {print $ON $i+1,"\t$names_orig[$i]\t$names_new[$i]\n";}
    close ($ON);
    open (my $OP,">",$fnum) or return ("[ERROR] build_pa_fast: could not open '$fnum' $!");
    print $OP "pos\tcluster_name_DIAMOND\tcluster_name_new\tPrevalence\tAvg_Copy_Number\n";
    for (my $i=0;$i<=$cid;$i++) {
        my $ns=$n_species[$i]||0;
        print $OP $i+1,"\t$names_orig[$i]\t$names_new[$i]\t",sprintf("%.3f",$ns/$nsp),"\t",sprintf("%.3f",$n_copies[$i]/$ns),"\n";
    }
    close ($OP);
    print "[INFO] Out_PhyleticPattern_Clusters_Names_Pos was created: '$fpos'\n";
    print "[INFO] Out_PhyleticPattern_FASTA_01 was created: '$ffas'\n";
    print "[INFO] Out_PhyleticPattern_Matrix_01 was created: '$f01'\n";
    print "[INFO] Out_PhyleticPattern_Matrix was created: '$fcnt'\n";
    print "[INFO] Out_PhyleticPattern_Clusters_Names_and_NumOfSpecies was created: '$fnum'\n";
    return ("OK",[$fpos,$ffas,$f01,$fcnt,$fnum]);
}

# ---------------------------------------------------------------- all-vs-all GCSj / GCSo
sub _pair_row {
    my ($S1,$S2,$v1,$s1,$v2,$s2)=@_;
    my $M11=unpack("%32b*",$v1 & $v2);
    my $M10=$s1-$M11; my $M01=$s2-$M11;
    my ($dJ,$J,$O_ref,$O);
    my $S1_size=$M10+$M11; my $S2_size=$M01+$M11;
    my $smallest=($S1_size<$S2_size)?$S1_size:$S2_size;
    my $smallest_name=$S2; if ($S2_size>$S1_size) {$smallest_name=$S1;}
    my $u=$M01+$M10+$M11;
    if ($u==0) {$dJ="NA";$J="NA";} else {$dJ=sprintf("%.3f",($M01+$M10)/$u); $J=sprintf("%.3f",$M11/$u);}
    if ($S1_size==0 or $smallest==0) {$O_ref="NA";$O="NA";}
    else {$O_ref=sprintf("%.5f",$M11/$S1_size); $O=sprintf("%.5f",$M11/$smallest);}
    return ("$S1,$S2,$J,$dJ,$O_ref,$O,$S1_size,$S2_size,$smallest,$smallest_name,$M10,$M01\n",$J,$O);
}

sub all_pairs_fast {
    my ($file,$out_prefix,$delimiter,$threads,$ignore_private)=@_;
    $ignore_private = (defined $ignore_private and uc($ignore_private) eq "Y") ? 1 : 0;
    $delimiter=","  if (!defined $delimiter);
    $threads=1      if (!defined $threads or $threads<1);
    my $out_pairs=$out_prefix.".all_pairs.csv";
    open (my $IN,"<",$file) || return ("[ERROR] all_pairs_fast: can't open '$file' $!");
    my $hdr=<$IN>;
    chomp ($hdr); my @og_ids=split(/\Q$delimiter\E/,$hdr); shift @og_ids;
    my (%vec,%size,@colsum); my $npos="NA"; my $changed=0;
    while (my $line=<$IN>) {
        chomp ($line);
        my $p=index($line,$delimiter); next if ($p<0);
        my $s=substr($line,0,$p); my $data=substr($line,$p+1);
        my @f=split(/\Q$delimiter\E/,$data,-1);
        if ($npos eq "NA") {$npos=scalar(@f);}
        elsif ($npos!=@f) {print "[ERROR] calc_jaccard_and_overlap_all_pairs: The number of position in species '$s' is differnt from the rest! [".scalar(@f)."!=$npos]\n";}
        foreach (@f) { if ($_ ne "0" and $_ ne "1") {$_=($_>0)?1:0; $changed=1;} }
        $vec{$s}=pack("b*",join("",@f)); $size{$s}=unpack("%32b*",$vec{$s});
        if ($ignore_private) { for (my $k=0;$k<@f;$k++) {$colsum[$k]+=$f[$k];} }
    }
    close ($IN);
    my $npos_analyzed=$npos;
    if ($ignore_private and $npos ne "NA") {
        # same rule as the legacy code: positions present in at most one genome (sum<=1) are removed from all genomes
        my @drop=grep {($colsum[$_]||0)<=1} (0..$npos-1);
        my $fpos=$out_prefix.".filtered_pos_below_1.txt";
        if (@drop) {
            open (my $FL,">",$fpos) || return ("[ERROR] all_pairs_fast: Can't open '$fpos' $!");
            print $FL "OG_id\tOG_Pos\tSum_Presence\n";
            foreach my $k (@drop) {print $FL join("\t",(defined $og_ids[$k]?$og_ids[$k]:""),$k,($colsum[$k]||0)),"\n";}
            close ($FL);
            my $mask=pack("b*",join("",map {(($colsum[$_]||0)>1)?1:0} (0..$npos-1)));
            foreach my $s (keys %vec) {$vec{$s}=$vec{$s} & $mask; $size{$s}=unpack("%32b*",$vec{$s});}
            $npos_analyzed=$npos-scalar(@drop);
            print "[INFO] calc_jaccard_and_overlap_all_pairs: Filtered ".scalar(@drop)." positions below or equal to the 'presence' cutoff [1]\n";
        }
    }
    print "[INFO] calc_jaccard_and_overlap_all_pairs: Notice, positions with copy number higher than 1 were treated as 1\n" if ($changed);
    print "[INFO] calc_jaccard_and_overlap_all_pairs: Total positions (genes) analyzed: $npos_analyzed\n";
    my @sp=sort(keys %vec); my $n=scalar(@sp);
    my (%GJ,%GO);
    open (my $OUT,">",$out_pairs) || return ("[ERROR] all_pairs_fast: Can't open '$out_pairs' $!");
    print $OUT "strain1_ref,strain2,J,dJ,O_ref,Overlap_coeeficient,S1_size,S2_size,smallest_genome_size,smallest_genome_name,M01,M10\n";
    my $emit=sub { my ($row,$J,$O,$S1,$S2)=@_; print $OUT $row; $GJ{$S1}{$S2}=$J; $GJ{$S2}{$S1}=$J; $GO{$S1}{$S2}=$O; $GO{$S2}{$S1}=$O; };
    if ($threads==1 or $n<50) {
        for (my $i=0;$i<$n-1;$i++) { for (my $j=$i+1;$j<$n;$j++) {
            my ($row,$J,$O)=_pair_row($sp[$i],$sp[$j],$vec{$sp[$i]},$size{$sp[$i]},$vec{$sp[$j]},$size{$sp[$j]});
            $emit->($row,$J,$O,$sp[$i],$sp[$j]);
        } }
    } else {
        # fork workers; worker w computes rows i with i % W == w (balances the triangle) into its own temp file
        my $W=($threads<$n) ? $threads : $n; my @pids; my @tmp;
        for my $w (0..$W-1) {
            my $tf="$out_prefix.tmp_pairs.$w"; push @tmp,$tf;
            my $pid=fork(); return ("[ERROR] all_pairs_fast: fork failed $!") if (!defined $pid);
            if ($pid==0) {
                open (my $TF,">",$tf) or exit 1;
                for (my $i=$w;$i<$n-1;$i+=$W) { for (my $j=$i+1;$j<$n;$j++) {
                    my ($row)=_pair_row($sp[$i],$sp[$j],$vec{$sp[$i]},$size{$sp[$i]},$vec{$sp[$j]},$size{$sp[$j]});
                    print $TF "$i\t$row";
                } }
                close ($TF); exit 0;
            }
            push @pids,$pid;
        }
        my $bad=0; foreach my $pid (@pids) {waitpid($pid,0); $bad++ if ($? != 0);}
        return ("[ERROR] all_pairs_fast: a worker failed") if ($bad);
        my (@fh,@pend);
        for my $w (0..$W-1) {open ($fh[$w],"<",$tmp[$w]) or return ("[ERROR] all_pairs_fast: can't read $tmp[$w]"); $pend[$w]=readline($fh[$w]);}
        for (my $i=0;$i<$n-1;$i++) {
            my $w=$i % $W;
            while (defined $pend[$w] and $pend[$w]=~/^(\d+)\t(.*)$/s and $1==$i) {
                my $row=$2; my @c=split(/,/,$row);
                $emit->($row,$c[2],$c[5],$c[0],$c[1]);
                $pend[$w]=readline($fh[$w]);
            }
        }
        close ($_) for (@fh); unlink @tmp;
    }
    close ($OUT);
    print "[INFO] All GCS distances were written to: '$out_pairs'\n";
    return ("OK",\%GJ,\%GO,$out_pairs);
}
1;
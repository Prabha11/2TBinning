#!/bin/perl 

use warnings;
use strict;

use Term::ProgressBar;
use Data::Dumper;
use Bio::SeqIO;
use Bio::Tools::SeqWords;
use Statistics::Descriptive;

my $name 		= $ARGV[0]; 
my $output_path = $ARGV[1]; 
my $input_path 	= $ARGV[2]; 
my $seqfile_n	= $ARGV[3];
my $length_min	= $ARGV[4];
my $length_max 	= $ARGV[5];
my $format = 'fasta';
my $input = $output_path.$name.'/'.$seqfile_n; #seq file
my $binneddataset ='sample_data/simBG/simBG_output.L2_BINS'; # get the binned file 

getBinnedContigarray($binneddataset);
compute_TNF($name,$length_min,$length_max);
compute_OFDEG_GC($name,$length_min,$length_max);

my @binnedidarray;

sub getBinnedContigarray{
	my ($binneddataset) = @_;

    open(my $fh, '<:encoding(UTF-8)', $binneddataset)
    or die "Could not open file '$binneddataset' $!";
    
    while (my $row = <$fh>) { 
        my @values = split('\t', $row);
        push (@binnedidarray, $values[0]);
    }
    close ($binneddataset);
}

sub compute_OFDEG_GC {

    my ($name, $length_min, $length_max) = @_;
    my $outputfilename = 'sample_data/simBG/simBG_unbinned_contigs.OFDEG';
    
    open(OUTPUT, '>'.$outputfilename);
	print OUTPUT 'id,ns,ofdeg,r,q,rms,gc,length',"\n";
	close(OUTPUT);
        
    # --- start with all sequences again ---    
    my $all_seqs = new Bio::SeqIO(-format=>'fasta',-file=>$input);
	my $count_seqs = 0;

	while (my $seqobj = $all_seqs->next_seq) {
	$count_seqs++;
	}
	
	print "Total Seqs: $count_seqs\n";
	my $p = Term::ProgressBar->new(+{name => 'Generating unbinned data View 1', count=>$count_seqs,ETA=>'linear'});
	$p->minor(0);
	my $add = 0;
	my $update = 0;
	$all_seqs = new Bio::SeqIO(-format=>'fasta',-file=>$input);
	while(my $seqobj = $all_seqs->next_seq) {
        my $id = $seqobj->id;
        my $isBinned = 0; # boolean variable

        for (@binnedidarray){
            if($_ eq $id){
                $isBinned = 1;
                last; # stopping if found
            }
        }
        if ($isBinned == 0) {
                      
			my $num_N = ($seqobj->seq =~ tr/Nn/Nn/); #normalize
			my $seq = \(uc $seqobj->seq); #make them uppercase			

			$add++;
			$update++;
			if ($update == 10) {
				$p->update($add); # update the progress bar
				$update = 0;
			}
			my $gc = GC_content_checked($seq); 
			my ($ofdeg,$r,$q,$rms) = ofdeg2009($seq,4,0.1,0.8,5);

			open(OUTPUT, '>>'.$outputfilename);
			print OUTPUT $seqobj->id,',',$num_N,',',$ofdeg,',',$r,',',$q,',',$rms,',',$gc,',',$seqobj->length,"\n";
			close(OUTPUT);
        }
			
	}
	print "done.\n";

}

sub compute_TNF {
	my ($dataset, $name, $length_min, $length_max) = @_;
	my $outputfilename = 'sample_data/simBG/simBG_unbinned_contigs.n4';
	my $kmers;

	open(OUTPUT, '>'.$outputfilename);
	$kmers = generate_Kmers(4);
	print OUTPUT 'id,length,ns,',join(',',@$kmers),"\n";
	close(OUTPUT);

	my $all_seqs = new Bio::SeqIO(-format=>'fasta',-file=>$input);
	my $count_seqs = 0;
	while (my $seqobj = $all_seqs->next_seq) {
		$count_seqs++;
	}
	my $p = Term::ProgressBar->new(+{name => 'Generating unbinned data View 2', count=>$count_seqs,ETA=>'linear'});
	$p->minor(0);
	my $add = 0;
	my $update = 0;
	$all_seqs = new Bio::SeqIO(-format=>'fasta',-file=>$input);
	my $num_errs = 0;

	while(my $seqobj = $all_seqs->next_seq) {
		my $id = $seqobj->id;
        my $isBinned = 0; # boolean variable

        for (@binnedidarray){
            if($_ eq $id){
                $isBinned = 1;
                last; # stopping if found
            }
        }
		if ($isBinned == 0) {
			my $num_N = ($seqobj->seq =~ tr/Nn/Nn/);
			my $seq = \(uc $seqobj->seq);
			$add++;
			$update++;
			if ($update == 10) {
				$p->update($add);
				$update = 0;
			}
			my $k4 = generate_Kmers(4);
			my $tetra = freq_Kmers_gen_symm($seq,4,1);
			if ($tetra == -1) {
				$num_errs++;
			} else {
				open(OUTPUT, '>>'.$outputfilename);
				my @vals = +();
				push @vals, $tetra->{$_} foreach @$k4;
				print OUTPUT $id,',',$seqobj->length,',',$num_N,',',join(',',@vals),"\n";
				close(OUTPUT);
			}
		}
	}
	print "Errors (",$num_errs/$count_seqs,"), done.\n";

}

sub ofdeg2009 {
	my ($seq,$k,$stepsize,$alpha,$coverage) = @_;
	my @err;
	my @ls;
	my $kmers = count_Kmers_gen($seq,$k);
	my $all_kmers = generate_k_mer_nucleotides($k);
	
	my $step = ($stepsize < 1) ? length($$seq)*$stepsize : $stepsize;
	my $lim = $alpha*length($$seq);
	my $Nc = 0;
	my $Ns = 0;
	for (my $i = $step; $i < $alpha*length($$seq); $i += $step) {
		my $stat = Statistics::Descriptive::Sparse->new();
		my $N = $coverage*(length($$seq)/$i);
		$Nc += $N;
		my ($errs,$alens) = sample_seq($seq,$i,$N,$k,$kmers,$all_kmers);
		#print join (",", @$alens);

		$stat->add_data($errs);
		push @err,$stat->mean();
		$stat = Statistics::Descriptive::Sparse->new();
		$stat->add_data($alens);
		push @ls,$stat->mean();
		$Ns++;
	} 
	my $fit = Statistics::Descriptive::Full->new();
	$fit->add_data(\@err);
	my %hash;
	@hash{'q','m','r','rms'} = $fit->least_squares_fit(@ls);
	return ($hash{m},$hash{r},$hash{q},$hash{rms});	
}

sub sample_seq {
	my ($seq,$l,$N,$k,$kmers,$all_kmers) = @_;
	my $err_lookup = ();
	my @errs;
	my @actual_lengths;
	for (my $i = 0; $i < $N; $i++) {
		my $idx = int(rand(length($$seq)-$l));
		my $subseq = substr($$seq,$idx,$l);
		my $countACGT = ($subseq =~ tr/ACGTacgt/ACGTacgt/);
		push @actual_lengths,$countACGT;
		if (exists($err_lookup->{substr($$seq,$idx,$l)})) {
			push @errs, $err_lookup->{substr($$seq,$idx,$l)};
		} else {
			my $ksub = count_Kmers_gen(\substr($$seq,$idx,$l),$k);
			push @errs,  diff_euclidean2($kmers,$ksub,$all_kmers);
			$err_lookup->{substr($$seq,$idx,$l)} = $errs[$#errs];
		}
	}
	return (\@errs,\@actual_lengths);
}

sub diff_euclidean2 {
	my ($h1, $h2,$kmers) = @_;
	my $e = 0;
	for my $k (@$kmers) {
		if (exists($h1->{$k}) && exists($h2->{$k})) {
			$e += ($h1->{$k} - $h2->{$k})**2;
		} elsif (exists($h1->{$k})) {
			$e += ($h1->{$k})**2;
		} elsif (exists($h2->{$k})) {
			$e += ($h2->{$k})**2;
		}
	}
	return sqrt($e);
}

sub count_Kmers_gen {
	my ($seq, $k) = @_;
	my $km_a = generate_k_mer_nucleotides($k);	
	my $kmers;
	map { $kmers->{$_} = 0 } @$km_a;
	for (my $i = 0; $i <= length($$seq)-$k; $i++) {
		$kmers->{substr($$seq, $i, $k)}++ if exists ($kmers->{substr($$seq,$i,$k)}); #return string start index i to k
	}
	#print Dumper $kmers;
	return ($kmers);
}

sub generate_k_mer_nucleotides {
	my ($k) = @_;
	my $kmers = ();
	my @bases = ('A', 'C', 'G', 'T');
	my @words = @bases;
	for (my $i=1; $i < $k; $i++) {
		my @newwords;
		foreach my $w (@words) {
			foreach my $b (@bases) {
				push (@newwords, $w.$b);
			}
		}
		undef @words;
		@words = @newwords;
	}

	return (\@words);
}

sub freq_Kmers_gen {
	my ($seq, $k) = @_;
	my $km_a = generate_Kmers($k);
	my $kmers;
	map { $kmers->{$_} = 0 } @$km_a;
	for (my $i = 0; $i <= length($$seq)-$k; $i++) {
		$kmers->{substr($$seq, $i, $k)}++ if exists ($kmers->{substr($$seq,$i,$k)});
	}	
	my $total = 0;
	$total += $kmers->{$_} for keys %$kmers;
	$kmers->{$_} /= $total for keys %$kmers;
	return ($kmers);
}

sub freq_Kmers_gen_symm {
	my ($seq, $k, $symm) = @_;
	my $km_a = generate_Kmers($k);
	my $kmers;
	map { $kmers->{$_} = 0 } @$km_a;
	for (my $i = 0; $i <= length($$seq)-$k; $i++) {
		$kmers->{substr($$seq, $i, $k)}++ if exists ($kmers->{substr($$seq,$i,$k)});
	}
	if ($symm == 1) {
		my $rev = reverse_complement($seq);	

		for (my $i = 0; $i <= length($$rev)-$k; $i++) {
			$kmers->{substr($$rev, $i, $k)}++ if exists ($kmers->{substr($$rev,$i,$k)});
		}
	}
	my $total = 0;
	$total += $kmers->{$_} for keys %$kmers;
	$kmers->{$_} /= $total for keys %$kmers;
	return ($kmers);
}

sub count_Kmers {
	my ($seq, $k) = @_;
	my $kmers = ();
	for (my $i = 0; $i <= length($$seq)-$k; $i++) {
		if (exists ($kmers->{substr($$seq,$i,$k)})) { 
			$kmers->{substr($$seq, $i, $k)}++;
		} else {
			$kmers->{substr($$seq,$i,$k)} = 1;
		}
	}
	return ($kmers);
}

sub GC_content_checked {
	my ($seq) = @_;
	my $count = ($$seq =~ tr/GCgc/GCgc/);
	my $countACGT = ($$seq =~ tr/ACGTacgt/ACGTacgt/);
	return ($count / $countACGT);
}

sub count_words_array_strand_unbiased {
	my ($seq,$k) = @_;
	my $hash = count_words($seq,$k);
	my $rev_seq = reverse_complement($seq);
	my $hash_rc = count_words($rev_seq,$k);
	my $kmers = generate_Kmers($k);
	my @array;
	my $total = 0;
	foreach my $mer (@$kmers) {
		if (exists($hash->{$mer})) {
			if (exists($hash_rc->{$mer})) {
				push @array, ($hash->{$mer}+$hash_rc->{$mer});
				$total += $hash->{$mer}+$hash_rc->{$mer};
			} else {
				push @array, $hash->{$mer};
				$total += $hash->{$mer};
			}
		} else {
			if (exists($hash_rc->{$mer})) {
				push @array, $hash_rc->{$mer};
				$total += $hash_rc->{$mer};
			} else {
				push @array, 0;
			}
		}
	}
	foreach (@array) {
		#$_ = $_ / $total;
		#$sum += $_;
	}
	return (\@array);
}

sub count_words_array_strand_biased {
	my ($seq,$k) = @_;
	my $hash = count_words($seq,$k);
	my $kmers = generate_Kmers($k);
	my @array;
	my $total = 0;
	foreach my $mer (@$kmers) {
		if (exists($hash->{$mer})) {
			push @array, $hash->{$mer};
			$total += $hash->{$mer};
			
		} else {
			push @array, 0;
		}
	}
	foreach (@array) {
		$_ = $_ / $total;
	}
	return (\@array);
}

sub reverse_complement {
	my ($dna) = @_;
	my $revcom = reverse $$dna;
	$revcom =~ tr/ACGTacgt/TGCAtgca/;
	return \$revcom;
}

sub count_words {
	my ($seq,$k) = @_;
	my $s = Bio::Seq->new(-seq=>$$seq,-alphabet=>'dna');
	my $c = Bio::Tools::SeqWords->new(-seq => $s);
	return ($c->count_overlap_words($k));
}

sub generate_Kmers {
	my ($k) = @_;
	my $kmers = ();
	my @bases = ('A', 'C', 'G', 'T');
	my @words = @bases;
	
	for (my $i=1; $i < $k; $i++) {
		my @newwords;
		foreach my $w (@words) {
			foreach my $b (@bases) {
				push (@newwords, $w.$b);
			}
		}
		undef @words;
		@words = @newwords;
	}
	return (\@words);
}

sub nf_markov_norm {
	
	my ($seq,$k,$m,$symm) = @_;
	my $ref = freq_Kmers_gen_symm($seq, $k, $symm);
	my $norm;
	$norm->{$m-1} = freq_Kmers_gen_symm($seq, $m-1, $symm);
	$norm->{$m} = freq_Kmers_gen_symm($seq, $m, $symm);
	
	my $normalised;
	
	for my $ref_n (keys %$ref) {
		my @nucs = split('',$ref_n);
		my @f1;
		my @f2;
		my @f3;
		push @f1,$nucs[$_] for 1..2;
		push @f2,$nucs[$_] for 0..2;
		push @f3,$nucs[$_] for 1..3;
		my $h1 = join('',@f1);
		my $h2 = join('',@f2);
		my $h3 = join('',@f3);
		
		if (($norm->{$m}->{$h2} * $norm->{$m}->{$h3}) == 0) {
			return -1;
		}
		my $norm_factor = $norm->{$m-1}->{$h1} / ($norm->{$m}->{$h2} * $norm->{$m}->{$h3});
		$normalised->{$ref_n} = $ref->{$ref_n} * $norm_factor;
	}
	return $normalised;
}

sub tetra_Teeling2004 {
	my ($seq) = @_;
	
	my $te = freq_Kmers_gen($seq,4);
	my $tr = freq_Kmers_gen($seq,3);
	my $di = freq_Kmers_gen($seq,2);
	
	my $z;
	
	for my $tetra (keys %$te) {
		if ($te->{$tetra} > 0) {
			my ($n1) = ($tetra =~ /([ACGT][ACGT][ACGT])[ACGT]/);
			my ($n2) = ($tetra =~ /[ACGT]([ACGT][ACGT][ACGT])/);
			my ($n3) = ($tetra =~ /[ACGT]([ACGT][ACGT])[ACGT]/);
			my $e = $tr->{$n1} * $tr->{$n2} / $di->{$n3};
			my $v = $e * (  ($di->{$n3} - $tr->{$n1}) * ($di->{$n3} - $tr->{$n2})  )/($di->{$n3}**2);
			$z->{$tetra} = ($te->{$tetra} - $e) / sqrt($v);
		} else {
			$z->{$tetra} = 0;
		}
	}
	return $z;

}

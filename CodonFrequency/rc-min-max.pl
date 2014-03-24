#!/usr/bin/perl
#Programmer: Rory Carmichael
#Purpose: Refactor the insanity of the original %min-max perl script

if ($#ARGV != 1 ) {
	print STDERR "USAGE: better-min-max.pl <qry_fasta_orfs> <codon_usage_table>\n";
	exit 1;
}

my $windowsize = 17;
my $curid;
my %seqs = ();
my %codon_freqs= ();
my %min = ();
my %max = ();
my %avg = ();
my $qry_file = $ARGV[0];
my $codon_file = $ARGV[1];

my %DNAtoAA = ('GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TGT' => 'C',
	       'TGC' => 'C', 'GAT' => 'D', 'GAC' => 'D', 'GAA' => 'E', 'GAG' => 'E',
	       'TTT' => 'F', 'TTC' => 'F', 'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G',
	       'GGG' => 'G', 'CAT' => 'H', 'CAC' => 'H', 'ATT' => 'I', 'ATC' => 'I',
	       'ATA' => 'I', 'AAA' => 'K', 'AAG' => 'K', 'TTG' => 'L', 'TTA' => 'L',
	       'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L', 'ATG' => 'M',
	       'AAT' => 'N', 'AAC' => 'N', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P',
	       'CCG' => 'P', 'CAA' => 'Q', 'CAG' => 'Q', 'CGT' => 'R', 'CGC' => 'R',
	       'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R', 'AGG' => 'R', 'TCT' => 'S',
	       'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S',
	       'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T', 'GTT' => 'V',
	       'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V', 'TGG' => 'W', 'TAT' => 'Y',
	       'TAC' => 'Y', 'TAA' => '*', 'TAG' => '*', 'TGA' => '*');

#HASH containing the number of synonymous codons for each amino acid, calculated because I am too lazy to enumerate them manually and that is why we invented computers
my %syncounts = ();
foreach my $key (keys(%DNAtoAA)) {
	$syncounts{$DNAtoAA{$key}}++;
}	

sub read_codon_table {
	my $tfile = shift;
#	print STDERR "Reading codon table from $tfile\n";
	open(IF, $tfile);
	foreach my $line (<IF>) {
		my @triplets = split(/\t+/,$line);
		foreach my $trip (@triplets) {
			my $codon = [split(/\s+/,$trip)]->[0];
			$freq = [ split(/\(/,$trip)]->[-1];
			$freq =~ s/\)//;
#			print "FREQ is $freq\n";
			$codon_freqs{$codon} = $freq;
		}
	}
	close(IF);
}

sub read_seqs {
	my $sfile = shift;
#	print STDERR "Reading sequences from $sfile\n";
	open(IF, $sfile);
	my $curid;
	foreach my $line (<IF>) {
		chomp $line;
		if ($line =~ /^>/) {
			$curid = $line;
		} else {
			if (defined($curid)) {
				$line = "\U$line";
				$seqs{$curid} .= $line;
			} else {
#				print STDERR "Sequence without id\n";
			}
		}
	}
	close(IF);
	while (my ($key,$val) = each(%seqs) ) {
		if (length($val) % 3  != 0 ) {
#			print STDERR "Removing bad sequence $key of length " . length($val) . "\n";
			undef($seqs{$key});
		}
	}
}
		
#print "Reading input\n";
&read_codon_table($codon_file);
&read_seqs($qry_file);
#print "Running Min-Max\n";

#Print the synonym counts for debugging purposes
#while (my ($key,$val) = each (%syncounts) ) {
#	print "$key\t$val\n";
#}

#print "Getting min max and average for each codon\n";
while (my ($key, $val) = each(%codon_freqs)) {
	my $AA = $DNAtoAA{$key};
#	print STDERR "$key $val\n";
	#Find out if this codon is a min for its amino acid
	if(!defined($min{$AA}) || $min{$AA} > $val) {
#		print STDERR "found new min\n";
		$min{$AA} = $val;
	}
	#Find out if this codon is a max for its amino acid
	if(!defined($max{$AA}) || $max{$AA} < $val) {
#		print STDERR "found new max\n";
		$max{$AA} = $val;
	}
	#Add this codon's contribution to the average
	if($syncounts{$AA} > 0) {
#		print STDERR "got new avg\n";
		$avg{$AA} += $val/$syncounts{$AA};
	} else {
#		print STDERR "We have no synonyms for $AA\n";
	}
}
#print "MINS:\n";
#foreach my $AA (keys(%syncounts)) {
#	print $AA . "\t" . $min{$AA} . "\n";
#}
#print "MAXS:\n";
#foreach my $AA (keys(%syncounts)) {
#	print $AA . "\t" . $max{$AA} . "\n";
#}
#print "AVGS:\n";
#foreach my $AA (keys(%syncounts)) {
#	print $AA . "\t" . $avg{$AA} . "\n";
#}
#
#print "Running Min-Max for Realz\n";
#now we have all the numbers we'll need to do the in window calculations... time to get our windows
while (my ($id,$seq) = each(%seqs)) {
	print $id . "\n";
	my $numcodons = length($seq)/3;
	for(my $i=0;$i<$numcodons;$i++) {
		if ($i*3+$windowsize*3 > length($seq)) {
			last;
		}
		my $window = substr($seq,$i*3,$windowsize*3);
		my $winval = 0;
		my $actfreq = 0;
		my $maxfreq = 0;
		my $minfreq = 0;
		my $avgfreq = 0;
		for (my $j=0;$j<$windowsize*3;$j+=3) {
			my $codon = substr($window,$j,3);
			my $AA = $DNAtoAA{$codon};
			$actfreq += $codon_freqs{$codon};
			$maxfreq += $max{$AA};
			$minfreq += $min{$AA};
			$avgfreq += $avg{$AA};
		}
		$actfreq /= $windowsize;
		$maxfreq /= $windowsize;
		$minfreq /= $windowsize;
		$avgfreq /= $windowsize;
		if ($actfreq == $avgfreq || $maxfreq == $avgfreq || $minfreq == $avgfreq) {
			#Just to be on the safe side, we watch out for naughty equality cases
			$winval = 0;
		} elsif($actfreq > $avgfreq) {
			#calculate %max
			$winval = ($actfreq-$avgfreq)/($maxfreq-$avgfreq)*100;
		} else {
			#calculate %min
			$winval = ($avgfreq-$actfreq)/($avgfreq-$minfreq)*-100;
		}
		print sprintf("%.2f", $winval) . ",";
	}
	print "\n";
}

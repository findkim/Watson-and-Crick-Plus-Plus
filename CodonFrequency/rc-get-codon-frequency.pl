#!/usr/bin/perl
#Programmer: Rory Carmichael
#Purpose: calculate codon frequencies of an orfeome for later use in min-max.pl

my $windowsize = 17;
my $curid;
my %seqs = ();
my %codon_freqs= ();
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

sub count_line_codons {
	my $curid = shift;
	my $numcodons = 0;
	if (length($seqs{$curid}) % 3 == 0) {
		for(my $i=0;$i<length($seqs{$curid});$i+=3) {
			my $codon = substr($seqs{$curid},$i,3);
			if (defined($DNAtoAA{$codon})) {
				$codon_freqs{$codon}++;
				$numcodons++;
			} else {
				print STDERR "Saw bad codon: $codon\n";
			}
		}
	} else {
		print STDERR "Bad seq, was length: " . length($seqs{$curid}) . "\n";
		undef($seqs{$curid});
	}
	return $numcodons;
}

my %min = ();
my %max = ();
my %avg = ();
my $numcodons = 0;
while ( my $line = <> ) {
	chomp $line;
	if ($line =~ /^>/) {
		if (defined($curid)) {
			$numcodons += &count_line_codons($curid);
		}
		$curid = $line;
	} else {
		if (defined($curid)) {
			$line = "\U$line";
			$seqs{$curid} .= $line;
		} else {
			print STDERR "Sequence without id\n";
		}
	}
}
$numcodons += &count_line_codons($curid);

#print the codon frequencies in the stupid and awful format that the japanese people made
my $i = 0;
foreach my $key (keys(%codon_freqs)) {
	print "$key " . $codon_freqs{$key}/$numcodons*1000 . "($codon_freqs{$key})"; 
	$i++;
	if($i == 4) { $i=0; print "\n"; } else { print "\t"; }
}
print "\n";

/*

	Created by Kim Ngo 14/3/23
	
	CodonFrequency.cpp
	
	Calculates frequency for each codon in the gene
	Create library of codon frequencies
	
*/

#include <iostream>
#include "Sequence.h"
#include "CodonFrequency.h"
#include <vector>
#include <iterator>
#include <string>
#include <fstream> //output file

using namespace std;

const int NUM_TYPE_OF_CODONS = 64;

CodonFrequency::CodonFrequency(vector<Sequence> seq) {

	//vector< pair<string,int> > temp = {
	codon = {
	/*
				 'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TGT' => 'C',
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
*/
		// Codon A
		{"GCT",0}, {"GCC",0}, {"GCA",0}, {"GCG",0},
	
		// Codon C
		{"TGT",0}, {"TGC",0},
	
		// Codon D
		{"GAT",0}, {"GAC",0},
	
		// Codon E
		{"GAA",0}, {"GAG",0},

		// Codon F
		{"TTT",0}, {"TTC",0},
		
		// Codon G
		{"GGT",0}, {"GGC",0}, {"GGA",0}, {"GGG",0},
		
		// Codon H
		{"CAT",0}, {"CAC",0},
		
		// Codon I
		{"ATT",0}, {"ATC",0}, {"ATA",0},
		
		// Codon K
		{"AAA",0}, {"AAG",0},
		
		// Codon L
		{"TTG",0}, {"TTA",0}, {"CTT",0}, {"CTC",0}, {"CTA",0}, {"CTG",0},
				
		// Codon M
		{"ATG",0},
		
		// Codon N
		{"AAT",0}, {"AAC",0},
		
		// Codon P
		{"CCT",0}, {"CCC",0}, {"CCA",0}, {"CCG",0},
		
		// Codon Q
		{"CAA",0}, {"CAG",0},

		// Codon R
		{"CGT",0}, {"CGC",0}, {"CGA",0}, {"CGG",0}, {"AGA",0}, {"AGG",0},
		
		// Codon S
		{"TCT",0}, {"TCC",0}, {"TCA",0}, {"TCG",0}, {"AGT",0}, {"AGC",0},
		
		// Codon T
		{"ACT",0}, {"ACC",0}, {"ACA",0}, {"ACG",0},
		
		// Codon V
		{"GTT",0}, {"GTC",0}, {"GTA",0}, {"GTG",0},
		
		// Codon W
		{"TGG",0},
		
		// Codon Y
		{"TAT",0}, {"TAC",0},
		
		// Stop Codon
		{"TAA",0}, {"TAG",0}, {"TGA",0},
	};

	incrOccurance(seq);
	
	for (int i = 0; i < seq.size(); i++) {	// vector of sequences
//		cout << "THE NUMBER OF CODON FOR THIS SEQUENCE = " << seq[i].getNumCodon() << endl;
		for (int j = 0; j < seq[i].getNumCodon(); j++) {
			codonFreqSeq.push_back(0.0);
		}
	}

	calcFreq(seq);

}


// Increments codon occurance for every triplet in the sequence
// Increments the number of codons in sequence with setter method
void CodonFrequency::incrOccurance(vector<Sequence> seq) {

	int count = 0;
	for (int k = 0; k < seq.size(); k++) {
		for (int i = 1; i < seq[k].getSeqLength(); i+=3) {
		// i starts at 1 to properly utilize i+=3: 1, 4, 7 ... 
		// instead of i starting at 0: 0, 3, 6 -- which skips position 3
	
			string triplet; // set of 3 characters from sequence
			triplet.append(seq[k].getSeq(),i-1,3);	// sets string to triplet of char at position i
		
//		cout << triplet << " ";perator<<

			// Increments # of occurances to corresponding string/codon
			for (int j = 0; j < codon.size(); j++) {
				if (triplet == codon[j].first) {
					codon[j].second++;
					count++;
//					cout << "The current count is " << count << "\tfor " << codon[j].first << endl;
				}
			}
		}
	}
	set_codonCount(count);	// Increments with setter; codonCount is private data member

}


// Setter for codonCount -- used for incrementing # of codons in sequence in incrOcc method
void CodonFrequency::set_codonCount(int c) {
	codonCount = c;
}


// Getter for codonCount -- used for calculating freq
int CodonFrequency::getCodonCount() {
	return (codonCount);
}


// Calculates freq -- #ofOcc/codonCount
void CodonFrequency::calcFreq(vector <Sequence> seq) {

for (int k = 0; k < seq.size(); k++) {
	for (int i = 1, count = 0; i < seq[k].getSeqLength(); i+=3, count++) {
	// i starts at 1 to properly utilize i+=3: 1, 4, 7 ... 
	// instead of i starting at 0: 0, 3, 6 -- which skips position 3
	
		string triplet; // set of 3 characters from sequence
		triplet.append(seq[k].getSeq(),i-1,3);	// sets string to triplet of char at position i

		// For every triplet in seq, it matches with one codon
		// Calculates freq for that codon, then sets into codonFreqSeq vector
		for (int j = 0; j < codon.size(); j++)
			if (triplet == codon[j].first) {
		
				codonFreqSeq[count] = (float) codon[j].second/getCodonCount()*1000;

//				cout << triplet << " " << codon[j].second << "/" << getCodonCount() << " = " << codon[j].second/getCodonCount() << endl;
			}
	}
}
}


// Prints # of occurances for each codon for the vector of sequences
void CodonFrequency::printCodonCount() {

	for (int i = 0, count = 0; i < codon.size(); i++, ++count) {
		if (count%4 == 0) cout << endl;
		cout << codon[i].first << " " << (float) codon[i].second/getCodonCount()*1000 << "(" << codon[i].second << ")  ";
	}
//	cout << endl;
//	cout << "The total number of codons is " << getCodonCount() << endl;
}

// Prints frequency for the sequence
void CodonFrequency::printFreq() {

	for (int i = 0; i < getCodonCount(); i++) {
		cout << codonFreqSeq[i] << " ";
	}

}

// Stores frequency and # of occurances for each codon for the vector of sequences in an output file
void CodonFrequency::outputFileCodonCount(ofstream &ofilename) {

	for (int i = 0, count = 0; i < codon.size(); i++, ++count) {
		if (count%4 == 0 && count != 0)
			ofilename << endl;
		ofilename << codon[i].first;
		ofilename << " ";
		ofilename << (float) codon[i].second/getCodonCount()*1000;
		ofilename << "(";
		ofilename << codon[i].second;
		ofilename << ")  ";
	}
}

// Stores each sequence as a series of frequencies in an output file
void CodonFrequency::outputfileFreq(ofstream &ofilename) {

	for (int i = 0; i < getCodonCount(); i++) {
		ofilename << codonFreqSeq[i] << " ";
	}

}

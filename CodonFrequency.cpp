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
#include <unordered_map>	// stringmap for unordered_multimap
#include <utility> // pair
#include <bitset>

using namespace std;

const int NUM_TYPE_OF_CODONS = 64;

CodonFrequency::CodonFrequency(vector<Sequence> seq) {


	cout << "Calculating codon frequency..." << endl;
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
		// Codon A Ala
		{"GCT",0}, {"GCC",0}, {"GCA",0}, {"GCG",0},
	
		// Codon C Tyr
		{"TGT",0}, {"TGC",0},
	
		// Codon D Asp
		{"GAT",0}, {"GAC",0},
	
		// Codon E Glu
		{"GAA",0}, {"GAG",0},

		// Codon F Phe
		{"TTT",0}, {"TTC",0},
		
		// Codon G Arg
		{"GGT",0}, {"GGC",0}, {"GGA",0}, {"GGG",0},
		
		// Codon H His
		{"CAT",0}, {"CAC",0},
		
		// Codon I Ile
		{"ATT",0}, {"ATC",0}, {"ATA",0},
		
		// Codon K Lys
		{"AAA",0}, {"AAG",0},
		
		// Codon L Leu
		{"TTG",0}, {"TTA",0}, {"CTT",0}, {"CTC",0}, {"CTA",0}, {"CTG",0},
				
		// Codon M Met
		{"ATG",0},
		
		// Codon N Asn
		{"AAT",0}, {"AAC",0},
		
		// Codon P Pro
		{"CCT",0}, {"CCC",0}, {"CCA",0}, {"CCG",0},
		
		// Codon Q Gln
		{"CAA",0}, {"CAG",0},

		// Codon R Arg
		{"CGT",0}, {"CGC",0}, {"CGA",0}, {"CGG",0}, {"AGA",0}, {"AGG",0},
		
		// Codon S Ser
		{"TCT",0}, {"TCC",0}, {"TCA",0}, {"TCG",0}, {"AGT",0}, {"AGC",0},
		
		// Codon T Thr
		{"ACT",0}, {"ACC",0}, {"ACA",0}, {"ACG",0},
		
		// Codon V Val
		{"GTT",0}, {"GTC",0}, {"GTA",0}, {"GTG",0},
		
		// Codon W Trp
		{"TGG",0},
		
		// Codon Y Try
		{"TAT",0}, {"TAC",0},
		
		// Stop Codon
		{"TAA",0}, {"TAG",0}, {"TGA",0}
	};

	incrOccurance(seq);
	
	for (int i = 0; i < seq.size(); i++) {	// vector of sequences
//		cout << "THE NUMBER OF CODON FOR THIS SEQUENCE = " << seq[i].getNumCodon() << endl;
		for (int j = 0; j < seq[i].getNumCodon(); j++) {
			codonFreqSeq.push_back(0.0);
		}
	}

//	calcFreq(seq);
	AAtoCodonMap = createMap(codonFreq);
	cout << "Finished calculations." << endl; 
}



// Increments codon occurance for every triplet in the sequence
// Increments the number of codons in sequence with setter method
void CodonFrequency::incrOccurance(vector<Sequence> seq) {

	int count = 0;
	for (int k = 0; k < seq.size(); k++) {
	
			string tempStr = seq[k].getSeq();
//			cout << tempStr << endl;

			for (string::iterator itr = tempStr.begin(), 
			end = tempStr.end(); itr != end; ++itr) {
				
				// A 00 -- G 01 -- C 10 -- T 11
				int binaryRep = 0;
				// First nucleotide
				//if (*itr == 'A') binaryRep += 0;					// 00xxxx
				if (*itr == 'G') binaryRep += 16;						// 01xxxx
				else if (*itr == 'C') binaryRep += 32;			// 10xxxx
				else if (*itr == 'T') binaryRep += (32+16);	// 11xxxx
				
				// Second nucleotide
				itr++;
				//if (*itr == 'A') binaryRep += 0;					// xx00xx
				if (*itr == 'G') binaryRep += 4;						// xx01xx
				else if (*itr == 'C') binaryRep += 8;				// xx10xx
				else if (*itr == 'T') binaryRep += (8+4);		// xx11xx
				
				// Third nucleotide
				itr++;
				//if (*itr == 'A') binaryRep += 0;					// xxxx00
				if (*itr == 'G') binaryRep += 1;						// xxxx01
				else if (*itr == 'C') binaryRep += 2;				// xxxx10
				else if (*itr == 'T') binaryRep += (2+1);		// xxxx11
				
//				cout << "The binary rep is " << binaryRep << endl;
				codonOcc[binaryRep]++;
				count++;
			}
			
/*	
		for (int i = 1; i < seq[k].getSeqLength(); i+=3) {
		// i starts at 1 to properly utilize i+=3: 1, 4, 7 ... 
		// instead of i starting at 0: 0, 3, 6 -- which skips position 3
	
			string triplet; // set of 3 characters from sequence
			triplet.append(seq[k].getSeq(),i-1,3);	// sets string to triplet of char at position i
			//		cout << triplet << " ";
			
			// Increments # of occurances to corresponding string/codon
			for (int j = 0; j < codon.size(); j++) {
				if (triplet == codon[j].first) {
					codon[j].second++;
					count++;
//					cout << "The current count is " << count << "\tfor " << codon[j].first << endl;
				}
			}
		}
*/	}
		set_codonCount(count);	// Increments with setter; codonCount is private data member

		for (int k = 0; k < NUM_TYPE_OF_CODONS; k++) {
			codonFreq[k] = (float) codonOcc[k] / count * 1000;
		}
}


unordered_multimap<string, pair<int, float> > CodonFrequency::createMap(float codonFreq[]) {

//	bitset<3> bits(string("010"));
//	cout << bits.to_ulong() << endl;
	
	unordered_multimap<string, pair<int, float> > AAtoCodonMap;
	
	for (int i = 0; i < NUM_TYPE_OF_CODONS; i++) {

		pair <int, float> temp (i, codonFreq[i]);
		cout << temp.first << " " << temp.second << endl;
		
			if (i <= 24 && i >= 27) {							// A
				pair <string, pair<int, float> > temp2 ("A", temp);
				AAtoCodonMap.insert (temp2);
			}
				//AAtoCodonMap["A"] = temp;
			else if (i <= 54 && i >= 55) {				// C
				pair <string, pair<int, float> > temp2 ("C", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 18 && i >= 19) {				// D
				pair <string, pair<int, float> > temp2 ("D", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 16 && i >= 17) {				// E
				pair <string, pair<int, float> > temp2 ("E", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 62 && i >= 63) {				// F
				pair <string, pair<int, float> > temp2 ("F", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 20 && i >= 23) {				// G
				pair <string, pair<int, float> > temp2 ("G", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 34 && i >= 35) {				// H
				pair <string, pair<int, float> > temp2 ("H", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 14 && i >= 15 || i == 12) {				// I
				pair <string, pair<int, float> > temp2 ("I", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 0 && i >= 1) {					// K
				pair <string, pair<int, float> > temp2 ("K", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 60 && i >= 61 || i <= 44 && i >= 47) {				// L
				pair <string, pair<int, float> > temp2 ("L", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i == 13) {										// M
				pair <string, pair<int, float> > temp2 ("M", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 2 && i >= 3) {				// N
				pair <string, pair<int, float> > temp2 ("N", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 40 && i >= 43) {				// P
				pair <string, pair<int, float> > temp2 ("P", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 31 && i >= 32) {				// Q
				pair <string, pair<int, float> > temp2 ("Q", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 36 && i >= 39 || i <= 4 && i >= 5) {				// R
				pair <string, pair<int, float> > temp2 ("R", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 56 && i >= 59 || i <= 6 && i >= 7) {				// S
				pair <string, pair<int, float> > temp2 ("S", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 8 && i >= 11) {					// T
				pair <string, pair<int, float> > temp2 ("T", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 28 && i >= 31) {				// V
				pair <string, pair<int, float> > temp2 ("V", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i == 53) {				// W
				pair <string, pair<int, float> > temp2 ("W", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 50 && i >= 51) {				// Y
				pair <string, pair<int, float> > temp2 ("Y", temp);
				AAtoCodonMap.insert (temp2);
			}
			else if (i <= 48 && i >= 49 || i == 52) {				// Stop
				pair <string, pair<int, float> > temp2 ("*", temp);
				AAtoCodonMap.insert (temp2);
			}
	}
	cout << "AAtoCodonMap contains: " << endl;
//	if (AAtoCodonMap.find("A"))
	cout << AAtoCodonMap.find("*")->first.size() << endl;
//	for (auto x: AAtoCodonMap)
	//for (int j = 0; j < AAtoCodonMap.find("*")->first.size(); j++)
	cout << AAtoCodonMap.find("*")->first << " " << AAtoCodonMap.find("*")->second.first << " " << AAtoCodonMap.find("*")->second.second << endl;
	cout << "In here" << endl;
	return AAtoCodonMap;
}


// Setter for codonCount -- used for incrementing # of codons in sequence in incrOcc method
void CodonFrequency::set_codonCount(int c) {
	codonCount = c;
}


// Getter for codonCount -- used for calculating freq
int CodonFrequency::getCodonCount() {
	return (codonCount);
}

/*
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
*/

// Prints # of occurances for each codon for the vector of sequences
void CodonFrequency::printCodonCount() {

	for (int i = 0, count = 0; i < codon.size(); i++, ++count) {
		if (count%4 == 0) cout << endl;
		cout << codon[i].first << " " << (float) codon[i].second/getCodonCount()*1000 << "(" << codon[i].second << ")  ";
	}
//	cout << endl;
//	cout << "The total number of codons is " << getCodonCount() << endl;
}

/*
// Prints frequency for the sequence
void CodonFrequency::printFreq() {

	for (int i = 0; i < getCodonCount(); i++) {
		cout << codonFreqSeq[i] << " ";
	}

}

void CodonFrequency::binary(int decimal) {
   int remainder;

   if(decimal <= 1) {
       std::cout << decimal;
       return;
   }
   remainder = decimal % 2;
   binary(decimal >> 1);    
   std::cout << remainder;
}
*/

// Stores frequency and # of occurances for each codon for the vector of sequences in an output file
void CodonFrequency::outputFileCodonCount(ofstream &ofilename) {

/*
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
*/
	for (int i = 0, count = 0; i < NUM_TYPE_OF_CODONS; i++, count++) {
		ofilename << codonFreq[i];
		ofilename << " (";
		ofilename << codonOcc[i];
		ofilename << ")";
		ofilename << endl;
	}
}

/*
// Stores each sequence as a series of frequencies in an output file
void CodonFrequency::outputfileFreq(ofstream &ofilename) {

	for (int i = 0; i < getCodonCount(); i++) {
		ofilename << codonFreqSeq[i] << " ";
	}

}
*/

/*

	Created by Kim Ngo 14/3/23
	
	CodonFrequency.cpp
	
	Calculates frequency for each codon in the gene
	Create library of codon frequencies
	Creates a library of maximum frequencies for each aa
	Creates a library of minimum frequencies for each aa
	
*/

#include <iostream>
#include "Sequence.h"
#include "CodonFrequency.h"
#include <vector>
#include <iterator>
#include <string>
#include <fstream> // output file
#include <map>	// multimap
#include <utility> // pair
#include <algorithm>
#include <deque>

using namespace std;

const int NUM_TYPE_OF_CODONS = 64;

CodonFrequency::CodonFrequency(char *filename, vector<Sequence> seq) {

	// initializes codonOcc array to 0
	for (int i = 0; i < NUM_TYPE_OF_CODONS; i++) {
		codonOcc[i] = 0;
	}

	cout << "Calculating codon frequency..." << endl;
/*
				 'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 
				 'TGT' => 'C', 'TGC' => 'C',
				 'GAT' => 'D', 'GAC' => 'D',
				 'GAA' => 'E', 'GAG' => 'E',
	       'TTT' => 'F', 'TTC' => 'F',
	       'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G', 
	       'CAT' => 'H', 'CAC' => 'H', 
	       'ATT' => 'I', 'ATC' => 'I', 'ATA' => 'I',
	       'AAA' => 'K', 'AAG' => 'K',
	       'TTG' => 'L', 'TTA' => 'L', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L', 
	       'ATG' => 'M',
	       'AAT' => 'N', 'AAC' => 'N',
	       'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P', 'CCG' => 'P',
	       'CAA' => 'Q', 'CAG' => 'Q',
	       'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R', 'AGG' => 'R', 
	       'TCT' => 'S', 'TCC' => 'S', 'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S',
	       'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
	       'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V', 'GTG' => 'V',
	       'TGG' => 'W',
	       'TAT' => 'Y', 'TAC' => 'Y',
	       'TAA' => '*', 'TAG' => '*', 'TGA' => '*');
*/
	// string of all codons
	string codonStr = "AAAAAGAACAATAGAAGGAGCAGTACAACGACCACTATAATGATCATTGAAGAGGACGATGGAGGGGGCGGTGCAGCGGCCGCTGTAGTGGTCGTTCAACAGCACCATCGACGGCGCCGTCCACCGCCCCCTCTACTGCTCCTTTAATAGTACTATTGATGGTGCTGTTCATCGTCCTCTTTATTGTTCTTT";
	codonToAAMap = createCodonToAAMap(codonStr, codonMap);

	vector<Sequence> Sequences = removeSeq(seq);
	calcFreq(Sequences);

	AAtoCodonMap = createMap(codonFreq);
	printMap(AAtoCodonMap);
	minMap = createMinMap(AAtoCodonMap);
	maxMap = createMaxMap(AAtoCodonMap);
	avgMap = createAvgMap(AAtoCodonMap);
	cout << "Finished calculations." << endl; 
	
	outputFileCF(filename);
}


// Does not initialize codonFreq for alignments
// Removes sequences that are divisible by 3
vector<Sequence> CodonFrequency::removeSeq(vector<Sequence> &Sequences){

	for (int i = 0; i < Sequences.size(); i++) {
  	if (Sequences[i].getSeqLength()%3 != 0) {
   		Sequences.erase (Sequences.begin()+i);
   		cout << Sequences[i].getSeqName() << " " << Sequences[i].getSeqDescription() << " has been removed from codon frequency calculations due to improper length." << endl;
 		}
 	}
 	return Sequences;
}


// Increments codon occurance for every triplet in the sequence
// Increments the number of codons in sequence with setter method
// Calculates codon frequency
void CodonFrequency::calcFreq(vector<Sequence> &seq) {

	int count = 0;
	
	for (int k = 0; k < seq.size(); k++) {

			string tempStr = seq[k].getSeq();
//			cout << tempStr << endl;

			for (string::iterator itr = tempStr.begin(); itr != tempStr.end(); ++itr) {
				
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

				codonOcc[binaryRep] = codonOcc[binaryRep] + 1;
//				cout << "The binary rep is " << binaryRep << "\tThe num of occ = " << codonOcc[binaryRep] << endl;
				count++;
			}
	}
	cout << endl;
	set_codonCount(count);	// Increments with setter; codonCount is private data member
	cout << "codon count = " << getCodonCount() << endl;
	for (int k = 0; k < NUM_TYPE_OF_CODONS; k++) {
		codonFreq[k] = (float) codonOcc[k] / count * 1000;
//		cout << k << "	" << codonFreq[k] << " = " << codonOcc[k] << " / " << count << endl;
	}
}


// Maps all codons with frequency to corresponding amino acid
multimap<char, pair<int, float> > CodonFrequency::createMap(float codonFreq[64]) {

	multimap<char, pair<int, float> > AAtoCodonMap;
	char AA;
	
	for (int i = 0; i < NUM_TYPE_OF_CODONS; i++) {

		pair <int, float> temp (i, codonFreq[i]);
//		cout << temp.first << " " << temp.second << endl;
			if (i >= 24 && i <= 27) {							// A
				AA = 'A';
			}
			else if (i >= 54 && i <= 55) {				// C
				AA = 'C';
			}
			else if (i >= 18 && i <= 19) {				// D
				AA = 'D';
			}
			else if (i >= 16 && i <= 17) {				// E
				AA = 'E';
			}
			else if (i >= 62 && i <= 63) {				// F
				AA = 'F';
			}
			else if (i >= 20 && i <= 23) {				// G
				AA = 'G';
			}
			else if (i >= 34 && i <= 35) {				// H
				AA = 'H';
			}
			else if (i == 12 || i == 14 || i == 15) {				// I
				AA = 'I';
			}
			else if (i >= 0 && i <= 1) {					// K
				AA = 'K';
			}
			else if (i >= 60 && i <= 61 || i >= 44 && i <= 47) {				// L
				AA = 'L';
			}
			else if (i == 13) {										// M
				AA = 'M';
			}
			else if (i >= 2 && i <= 3) {				// N
				AA = 'N';
			}
			else if (i >= 40 && i <= 43) {				// P
				AA = 'P';
			}
			else if (i >= 32 && i <= 33) {				// Q
				AA = 'Q';
			}
			else if (i >= 36 && i <= 39 || i >= 4 && i <= 5) {				// R
				AA = 'R';
			}
			else if (i >= 56 && i <= 59 || i >= 6 && i <= 7) {				// S
				AA = 'S';
			}
			else if (i >= 8 && i <= 11) {					// T
				AA = 'T';
			}
			else if (i >= 28 && i <= 31) {				// V
				AA = 'V';
			}
			else if (i == 53) {				// W
				AA = 'W';
			}
			else if (i >= 50 && i <= 51) {				// Y
				AA = 'Y';
			}
			else if (i >= 48 && i <= 49 || i == 52) {		// Stop codons
				AA = 'Z';
			}
//			cout << AA << " " << codonFreq[i] << endl;
			pair<char, pair<int, float> > temp2 (AA, temp);
			AAtoCodonMap.insert (temp2);
	}
	return AAtoCodonMap;
}


// Prints each codon that codes for the amino acid
void CodonFrequency::printMap(multimap<char, pair<int, float> > AAtoCodonMap) {

	for (char AA = 'A'; AA <= 'Z'; AA++) {
		pair< multimap<char, pair<int, float> >::iterator, multimap<char, pair<int, float> >::iterator> ret;
		
		if (AAtoCodonMap.count(AA) <= 0) continue;
		
		ret = AAtoCodonMap.equal_range(AA);
		cout << AA << " =>";

		for (multimap<char, pair<int, float> >::iterator it = ret.first; it!=ret.second; ++it) {
			cout << ' ' << it->second.first;
		}
		cout << endl;
	}
}


// Returns the Codon with lowest frequency for that amino acid
float CodonFrequency::findMin(multimap<char, pair<int, float> > AAtoCodonMap, char AA) {

	float min = 100;
	
	pair< multimap<char, pair<int, float> >::iterator, multimap<char, pair<int, float> >::iterator> ret;
	
	ret = AAtoCodonMap.equal_range(AA);

	for (multimap<char, pair<int, float> >::iterator it = ret.first; it != ret.second; ++it) {

		if (min > it->second.second)
			min = it->second.second;
	}
	return min;
}


// Returns the Codon with highest frequency for that amino acid
float CodonFrequency::findMax(multimap<char, pair<int, float> > AAtoCodonMap, char AA) {

	float max = 0;
	
	pair< multimap<char, pair<int, float> >::iterator, 
		multimap<char, pair<int, float> >::iterator> ret;
	
	ret = AAtoCodonMap.equal_range(AA);

	for (multimap<char, pair<int, float> >::iterator it = ret.first; it!=ret.second; ++it) {

		if (max < it->second.second)
			max = it->second.second;
	}
	return max;
}


// Maps lowest frequency of codons to each amino acid
vector<float> CodonFrequency::createMinMap(multimap<char, pair<int, float> > AAtoCodonMap) {

	cout << "------minFreq------" << endl;

	float minFreq;
	vector<float> minMap;

	// Loops through all amino acids
	// If char rep of AA exists
	// Calculate smallest frequency to minMap
	for (char AA = 'A'; AA <= 'Z'; AA++) {
	
		int addressAA = AA - 65;

		if (AA == 'B' || AA == 'J' || AA == 'O' || AA == 'U' || AA == 'X') { 
			minMap.push_back(-1);
//			continue;
		}

		else if (AAtoCodonMap.count(AA) > 1) {
			
			minFreq = findMin(AAtoCodonMap, AA);
			minMap.push_back(minFreq);
//			continue;
			
		} else if (AAtoCodonMap.count(AA) == 1) {
			minMap.push_back(AAtoCodonMap.find(AA)->second.second);
		}
//		to print values: comment out the continues in if statements
		cout << addressAA << " " << AA << " " << minMap[addressAA] << endl;
	}
	return minMap;
}


// Maps highest frequency of codons to each amino acid
vector<float> CodonFrequency::createMaxMap(multimap<char, pair<int, float> > AAtoCodonMap){

	cout << "------maxFreq------" << endl;

	vector<float> maxMap;
	float maxFreq;

	// Loops through all amino acids
	// If char rep of AA exists
	// Calculate highest frequency to maxMap
	for (char AA = 'A'; AA <= 'Z'; AA++) {
		
		int addressAA = AA - 65;
		
		if (AA == 'B' || AA == 'J' || AA == 'O' || AA == 'U' || AA == 'X') {
		maxMap.push_back(-1);
//		continue;
		}

		else if (AAtoCodonMap.count(AA) > 1) {
			
			maxFreq = findMax(AAtoCodonMap, AA);
			maxMap.push_back(maxFreq);
//			continue;
			
		} else if (AAtoCodonMap.count(AA) == 1) {
			maxMap.push_back(AAtoCodonMap.find(AA)->second.second);
		}
//		to print values: comment out the continues in if statements
		cout << addressAA << " " << AA << " " << maxMap[addressAA] << endl;
	}
	return maxMap;
}



// Maps average frequency of codons to each amino acid
vector<float> CodonFrequency::createAvgMap(multimap<char, pair<int, float> > AAtoCodonMap) {
	
	cout << "------avgFreq------" << endl;
	
	int addressAA;
	int codonPerAA;	// Numbers of codons that code for an amino acid
	float sumFreq;	// Sum of codon frequencies for an amino acid
	float avgFreq;	// sumFreq/codonPerAA to calculate average frequency for an amino acid
	vector<float> avgMap;
	
	for (char AA = 'A'; AA <= 'Z'; AA++) {

		addressAA = AA-65;
		codonPerAA = 0;
		sumFreq = 0;
		avgFreq = 0;

		if (AA == 'B' || AA == 'J' || AA == 'O' || AA == 'U' || AA == 'X') { 
			avgMap.push_back(-1);
//			continue;
		}

		else if (AAtoCodonMap.count(AA) > 1) {
		
			pair< multimap<char, pair<int, float> >::iterator, 
			multimap<char, pair<int, float> >::iterator> ret;
	
			ret = AAtoCodonMap.equal_range(AA);
			
			for (multimap<char, pair<int, float> >::iterator it = ret.first; it != ret.second; ++it) {
//				cout << AA << " " << it->second.second << endl;
				sumFreq += it->second.second;
				codonPerAA++;
			}

			avgFreq = sumFreq/codonPerAA;
//			cout << sumFreq << " / " << codonPerAA << " = " << avgFreq << endl;
			avgMap.push_back(avgFreq);
//			continue;
			
		} else if (AAtoCodonMap.count(AA) == 1) {
			avgMap.push_back(AAtoCodonMap.find(AA)->second.second);
//			continue;
		}
		cout << addressAA << " " << AA << " " << avgMap[addressAA] << endl;
	}
	return avgMap;
}


vector<float> CodonFrequency::getMinMap() {
	return minMap;
}
vector<float> CodonFrequency::getMaxMap() {
	return maxMap;
}
vector<float> CodonFrequency::getAvgMap() {
	return avgMap;
}
float * CodonFrequency::getCodonFreq() {
	return codonFreq;
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
// Prints # of occurances for each codon for the vector of sequences
void CodonFrequency::printCodonCount() {

	for (int i = 0, count = 0; i < codon.size(); i++, ++count) {
		if (count%4 == 0) cout << endl;
		cout << codon[i].first << " " << (float) codon[i].second/getCodonCount()*1000 << "(" << codon[i].second << ")  ";
	}
//	cout << endl;
//	cout << "The total number of codons is " << getCodonCount() << endl;
}
*/

// Converts decimal number to binary
// Useful for converting decimal rep of codon to binary rep
// Returns a 6 bit binary string
string CodonFrequency::decimalToBinary(int n) {

	string binaryStr;
	deque <char> binary;

	while(n>1) {

		binary.push_back((n%2)+'0');
		n = n/2;
	}
	
	if (n == 1) {
		binary.push_back('1');
	}
	
	while (binary.size() < 6) binary.push_back('0');

	for (deque<char>::iterator it = binary.begin(); it!=binary.end(); ++it) {
		binaryStr.insert(binaryStr.begin(), *it);
	}
	return binaryStr;
}



// Converts a 6 bit binary string to codon
string CodonFrequency::binaryToCodon(string binaryStr) {
	
	string codon;

	string::iterator itr = binaryStr.begin(), end = binaryStr.end();
	for (string::iterator itr = binaryStr.begin(); itr != end; ++itr) {

		char a = *itr;
		itr++;
		char b = *itr;

//		cout << "A is " << a << " B is " << b << endl;
		if (a == '0' && b == '0') codon.append("A");
		else if (a == '0' && b == '1') codon.append("G");
		else if (a == '1' && b == '0') codon.append("C");
		else if (a == '1' && b == '1') codon.append("T");
	}
//	cout << codon << endl;
	return codon;
}


// Converts codon string of 3 nucleotides to binary representation of codon
// A = 00; G = 01; C = 10; T = 11;
int CodonFrequency::codonStrToBinaryRep(string codonStr) {

	int binaryRep = 0;

	for (string::iterator itr = codonStr.begin(), end = codonStr.end(); itr != end; ++itr) {

		// First nucleotide
		if (*itr == 'A') binaryRep += 0;						// 00xxxx
		else if (*itr == 'G') binaryRep += 16;			// 01xxxx
		else if (*itr == 'C') binaryRep += 32;			// 10xxxx
		else if (*itr == 'T') binaryRep += (32+16);	// 11xxxx
//		cout << binaryRep << endl;
		
		// Second nucleotide
		itr++;
		if (*itr == 'A') binaryRep += 0;						// xx00xx
		else if (*itr == 'G') binaryRep += 4;				// xx01xx
		else if (*itr == 'C') binaryRep += 8;				// xx10xx
		else if (*itr == 'T') binaryRep += (8+4);		// xx11xx
//		cout << binaryRep << endl;
		// Third nucleotide
		itr++;
		if (*itr == 'A') binaryRep += 0;						// xxxx00
		else if (*itr == 'G') binaryRep += 1;				// xxxx01
		else if (*itr == 'C') binaryRep += 2;				// xxxx10
		else if (*itr == 'T') binaryRep += (2+1);		// xxxx11
//		cout << binaryRep << endl;
	}
	return binaryRep;
}


// Converts codon to amino acid number representation 0-25
int * CodonFrequency::createCodonToAAMap(string codonStr, int CodontoAA[64]) {

	for (int j = 0; j < codonStr.size(); j+=3) {

		string triplet = codonStr.substr(j,3);

		int binaryRep = codonStrToBinaryRep(triplet);
//		cout << triplet << " " << binaryRep << endl;
		int i = binaryRep;

		if (i >= 24 && i <= 27) {							// A
			CodontoAA[binaryRep] = 'A'-65;
		}
		else if (i >= 54 && i <= 55) {				// C
			CodontoAA[binaryRep] = 'C'-65;
		}
		else if (i >= 18 && i <= 19) {				// D
			CodontoAA[binaryRep] = 'D'-65;
		}
		else if (i >= 16 && i <= 17) {				// E
			CodontoAA[binaryRep] = 'E'-65;
		}
		else if (i >= 62 && i <= 63) {				// F
			CodontoAA[binaryRep] = 'F'-65;
		}
		else if (i >= 20 && i <= 23) {				// G
			CodontoAA[binaryRep] = 'G'-65;
		}
		else if (i >= 34 && i <= 35) {				// H
			CodontoAA[binaryRep] = 'H'-65;
		}
		else if (i == 12 || i == 14 || i == 15) {				// I
			CodontoAA[binaryRep] = 'I'-65;
		}
		else if (i >= 0 && i <= 1) {					// K
			CodontoAA[binaryRep] = 'K'-65;
		}
		else if (i >= 60 && i <= 61 || i >= 44 && i <= 47) {				// L
			CodontoAA[binaryRep] = 'L'-65;
		}
		else if (i == 13) {										// M
			CodontoAA[binaryRep] = 'M'-65;
		}
		else if (i >= 2 && i <= 3) {				// N
			CodontoAA[binaryRep] = 'N'-65;
		}
		else if (i >= 40 && i <= 43) {				// P
			CodontoAA[binaryRep] = 'P'-65;
		}
		else if (i >= 32 && i <= 33) {				// Q
			CodontoAA[binaryRep] = 'Q'-65;
		}
		else if (i >= 36 && i <= 39 || i >= 4 && i <= 5) {				// R
			CodontoAA[binaryRep] = 'R'-65;
		}
		else if (i >= 56 && i <= 59 || i >= 6 && i <= 7) {				// S
			CodontoAA[binaryRep] = 'S'-65;
		}
		else if (i >= 8 && i <= 11) {					// T
			CodontoAA[binaryRep] = 'T'-65;
		}
		else if (i >= 28 && i <= 31) {				// V
			CodontoAA[binaryRep] = 'V'-65;
		}
		else if (i == 53) {				// W
			CodontoAA[binaryRep] = 'W'-65;
		}
		else if (i >= 50 && i <= 51) {				// Y
			CodontoAA[binaryRep] = 'Y'-65;
		}
		else if (i >= 48 && i <= 49 || i == 52) {		// Stop codons
			CodontoAA[binaryRep] = 'Z'-65;
		}
		else {
			CodontoAA[binaryRep] = -1;
		}
	}
//	for (int i = 0; i < 64; i++) {
//		cout << i << " " << CodontoAA[i] << endl;
//	}
	return CodontoAA;
}


int * CodonFrequency::getCodonToAAMap() {
	return codonToAAMap;
}


// Creates an output file with the inputfile name with .cf appended
// Output file contains codon, codon frequency, and # of occurances
void CodonFrequency::outputFileCF(char *file) {

  string ofilename(file);
	ofilename.append(".cf");

  ofstream ofile;
  ofile.open (ofilename.c_str());

  if (ofile.is_open()) {
  
		for (int i = 0, count = 0; i < NUM_TYPE_OF_CODONS; i++, count++) {
	
			string binaryStr = decimalToBinary(i);
			string codon = binaryToCodon(binaryStr);

			ofile << codon;
			ofile << " ";
			ofile << codonFreq[i];
			ofile << " (";
			ofile << codonOcc[i];
			ofile << ")";
			ofile << endl;
		}

	  ofile.close();
		cout << ofilename << " has been created." << endl;

	} else cout << "Unable to open " << ofilename << endl;
}

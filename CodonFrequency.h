/*

	Created by Kim Ngo 14/3/23
	
	CodonFrequency.cpp
	
	Calculates frequency for each codon in the gene
	Create library of codon frequencies
	Creates a library of maximum frequencies for each aa
	Creates a library of minimum frequencies for each aa
	
*/

#ifndef CODONFREQUENCY_H
#define CODONFREQUENCY_H

#include "Sequence.h"
#include <iostream>
#include <vector>
#include <iterator>
#include <string>
#include <utility> // pair
#include <map>	// stringmap for multimap

using namespace std;

class CodonFrequency {

	public:
		CodonFrequency(vector<Sequence>);			
			// Initializes codon vector with codon triplets and # of occurances
		//~CodonFrequency();
		int getCodonCount();
			// Getter for codonCount
		void set_codonCount(int);
			// Setter for codonCount
		void calcFreq(vector<Sequence>);	
			// Increments codon occurance & counts # of codons
			// Calculates freq -- #ofOcc/codonCount
		void printCodonCount();
			// Prints number of occurances and count
		void printFreq();
			// Prints codon frequency for a seq
		void outputFileCodonCount(ofstream &);	
			// Spits Codon|Codon Freq|count into output file name string
		string decimalToBinary(int);
			// Converts decimal number to 6 bit string
			// Used for converting decimal codon rep to binary codon rep
		string binaryToCodon(string);
			// Converts 6 bit binary string to codon string
			// A = 00; G = 01; C = 10; T = 11;
	
	private:
		vector< pair<string, int> > codon;
		int codonOcc [64];
			// Stores # of occurance for each codon
		float codonFreq [64];
			// Stores frequency for each codon
		int codonCount; 	
			// Total number of codons in sequences used to calculate frequency
		multimap<char, pair<int, float> > AAtoCodonMap;
			// Maps all codons with frequency to corresponding amino acid
		multimap<char, pair<int, float> > createMap(float []);
			// Creates mapping from codon & freq to amino acid
		void printMap(multimap<char, pair<int, float> > );
			// Prints each codon that codes for the amino acid
		float findMin(multimap<char, pair<int, float> >, char);
			// Returns codon with lowest frequency for that amino acid
		float findMax(multimap<char, pair<int, float> >, char);
			// Returns codon with highest frequency for that amino acid
		map<char, float> createMinMap(multimap<char, pair<int, float> >);
			// Maps lowest frequency of codons to each amino acid
		map<char, float> createMaxMap(multimap<char, pair<int, float> >);
			// Maps highest frequency of codons to each amino acid
		map<char, float> minMap;
		map<char, float> maxMap;

};

#endif

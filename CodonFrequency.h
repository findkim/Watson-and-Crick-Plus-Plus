/*

	Created by Kim Ngo 14/3/23
	
	CodonFrequency.cpp
	
	Calculates frequency for each codon in the gene
	Create library of codon frequencies
	
*/

#ifndef CODONFREQUENCY_H
#define CODONFREQUENCY_H

#include <iostream>
#include <vector>
#include <iterator>
#include <string>
#include <utility> // pair

using namespace std;

class CodonFrequency {

	public:
		CodonFrequency(string);			
			// Initializes codon vector with codon triplets and # of occurances
		float getCodonCount();				// Getter for codonCount
		void set_codonCount(int);		// Setter for codonCount
		void incrOccurance(string);	// Increments codon occurance & counts # of codons
		void calcFreq(string);						// Calculates freq -- #ofOcc/codonCount
		void printFreq();								// Prints codon frequency
	
	private:
		vector< pair<string, int> > codon;
			// Stores codon and # of occurances
		vector <float> codonFreqSeq;
			// Vector of codon frequencies that correspond with the seq using codon vector
		int codonCount; 	
			// Total number of codons in sequence used to calculate frequency

};

#endif

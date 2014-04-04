/*

	Created by Kim Ngo 14/3/23
	
	CodonFrequency.cpp
	
	Calculates frequency for each codon in the gene
	Create library of codon frequencies
	
*/

#ifndef CODONFREQUENCY_H
#define CODONFREQUENCY_H

#include "Sequence.h"
#include <iostream>
#include <vector>
#include <iterator>
#include <string>
#include <utility> // pair

using namespace std;

class CodonFrequency {

	public:
		CodonFrequency(vector<Sequence>);			
			// Initializes codon vector with codon triplets and # of occurances
		//~CodonFrequency();
		int getCodonCount();				// Getter for codonCount
		void set_codonCount(int);		// Setter for codonCount
		void incrOccurance(vector<Sequence>);	// Increments codon occurance & counts # of codons
		void calcFreq(vector<Sequence>);						// Calculates freq -- #ofOcc/codonCount
		void printCodonCount();			// Prints number of occurances and count
		void printFreq();						// Prints codon frequency for a seq
		void outputFileCodonCount(ofstream &);	
			// Spits Codon|Codon Freq|count into output file name string
		void outputfileFreq(ofstream &);
			// Spits each sequence as a series of frequencies in an output file
		void binary(int);
	
	private:
		vector< pair<string, int> > codon;
		int codonOcc [64];
		float codonFreq [64];
			// Stores codon and # of occurances
		vector <float> codonFreqSeq;
			// Vector of codon frequencies that correspond with the seq using codon vector
		int codonCount; 	
			// Total number of codons in sequence used to calculate frequency

};

#endif

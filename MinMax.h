/*

	Created by Kim Ngo 14/4/7
	
	MinMax.h
	
	Uses CodonFrequency class for min, max, and actual freq to calculate minmax using Clark & Clarke method
	Sliding window of 17 codons within each orfeome.
	Using a 17 codon window to identify clusters of rare codons instead of individual rare codons

*/

#ifndef MINMAX_H
#define MINMAX_H

#include "Sequence.h"
#include "CodonFrequency.h"
#include <map>
#include <vector>

using namespace std;

class MinMax : public CodonFrequency {

	public:
			MinMax(vector<Sequence>);
			
	private:
		vector< vector<float> > minMaxSequences;
			// Vector of %Min, %Max calculations in order of the sequence.
		vector<float> calcMinMax(Sequence, vector<float>, vector<float>, vector<float>, float *, int *);
			// a sequence, minMap, maxMap, avgMap, codonFreq, codonToAAMap
			// Calculates min, max, avg, and actual average frequencies for a codon window
};

#endif



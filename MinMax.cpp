/*

	Created by Kim Ngo 14/4/7
	
	MinMax.h
	
	Uses CodonFrequency class for min, max, and actual freq to calculate minmax using Clark & Clarke method
	Sliding window of 17 codons within each orfeome.
	Using a 17 codon window to identify clusters of rare codons instead of individual rare codons

*/

#include "MinMax.h"
#include "Sequence.h"
#include "CodonFrequency.h"
#include <vector>

const int WINDOWSIZE = 17;

using namespace std;

MinMax::MinMax(vector<Sequence> seq) : CodonFrequency(seq) {
	
	vector<float> minMap = getMinMap();
	vector<float> maxMap = getMaxMap();
	vector<float> avgMap = getAvgMap();
	float *codonFreq = getCodonFreq();
	int *codonToAAMap = getCodonToAAMap();
	
	for (int i = 0; i < seq.size(); i++) {
	
		vector<float> minMaxSeq = calcMinMax(seq[i], minMap, maxMap, avgMap, codonFreq, codonToAAMap);
//		minMaxSequences[i].push_back(minMaxSeq);
//		minMaxSequences.push_back(vector<float>());
//		minMaxSequences.back().push_back(minMaxSeq);
	}
}


// Calculates the max or min frequency for a 17 codon window
vector<float> calcMinMax(Sequence seq, vector<float> minMap, vector<float> maxMap, vector<float> avgMap, float *codonFreq, int *codonToAAMap) {

	int numCodonRep, AA;
	float minFreqWindowSum, maxFreqWindowSum, avgFreqWindowSum, actualFreqWindowSum;
	float minFreqWindowAvg, maxFreqWindowAvg, avgFreqWindowAvg, actualFreqWindowAvg;
	float percentMax, percentMin;
	vector<float> MinMax;

	for (int i = 0; i < seq.getSeqLength()-WINDOWSIZE; i++) {
	
		minFreqWindowSum = 0;
		maxFreqWindowSum = 0;
		avgFreqWindowSum = 0;
		actualFreqWindowSum = 0;
	
		for (int j = 0; j < WINDOWSIZE*3; j+=3) {
			
			string seqStr = seq.getSeq();
			string triplet = seqStr.substr(i+j,3);
			numCodonRep = CodonFrequency::codonStrToBinaryRep(triplet);
			AA = codonToAAMap[numCodonRep];
			
			minFreqWindowSum += minMap[AA];
			maxFreqWindowSum += maxMap[AA];
			avgFreqWindowSum += avgMap[AA];
			actualFreqWindowSum += codonFreq[numCodonRep];
		}
		
		minFreqWindowAvg = minFreqWindowSum / WINDOWSIZE;
		maxFreqWindowAvg = maxFreqWindowSum / WINDOWSIZE;
		avgFreqWindowAvg = avgFreqWindowSum / WINDOWSIZE;
		actualFreqWindowAvg = actualFreqWindowSum / WINDOWSIZE;
		
		percentMax = (actualFreqWindowAvg - avgFreqWindowAvg) / (maxFreqWindowAvg - avgFreqWindowAvg) * 100;
		
		if (percentMax < 0) {
			percentMin = (-1) * (avgFreqWindowAvg - actualFreqWindowAvg) / (avgFreqWindowAvg - minFreqWindowAvg) * 100;
			MinMax.push_back(percentMin);
			
		} else MinMax.push_back(percentMax);
	}
	return MinMax;
}

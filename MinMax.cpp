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
#include <fstream> // output file
#include <utility>
#include <string>

const int WINDOWSIZE = 17;

using namespace std;

MinMax::MinMax(char *filename, vector<Sequence> seq)
	: CodonFrequency(filename, seq) {
	
	vector<float> minMap = getMinMap();
	vector<float> maxMap = getMaxMap();
	vector<float> avgMap = getAvgMap();
	float *codonFreq = getCodonFreq();
	int *codonToAAMap = getCodonToAAMap();
	
	for (int i = 0; i < seq.size(); i++) {

		vector<float> minMaxSeq = 
		calcMinMax(seq[i], minMap, maxMap, avgMap, codonFreq, codonToAAMap);

		string nameDescription = seq[i].getSeqName();
		nameDescription.append("|");
		nameDescription.append(seq[i].getSeqDescription());

		pair<string, vector<float> > temp(nameDescription, minMaxSeq);
		minMaxSequences.push_back(temp);
	}
	outputFileMM(filename, minMaxSequences);
}


// Calculates the max or min frequency for a 17 codon window
vector<float> MinMax::calcMinMax(Sequence &seq, vector<float> &minMap, vector<float> &maxMap, vector<float> &avgMap, float *codonFreq, int *codonToAAMap) {

	int numCodonRep, AA;
	float minFreqWindowSum, maxFreqWindowSum, avgFreqWindowSum, actualFreqWindowSum;
	float minFreqWindowAvg, maxFreqWindowAvg, avgFreqWindowAvg, actualFreqWindowAvg;
	float percentMax, percentMin;
	vector<float> MinMaxValues;
	
	string seqStr = seq.getSeq();

	for (int i = 0; i < seq.getSeqLength()-WINDOWSIZE*3; i+=3) {
	
		minFreqWindowSum = 0;
		maxFreqWindowSum = 0;
		avgFreqWindowSum = 0;
		actualFreqWindowSum = 0;
	
		for (int j = 0; j < WINDOWSIZE*3; j+=3) {
		
			string triplet = seqStr.substr(i+j,3);
//			cout << "innner for loop" << endl;
			numCodonRep = CodonFrequency::codonStrToBinaryRep(triplet);
//			cout << seqStr << endl;
//			cout << "numcodonRep = " << numCodonRep << endl;
			AA = codonToAAMap[numCodonRep];
//			cout << "AA " << AA << " " << codonToAAMap[numCodonRep] << endl;

			minFreqWindowSum += minMap[AA];
//			cout << triplet << " " << minMap[AA] << " " << minFreqWindowSum << endl;
			maxFreqWindowSum += maxMap[AA];
			avgFreqWindowSum += avgMap[AA];
			actualFreqWindowSum += codonFreq[numCodonRep];
		}
		
		minFreqWindowAvg = minFreqWindowSum / WINDOWSIZE;
		maxFreqWindowAvg = maxFreqWindowSum / WINDOWSIZE;
		avgFreqWindowAvg = avgFreqWindowSum / WINDOWSIZE;
		actualFreqWindowAvg = actualFreqWindowSum / WINDOWSIZE;
		
//		cout << "minFreqWindowAvg " << minFreqWindowAvg << endl;
		
		percentMax = (actualFreqWindowAvg - avgFreqWindowAvg) / (maxFreqWindowAvg - avgFreqWindowAvg) * 100;
		
		if (percentMax < 0) {
			percentMin = (-1) * (avgFreqWindowAvg - actualFreqWindowAvg) / (avgFreqWindowAvg - minFreqWindowAvg) * 100;
			MinMaxValues.push_back(percentMin);
			
		} else MinMaxValues.push_back(percentMax);

	}
	return MinMaxValues;
}


// Stores MinMax frequency for each codon for the vector of sequences in an output file
void MinMax::outputFileMM(char *file, vector< pair< string, vector<float> > > minMaxSequences) {

	string filename(file);
	filename.append(".mm");	

	ofstream ofile;
	ofile.open (filename.c_str());
	
	if (ofile.is_open()) {
	
		for (int i = 0; i < minMaxSequences.size(); i++) {

			ofile << minMaxSequences[i].first;
			ofile << endl;

			for (int j = 0; j < minMaxSequences[i].second.size(); j++) {
			
				ofile << minMaxSequences[i].second[j];
				ofile << ",";
			}
			ofile << endl;
		}
		ofile.close();
		cout << filename << " has been created." << endl;
		
	} else cout << "Unable to open " << filename << endl;
}

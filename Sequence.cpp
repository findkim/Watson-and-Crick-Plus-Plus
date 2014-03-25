//
//  Sequence.cpp
//  
//
//  Created by Xuanyi Li on 14-3-8.
//
//

#include "Sequence.h"
#include "CodonFrequency.h"	// Composition, each seq can calculate codon freq
#include <iostream>
#include <string>
#include <vector>

using namespace std;

Sequence :: Sequence(string name, string description, string sequence){
    seqName = name;
    seqDescription = description;
    seq = sequence;
    seqLength = seq.size();
    
   	// Does not initialize codonFreq for alignments
    if (seq.find("-"))
    	codonFreq = new CodonFrequency (seq);
	  else codonFreq = NULL;
}

string Sequence::getSeqName(){
	return seqName;
}


string Sequence::getSeqDescription(){
	return seqDescription;
}


int Sequence::getSeqLength(){
	return seqLength;
}


// Prints reference number, description, seq length
void Sequence::print(){

	cout << "The sequence name is: " << getSeqName() << endl;
	cout << "The sequence description: " << getSeqDescription() << endl;
	cout << "The sequence length is " << getSeqLength() << endl;
	
	// Does not print for alignments
	if (codonFreq) {
		codonFreq->printFreq();
		cout << endl;
	}
}


// Prints sequence
void Sequence::printSeq(){
    cout << seq << endl;

	cout << "The codon count is " << codonFreq->getCodonCount() << endl;
	codonFreq->printFreq();

	/*ostream_iterator<char> output (cout, "");
	copy (seq.begin(), seq.end(), output);
	cout << endl; */
}

char Sequence :: operator[](int i){
    if (i<0 || i>seq.size()) {
		if (i < 0 || i >= getSeqLength()) {
    		throw "Subscript out of range";
        }
    }
    return seq[i];
}

void Sequence :: addGap(int i){
    seq.insert(i+1,"-");
}

//
//  Sequence.cpp
//  
//
//  Created by Xuanyi Li on 14-3-8.
//
//

#include "Sequence.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <list>			// STL lists
#include <vector>
#include <iterator>	// ostream iterators

const int FASTA_DESCRIPTION_LINE_NUM_ARG = 5;

using namespace std;

Sequence :: Sequence(string name, string description, string sequence){
    seqName = name;
    seqDescription = description;
    seq = sequence;
    seqLength = seq.size();
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
}


// Prints sequence
void Sequence::printSeq(){
    cout << seq << endl;

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

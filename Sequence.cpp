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
#include <vector>

using namespace std;

Sequence :: Sequence(string name, string description, string sequence){
    seqName = name;
    seqDescription = description;
    seq = sequence;
    seqLength = seq.size();
}

string Sequence::getSeq(){
	return seq;
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


int Sequence::getNumCodon(){
	if (getSeqLength()%3 == 0)
		return getSeqLength()/3;
	else return 0;
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

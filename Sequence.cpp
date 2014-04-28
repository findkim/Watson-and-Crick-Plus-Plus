//
//  Sequence.cpp
//  
//
//  Created by Xuanyi Li on 14-3-8.

//  Store name, description and actual sequence (get and set functions)
//  add or delete gaps in the sequence

#include "Sequence.h"
#include <iostream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;

Sequence :: Sequence(string name, string description, string sequence){
    seqName = name;
    seqDescription = description;
    seq = sequence;
    seqLength = seq.size();
}

Sequence :: Sequence(){};

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
	return seq.size();
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

string Sequence :: operator[](int i){
    // check out of range
    if (i<0 || i>seq.size()) {
		if (i < 0 || i >= getSeqLength()) {
    		throw "Subscript out of range";
        }
    }
    // convert a char to a string
    stringstream ss;
    string s;
    ss << seq[i];
    ss >> s;
    return s;
}
// add gap before i
void Sequence :: addGap(int i){
    seq.insert(i,"-");
}
void Sequence :: setSeq(string s){
    seq = s;
}
// remove the gap added to the front (check if there is a gap in the alignment class)
void Sequence :: removeGapfront(){
    seq.erase(seq.begin());
}
void Sequence :: addGaptoEnd(int g){
    for (int i = 0; i<g; i++) {
        seq+="-";
    }
}

//Add a Domain to the member vector Domains
void Sequence :: addDomain(Domain addDomain)
{
	Domains.push_back(addDomain);
}




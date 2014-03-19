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
Sequence :: Sequence(char *filename){

	int lengthCounter = 0;
	int i = 0;

	cout << filename << endl;		// Prints filename
	ifstream file;
	file.open(filename);	// Opens file
	
	while (!file.eof()) { // Reads until reaches end of file
	
		string line;
		char value;
		
		getline(file, line);									// First line of FASTA file		
		istringstream ss(line);
		string token;
		
		while(getline(ss, token, '|')) {			// Parses first line of FASTA file
			if (i == 3) seqName = token;				// 3rd token is ref number
			if (i == 4) seqDescription = token;	// 4th token is seq description
			i++;
		}

		while (getline(file, line)) {					// The remaining lines is the sequence
			
			stringstream split(line);	// Parses line into chars
			while (split >> value) {
				seq+=value;		// Adds each char value to list
				lengthCounter++;				// Counts length of nucleotides/amino acids
			}
		}
	}

	seqLength = lengthCounter;

	file.close();
	cout << "Done reading " << filename << endl;
    
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

	cout << "The reference number is " << getSeqName() << endl;
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
            throw out_of_range("Subscript out of range");
    		throw "Subscript out of range";
        }
    }
    return seq[i];
}

void Sequence :: addGap(int i){
    seq.insert(i+1,"-");
}

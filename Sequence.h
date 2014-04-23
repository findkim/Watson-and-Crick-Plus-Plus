//
//  Sequence.h
//  
//
//  Created by Xuanyi Li on 14-3-8.
//
//

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <iostream>
#include <string>

using namespace std;

class Sequence{
public:
    Sequence(string, string, string);			// Constructor with parameter: name, description and sequence.
    void addGap(int); 		// add gap after int
    char operator[](int);
    void print();					// Prints name, description, legnth.
    void printSeq();			// Prints sequence
    string getSeqName();	// Accessor functions for name, seq, description, length
    string getSeq();
    string getSeqDescription();
    int getSeqLength();
    
private:
    string seqName;				// Reference number for DNA/Protein sequence
    string seqDescription;// Description of seqeuence
    string seq;				// Sequence of nucleotides or amino acids
    int seqLength;				// Length of string--converted string to int
};
#endif /* defined(____Sequence__) */

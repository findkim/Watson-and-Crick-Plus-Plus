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
#include <list>
#include <string>

using namespace std;

class Sequence{

	public:
    Sequence(char *);			// Constructor takes name of FASTA file
    void addGap(int); 		// add gap after int
    char operator[](int);
    void print();					// Prints name, description, legnth.
    void printSeq();			// Prints sequence

    string getSeqName();	// Accessor functions for name, seq, description, length
    string getSeq();
    string getSeqDescription();
    int getSeqLength();
    char *location_ptr;			// Stores location on sequence

	private:
    string seqName;				// Reference number for DNA/Protein sequence
    string seqDescription;// Description of seqeuence
    list<char> seq;				// Sequence of nucleotides or amino acids
    int seqLength;				// Length of string--converted string to int
};
#endif /* defined(____Sequence__) */

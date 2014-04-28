//
//  Sequence.h
//  
//
//  Created by Xuanyi Li on 14-3-8.

//  Store name, description and actual sequence (get and set functions)
//  add or delete gaps in the sequence

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "Domain.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class Sequence{
public:
    Sequence(string, string, string);			// Constructor with parameter: name, description and sequence.
    Sequence(); // default constructor
    void addGap(int); 		// add gap before int
    string operator[](int);
    void print();					// Prints name, description, legnth.
    void printSeq();			// Prints sequence
    string getSeqName();	// Accessor functions for name, seq, description, length
    string getSeq();
    string getSeqDescription();
    int getSeqLength();
    int getNumCodon();			// Calculates the number of codons if seq is of proper length
    char *location_ptr;			// Stores location on sequence
    void addDomain(Domain);
    void setSeq(string);
    void removeGapfront(); // remove the gap added to the beginning of sequence
    void addGaptoEnd(int); // add number of gaps to the end of the sequence

    
private:
    string seqName;				// Reference number for DNA/Protein sequence
    string seqDescription;// Description of seqeuence
    string seq;				// Sequence of nucleotides or amino acids
    int seqLength;				// Length of string--converted string to int
    int numCodon;					// Numbers of codons if of proper length (mult of 3)
    vector<Domain> Domains;		//domains present on the sequence	

};
#endif /* defined(____Sequence__) */

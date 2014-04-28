//
//  ExtractSequence.h
//  
//
//  Created by Xuanyi Li on 14-3-18.

//  Read in fasta file and store the sequences in a vector of sequences
//  get functions; add gaps to all sequences in the vector.

#ifndef ____ExtractSequence__
#define ____ExtractSequence__

#include <iostream>
#include <string>
#include <vector>
#include "Sequence.h"
#include "MinMax.h"

using namespace std;

class ExtractSequence{

	public:
    ExtractSequence(char *); // constructor takes in filename
    ExtractSequence(); // default constructor
    vector<Sequence> getVectorOfSequences();	// Returns the vector of sequences
    Sequence getSequence(string);// searching for Sequence by name
    void getHeader(string); // get name and description from header line
    int getSize(); // get the size of the Sequences vector
    void printSequences();
    Sequence operator[](int);
    void addSequence(Sequence);
    void remove1Seq(int); // remove 1 sequence from the vector
    void addGapstoAllEnds(int); // add gaps to the ends of all sequences
    void addGapstoAll(int); // add one gap to all sequences at a position
private:
    vector<Sequence> Sequences; // vector of Sequence object
    vector<string> headers; // store name and description.

    void addDomains(string, vector<Sequence>); //Takes a text file and adds Domains to their respective Sequences	

};
#endif /* defined(____ExtractSequence__) */

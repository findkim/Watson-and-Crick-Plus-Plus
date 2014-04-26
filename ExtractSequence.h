//
//  ExtractSequence.h
//  
//
//  Created by Xuanyi Li on 14-3-18.
//
//

#ifndef ____ExtractSequence__
#define ____ExtractSequence__

#include <iostream>
#include <string>
#include <vector>
#include "Sequence.h"
#include "CodonFrequency.h"

using namespace std;

class ExtractSequence{
public:
    ExtractSequence(char *); // constructor takes in filename
    ExtractSequence(); // default constructor
    void getHeader(string); // get name and description from header line
    Sequence getSequence(string);// searching for Sequence by name
    int getSize(); // get the size of the Sequences vector
    ~ExtractSequence();
    void printSequences();
    //void outputfile(char *);
    //vector<Sequence> removeSeq(vector < Sequence >);
    Sequence operator[](int);
    void addSequence(Sequence);
    void remove1Seq(int); // remove 1 sequence from the vector
    void addGapstoAllEnds(int); // add gaps to the ends of all sequences
    void addGapstoAll(int); // add one gap to all sequences at a position
private:
    vector<Sequence> Sequences; // vector of Sequence object
    vector<string> headers; // store name and description.
    //CodonFrequency *codonFreq;	// Composition; calculates codon frequency for that seq (not for alignments)
};
#endif /* defined(____ExtractSequence__) */

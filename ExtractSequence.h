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
#include "MinMax.h"

using namespace std;

class ExtractSequence{
public:
    ExtractSequence(char *); // constructor takes in filename
    void getHeader(string); // get name and description from header line
    Sequence getSequence(string);// searching for Sequence by name
    int getSize(); // get the size of the Sequences vector
    ~ExtractSequence();
    void printSequences();
    void outputfileCF(char *);
    void outputfileMM(char *);
    vector<Sequence> removeSeq(vector < Sequence >);
    Sequence operator[](int);
    
private:
    vector<Sequence> Sequences; // vector of Sequence object
    vector<string> headers; // store name and description.
    MinMax *minMax;	// Composition: calculates minmax values from CodonFreq
};
#endif /* defined(____ExtractSequence__) */

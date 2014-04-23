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

using namespace std;

class ExtractSequence{
public:
    ExtractSequence(char *); // constructor takes in filename
    void printSequences(int); // print sequences: name description and sequence
    void getHeader(string); // get name and description from header line
    Sequence getSequence(string);// searching for Sequence by name
    int getSize(); // get the size of the Sequences vector
    Sequence operator[](int);
    
private:
    vector<Sequence> Sequences; // vector of Sequence object
    vector<string> headers; // store name and description.
};
#endif /* defined(____ExtractSequence__) */

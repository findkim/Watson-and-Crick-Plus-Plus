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
    ExtractSequence(char *);
    void printSequences();
    void getHeader(string);
    /*void calAlignment();
    void displayAlignment(); */
private:
    vector<Sequence> Sequences; // vector of Sequence object
    vector<string> headers; // store name and description.
};
#endif /* defined(____ExtractSequence__) */

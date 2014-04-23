//
//  AlignmentCIGAR.h
//  
//
//  Created by Xuanyi Li on 14-3-21.
//
//  

#ifndef ____AlignmentCIGAR__
#define ____AlignmentCIGAR__

#include <iostream>
#include "ExtractSequence.h"
#include "Sequence.h"
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

class AlignmentCIGAR{
public:
    AlignmentCIGAR(char *filename);
    int checkAligned(); // check whether the vector of sequence is aligned already
    string cigarOneSeq(int); // Compare one Sequence to the reference and return this cigar representation.
    void setCigar(); // transform all the sequences in ExtractSequence into cigar string
    void printCigar(); // print the cigar string for the alignment
private:
    ExtractSequence sequences; // stored the input sequences
    string ref; // reference sequence for alignment
    vector<string> cigar; // store the alignment format before output to file
};
#endif /* defined(____AlignmentCIGAR__) */

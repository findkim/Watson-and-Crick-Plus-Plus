//
//  PairAlignment.h
//  
//
//  Created by 李 萱伊 on 14-4-8.
//
//

#ifndef ____PairAlignment__
#define ____PairAlignment__

#include <iostream>
#include "ExtractSequence.h"
#include <vector>
#include "Sequence.h"
#include <string>

class PairAlignment{
public:
    PairAlignment(ExtractSequence);
    void computeTables(); // updating the score and direction table
    void updateAlign(); // align the two sequence by adding gap to align 1 and 2
    void printAlignment(); // print the pairwise alignment TESTING
    
private:
    string align1;
    string align2;
    string seq1;
    string seq2;
    vector<vector<int> > score; // +2 for match; -1 for mismatch; -2 for gap
    vector<vector<int> > direction; // 0 for diagonal; 1 for horizontal; 2 for vertical
    int max (int,int,int); // helper function for finding the max number
};
#endif /* defined(____PairAlignment__) */

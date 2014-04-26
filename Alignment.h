//
//  Alignment.h
//  
//
//  Created by 李 萱伊 on 14-4-8.
//
//

#ifndef ____Alignment__
#define ____Alignment__

#include <iostream>
#include "ExtractSequence.h"
#include <vector>
#include "Sequence.h"
#include <string>

class Alignment{
public:
    Alignment(char *);
    void computeTables(int); // updating the score and direction table (align one sequence to the center star)
    void updateAlign(int); // align the two sequence by adding gap to align 1 and 2
    void printAlignment(); // print all aligned sequences
    int compute_distance(int,int); // compute the distance between two sequences
    void write_distance();// fill the 2d vector distance
    int findCenterStar();
    void alignAll(); // multiple sequence alignment
    void printTables(); //TEST
    
private:
    ExtractSequence seqs; // store string representation of sequences to be aligned
    ExtractSequence aligned; // store aligned sequences as Sequence objects with names and descriptions.
    vector<vector<int> > score; // +2 for match; -1 for mismatch; -2 for gap
    vector<vector<int> > direction; // 0 for diagonal; 1 for horizontal; 2 for vertical
    vector<vector<int> > distance;
    int max (vector<int>); // helper function for finding the max number
    int min (vector<int>); // helper function for finding the min number
    Sequence centerStar;
    string starAlign; // temperarilly storing aligned centerstar
};
#endif /* defined(____PairAlignment__) */

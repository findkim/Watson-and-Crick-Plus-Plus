//
//  NeighborJoining.h
//  
//
//  Created by 李 萱伊 on 14-4-2.
//
//

#ifndef ____NeighborJoining__
#define ____NeighborJoining__

#include <iostream>
#include "ExtractSequence.h"
#include "Sequence.h"
#include <fstream>
#include <vector>
using namespace std;

class NeighborJoining{
public:
    NeighborJoining(char *filename);
    int compute_distance(int,int);// compute the distance between two sequences
    void write_distance();// fill the 2d vector distance
    void print_distance();
private:
    ExtractSequence sequences;
    vector<vector<int> > distance;
};
#endif /* defined(____NeighborJoining__) */

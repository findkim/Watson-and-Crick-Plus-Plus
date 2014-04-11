//
//  NeighborJoining.cpp
//  
//
//  Created by 李 萱伊 on 14-4-2.
//
//

#include "NeighborJoining.h"
#include <iostream>
#include "ExtractSequence.h"
#include "Sequence.h"
#include <fstream>
#include <vector>
using namespace std;

NeighborJoining :: NeighborJoining(char *filename):sequences(filename){};
int NeighborJoining :: compute_distance(int i, int j){
    int r = sequences[i].getSeqLength();
    int c = sequences[j].getSeqLength();
    int score = 0;
    int k = 0;
    while (k<r && k<c) {
        if (sequences[i][k] != sequences[j][k]) {
            score ++;
        }
        k ++;
    }
    return score;
}
void NeighborJoining :: write_distance(){
    int size = sequences.getSize();
    for (int i = 0; i<size; i++) {
        vector<int> rdistance;
        for (int j = 0; j<size; j++) {
            rdistance.push_back(compute_distance(i,j));
        }
        distance.push_back(rdistance);
    }
}
void NeighborJoining :: print_distance(){
    for (int i = 0; i<distance.size(); i++) {
        for (int j = 0; j<distance[i].size(); j++) {
            cout << distance[i][j] << " ";
        }
        cout << endl;
    }
}
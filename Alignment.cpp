//
//  Alignment.cpp
//  
//
//  Created by 李 萱伊 on 14-4-8.
//
//

#include "Alignment.h"
#include "ExtractSequence.h"
#include <iostream>
#include <vector>
#include "Sequence.h"
#include <string>

using namespace std;

Alignment :: Alignment(ExtractSequence sequences){
    for (int i = 0; i<sequences.getSize(); i++) {
        seqs.push_back(sequences[i].getSeq());
    }
}
void Alignment :: computeTables(int i){
    int r = seqs[i].size();
    int c = seqs[centerStar].size();
    // updating first rows in score and direction
    vector<int> rscore0;
    vector<int> rdirection0;
    for (int k = 0; k<=c; k++) {
        rscore0.push_back(-2*k); // all gaps
        rdirection0.push_back(1); // first row all horizontal
    }
    score.push_back(rscore0);
    direction.push_back(rdirection0);
    // updating the rest of the tables
    for (int i = 1; i<=r; i++) {
        vector<int> rscore;
        vector<int> rdirection;
        rscore.push_back(-2*i); // all gaps
        rdirection.push_back(2); // first column all vertical
        char base1 = align1[i-1]; // the base in the first sequence to be compaired
        for (int j = 1; j<=c; j++) {
            int diagonal;
            int horizontal = rscore[j-1] - 2; // mismatch from left
            int vertical = score[i-1][j] - 2; // mismatch from above
            // check diagonal
            if (base1 == align2[j-1]) {
                diagonal = score[i-1][j-1] + 2;
            }
            else{
                diagonal = score[i-1][j-1] - 1;
            }
            vector<int> scores;
            scores.push_back(diagonal);
            scores.push_back(horizontal);
            scores.push_back(vertical);
            int maxscore = max(scores);
            if (maxscore == diagonal) {
                rdirection.push_back(1);
            }
            else if (maxscore == vertical){
                rdirection.push_back(2);
            }
            else{
                rdirection.push_back(0);
            }
            rscore.push_back(maxscore);
        }
        score.push_back(rscore);
        direction.push_back(rdirection);
    }
}
void Alignment :: updateAlign(int i){
    int r = seqs[i].size();
    int c = seqs[centerStar].size();
    // table is r+1 * c+1
    seqs[i].insert(0,"-");
    align1.insert(0,1,seq1[r]);
    align2.insert(0,1,seq2[c]);
    int i = r;
    int j = c;
    int d = direction[i][j];
    while (j>1 || i>1) {
        int newd,newi,newj;
        // diagonal, match/mismatch
        if (d==0) {
            align1.insert(0,1,seq1[i-1]);
            align2.insert(0,1,seq2[j-1]);
            newd = direction[i-1][j-1];
            newi = i - 1;
            newj = j - 1;
        }
        // horizontal, gap in align1
        else if (d == 1){
            align1.insert(0,"-");
            align2.insert(0,1,seq2[j-1]);
            newd = direction[i][j-1];
            newi = i;
            newj = j - 1;
        }
        // vertical, gap in align2
        else{
            align2.insert(0,"-");
            align1.insert(0,1,seq1[i-1]);
            newd = direction[i-1][j];
            newi = i - 1;
            newj = j;
        }
        d = newd;
        i = newi;
        j = newj;
    }
}
void Alignment :: printAlignment(){
    cout << "ALIGN1: "<< endl;
    cout << align1 << endl;
    cout << "ALIGN2: "<< endl;
    cout << align2 << endl;
}
int Alignment :: max(vector<int> x){
    int max = 0;
    for (int i = 0; i<x.size(); i++) {
        if (x[i]>x[max]) {
            max = i;
        }
    }
    return max;
}
int Alignment :: compute_distance(int i, int j){
    int r = seqs[i].length();
    int c = seqs[j].length();
    int score = 0;
    int k = 0;
    while (k<r && k<c) {
        if (seqs[i][k] != seqs[j][k]) {
            score ++;
        }
        k ++;
    }
    return score;
}
void Alignment :: write_distance(){
    int size = seqs.size();
    for (int i = 0; i<size; i++) {
        vector<int> rdistance;
        for (int j = 0; j<size; j++) {
            rdistance.push_back(compute_distance(i,j));
        }
        distance.push_back(rdistance);
    }
}
void Alignment :: findCenterStar(){
    vector<int> sumofdistance; // sum of distance from one string to all other strings
    for (int i = 0; i<distance.size(); i++) {
        int sum = 0;
        for (int j = 0; j<distance[i].size(); j++) {
            sum += distance[i][j];
        }
        sumofdistance.push_back(sum);
    }
    centerStar = max(sumofdistance);
}

//
//  PairAlignment.cpp
//  
//
//  Created by 李 萱伊 on 14-4-8.
//
//

#include "PairAlignment.h"
#include "ExtractSequence.h"
#include <iostream>
#include <vector>
#include "Sequence.h"
#include <string>

using namespace std;

PairAlignment :: PairAlignment(ExtractSequence seqs){
    seq1 = seqs[0].getSeq();
    seq2 = seqs[1].getSeq();
}
void PairAlignment :: computeTables(){
    int r = seq1.size();
    int c = seq2.size();
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
            int maxscore = max(diagonal,horizontal,vertical);
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
void PairAlignment :: updateAlign(){
    int r = seq1.size();
    int c = seq2.size();
    // table is r+1 * c+1
    seq1.insert(0,"-");
    seq2.insert(0,"-");
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
void PairAlignment :: printAlignment(){
    cout << "ALIGN1: "<< endl;
    cout << align1 << endl;
    cout << "ALIGN2: "<< endl;
    cout << align2 << endl;
}
int PairAlignment :: max(int a, int b, int c){
    if (a>=b && a>=c) {
        return a;
    }
    else if (b>=a && b>=c){
        return b;
    }
    else{
        return c;
    }
}
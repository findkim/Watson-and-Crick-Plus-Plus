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

Alignment :: Alignment(char *filename):seqs(filename){
    write_distance();
    int s = findCenterStar();
    centerStar = seqs[s]; // store the Sequence object for center star
    seqs.remove1Seq(s); // remove the center star sequence from seqs
    starAlign = centerStar.getSeq();
}
void Alignment :: computeTables(int t){
    direction.clear();
    score.clear();
    int r = seqs[t].getSeqLength();
    int c = centerStar.getSeqLength();
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
        string base1 = seqs[t][i-1]; // the base in the first sequence to be compaired
        for (int j = 1; j<=c; j++) {
            int diagonal = 0;
            int horizontal = rscore[j-1] - 2; // mismatch from left
            int vertical = score[i-1][j] - 2; // mismatch from above
            // check diagonal
            if (base1 == centerStar[j-1]) {
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
            if (scores[maxscore] == diagonal) {
                rdirection.push_back(0);
            }
            else if (scores[maxscore] == vertical){
                rdirection.push_back(2);
            }
            else{
                rdirection.push_back(1);
            }
            rscore.push_back(scores[maxscore]);
        }
        score.push_back(rscore);
        direction.push_back(rdirection);
    }
}
void Alignment :: updateAlign(int t){
    /*if (centerStar[0]=="-") {
        centerStar.removeGapfront();
    } */
    int r = seqs[t].getSeqLength();
    int c = centerStar.getSeqLength();
    // table is r+1 * c+1
    Sequence se = seqs[t];
    se.addGap(0);
    centerStar.addGap(0);
    string a;
    string b;
    a.insert(0,se[r]);
    b.insert(0,centerStar[c]);
    int i = r;
    int j = c;
    int d = direction[i][j];
    while (j>1 || i>1) {
        int newd,newi,newj;
        // diagonal, match/mismatch
        if (d==0) {
            a.insert(0,se[i-1]);
            b.insert(0,centerStar[j-1]);
            newd = direction[i-1][j-1];
            newi = i - 1;
            newj = j - 1;
        }
        // horizontal, gap in align1
        else if (d == 1){
            a.insert(0,"-");
            b.insert(0,centerStar[j-1]);
            newd = direction[i][j-1];
            newi = i;
            newj = j - 1;
        }
        // vertical, gap in align2
        else{
            b.insert(0,"-");
            a.insert(0,se[i-1]);
            newd = direction[i-1][j];
            newi = i - 1;
            newj = j;
        }
        d = newd;
        i = newi;
        j = newj;
    }
    // PROBLEM
    for (int i = 0; i<b.size() && i < starAlign.size(); i++) {
        // new gap in the aligned
        if (b[i]=='-' && starAlign[i]!='-') {
            starAlign.insert(i,"-");
            aligned.addGapstoAll(i);
        }
        // add existed gap to the newly aligned
        if (starAlign[i]=='-' && b[i]!='-') {
            b.insert(i,"-");
            a.insert(i,"-");
        }
    }
    if (b.size()>starAlign.size()) {
        int num = b.size()-starAlign.size();
        for (int i = 0; i<num; i++) {
            starAlign+="-";
        }
        aligned.addGapstoAllEnds(num);
    }
    else{
        int num = starAlign.size() - b.size();
        for (int i = 0; i<num; i++) {
            a+="-";
        }
    }
    Sequence s = seqs[t];
    s.setSeq(a);
    aligned.addSequence(s);
    centerStar.setSeq(b);
}
void Alignment :: printAlignment(){
    // center star
    cout << centerStar.getSeqName() << "|" <<centerStar.getSeqDescription() << "    ";
    centerStar.printSeq();
    cout << endl;
    for (int i = 0; i<aligned.getSize(); i++) {
        cout << aligned[i].getSeqName() << "|" <<aligned[i].getSeqDescription() << "    ";
        aligned[i].printSeq();
        cout << endl;
    }
    /*
    for (int i = 0; i<aligned.getSize(); i++) {
        cout << aligned[i].getSeqLength() << "    ";
    } */
    cout << endl;
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

int Alignment :: min(vector<int> x){
    int min = 0;
    for (int i = 0; i<x.size(); i++) {
        if (x[i]<x[min]) {
            min = i;
        }
    }
    return min;
}

int Alignment :: compute_distance(int i, int j){
    int r = seqs[i].getSeqLength();
    int c = seqs[j].getSeqLength();
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
    int size = seqs.getSize();
    for (int i = 0; i<size; i++) {
        vector<int> rdistance;
        for (int j = 0; j<size; j++) {
            rdistance.push_back(compute_distance(i,j));
        }
        distance.push_back(rdistance);
    }
}
int Alignment :: findCenterStar(){
    vector<int> sumofdistance; // sum of distance from one string to all other strings
    for (int i = 0; i<distance.size(); i++) {
        int sum = 0;
        for (int j = 0; j<distance[i].size(); j++) {
            sum += distance[i][j];
        }
        sumofdistance.push_back(sum);
    }
    int star = min(sumofdistance);
    return star;
    
}

void Alignment :: alignAll(){
    for (int i = 0; i<seqs.getSize(); i++) {
        computeTables(i);
        updateAlign(i);
    }
    //centerStar.setSeq(starAlign);
    
}
void Alignment :: printTables(){
    cout << "SCORE" << endl;
    for (int i = 0; i<score.size(); i++) {
        for (int j = 0; j<score[i].size(); j++) {
            cout << score[i][j] << " ";
        }
        cout << endl<<endl;
    }
    cout << "DIRECTION" << endl;
    for (int i = 0; i<direction.size(); i++) {
        for (int j = 0; j<direction[i].size(); j++) {
            cout << direction[i][j] << " ";
        }
        cout << endl;
    }
}
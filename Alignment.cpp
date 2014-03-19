//
//  Alignment.cpp
//  
//
//  Created by Xuanyi Li on 14-3-18.
//
//

#include "Alignment.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
using namespace std;

// I basically just copied your loadfile code... Alignment files contains multiple sequences.
Alignment :: Alignment(char *filename){
    ifstream file;
    file.open(filename);
    string sequence;
    string token,name,description,line;
    int i=0;
    // partition before and after '|'
    getline(file,line);
    istringstream ss(line);
    while (getline(ss,token,'|')) {
        if (i==0) {
            name = token;
        }
        else{
            description = token;
        }
    }
    while (!file.eof()) {
        getline(file,line);
        if (line[0]=='>') {
            Sequence seq (name,description,sequence);
            AlignSeqs.push_back(seq);
            sequence.clear();
            istringstream ss(line);
            int i=0;
            // partition before and after '|'
            while (getline(ss,token,'|')) {
                if (i==0) {
                    name = token;
                }
                else{
                    description = token;
                }
            }
        }
        else{
            sequence += line;
        }
    }
    Sequence seq (name,description,sequence);
    AlignSeqs.push_back(seq);
}
void Alignment :: printAlignment(){
    for (int i = 0; i<AlignSeqs.size(); i++) {
        AlignSeqs[i].print();
        AlignSeqs[i].printSeq();
        cout << "-----------------------" << endl;
    }
    cout << "Number of sequences: " << AlignSeqs.size() << endl;
}

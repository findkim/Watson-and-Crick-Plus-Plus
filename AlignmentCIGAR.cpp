//
//  AlignmentCIGAR.cpp
//  
//
//  Created by 李 萱伊 on 14-3-21.

//  transform an alignment into cigar string
//

#include "AlignmentCIGAR.h"
#include <iostream>
#include "ExtractSequence.h"
#include "Sequence.h"
#include <sstream>
#include <fstream>
#include <vector>
using namespace std;

AlignmentCIGAR :: AlignmentCIGAR(char *filename):sequences(filename){
    ref = sequences[0].getSeq(); // always set the first Sequence as reference
}
int AlignmentCIGAR :: checkAligned(){
    if (sequences.getSize()==1) {
        cout << "WARNING: File contains only one sequence." << endl;
        return 1;
    }
    int l = sequences[0].getSeqLength(); // l is the length of the first Sequence in the vector
    // for loop checking if all the Sequence objects have the same length
    for (int i = 1; i<sequences.getSize(); i++) {
        if (sequences[i].getSeqLength() != l) {
            return 0;
        }
    }
    return 1;
}

string AlignmentCIGAR :: cigarOneSeq(int i){
    string c;
    string temp; // for storing cigar segments
    int match = 1, insertion = 1, deletion = 1, soft = 1;
    int change = 0; // keep track of what was incremented last time: 1: match  2: insertion  3: deletion 4: soft clipping
    if (i == 0) {
        cout << "WARNING: Sequence is the reference." << endl;
        return c;
    }
    // append to cigar string when the state changes.
    for (int j = 0; j<ref.size(); j++) {
        ostringstream oss;
        // both sequence and reference have nucleotide
        if (sequences[i][j]!="-" && ref[j]!='-') {
            if (change == 1) {
                match ++; // increment
            }
            // all other cases, add to cigar string and change to match state
            else if (change == 2){
                oss << temp << insertion;
                c.append(oss.str());
                c.append("I");
                insertion = 1;
            }
            else if (change == 3){
                oss << temp << deletion;
                c.append(oss.str());
                c.append("D");
                deletion = 1;
            }
            else if (change == 4){
                oss << temp << soft;
                c.append(oss.str());
                c.append("S");
                soft = 1;
            }
            change = 1;
        }
        // sequence has nucleotide and reference has gap
        if(sequences[i][j]!="-" && ref[j]=='-'){
            if (change == 1) {
                oss << temp << match;
                c.append(oss.str());
                c.append("M");
                match = 1;
            }
            else if (change == 2){
                insertion ++; // increment
            }
            else if (change == 3){
                oss << temp << deletion;
                c.append(oss.str());
                c.append("D");
                deletion = 1;
            }
            else if (change == 4){
                oss << temp << soft;
                c.append(oss.str());
                c.append("S");
                soft = 1;
            }
            change = 2;
        }
        // sequence has gap and reference has nucleotide
        if (sequences[i][j]=="-" && ref[j]!='-') {
            if (change == 1) {
                oss << temp << match;
                c.append(oss.str());
                c.append("M");
                match = 1;
            }
            else if (change == 2){
                oss << temp << insertion;
                c.append(oss.str());
                c.append("I");
                insertion = 1;
            }
            else if (change == 3){
                deletion ++; // increment
            }
            else if (change == 4){
                oss << temp << soft;
                c.append(oss.str());
                c.append("S");
                soft = 1;
            }
            change = 3;
        }
        // both sequence and reference has gaps
        if (sequences[i][j]=="-" && ref[j]=='-') {
            if (change == 1) {
                oss << temp << match;
                c.append(oss.str());
                c.append("M");
                match = 1;
            }
            else if (change == 2){
                oss << temp << insertion;
                c.append(oss.str());
                c.append("I");
                insertion = 1;
            }
            else if (change == 3){
                oss << temp << deletion;
                c.append(oss.str());
                c.append("D");
                deletion = 1;
            }
            else if (change == 4){
                soft ++; // increment
            }
            change = 4;
        }
        temp.clear();
    }
    // the last case
    ostringstream oss;
    if (change == 1) {
        oss << temp << match;
        c.append(oss.str());
        c.append("M");
    }
    if (change == 2) {
        oss << temp << insertion;
        c.append(oss.str());
        c.append("I");
    }
    if (change == 3) {
        oss << temp << deletion;
        c.append(oss.str());
        c.append("D");
    }
    if (change == 4) {
        oss << temp << soft;
        c.append(oss.str());
        c.append("S");
    }
    c.append("\n");
    return c;
}
void AlignmentCIGAR :: setCigar(){
    for (int i = 1; i<sequences.getSize(); i++) {
        string c = cigarOneSeq(i);
        cigar.push_back(c);
    }
}
void AlignmentCIGAR :: printCigar(){
    cout << "reference name: " << sequences[0].getSeqName()<< endl;
    for (int i = 1; i<=cigar.size(); i++) {
        cout << sequences[i].getSeqName()<<"    "<<cigar[i-1]<<endl; // sequence name and cigar string tab seperated.
    }
}
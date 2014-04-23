//
//  ExtractSequence.cpp
//  
//
//  Created by Xuanyi Li on 14-3-18.
//
//

#include "ExtractSequence.h"
#include "Sequence.h"
#include "CodonFrequency.h"
#include "MinMax.h"
#include <iostream>
#include <string>
#include "string.h"
#include <sstream>
#include <fstream>	// outfile
#include <vector>
using namespace std;

ExtractSequence :: ExtractSequence(char *filename){
    ifstream file;
    file.open(filename);
    string sequence;
    string line;
    int i=0;
    // first line: always a header line
    getline(file,line);
    getHeader(line);
    while (!file.eof()) {
        getline(file,line);
        if (line[0]=='>') { // header line starts with ">"
            Sequence seq (headers[0],headers[1],sequence); // create a Sequence object with name description and sequence
            Sequences.push_back(seq); // add the Sequence object to the Sequences vector.
            sequence.clear();// clear the sequence string to start storing a new sequence
            headers.clear(); // clear header vector to start storing new name and description
            getHeader(line); // get new name and description

        } else {
            sequence += line;
        }
    }
    // store the last sequence in the file.
    Sequence seq (headers[0],headers[1],sequence);
    Sequences.push_back(seq);
    // If file is not file of alignments
    // Removes sequences that aren't proper length
    if (sequence.find("-") && sequence.find("M")) {
    	codonFreq = new CodonFrequency (filename, removeSeq(Sequences));
	   	minMax = new MinMax (filename, Sequences);
//      outputfileCF(filename);
      
	  } else {
	  	codonFreq = NULL;
	  	minMax = NULL;
	  }
	  
    file.close();
}

ExtractSequence :: ~ExtractSequence() {
	if (codonFreq)
		delete codonFreq;
	if (minMax)
		delete minMax;
}


// Does not initialize codonFreq for alignments
// Removes sequences that are divisible by 3
vector < Sequence > ExtractSequence::removeSeq(vector<Sequence> Sequences){

	for (int i = 0; i < Sequences.size(); i++) {
  	if (Sequences[i].getSeqLength()%3 != 0) {
   		Sequences.erase (Sequences.begin()+i);
   		cout << Sequences[i].getSeqName() << " " << Sequences[i].getSeqDescription() << " has been removed from codon frequency calculations due to improper length." << endl;
 		}
 	}
 	return Sequences;
}

/*
// Creates an output file with the inputfile name with .cf appended
// Output file contains codon frequencies from sequences
void ExtractSequence::outputfileCF(char *filename){

  string ofilename(filename);
	ofilename.append(".cf");

  ofstream ofile;
  ofile.open (ofilename.c_str());

  if (ofile.is_open()) {

  	codonFreq->outputFileCodonCount(ofile);
	  ofile.close();
	  
		cout << ofilename << " has been created." << endl;
	} else cout << "Unable to open " << ofilename << endl;
	
}
*/

void ExtractSequence :: printSequences(){

		// Does not print for alignments
/*		if (codonFreq) {
			codonFreq->printCodonCount();
			cout << endl;
			cout << "-----------------------" << endl;
		
		} else {
*/
		  for (int i = 0; i<Sequences.size(); i++) {
		      Sequences[i].print();
		      Sequences[i].printSeq();
		      cout << "-----------------------" << endl;
		  }
//		}
//    cout << "Number of sequences: " << Sequences.size() << endl;
}


// get name and description from header line
void ExtractSequence :: getHeader(string line){
    int i=0;
    int fasta = 0; // if it is a fasta like header -> 0
    string token; // store partitions in the string
    istringstream ss(line);
    // iterate through partitions of "|"
    while (getline(ss,token,'|')) {
        if (i==0) {
            if (token == ">gi") {  // if header start with ">gi": than it's a fasta style header, continue.
                i++;
                fasta = 1;
            }
            else{
                headers.push_back(token); // if it is not a fasta style header, store its first token as name, which is the first element of vector header
            }
        }
        if (i==1){
            if (fasta != 1) {  // continute b/c it's fasta header
                headers.push_back(token); // if it is not a fasta header, store token as second element of vector header.
            }
        }
        // for fasta header, i == 4 is the name and i==5 is the description
        if (i == 4) headers.push_back(token);				// 3rd token is ref number
        if (i == 5) headers.push_back(token);	// 4th token is seq description
        i++;
    }
}
Sequence ExtractSequence :: getSequence (string name){
    for (int i = 0; i<Sequences.size(); i++) {
        if (Sequences[i].getSeqName()==name) {
            return Sequences[i];
        }
    }
    // throw exception????
}
int ExtractSequence :: getSize(){
    return Sequences.size();
}
Sequence ExtractSequence :: operator[](int i){
    return Sequences[i];
}


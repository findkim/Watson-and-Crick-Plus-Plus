//
//  ExtractSequence.cpp
//  
//
//  Created by Xuanyi Li on 14-3-18.

//  Read in fasta file and store the sequences in a vector of sequences
//  get functions; add gaps to all sequences in the vector.

#include "ExtractSequence.h"
#include "Sequence.h"
#include "Domain.h"
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

    /* Sean here. Don't want to mess anything up for the end of the project, so adding the functionality
	where it should be but commented out

    addDomains(filename, Sequences);

    filename is a string, which is the filename of the txt file containing the info on the Domains. I'm assuming
    it is formatted the way the file Aaron sent me was. 
	
    */  

    file.close();
}
ExtractSequence :: ExtractSequence(){}

// Returns vector of sequences
vector<Sequence> ExtractSequence::getVectorOfSequences() {
	return Sequences;
}


void ExtractSequence :: printSequences() {

  for (int i = 0; i<Sequences.size(); i++) {
    Sequences[i].print();
    Sequences[i].printSeq();
    cout << "-----------------------" << endl;
  }
}


// get name and description from header line
void ExtractSequence :: getHeader(string line) {
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
}


int ExtractSequence :: getSize(){
    return Sequences.size();
}


Sequence ExtractSequence :: operator[](int i){
    return Sequences[i];
}

//To (most likely) be used in the constructor, takes a filename for the domains, and 
//creates a vector of domains, then adds them to a sequence if their ID's match
void ExtractSequence::addDomains(string filename,vector<Sequence> Seq_Vec)
{
	vector<Domain> domains2Add;

	ifstream file;
	file.open(filename.c_str());

	while(!file.eof())
	{
		string buffer;

		getline(file,buffer);
		
		char ID[30];
		char Type[10];
		int start,end;

		sscanf(buffer.c_str(),"%s %*[^:]:%s %*s %d-%d",ID,Type,&start,&end);
		
		string ID_str(ID);
		string Type_str(Type);
		Domain tempDom(ID,Type,start,end);
		
		domains2Add.push_back(tempDom);
	}

	//Not the most efficient way to do this, but I'll leave it up for future optimization if necessary.
	//Better way would be to use a map with ID as key, but that would require me changing things I'm not
	//comfortable with changing. -Sean
	for(vector<Sequence>::iterator seq_it= Seq_Vec.begin(); seq_it != Seq_Vec.end(); ++seq_it)
	{
		for(vector<Domain>::iterator dom_it= domains2Add.begin(); dom_it != domains2Add.end() ; ++dom_it)
		{
			if((*seq_it).getSeqName()==(*dom_it).getID())
			{
				(*seq_it).addDomain(*dom_it);
			}
		}
	}
	
}

void ExtractSequence :: addSequence(Sequence s){
    Sequences.push_back(s);
}
void ExtractSequence :: remove1Seq(int s){
    Sequences.erase(Sequences.begin()+s);
}
void ExtractSequence :: addGapstoAllEnds(int num){
    for (int i = 0; i<Sequences.size(); i++) {
        Sequences[i].addGaptoEnd(num);
    }
}
void ExtractSequence :: addGapstoAll(int num){
    for (int i = 0; i<Sequences.size(); i++) {
        Sequences[i].addGap(num);
    }
}


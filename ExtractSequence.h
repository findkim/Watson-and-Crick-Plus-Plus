//
//  ExtractSequence.h
//  
//
//  Created by Xuanyi Li on 14-3-18.
//
//

#ifndef ____ExtractSequence__
#define ____ExtractSequence__

#include <iostream>
#include <string>
#include <vector>
#include "Sequence.h"
#include "MinMax.h"

using namespace std;

class ExtractSequence{

	public:
    ExtractSequence(char *); // constructor takes in filename
    vector<Sequence> getVectorOfSequences();	// Returns the vector of sequences
		Sequence getSequence(string);// searching for Sequence by name
    void getHeader(string); // get name and description from header line
    int getSize(); // get the size of the Sequences vector
    void printSequences();
    Sequence operator[](int);
    
	private:
    vector<Sequence> Sequences; // vector of Sequence object
    vector<string> headers; // store name and description.

};
#endif /* defined(____ExtractSequence__) */

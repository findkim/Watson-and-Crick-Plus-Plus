//
//  Main.cpp
//  
//	Driver program for Sequence class
//	Reads file; parses file; prints name, description, and sequence
//
//  Created by Kim Ngo on 14-3-15.
//
//

#include "Sequence.h"
#include "ExtractSequence.h"
#include "AlignmentCIGAR.h"
#include "NeighborJoining.h"
#include "PairAlignment.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

int main() {
    // test loading file with multiple sequences
    ExtractSequence testSequences("ultimateORFS18608.fasta");
    PairAlignment testPair(testSequences);
    testPair.computeTables();
    testPair.updateAlign();
    testPair.printAlignment();
    //testSequences.printSequences(0);
    //cout << "++++++++++++++++++++++++++++++++++++" << endl;
    //testSequences.printSequences(1);
    //cout << "++++++++++++++++++++++++++++++++++++" << endl;
    //AlignmentCIGAR testCigar("nucSampleAlignment.fa");
    //string c = testCigar.cigarOneSeq(1);
    //cout << c << endl;
    //testCigar.setCigar();
    //testCigar.printCigar();
    //NeighborJoining testNJ("nucSampleAlignment.fa");
    //testNJ.write_distance();
    //testNJ.print_distance();
	return 0;
}

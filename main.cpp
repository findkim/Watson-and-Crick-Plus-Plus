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
#include "Alignment.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

int main() {
    // test loading file with multiple sequences
    //ExtractSequence testSequences("ultimateORFS18608.fasta");
    //Alignment testPair(testSequences);
    //testPair.computeTables();
    //testPair.updateAlign();
    //testPair.printAlignment();
    //testSequences.printSequences(0);
    //cout << "++++++++++++++++++++++++++++++++++++" << endl;
    //testSequences.printSequences(1);
    //cout << "++++++++++++++++++++++++++++++++++++" << endl;
//    AlignmentCIGAR testCigar("ultimateORFS18608.fasta");
    //string c = testCigar.cigarOneSeq(1);
    //cout << c << endl;
//    testCigar.setCigar();
//    testCigar.printCigar();
    //NeighborJoining testNJ("ultimateORFS18608.fasta");
    //testNJ.write_distance();
    //testNJ.print_distance();
//	ExtractSequence test("NM_000927.4.fas");
//	test.printSequences();
//    ExtractSequence testSequences("nucSampleAlignment.fa");
//    testSequences.printSequences();
  
  // testing to verify the frequency calculations
//  ExtractSequence codFreqTest("Ecol_test.fasta.short");
  ExtractSequence codFreqTest("Ecol_test.fasta");
//  ExtractSequence codFreqTest("Ypes.fasta");
//  ExtractSequence codFreqTest("Ecol_test.fasta.short.really");
//  codFreqTest.printSequences();

	return 0;
}

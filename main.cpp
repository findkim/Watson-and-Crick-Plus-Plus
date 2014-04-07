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
#include "CodonFrequency.h"
#include <iostream>

using namespace std;

int main() {
//	ExtractSequence test("NM_000927.4.fas");
//	test.printSequences();
//    ExtractSequence testSequences("nucSampleAlignment.fa");
//    testSequences.printSequences();
  
  // testing to verify the frequency calculations
//  ExtractSequence codFreqTest("Ecol_test.fasta.short");
//  ExtractSequence codFreqTest("Ecol_test.fasta");
//  ExtractSequence codFreqTest("Ypes.fasta");
  ExtractSequence codFreqTest("Ecol_test.fasta.short.really");
//  codFreqTest.printSequences();

	return 0;
}

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
#include "MinMax.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

int main() {
    // Alignment
    Alignment testAlign("ultimateORFS3710.fasta");
    testAlign.alignAll();
    testAlign.printAlignment();
    // CIGAR string
    AlignmentCIGAR testCigar("ultimateORFS3710_aligned.fasta");
    testCigar.setCigar();
    testCigar.printCigar();
  ExtractSequence codFreq("Ypes.fasta");
  MinMax calcMinMax("Ypes.fasta", codFreq.getVectorOfSequences());
	return 0;
}

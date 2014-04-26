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
    AlignmentCIGAR testCigar("ultimateORFS18608.fasta");
    //string c = testCigar.cigarOneSeq(1);
    //cout << c << endl;
    testCigar.setCigar();
    testCigar.printCigar();
    Alignment testAlign("ultimateORFS3710.fasta");
    testAlign.alignAll();
    testAlign.printAlignment();
    //testAlign.printTables(); // TEST
    
    /*string a1 = "ATCT--TGA";
    string a = "-A-C-T-G-";
    string b = "A--CT-G";
    string b1 = "TC-G-TA";
    cout << a << endl;
    cout << b << endl;
    cout << a1 << endl;
    cout << b1 << endl;
    for (int i = 0; i<a.size() && i<b.size(); i++) {
        if (a[i]=='-' && b[i]!='-') {
            b.insert(i,"-");
            b1.insert(i,"-");
        }
        if (a[i]!='-' && b[i]=='-') {
            a.insert(i,"-");
            a1.insert(i,"-");
        }
    }
    if (a.size()>b.size()) {
        int num = a.size()-b.size();
        for (int i = 0; i<num; i++) {
            b+="-";
            b1+="-";
        }
    }
    else{
        int num = b.size()-a.size();
        for (int i = 0; i<num; i++) {
            a+="-";
            a1+="-";
        }
    }
    cout << a << endl;
    cout << b << endl;
    cout << a1 << endl;
    cout << b1 << endl; */

    //NeighborJoining testNJ("ultimateORFS18608.fasta");
    //testNJ.write_distance();
    //testNJ.print_distance();
//	ExtractSequence test("NM_000927.4.fas");
//	test.printSequences();
//    ExtractSequence testSequences("nucSampleAlignment.fa");
//    testSequences.printSequences();
  
  // testing to verify the frequency calculations
//  ExtractSequence codFreqTest("Ecol_test.fasta.short");
//  ExtractSequence codFreqTest("Ecol_test.fasta");
  //ExtractSequence codFreqTest("Ypes.fasta");
//  ExtractSequence codFreqTest("Ecol_test.fasta.short.really");
//  codFreqTest.printSequences();

	return 0;
}

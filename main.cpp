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
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

using namespace std;

int main() {

	Sequence test("NM_000927.4.fas");
	test.print();
	test.printSeq();

	return 0;
}

//
//  Alignment.h
//  
//
//  Created by Xuanyi Li on 14-3-18.
//
//

#ifndef ____Alignment__
#define ____Alignment__

#include <iostream>
#include <string>
#include <vector>
#include "Sequence.h"
using namespace std;

class Alignment{
public:
    Alignment(char *);
    void printAlignment();
    /*void calAlignment();
    void displayAlignment(); */
private:
    vector<Sequence> AlignSeqs; // vector of Sequence object
};
#endif /* defined(____Alignment__) */

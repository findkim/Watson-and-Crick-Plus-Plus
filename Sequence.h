//
//  Sequence.h
//  
//
//  Created by Xuanyi Li on 14-3-8.
//
//

#ifndef ____Sequence__
#define ____Sequence__

#include <iostream>
#include <string>
using namespace std;

class Sequence{
public:
    Sequence();
    void addGap(int); // add gap after int
    char operator[](int);
private:
    string seq;
};
#endif /* defined(____Sequence__) */

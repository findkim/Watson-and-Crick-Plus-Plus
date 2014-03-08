//
//  Sequence.cpp
//  
//
//  Created by Xuanyi Li on 14-3-8.
//
//

#include "Sequence.h"
#include <iostream>
#include <string>
using namespace std;

Sequence :: Sequence(){
    
}
void Sequence :: addGap(int i){
    seq.insert(i+1,'-');
}
char Sequence :: operator[](int i){
    if (i<0 || i>seq.size()) {
        throw out_of_range("Subscript out of range");
    }
    return seq[i];
}
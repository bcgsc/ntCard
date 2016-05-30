#include "../src/hip_estimator.h"
#include <iostream>
#include <string>
#include <fstream>
#include <cstdio>

using namespace std;

void printBin(uint64_t n) {
    int count=0,count1=0;
    while (++count<=64) {
        if (n & ((uint64_t)1<<63)) {
            printf("1");
            ++count1;
        }
        else
            printf("0");

        n <<= 1;
    }
    printf("\n");
}

void test() {
    uint64_t hVal= 0x00F000000000FF0;
    printBin(hVal);
    cout << __builtin_clzll(hVal) << "\n";
    cout << __builtin_ctzll(hVal) << "\n";
    exit(0);
}

int main(int argc, const char *argv[]) {
    test();
    hip_estimator::hip_estimator<std::string> h;
    std::ifstream in(argv[1]);
    int k=atoi(argv[2]);
    string hseq,seq;
    while(getline(in,hseq)) {
        getline(in,seq);
        getline(in,hseq);
        getline(in,hseq);
        for(unsigned i=0; i<seq.length()-k+1; i++)
            h.insert(seq.substr(i,k));
    }
    std::cout << h.count() << std::endl;
    return 0;
}


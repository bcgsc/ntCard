#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
//#include "hip_estimator.h"
#include "nthash.hpp"

using namespace std;

namespace opt {
    unsigned threads=1;
    unsigned kmerLen=64;
    const unsigned nhash=4;
}

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

/*void hipTest(const char *inName, const int opt::kmerLen){
	hip_estimator::hip_estimator<std::string> h;
    std::ifstream in(inName);
	string hseq,seq;
    while(getline(in,hseq)) {
        getline(in,seq);
        getline(in,hseq);
        getline(in,hseq);
        for(unsigned i=0; i<seq.length()-opt::kmerLen+1;i++)
            h.insert(seq.substr(i,opt::kmerLen));
    }
    in.close();
    std::cout << h.count() << std::endl;
}*/

void nthTest(const char *inName) {
    std::ifstream in(inName);
    #pragma omp parallel
    {
        uint64_t mVec[opt::nhash];
        for (unsigned j=0; j<opt::nhash; j++) mVec[j]=0;
        uint64_t hVec[opt::nhash];
        bool good = true;
        for(string seq, hseq; good;) {
            #pragma omp critical(in)
            {
                good = getline(in, hseq);
                good = getline(in, seq);
                good = getline(in, hseq);
                good = getline(in, hseq);
            }
            if(good) {
                string kmer = seq.substr(0, opt::kmerLen);
                NTM64(kmer.c_str(), opt::kmerLen, opt::nhash, hVec); // initial hash vector
                if(hVec[0])
                for (unsigned j=0; j<opt::nhash; j++) {
                    //if(!hVec[j]) continue;
                    int run0 = __builtin_clzll(hVec[j]);
                    if(run0 > mVec[j]) mVec[j] = run0;
                }
                for (size_t i = 0; i < seq.length() - opt::kmerLen; i++) {
                    NTM64(seq[i], seq[i+opt::kmerLen], opt::kmerLen, opt::nhash, hVec); // consecutive hash vectors
                    if(hVec[0])
                    for (unsigned j=0; j<opt::nhash; j++) {
                        //if(!hVec[j]) continue;
                        int run0 = __builtin_clzll(hVec[j]);
                        if(run0 > mVec[j]) mVec[j] = run0;
                    }
                }
            }
        }
        double avg=0.0;
        for (unsigned j=0; j<opt::nhash; j++)
            avg+=mVec[j];
        std::cout << std::fixed << std::setprecision(11) << pow(2,avg/opt::nhash) << "\n";
        for (unsigned j=0; j<opt::nhash; j++)
            printf("%d\t",mVec[j]);//cerr << mVec[j] << "\n";
    }

    /*string hseq,seq;
    uint64_t hVec[opt::nhash]={0,0,0,0};
    int mVec[opt::nhash]={0,0,0,0};
    while(getline(in,hseq)) {
        getline(in,seq);
        getline(in,hseq);
        getline(in,hseq);

        string kmer = seq.substr(0, opt::kmerLen);
        NTM64(kmer.c_str(), opt::kmerLen, opt::nhash, hVec); // initial hash vector
        for (unsigned j=0; j<opt::nhash; j++) {
            if(!hVec[j]) continue;
            int run0 = __builtin_clzll(hVec[j]);
            // GreaterThanCAS
            while(true) {
                int oVal = mVec[j];
                if(run0<=oVal) break;
                if(__sync_bool_compare_and_swap(mVec+j,oVal,run0)) break;
            }
            
            //atomic
            //if(run0 > mVec[j]) mVec[j] = run0;
        }
        for (size_t i = 0; i < seq.length() - opt::kmerLen; i++) {
            NTM64(seq[i], seq[i+opt::kmerLen], opt::kmerLen, opt::nhash, hVec); // consecutive hash vectors
            for (unsigned j=0; j<opt::nhash; j++) {
                if(!hVec[j]) continue;
                int run0 = __builtin_clzll(hVec[j]);
                
                //if(run0 > mVec[j]) mVec[j] = run0;
                // GreaterThanCAS
                while(true) {
                    int oVal = mVec[j];
                    if(run0<=oVal) break;
                    if(__sync_bool_compare_and_swap(mVec+j,oVal,run0)) break;
                }
                
            }
        }
    }*/
    in.close();
    //int sum=0;
    //for (unsigned j=0; j<opt::nhash; j++)
        //sum+=mVec[j];
    //cout << sum/4.0 << "\n";
    //return sum;

}

int main(int argc, const char *argv[]){
    ///test();
    // hipTest(argv[1], atoi(argv[2]));
    nthTest(argv[1]);
    return 0;
}

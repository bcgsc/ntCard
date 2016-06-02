#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "nthash.hpp"

using namespace std;

namespace opt {
    unsigned threads=1;
    unsigned kmerLen=64;
    const unsigned nhash=7;
    const unsigned nbuck=512;
    const unsigned nbits=9;
}

int main(int argc, const char *argv[]){
    std::ifstream in(argv[1]);
    uint64_t tVec[opt::nbuck];
	for (unsigned j=0; j<opt::nbuck; j++)
		tVec[j]=0;
	
	size_t buffCard=0,buffSize=20000000;
	//vector<string> readBuffer(buffSize);
	string *readBuffer = new string [buffSize];
	bool good=true;
	#pragma omp parallel 
    {
		uint64_t mVec[opt::nbuck];
		for (unsigned j=0; j<opt::nbuck; j++)
			mVec[j]=0;
		
		
		while(good){
			
			#pragma omp single
			{			
				string seq, hseq;
				for(buffCard=0; buffCard<buffSize; buffCard++) {
					good = getline(in, hseq);
					good = getline(in, seq);
					good = getline(in, hseq);
					good = getline(in, hseq);
					if(good) readBuffer[buffCard]=seq;
					else break;					
				}			
				//std::cerr << buffCard << "\t";
			}
			int rInd;
			#pragma omp for
			for(rInd=0; rInd<buffCard; rInd++) {
				string seq = readBuffer[rInd];
				string kmer = seq.substr(0, opt::kmerLen);
				uint64_t hVal = NTP64(kmer.c_str(), opt::kmerLen);
				if(hVal>>opt::nbits) {
					int run0 = __builtin_clzll(hVal&(~(uint64_t)opt::nbuck-1));
					if(run0 > mVec[hVal&(opt::nbuck-1)]) mVec[hVal&(opt::nbuck-1)]=run0;
				}
				for (size_t i = 0; i < seq.length() - opt::kmerLen; i++) {
					hVal = NTP64(hVal, seq[i], seq[i+opt::kmerLen], opt::kmerLen);
					if(hVal&(~(uint64_t)opt::nbuck-1)) {
						int run0 = __builtin_clzll(hVal&(~(uint64_t)opt::nbuck-1));;
						if(run0 > mVec[hVal&(opt::nbuck-1)]) mVec[hVal&(opt::nbuck-1)]=run0;
					}
					
				}
			}
		}
		
		//bool good = true;
		/*for(string seq, hseq; good;) {
			#pragma omp critical(in)
			{
				good = getline(in, hseq);
				good = getline(in, seq);
				good = getline(in, hseq);
				good = getline(in, hseq);
			}*/
			/*if(good) {
				string kmer = seq.substr(0, opt::kmerLen);
				uint64_t hVal = NTP64(kmer.c_str(), opt::kmerLen);
				if(hVal>>opt::nbits) {
					int run0 = __builtin_clzll(hVal&(~(uint64_t)opt::nbuck-1));
					if(run0 > mVec[hVal&(opt::nbuck-1)]) mVec[hVal&(opt::nbuck-1)]=run0;
				}
				for (size_t i = 0; i < seq.length() - opt::kmerLen; i++) {
					hVal = NTP64(hVal, seq[i], seq[i+opt::kmerLen], opt::kmerLen);
					if(hVal&(~(uint64_t)opt::nbuck-1)) {
						int run0 = __builtin_clzll(hVal&(~(uint64_t)opt::nbuck-1));;
						if(run0 > mVec[hVal&(opt::nbuck-1)]) mVec[hVal&(opt::nbuck-1)]=run0;
					}
					
				}
			}*/
		//}
		#pragma omp critical
		{
			for (unsigned j=0; j<opt::nbuck; j++)
				if(tVec[j] < mVec[j]) tVec[j] = mVec[j];
		}
	}
	delete [] readBuffer;
    in.close();
	
	double pEst = 0.0, zEst = 0.0, eEst = 0.0, alpha = 0.0;

	if(opt::nbuck==16)	
		alpha = 2*0.673;
	else if (opt::nbuck == 32)
		alpha = 2*0.697;
	else if(opt::nbuck==64)
		alpha = 2*0.709;
	else
		alpha = 2*0.7213/(1 + 1.079/opt::nbuck);
	
	for (unsigned j=0; j<opt::nbuck-1; j++) 
		pEst += 1.0/((uint64_t)1<<tVec[j]);		
	zEst = 1.0/pEst;
	eEst = alpha * opt::nbuck * opt::nbuck * zEst;
	
	std::cout << (unsigned long long) eEst << "\n";
    
    return 0;
}

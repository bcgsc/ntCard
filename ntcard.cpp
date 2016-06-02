#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "nthash.hpp"

using namespace std;

namespace opt {
    unsigned threads=1;
    unsigned kmerLen=64;
    const unsigned nhash=7;
    const unsigned nbuck=512;
    const unsigned nbits=9;
}

void qSort(uint64_t a[], int l, int r) {
    long i = l-1, j = r, v = a[r];
    int tmp;
    if (r <= l) return;
    for (;;) {
        while (a[++i] > v);
        while (v > a[--j]) if (j == l) break;
        if (i >= j) break;
        tmp = a[i];
        a[i] = a[j];
        a[j]=tmp;
    }
    tmp = a[i];
    a[i] = a[r];
    a[r]=tmp;
    qSort(a, l, i-1);
    qSort(a, i+1, r);
}

void nthTest(const char *inName) {
    std::ifstream in(inName);
    uint64_t mVec[opt::nhash];
    for (unsigned j=0; j<opt::nhash; j++) mVec[j]=0;
    uint64_t hVec[opt::nhash];
    bool good = true;
    for(string seq, hseq; good;) {
        good = getline(in, hseq);
        good = getline(in, seq);
        good = getline(in, hseq);
        good = getline(in, hseq);
        if(good) {
            string kmer = seq.substr(0, opt::kmerLen);
            NTM64(kmer.c_str(), opt::kmerLen, opt::nhash, hVec); // initial hash vector
            if(hVec[0])
            for (unsigned j=0; j<opt::nhash; j++) {
                //int run0 = __builtin_clzll(hVec[j]);
                int run0 = std::max(__builtin_clzll(hVec[j]),__builtin_ctzll(hVec[j]));
                if(run0 > mVec[j]) mVec[j] = run0;
            }
            for (size_t i = 0; i < seq.length() - opt::kmerLen; i++) {
                NTM64(seq[i], seq[i+opt::kmerLen], opt::kmerLen, opt::nhash, hVec); // consecutive hash vectors
                if(hVec[0])
                for (unsigned j=0; j<opt::nhash; j++) {
                    //int run0 = __builtin_clzll(hVec[j]);
                    int run0 = std::max(__builtin_clzll(hVec[j]),__builtin_ctzll(hVec[j]));
                    if(run0 > mVec[j]) mVec[j] = run0;
                }
            }
        }
    }
    in.close();
    double avg=0.0;
    for (unsigned j=0; j<opt::nhash; j++)
        avg+=mVec[j];
    std::cout << std::fixed << std::setprecision(11) << pow(2,avg/opt::nhash) << "\n";
    for (unsigned j=0; j<opt::nhash; j++)
        cout << mVec[j] << "\n";
}

void nthTestBuff(const char *inName) {
    std::ifstream in(inName);
    uint64_t mVec[opt::nhash][opt::nbuck];
    for (unsigned i=0; i<opt::nhash; i++)
		for (unsigned j=0; j<opt::nbuck; j++)
			mVec[i][j]=0;
    uint64_t hVec[opt::nhash];
    bool good = true;
    for(string seq, hseq; good;) {
        good = getline(in, hseq);
        good = getline(in, seq);
        good = getline(in, hseq);
        good = getline(in, hseq);
        if(good) {
            string kmer = seq.substr(0, opt::kmerLen);
            NTM64(kmer.c_str(), opt::kmerLen, opt::nhash, hVec); // initial hash vector
            if(hVec[0])
            for (unsigned j=0; j<opt::nhash; j++) {
                int run0 = std::max(__builtin_clzll(hVec[j]),__builtin_ctzll(hVec[j]));
                if(run0 > mVec[j][hVec[j]&(opt::nbuck-1)]) mVec[j][hVec[j]&(opt::nbuck-1)] = run0;
            }
            for (size_t i = 0; i < seq.length() - opt::kmerLen; i++) {
                NTM64(seq[i], seq[i+opt::kmerLen], opt::kmerLen, opt::nhash, hVec); // consecutive hash vectors
                if(hVec[0])
                for (unsigned j=0; j<opt::nhash; j++) {
                    int run0 = std::max(__builtin_clzll(hVec[j]),__builtin_ctzll(hVec[j]));
                    if(run0 > mVec[j][hVec[j]&(opt::nbuck-1)]) mVec[j][hVec[j]&(opt::nbuck-1)] = run0;
                }
            }
        }
    }
    in.close();
    //double avg=0.0;
    //for (unsigned j=0; j<opt::nhash; j++)
        //avg+=mVec[j];
    //std::cout << std::fixed << std::setprecision(11) << pow(2,avg/opt::nhash) << "\n";
    //for (unsigned j=0; j<opt::nhash; j++)
        //cout << mVec[j] << "\n";
	cerr << "Done.\n";
	
	for (unsigned i=0; i<opt::nhash; i++)
		qSort(mVec[i],0,opt::nbuck-1);
	for (unsigned i=0; i<opt::nhash; i++) {
		for (unsigned j=0; j<4; j++)
			cout << mVec[i][j] << " "; 
		cout << "\n";
	}
}

void ntCardBuff(const char *inName) {
    std::ifstream in(inName);
    uint64_t tVec[opt::nbuck];
	for (unsigned j=0; j<opt::nbuck; j++)
		tVec[j]=0;
    #pragma omp parallel 
    {
	uint64_t mVec[opt::nbuck];
	for (unsigned j=0; j<opt::nbuck; j++)
		mVec[j]=0;
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
            uint64_t hVal = NTP64(kmer.c_str(), opt::kmerLen); // initial hash vector
            if(hVal>>opt::nbits) {
				int run0 = __builtin_clzll(hVal&(~(uint64_t)opt::nbuck-1));
				if(run0 > mVec[hVal&(opt::nbuck-1)]) mVec[hVal&(opt::nbuck-1)]=run0;
			}
            for (size_t i = 0; i < seq.length() - opt::kmerLen; i++) {
                hVal = NTP64(hVal, seq[i], seq[i+opt::kmerLen], opt::kmerLen); // consecutive hash vectors
				if(hVal&(~(uint64_t)opt::nbuck-1)) {
					int run0 = __builtin_clzll(hVal&(~(uint64_t)opt::nbuck-1));;
					if(run0 > mVec[hVal&(opt::nbuck-1)]) mVec[hVal&(opt::nbuck-1)]=run0;
				}
                
            }
        }
    }
    #pragma omp critical(out)
    {
	for (unsigned j=0; j<opt::nbuck; j++)
		if(tVec[j] < mVec[j]) tVec[j] = mVec[j];
	}
	}
    in.close();
    //double avg=0.0;
    //for (unsigned j=0; j<opt::nhash; j++)
        //avg+=mVec[j];
    //std::cout << std::fixed << std::setprecision(11) << pow(2,avg/opt::nhash) << "\n";
    //for (unsigned j=0; j<opt::nhash; j++)
        //cout << mVec[j] << "\n";
	
	//qSort(mVec,0,15);
	
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
	
	//std::cout << std::fixed << std::setprecision(10) << eEst << "\n";
	std::cout << (unsigned long long) eEst << "\n";
	
	//for (unsigned j=0; j<opt::nbuck-1; j++)
		//cout << mVec[j] << " "; 
	cout << "\n";
}

int main(int argc, const char *argv[]){
    //ntCardBuff(argv[1]);
    
    std::ifstream in(argv[1]);
    uint64_t tVec[opt::nbuck];
	for (unsigned j=0; j<opt::nbuck; j++)
		tVec[j]=0;
    #pragma omp parallel 
    {
	uint64_t mVec[opt::nbuck];
	for (unsigned j=0; j<opt::nbuck; j++)
		mVec[j]=0;
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
            uint64_t hVal = NTP64(kmer.c_str(), opt::kmerLen); // initial hash vector
            if(hVal>>opt::nbits) {
				int run0 = __builtin_clzll(hVal&(~(uint64_t)opt::nbuck-1));
				if(run0 > mVec[hVal&(opt::nbuck-1)]) mVec[hVal&(opt::nbuck-1)]=run0;
			}
            for (size_t i = 0; i < seq.length() - opt::kmerLen; i++) {
                hVal = NTP64(hVal, seq[i], seq[i+opt::kmerLen], opt::kmerLen); // consecutive hash vectors
				if(hVal&(~(uint64_t)opt::nbuck-1)) {
					int run0 = __builtin_clzll(hVal&(~(uint64_t)opt::nbuck-1));;
					if(run0 > mVec[hVal&(opt::nbuck-1)]) mVec[hVal&(opt::nbuck-1)]=run0;
				}
                
            }
        }
    }
    #pragma omp critical(out)
    {
	for (unsigned j=0; j<opt::nbuck; j++)
		if(tVec[j] < mVec[j]) tVec[j] = mVec[j];
	}
	}
    in.close();
    //double avg=0.0;
    //for (unsigned j=0; j<opt::nhash; j++)
        //avg+=mVec[j];
    //std::cout << std::fixed << std::setprecision(11) << pow(2,avg/opt::nhash) << "\n";
    //for (unsigned j=0; j<opt::nhash; j++)
        //cout << mVec[j] << "\n";
	
	//qSort(mVec,0,15);
	
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
	
	//std::cout << std::fixed << std::setprecision(10) << eEst << "\n";
	std::cout << (unsigned long long) eEst << "\n";
	
	//for (unsigned j=0; j<opt::nbuck-1; j++)
		//cout << mVec[j] << " "; 
	cout << "\n";
    
    return 0;
}

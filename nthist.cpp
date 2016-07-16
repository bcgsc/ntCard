#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <getopt.h>
#include <cassert>
#include <cmath>

#include "nthash.hpp"
//#include "murmur.hpp"
#include "Uncompress.h"

#ifdef _OPENMP
# include <omp.h>
#endif

#define PROGRAM "nthist"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.0.0 \n"
    "Written by Hamid Mohamadi.\n"
    "Copyright 2016 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... FILES...\n"
    "Estimates the number of k-mers in FILES(>=1).\n"
    "Accepatble file formats: fastq, fasta, sam, bam, gz, bz, zip.\n"
    "\n"
    " Options:\n"
    "\n"
    "  -t, --threads=N	use N parallel threads [1]\n"
    "  -k, --kmer=N	the length of kmer [64]\n"
    "  -c, --canonical	get the estimate for cannonical form\n"
    "      --help	display this help and exit\n"
    "      --version	output version information and exit\n"
    "\n"
    "Report bugs to hmohamadi@bcgsc.ca.\n";


using namespace std;

namespace opt {
unsigned nThrd=1;
unsigned nHash=7;
unsigned kmLen=64;
unsigned nBuck=65536;
unsigned nBits=16;
unsigned rBuck=4194304;
unsigned rBits=22	;
unsigned sBits=16;
size_t totKmer=0;
bool canon=false;
bool samH=true;
}

static const char shortopts[] = "t:k:b:s:hc";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads",	required_argument, NULL, 't' },
    { "kmer",	required_argument, NULL, 'k' },
    { "bit",	required_argument, NULL, 'b' },
    { "sit",	required_argument, NULL, 's' },
    { "hash",	required_argument, NULL, 'h' },
    { "canonical",	no_argument, NULL, 'c' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

unsigned getftype(std::ifstream &in, std::string &samSeq) {
    std::string hseq;
    getline(in,hseq);
    if(hseq[0]=='>') {
        return 1;
    }
    if(hseq[0]=='@') {
        if( (hseq[1]=='H'&& hseq[2]=='D') ||
                (hseq[1]=='S'&& hseq[2]=='Q') ||
                (hseq[1]=='R'&& hseq[2]=='G') ||
                (hseq[1]=='P'&& hseq[2]=='G') ||
                (hseq[1]=='C'&& hseq[2]=='O') ) {
            return 2;
        }
        else
            return 0;
    }
    opt::samH=false;
    samSeq=hseq;
    return 2;
}

inline void ntComp(const uint64_t hVal, uint8_t *mVec, uint8_t *m_Counter) {
    if(hVal&(~((uint64_t)opt::nBuck-1))) {
        uint8_t run0 = __builtin_clzll(hVal&(~((uint64_t)opt::nBuck-1)));
        if(run0==opt::sBits) {
            size_t shVal=hVal&(opt::rBuck-1);
            if(m_Counter[shVal]<255) ++m_Counter[shVal];
        }
        if(run0 > mVec[hVal&(opt::nBuck-1)]) mVec[hVal&(opt::nBuck-1)]=run0;
    }
}

void getEfq(size_t &parkCount, uint8_t *mVec, uint8_t *m_Counter, std::ifstream &in) {
    bool good = true;
    for(string seq, hseq; good;) {
        good = getline(in, seq);
        good = getline(in, hseq);
        good = getline(in, hseq);

        if(good && seq.length()>=opt::kmLen) {
            //parkCount+=seq.length()-opt::kmLen+1;

            string kmer;
            uint64_t hVal, fhVal=0, rhVal=0;

            /*bool hGood=false;
            //unsigned posN=1;
            size_t seqIndex=0;
            while(!hGood && seqIndex<seq.length()-opt::kmLen+1) {
                kmer = seq.substr(seqIndex, opt::kmLen);
                hGood = NTPC64(kmer.c_str(), opt::kmLen, fhVal, rhVal, hVal);
                seqIndex++;
            }
            if(hGood) ntComp(hVal,mVec,m_Counter);
            //for (unsigned i = seqIndex-1; i < seq.length() - opt::kmLen; i++) {*/
            unsigned seqIndex=0;
            while(seqIndex<seq.length()-opt::kmLen +1) {
				//bool rollFlag = true;
                bool rollFlag = (seqIndex==0)? false : true;
                
                if(seedTab[seq[seqIndex+opt::kmLen-1]]==seedN) {                    
                    seqIndex+=opt::kmLen;
                    rollFlag=false;
                }
                if(!rollFlag) {
                    hGood=false;                    
                    while(!hGood && seqIndex<seq.length()-opt::kmLen+1) {                    
                        kmer = seq.substr(seqIndex, opt::kmLen);
                        hGood = NTPC64(kmer.c_str(), opt::kmLen, fhVal, rhVal, hVal);
						seqIndex++;
                    }
                    if(hGood) ntComp(hVal,mVec,m_Counter);
                }
                else {
                    hVal=NTPC64(fhVal, rhVal, seq[seqIndex-1], seq[seqIndex-1+opt::kmLen], opt::kmLen);
                    seqIndex++;
                    ntComp(hVal,mVec,m_Counter);
                }
            }           
        }

        good = getline(in, hseq);
    }
}

/*void getEfq(size_t &parkCount, uint8_t *mVec, uint8_t *m_Counter, std::ifstream &in) {
    bool good = true;
    for(string seq, hseq; good;) {
        good = getline(in, seq);
        good = getline(in, hseq);
        good = getline(in, hseq);

        if(good && seq.length()>=opt::kmLen) {
            parkCount+=seq.length()-opt::kmLen+1;
            uint64_t hVal, fhVal=0, rhVal=0;
            string kmer = seq.substr(0, opt::kmLen);
            if(!opt::canon) hVal=NTP64(kmer.c_str(), opt::kmLen);
            else hVal=NTPC64(kmer.c_str(), opt::kmLen, fhVal, rhVal);
            ntComp(hVal,mVec,m_Counter);
            for (size_t i = 0; i < seq.length() - opt::kmLen; i++) {
                if(!opt::canon) hVal = NTP64(hVal, seq[i], seq[i+opt::kmLen], opt::kmLen);
                else hVal=NTPC64(fhVal, rhVal, seq[i], seq[i+opt::kmLen], opt::kmLen);
                ntComp(hVal,mVec,m_Counter);
            }
        }

        good = getline(in, hseq);
    }
}*/

/*void getEfq(size_t &parkCount, uint8_t *mVec, uint8_t *m_Counter, std::ifstream &in) {
    bool good = true;
    for(string seq, hseq; good;) {
        good = getline(in, seq);
        good = getline(in, hseq);
        good = getline(in, hseq);

        if(good && seq.length()>=opt::kmLen) {
            parkCount+=seq.length()-opt::kmLen+1;
            for (size_t i = 0; i < seq.length() - opt::kmLen +1 ; i++) {
				string kmer = seq.substr(i, opt::kmLen);
				if(opt::canon) getCanon(kmer);
				uint64_t hVal=MurmurHash64A(kmer.c_str(),opt::kmLen,0);
				ntComp(hVal,mVec,m_Counter);
            }
        }

        good = getline(in, hseq);
    }
}*/


void getEfa(size_t &parkCount, uint8_t *mVec, uint8_t *m_Counter, std::ifstream &in) {
    bool good = true;
    for(string seq, hseq; good;) {
        good = getline(in, seq);
        if(good && seq.length()>=opt::kmLen) {
            parkCount+=seq.length()-opt::kmLen+1;
            string kmer = seq.substr(0, opt::kmLen);
            uint64_t hVal, fhVal=0, rhVal=0;
            if(!opt::canon) hVal=NTP64(kmer.c_str(), opt::kmLen);
            else hVal=NTPC64(kmer.c_str(), opt::kmLen, fhVal, rhVal);
            ntComp(hVal,mVec,m_Counter);
            for (size_t i = 0; i < seq.length() - opt::kmLen; i++) {
                if(!opt::canon) hVal = NTP64(hVal, seq[i], seq[i+opt::kmLen], opt::kmLen);
                else hVal=NTPC64(fhVal, rhVal, seq[i], seq[i+opt::kmLen], opt::kmLen);
                ntComp(hVal,mVec,m_Counter);
            }
        }
        good = getline(in, hseq);
    }
}

void getEsm(size_t &parkCount, uint8_t *mVec, uint8_t *m_Counter, std::ifstream &in, const std::string &samSeq) {
    std::string samLine,seq;
    std::string s1,s2,s3,s4,s5,s6,s7,s8,s9,s11;
    if(opt::samH) {
        while(getline(in,samLine))
            if (samLine[0]!='@') break;
    }
    else
        samLine=samSeq;
    do {
        std::istringstream iss(samLine);
        iss>>s1>>s2>>s3>>s4>>s5>>s6>>s7>>s8>>s9>>seq>>s11;
        if(seq.length()>=opt::kmLen) {
            parkCount+=seq.length()-opt::kmLen+1;
            string kmer = seq.substr(0, opt::kmLen);
            uint64_t hVal, fhVal=0, rhVal=0;
            if(!opt::canon) hVal=NTP64(kmer.c_str(), opt::kmLen);
            else hVal=NTPC64(kmer.c_str(), opt::kmLen, fhVal, rhVal);
            ntComp(hVal,mVec,m_Counter);
            for (size_t i = 0; i < seq.length() - opt::kmLen; i++) {
                if(!opt::canon) hVal = NTP64(hVal, seq[i], seq[i+opt::kmLen], opt::kmLen);
                else hVal=NTPC64(fhVal, rhVal, seq[i], seq[i+opt::kmLen], opt::kmLen);
                ntComp(hVal,mVec,m_Counter);
            }
        }
    } while(getline(in,samLine));
}

inline unsigned popCnt(unsigned char x) {
    return ((0x876543210 >>
             (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >>
            ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2))
           & 0xf;
}

int main(int argc, char** argv) {
    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 't':
            arg >> opt::nThrd;
            break;
        case 'b':
            arg >> opt::nBits;
            break;
        case 's':
            arg >> opt::rBits;
            break;
        case 'h':
            arg >> opt::nHash;
            break;
        case 'k':
            arg >> opt::kmLen;
            break;
        case 'c':
            opt::canon=true;
            break;
        case OPT_HELP:
            std::cerr << USAGE_MESSAGE;
            exit(EXIT_SUCCESS);
        case OPT_VERSION:
            std::cerr << VERSION_MESSAGE;
            exit(EXIT_SUCCESS);
        }
        if (optarg != NULL && !arg.eof()) {
            std::cerr << PROGRAM ": invalid option: `-"
                      << (char)c << optarg << "'\n";
            exit(EXIT_FAILURE);
        }
    }
    if (argc - optind < 1) {
        std::cerr << PROGRAM ": missing arguments\n";
        die = true;
    }
    if (die) {
        std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }
    vector<string> inFiles;
    for (int i = optind; i < argc; ++i) {
        string file(argv[i]);
        if(file[0]=='@') {
            string inName;
            ifstream inList(file.substr(1,file.length()).c_str());
            while(getline(inList,inName))
                inFiles.push_back(inName);
        }
        else
            inFiles.push_back(file);
    }

    opt::nBuck = ((unsigned)1) << opt::nBits;

    opt::rBuck = ((unsigned)1) << opt::rBits;

    uint8_t *tVec = new uint8_t [opt::nBuck];
    for (unsigned j=0; j<opt::nBuck; j++) tVec[j]=0;

    uint8_t *t_Counter = new uint8_t [opt::rBuck];
    for(size_t i=0; i<opt::rBuck; i++) t_Counter[i]=0;

#ifdef _OPENMP
    omp_set_num_threads(opt::nThrd);
#endif

    #pragma omp parallel
    {
        size_t parkCount=0;
        uint8_t *mVec = new uint8_t [opt::nBuck];
        for (unsigned j=0; j<opt::nBuck; j++) mVec[j]=0;

        uint8_t *m_Counter = new uint8_t [opt::rBuck];
        for(size_t i=0; i<opt::rBuck; i++) m_Counter[i]=0;

        #pragma omp for schedule(dynamic) nowait
        for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
            std::ifstream in(inFiles[inFiles.size()-file_i-1].c_str());
            std::string samSeq;
            unsigned ftype = getftype(in,samSeq);
            if(ftype==0)
                getEfq(parkCount, mVec, m_Counter, in);
            else if (ftype==1)
                getEfa(parkCount, mVec, m_Counter, in);
            else if (ftype==2)
                getEsm(parkCount, mVec, m_Counter, in, samSeq);
            in.close();
        }
        #pragma omp critical(vmrg)
        {
            for (unsigned j=0; j<opt::nBuck; j++)
                if(tVec[j]<mVec[j])
                    tVec[j]=mVec[j];
        }
        delete [] mVec;

        #pragma omp critical(cmrg)
        {
            for(size_t i=0; i<opt::rBuck; i++) {
                if(t_Counter[i]+m_Counter[i]<255)
                    t_Counter[i] += m_Counter[i];
            }
        }
        delete [] m_Counter;

        #pragma omp atomic
        opt::totKmer+=parkCount;
    }

    double pEst = 0.0, zEst = 0.0, eEst = 0.0, alpha = 0.0;
    alpha = 1.4426/(1 + 1.079/opt::nBuck);
    if(opt::canon) alpha/=2;

    for (unsigned j=0; j<opt::nBuck; j++)
        pEst += 1.0/((uint64_t)1<<tVec[j]);
    zEst = 1.0/pEst;
    eEst = alpha * opt::nBuck * opt::nBuck * zEst;

    delete [] tVec;

    unsigned singlton=0,totton=0;
    unsigned kHist[256];
    for(size_t i=0; i<256; i++) kHist[i]=0;

    for(size_t i=0; i<opt::rBuck; i++) {
        ++kHist[t_Counter[i]];
        if(t_Counter[i]==1) ++singlton;
        if(t_Counter[i]) ++totton;
    }

    delete [] t_Counter;

    double eRatio = 1.0*singlton/((opt::rBuck-totton)*(log(opt::rBuck)-log(opt::rBuck-totton)));
    double sEst = eRatio*eEst;

    std::cout << "F0, HLL Exp# of distnt kmers(k=" << opt::kmLen << "): " << (unsigned long long) eEst << "\n";
    std::cout << "f1, HLL Exp# of unique kmers(k=" << opt::kmLen << "): " << (unsigned long long) sEst << "\n";
    std::cout << "F1, HLL Exct# of total kmers(k=" << opt::kmLen << "): " << opt::totKmer << "\n\n";

    /*for(size_t i=1; i<32; i++) {
    	double eRatio = 1.0*kHist[i]/((opt::rBuck-totton)*(log(opt::rBuck)-log(opt::rBuck-totton)));
    	double sEst = eRatio*eEst;
    	std::cout << i << ": " << (unsigned long long)sEst << "\n";
    }*/

    //size_t F0= totton * ((size_t)1<<16);
    //std::cout << "F0: " << F0 << "\n";

    for(size_t i=0; i<33; i++) {
        //double eRatio = 1.0*kHist[i]/((opt::rBuck-totton)*(log(opt::rBuck)-log(opt::rBuck-totton)));
        //double sEst = eRatio*F0;
        //std::cout << i << ": " << (unsigned long long)sEst << "\n";
        std::cout << i << ": " << kHist[i] << "\n";
    }
    std::cout << "F0_sample or X(!=0): " << totton << "\n";

    size_t X0=opt::rBuck-totton;

    double F0_INANC= (opt::rBits*log(2)-log(X0)) * 1.0* ((size_t)1<<(opt::sBits+opt::rBits));
    double f1_INANC = eRatio*F0_INANC;
    //double F0_kmerstr = 1.0*((size_t)1<<opt::sBits)*(log(X0)-opt::rBits*log(2))/(log(opt::rBuck-1)-opt::rBits*log(2));
    //std::cout << "F0_kmerstr: " << (unsigned long long)F0_kmerstr << "\n";
    std::cout << "F0: " << (unsigned long long)F0_INANC << "\n";
    std::cout << "f1: " << (unsigned long long)f1_INANC << "\n";



    return 0;
}

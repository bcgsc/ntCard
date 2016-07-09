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
#include "Uncompress.h"

#ifdef _OPENMP
# include <omp.h>
#endif

#define PROGRAM "ntcardm"

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
    "  -k, --kmer=N1,N2,N3,...	the lengths of kmers\n"
    "  -c, --canonical	get the estimate for cannonical form\n"
    "      --help	display this help and exit\n"
    "      --version	output version information and exit\n"
    "\n"
    "Report bugs to hmohamadi@bcgsc.ca.\n";


using namespace std;

namespace opt {
unsigned nThrd=1;
unsigned nBuck=65536;
unsigned nBits=16;
unsigned sBuck=4194304;
unsigned sBits=22;
unsigned nk=0;
bool canon=false;
bool samH=true;
string kmerSet;
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

inline void ntComp(const uint64_t hVal, uint8_t *mVec, uint8_t *m1_Counter, uint8_t *m2_Counter) {
    if(hVal&(~((uint64_t)opt::nBuck-1))) {
        uint8_t run0 = __builtin_clzll(hVal&(~((uint64_t)opt::nBuck-1)));
        if(run0==16) {
			size_t shVal=hVal&(opt::sBuck-1);
			if((m2_Counter[shVal / 8] & (1 << (7 - shVal % 8))) == 0) {
                if((m1_Counter[shVal / 8] & (1 << (7 - shVal % 8))) == 0)
                    m1_Counter[shVal/8] |= (1 << (7 - shVal % 8));
                else
                    m2_Counter[shVal/8] |= (1 << (7 - shVal % 8));
            }
        }
        if(run0 > mVec[hVal&(opt::nBuck-1)]) mVec[hVal&(opt::nBuck-1)]=run0;
    }
}

void getEfq(size_t *parkVec, uint8_t *mVec, uint8_t *m1_Counter, uint8_t *m2_Counter, std::ifstream &in, const std::vector<unsigned> &kList) {
    bool good = true;
    for(string seq, hseq; good;) {
        good = getline(in, seq);
        good = getline(in, hseq);
        good = getline(in, hseq);

        for(unsigned k_index=0; k_index<opt::nk; k_index++) {
            unsigned kmLen=kList[k_index];
            if(good && seq.length()>=kmLen) {
                parkVec[k_index]+=seq.length()-kmLen+1;
                uint64_t hVal, fhVal=0, rhVal=0;
                string kmer = seq.substr(0, kmLen);
                if(!opt::canon) hVal=NTP64(kmer.c_str(), kmLen);
                else hVal=NTPC64(kmer.c_str(), kmLen, fhVal, rhVal);
                ntComp(hVal,mVec+k_index*opt::nBuck,m1_Counter+k_index*((opt::sBuck + 7)/8),m2_Counter+k_index*((opt::sBuck + 7)/8));
                for (size_t i = 0; i < seq.length() - kmLen; i++) {
                    if(!opt::canon) hVal = NTP64(hVal, seq[i], seq[i+kmLen], kmLen);
                    else hVal=NTPC64(fhVal, rhVal, seq[i], seq[i+kmLen], kmLen);
                    ntComp(hVal,mVec+k_index*opt::nBuck,m1_Counter+k_index*((opt::sBuck + 7)/8),m2_Counter+k_index*((opt::sBuck + 7)/8));
                }
            }
        }

        good = getline(in, hseq);
    }
}

void getEfa(size_t *parkVec, uint8_t *mVec, uint8_t *m1_Counter, uint8_t *m2_Counter, std::ifstream &in, const std::vector<unsigned> &kList) {
    bool good = true;
    for(string seq, hseq; good;) {
        good = getline(in, seq);

        for(unsigned k_index=0; k_index<opt::nk; k_index++) {
            unsigned kmLen=kList[k_index];
            if(good && seq.length()>=kmLen) {
                parkVec[k_index]+=seq.length()-kmLen+1;
                uint64_t hVal, fhVal=0, rhVal=0;
                string kmer = seq.substr(0, kmLen);
                if(!opt::canon) hVal=NTP64(kmer.c_str(), kmLen);
                else hVal=NTPC64(kmer.c_str(), kmLen, fhVal, rhVal);
                ntComp(hVal,mVec+k_index*opt::nBuck,m1_Counter+k_index*((opt::sBuck + 7)/8),m2_Counter+k_index*((opt::sBuck + 7)/8));
                for (size_t i = 0; i < seq.length() - kmLen; i++) {
                    if(!opt::canon) hVal = NTP64(hVal, seq[i], seq[i+kmLen], kmLen);
                    else hVal=NTPC64(fhVal, rhVal, seq[i], seq[i+kmLen], kmLen);
                    ntComp(hVal,mVec+k_index*opt::nBuck,m1_Counter+k_index*((opt::sBuck + 7)/8),m2_Counter+k_index*((opt::sBuck + 7)/8));
                }
            }
        }

        good = getline(in, hseq);
    }
}

void getEsm(size_t *parkVec, uint8_t *mVec, uint8_t *m1_Counter, uint8_t *m2_Counter, std::ifstream &in, const std::vector<unsigned> &kList, const std::string &samSeq) {
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
        for(unsigned k_index=0; k_index<opt::nk; k_index++) {
            unsigned kmLen=kList[k_index];
            if(seq.length()>=kmLen) {
                parkVec[k_index]+=seq.length()-kmLen+1;
                uint64_t hVal, fhVal=0, rhVal=0;
                string kmer = seq.substr(0, kmLen);
                if(!opt::canon) hVal=NTP64(kmer.c_str(), kmLen);
                else hVal=NTPC64(kmer.c_str(), kmLen, fhVal, rhVal);
                ntComp(hVal,mVec+k_index*opt::nBuck,m1_Counter+k_index*((opt::sBuck + 7)/8),m2_Counter+k_index*((opt::sBuck + 7)/8));
                for (size_t i = 0; i < seq.length() - kmLen; i++) {
                    if(!opt::canon) hVal = NTP64(hVal, seq[i], seq[i+kmLen], kmLen);
                    else hVal=NTPC64(fhVal, rhVal, seq[i], seq[i+kmLen], kmLen);
                    ntComp(hVal,mVec+k_index*opt::nBuck,m1_Counter+k_index*((opt::sBuck + 7)/8),m2_Counter+k_index*((opt::sBuck + 7)/8));
                }
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
            arg >> opt::sBits;
            break;
        case 'k':
            arg >> opt::kmerSet;
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

    std::istringstream ksstrm(opt::kmerSet);
    vector<unsigned> kList;
    unsigned kmerVal;
    while(ksstrm>>kmerVal) {
        kList.push_back(kmerVal);
        if(ksstrm.peek()==',')ksstrm.ignore();
        opt::nk++;
    }

    opt::nBuck = ((unsigned)1) << opt::nBits;
    opt::sBuck = ((unsigned)1) << opt::sBits;

    size_t *totkVec= new size_t [opt::nk];
    for(unsigned i=0; i<opt::nk; i++) totkVec[i]=0;

    uint8_t *tVec = new uint8_t [opt::nk*opt::nBuck];
    for(unsigned i=0; i<opt::nk; i++)
        for(unsigned j=0; j<opt::nBuck; j++)
            tVec[i*opt::nBuck+j]=0;

    uint8_t *t1_Counter = new uint8_t [opt::nk*((opt::sBuck + 7)/8)];
    uint8_t *t2_Counter = new uint8_t [opt::nk*((opt::sBuck + 7)/8)];
    for(unsigned i=0; i<opt::nk; i++)
        for(size_t j=0; j<(opt::sBuck + 7)/8; j++)
            t1_Counter[i*((opt::sBuck + 7)/8)+j]=t2_Counter[i*((opt::sBuck + 7)/8)+j]=0;


#ifdef _OPENMP
    omp_set_num_threads(opt::nThrd);
#endif

    #pragma omp parallel
    {
        size_t *parkVec= new size_t [opt::nk];
        for(unsigned i=0; i<opt::nk; i++) parkVec[i]=0;
        uint8_t *mVec = new uint8_t [opt::nk*opt::nBuck];
        for(unsigned i=0; i<opt::nk; i++)
            for(unsigned j=0; j<opt::nBuck; j++)
                mVec[i*opt::nBuck+j]=0;

        uint8_t *m1_Counter = new uint8_t [opt::nk*((opt::sBuck + 7)/8)];
        uint8_t *m2_Counter = new uint8_t [opt::nk*((opt::sBuck + 7)/8)];
        for(unsigned i=0; i<opt::nk; i++)
            for(size_t j=0; j<(opt::sBuck + 7)/8; j++)
                m1_Counter[i*((opt::sBuck + 7)/8)+j]=m2_Counter[i*((opt::sBuck + 7)/8)+j]=0;

        #pragma omp for schedule(dynamic) nowait
        for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
            std::ifstream in(inFiles[inFiles.size()-file_i-1].c_str());
            std::string samSeq;
            unsigned ftype = getftype(in,samSeq);

            if(ftype==0)
                getEfq(parkVec, mVec, m1_Counter, m2_Counter, in, kList);
            else if (ftype==1)
                getEfa(parkVec, mVec, m1_Counter, m2_Counter, in, kList);
            else if (ftype==2)
                getEsm(parkVec, mVec, m1_Counter, m2_Counter, in, kList, samSeq);
            in.close();
        }
        #pragma omp critical(vmrg)
        {
            for(unsigned i=0; i<opt::nk; i++)
                for (unsigned j=0; j<opt::nBuck; j++)
                    if(tVec[i*opt::nBuck+j]<mVec[i*opt::nBuck+j])
                        tVec[i*opt::nBuck+j]=mVec[i*opt::nBuck+j];
        }
        delete [] mVec;

        #pragma omp critical(cmrg)
        {
            for(unsigned i=0; i<opt::nk; i++)
                for(size_t j=0; j<(opt::sBuck + 7)/8; j++) {
                    t2_Counter[i*((opt::sBuck + 7)/8)+j] |= m2_Counter[i*((opt::sBuck + 7)/8)+j] | (m1_Counter[i*((opt::sBuck + 7)/8)+j]&t1_Counter[i*((opt::sBuck + 7)/8)+j]);
                    t1_Counter[i*((opt::sBuck + 7)/8)+j] |= m1_Counter[i*((opt::sBuck + 7)/8)+j];
                }
        }
        delete [] m1_Counter;
        delete [] m2_Counter;

        #pragma omp critical(tmrg)
        {
            for(unsigned i=0; i<opt::nk; i++) {
                totkVec[i]+=parkVec[i];
            }
        }
        delete [] parkVec;
    }

    std::cout << "Q\tk\tF0\tf1\tF1\n";
    for(unsigned i=0; i<opt::nk; i++)
    {
        double pEst = 0.0, zEst = 0.0, eEst = 0.0, alpha = 0.0;
        alpha = 1.4426/(1 + 1.079/opt::nBuck);
        if(opt::canon) alpha/=2;

        for (unsigned j=0; j<opt::nBuck; j++)
            pEst += 1.0/((uint64_t)1<<tVec[i*opt::nBuck+j]);
        zEst = 1.0/pEst;
        eEst = alpha * opt::nBuck * opt::nBuck * zEst;

        unsigned singlton=0,totton=0;
        for(size_t j=0; j<(opt::sBuck + 7)/8; j++) {
            singlton+=popCnt(t1_Counter[i*((opt::sBuck + 7)/8)+j]&(~t2_Counter[i*((opt::sBuck + 7)/8)+j]));
            totton+=popCnt(t1_Counter[i*((opt::sBuck + 7)/8)+j]|t2_Counter[i*((opt::sBuck + 7)/8)+j]);
        }

        double eRatio = 1.0*singlton/((opt::sBuck-totton)*(log(opt::sBuck)-log(opt::sBuck-totton)));
        double sEst = eRatio*eEst;

        std::cout<< "0\t"<< kList[i] << "\t" <<(unsigned long long) eEst << "\t" << (unsigned long long) sEst << "\t" << totkVec[i] << "\n";
    }
    delete [] tVec;
    delete [] t1_Counter;
    delete [] t2_Counter;
    delete [] totkVec;

    return 0;
}

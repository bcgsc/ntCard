#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <getopt.h>
#include <cassert>

#include "nthash.hpp"
#include "Uncompress.h"

#define PROGRAM "ntcard"

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
bool canon=false;
bool samH=true;
}

static const char shortopts[] = "t:k:b:hc";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads",	required_argument, NULL, 't' },
    { "kmer",	required_argument, NULL, 'k' },
    { "bit",	required_argument, NULL, 'b' },
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

void getEfq(uint8_t *mVec, std::ifstream &in) {
    bool good = true;
    for(string seq, hseq; good;) {
        good = getline(in, seq);
        good = getline(in, hseq);
        good = getline(in, hseq);
        if(good && seq.length()>=opt::kmLen) {
            string kmer = seq.substr(0, opt::kmLen);
            uint64_t hVal, fhVal=0, rhVal=0;
            if(!opt::canon) hVal=NTP64(kmer.c_str(), opt::kmLen);
            else hVal=NTPC64(kmer.c_str(), opt::kmLen, fhVal, rhVal);
            if(hVal&(~((uint64_t)opt::nBuck-1))) {
                uint8_t run0 = __builtin_clzll(hVal&(~((uint64_t)opt::nBuck-1)));
                if(run0 > mVec[hVal&(opt::nBuck-1)]) mVec[hVal&(opt::nBuck-1)]=run0;
            }
            for (size_t i = 0; i < seq.length() - opt::kmLen; i++) {
                if(!opt::canon) hVal = NTP64(hVal, seq[i], seq[i+opt::kmLen], opt::kmLen);
                else hVal=NTPC64(fhVal, rhVal, seq[i], seq[i+opt::kmLen], opt::kmLen);
                if(hVal&(~((uint64_t)opt::nBuck-1))) {
                    uint8_t run0 = __builtin_clzll(hVal&(~((uint64_t)opt::nBuck-1)));
                    if(run0 > mVec[hVal&(opt::nBuck-1)]) mVec[hVal&(opt::nBuck-1)]=run0;
                }
            }
        }
        good = getline(in, hseq);
    }
}

void getEfa(uint8_t *mVec, std::ifstream &in) {
    bool good = true;
    for(string seq, hseq; good;) {
        good = getline(in, seq);
        if(good && seq.length()>=opt::kmLen) {
            string kmer = seq.substr(0, opt::kmLen);
            uint64_t hVal, fhVal=0, rhVal=0;
            if(!opt::canon) hVal=NTP64(kmer.c_str(), opt::kmLen);
            else hVal=NTPC64(kmer.c_str(), opt::kmLen, fhVal, rhVal);
            if(hVal&(~((uint64_t)opt::nBuck-1))) {
                uint8_t run0 = __builtin_clzll(hVal&(~((uint64_t)opt::nBuck-1)));
                if(run0 > mVec[hVal&(opt::nBuck-1)]) mVec[hVal&(opt::nBuck-1)]=run0;
            }
            for (size_t i = 0; i < seq.length() - opt::kmLen; i++) {
                if(!opt::canon) hVal = NTP64(hVal, seq[i], seq[i+opt::kmLen], opt::kmLen);
                else hVal=NTPC64(fhVal, rhVal, seq[i], seq[i+opt::kmLen], opt::kmLen);
                if(hVal&(~((uint64_t)opt::nBuck-1))) {
                    uint8_t run0 = __builtin_clzll(hVal&(~((uint64_t)opt::nBuck-1)));
                    if(run0 > mVec[hVal&(opt::nBuck-1)]) mVec[hVal&(opt::nBuck-1)]=run0;
                }
            }
        }
        good = getline(in, hseq);
    }
}

void getEsm(uint8_t *mVec, std::ifstream &in, const std::string &samSeq) {
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
            string kmer = seq.substr(0, opt::kmLen);
            uint64_t hVal, fhVal=0, rhVal=0;
            if(!opt::canon) hVal=NTP64(kmer.c_str(), opt::kmLen);
            else hVal=NTPC64(kmer.c_str(), opt::kmLen, fhVal, rhVal);
            if(hVal&(~((uint64_t)opt::nBuck-1))) {
                uint8_t run0 = __builtin_clzll(hVal&(~((uint64_t)opt::nBuck-1)));
                if(run0 > mVec[hVal&(opt::nBuck-1)]) mVec[hVal&(opt::nBuck-1)]=run0;
            }
            for (size_t i = 0; i < seq.length() - opt::kmLen; i++) {
                if(!opt::canon) hVal = NTP64(hVal, seq[i], seq[i+opt::kmLen], opt::kmLen);
                else hVal=NTPC64(fhVal, rhVal, seq[i], seq[i+opt::kmLen], opt::kmLen);
                if(hVal&(~((uint64_t)opt::nBuck-1))) {
                    uint8_t run0 = __builtin_clzll(hVal&(~((uint64_t)opt::nBuck-1)));
                    if(run0 > mVec[hVal&(opt::nBuck-1)]) mVec[hVal&(opt::nBuck-1)]=run0;
                }
            }
        }
    } while(getline(in,samLine));
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
        inFiles.push_back(file);
    }

    opt::nBuck = ((unsigned)1) << opt::nBits;
    uint8_t *mVec = new uint8_t [opt::nBuck];
    for (unsigned j=0; j<opt::nBuck; j++) mVec[j]=0;

    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
        std::ifstream in(inFiles[file_i].c_str());
        std::string samSeq;
        unsigned ftype = getftype(in,samSeq);
        if(ftype==0)
            getEfq(mVec, in);
        else if (ftype==1)
            getEfa(mVec, in);
        else if (ftype==2)
            getEsm(mVec, in, samSeq);
        in.close();
    }

    double pEst = 0.0, zEst = 0.0, eEst = 0.0, alpha = 0.0;
    alpha = 1.4426/(1 + 1.079/opt::nBuck);
    for (unsigned j=0; j<opt::nBuck-1; j++)
        pEst += 1.0/((uint64_t)1<<mVec[j]);
    zEst = 1.0/pEst;
    eEst = alpha * opt::nBuck * opt::nBuck * zEst;
    delete [] mVec;

    std::cout << "Exp# of kmers(k=" << opt::kmLen << "): " << (unsigned long long) eEst << "\n";
    return 0;
}

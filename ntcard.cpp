#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <getopt.h>

#include "nthash.hpp"
#include "Uncompress.h"


#define PROGRAM "ntcard"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.0.0 \n"
    "Written by Hamid Mohamadi.\n"
    "Copyright 2016 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... FILES...\n"
    "Estimates the number of k-mers in FILES (>=1) seprated by space.\n"
    "\n"
    " Options:\n"
    "\n"
    "  -t, --threads=N	use N parallel threads [1]\n"
    "  -k, --kmer=N		the length of kmer [64]\n"
    "  -f, --fastq	input files in FASTA format (default:FASTQ)\n"
    "  -c, --canonical	get the estimate for cannonical form\n"
    "      --help	display this help and exit\n"
    "      --version	output version information and exit\n"
    "\n"
    "Report bugs to hmohamadi@bcgsc.ca.\n";


using namespace std;

namespace opt {
unsigned nThrd=1;
unsigned nHash=7;
unsigned kmLen=128;
unsigned nBuck=65536;
unsigned nBits=16;
static bool fasta=false;
static bool canon=false;
}

static const char shortopts[] = "t:k:b:h:fc";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads",	required_argument, NULL, 't' },
    { "kmer",	required_argument, NULL, 'k' },
    { "bit",	required_argument, NULL, 'b' },
    { "hash",	required_argument, NULL, 'h' },
    { "fasta",	no_argument, NULL, 'f' },
    { "canonical",	no_argument, NULL, 'c' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

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
        case 'f':
            opt::fasta=true;
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
        bool good = true;
        for(string seq, hseq; good;) {
            good = getline(in, hseq);
            good = getline(in, seq);
            if(!opt::fasta) {
                good = getline(in, hseq);
                good = getline(in, hseq);
            }
            if(good) {
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
        }
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

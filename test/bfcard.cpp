/*
 *
 * nttest.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */


#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <unistd.h>
#include <getopt.h>
#include "seqgen.hpp"
#include "BloomFilter.hpp"
#include "nthash.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "bfCard"

static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.0.0 \n"
    "Written by Hamid Mohamadi.\n"
    "Copyright 2016 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... QUERY\n"
    "Report bugs to hmohamadi@bcgsc.ca.\n";

namespace opt {
unsigned threads=1;
unsigned kmerLen=64;
unsigned ibits = 8;
unsigned nhash=1;
unsigned nz;
size_t sgene;
unsigned method;
bool fastq = false;
bool inpFlag = false;
bool uniformity = false;
}

using namespace std;

static const char shortopts[] = "k:b:h:j:g:m:a:i:u:f:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads",	required_argument, NULL, 'j' },
    { "kmer",	required_argument, NULL, 'k' },
    { "bit",	required_argument, NULL, 'b' },
    { "hash",	required_argument, NULL, 'h' },
    { "tlen",	required_argument, NULL, 'g' },
    { "input",	no_argument, NULL, 'i' },
    { "uniformity",	no_argument, NULL, 'u' },
    { "fastq",	no_argument, NULL, 'f' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

static const string itm[]= {"nthash","city","murmur","xxhash","ntbase"};

static const unsigned char b2r[256] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //0
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //1
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //2
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //3
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', //4   'A' 'C' 'G'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', //5   'T'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', //6   'a' 'c' 'g'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', //7   't'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //8
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //9
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //10
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //11
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //12
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //13
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //14
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //15
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};

void getCanon(std::string &bMer) {
    int p=0, hLen=(opt::kmerLen-1)/2;
    while (bMer[p] == b2r[(unsigned char)bMer[opt::kmerLen-1-p]]) {
        ++p;
        if(p>=hLen) break;
    }
    if (bMer[p] > b2r[(unsigned char)bMer[opt::kmerLen-1-p]]) {
        for (int lIndex = p, rIndex = opt::kmerLen-1-p; lIndex<=rIndex; ++lIndex,--rIndex) {
            char tmp = b2r[(unsigned char)bMer[rIndex]];
            bMer[rIndex] = b2r[(unsigned char)bMer[lIndex]];
            bMer[lIndex] = tmp;
        }
    }
}

void loadSeq(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        myFilter.insertF(kmer.c_str());
    }
}

void loadSeqr(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;

    string kmer = seq.substr(0,opt::kmerLen);
    uint64_t fhVal;
    myFilter.insertF(kmer.c_str(), fhVal);
    for (size_t i = 1; i < seq.size() - opt::kmerLen + 1; i++) {
        myFilter.insertF(fhVal, seq[i-1], seq[i+opt::kmerLen-1]);
    }
}

void loadSeqm(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        myFilter.insertMur(kmer.c_str());
    }
}

void loadSeqc(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        myFilter.insertCit(kmer.c_str());
    }
}

void loadSeqx(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        myFilter.insertXxh(kmer.c_str());
    }
}

void loadBf(BloomFilter &myFilter, const char* faqFile) {
    ifstream uFile(faqFile);
    bool good = true;
    #pragma omp parallel
    for(string line, hline; good;) {
        #pragma omp critical(uFile)
        {
            good = getline(uFile, hline);
            good = getline(uFile, line);
			good = getline(uFile, hline);
			good = getline(uFile, hline);
        }
        if(good)
			loadSeqr(myFilter, line);
    }
    uFile.close();
}

void nthashBF(const char *geneName) {
#ifdef _OPENMP
    omp_set_num_threads(opt::threads);
#endif
	opt::method = 0;
    std::cerr<<"#threads="<<opt::threads << "\n";
	std::cerr<<"method="<<itm[opt::method]<<"\n";
	std::cerr<<"kmerl="<<opt::kmerLen<<"\n";
	std::cerr<<"nhash="<<opt::nhash<<"\n";
	double sTime = omp_get_wtime();
	BloomFilter myFilter(opt::ibits*opt::sgene , opt::nhash, opt::kmerLen);
	loadBf(myFilter, geneName);
	cerr << "|popBF|=" << myFilter.getPop() << "\n";
	cerr << "load_time=" <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
}

int main(int argc, char** argv) {
    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
        case '?':
            die = true;
            break;
        case 'j':
            arg >> opt::threads;
            break;
        case 'b':
            arg >> opt::ibits;
            break;
        case 'g':
            arg >> opt::sgene;
            break;
        case 'h':
            arg >> opt::nhash;
            break;
        case 'k':
            arg >> opt::kmerLen;
            break;
        case 'a':
            arg >> opt::method;
            break;
        case 'i':
            opt::inpFlag=true;
            break;
        case 'f':
            opt::fastq=true;
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
    if (argc - optind != 1) {
        std::cerr << PROGRAM ": missing arguments\n";
        die = true;
    }

    if (die) {
        std::cerr << "Try `" << PROGRAM
                  << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }

	std::cerr<<"bit/i="<<opt::ibits<<"\n";
	std::cerr<<"sgene="<<opt::sgene<<"\n\n";

	const char *geneName(argv[argc-1]);
	nthashBF(geneName);

    return 0;
}

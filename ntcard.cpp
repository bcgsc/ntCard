#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <getopt.h>
#include <cassert>
#include <cmath>

#include "ntHashIterator.hpp"
#include "Uncompress.h"

#ifdef _OPENMP
# include <omp.h>
#endif

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
    "  -c, --cov=N	the maximum coverage of kmer in output [64]\n"

    "      --help	display this help and exit\n"
    "      --version	output version information and exit\n"
    "\n"
    "Report bugs to hmohamadi@bcgsc.ca.\n";


using namespace std;

namespace opt {
unsigned nThrd=1;
unsigned kmLen=64;
unsigned rBuck=33554432;
unsigned rBits=27;
unsigned sBits=11;
unsigned sMask=1023;
unsigned covMax=64;
unsigned nSamp=2;
bool samH=true;
}

static const char shortopts[] = "t:s:r:k:c:f:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads",	required_argument, NULL, 't' },
    { "kmer",	required_argument, NULL, 'k' },
    { "cov",	required_argument, NULL, 'c' },
    { "rbit",	required_argument, NULL, 'r' },
    { "sbit",	required_argument, NULL, 's' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

size_t getInf(const char* inFile) {
    std::ifstream in(inFile, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

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

inline void ntComp(const uint64_t hVal, uint16_t *t_Counter) {
    uint64_t indBit=opt::nSamp;
    if(hVal>>(63-opt::sBits) == 1) indBit=0;
    if(hVal>>(64-opt::sBits) == opt::sMask) indBit=1;
    if(indBit < opt::nSamp) {
        size_t shVal=hVal&(opt::rBuck-1);
        #pragma omp atomic
        ++t_Counter[indBit*opt::rBuck+shVal];
    }

}

inline void ntRead(const string &seq, uint16_t *t_Counter, size_t &totKmer) {
    ntHashIterator itr(seq, opt::kmLen);
    while (itr != itr.end()) {
        ntComp((*itr),t_Counter);
        ++itr;
        ++totKmer;
    }
}

void getEfq(std::ifstream &in, uint16_t *t_Counter, size_t &totKmer) {
    bool good = true;
    for(string seq, hseq; good;) {
        good = getline(in, seq);
        good = getline(in, hseq);
        good = getline(in, hseq);
        if(good && seq.length()>=opt::kmLen)
            ntRead(seq, t_Counter, totKmer);
        good = getline(in, hseq);
    }
}

void getEfa(std::ifstream &in, uint16_t *t_Counter, size_t &totKmer) {
    bool good = true;
    for(string seq, hseq; good;) {
        string line;
        good = getline(in, seq);
        while(good&&seq[0]!='>') {
            line+=seq;
            good = getline(in, seq);
        }
        if(line.length()>=opt::kmLen)
            ntRead(line, t_Counter, totKmer);
    }
}

void getEsm(std::ifstream &in, const std::string &samSeq, uint16_t *t_Counter, size_t &totKmer) {
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
        if(seq.length()>=opt::kmLen)
            ntRead(seq, t_Counter, totKmer);
    } while(getline(in,samLine));
}

int main(int argc, char** argv) {

    double sTime = omp_get_wtime();

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
        case 's':
            arg >> opt::sBits;
            break;
        case 'r':
            arg >> opt::rBits;
            break;
        case 'c':
            arg >> opt::covMax;
            break;
        case 'k':
            arg >> opt::kmLen;
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

    size_t totalSize=0;
    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i)
        totalSize += getInf(inFiles[file_i].c_str());
    if(totalSize<50000000000) opt::sBits=7;

    size_t totalKmers=0;
    opt::rBuck = ((unsigned)1) << opt::rBits;
    opt::sMask = (((unsigned)1) << (opt::sBits-1))-1;
    uint16_t *t_Counter = new uint16_t [opt::nSamp*opt::rBuck]();

#ifdef _OPENMP
    omp_set_num_threads(opt::nThrd);
#endif

    #pragma omp parallel for schedule(dynamic)
    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
        size_t totKmer=0;
        std::ifstream in(inFiles[file_i].c_str());
        std::string samSeq;
        unsigned ftype = getftype(in,samSeq);
        if(ftype==0)
            getEfq(in, t_Counter, totKmer);
        else if (ftype==1)
            getEfa(in, t_Counter, totKmer);
        else if (ftype==2)
            getEsm(in, samSeq, t_Counter, totKmer);
        in.close();
        #pragma omp atomic
        totalKmers+=totKmer;
    }

    unsigned x[opt::nSamp][65536];
    for(size_t i=0; i<opt::nSamp; i++)
        for(size_t j=0; j<65536; j++)
            x[i][j]=0;

    for(size_t i=0; i<opt::nSamp; i++)
        for(size_t j=0; j<opt::rBuck; j++)
            ++x[i][t_Counter[i*opt::rBuck+j]];

    delete [] t_Counter;

    double xMean[65536];
    for(size_t i=0; i<65536; i++) xMean[i]=0.0;
    for(size_t i=0; i<65536; i++) {
        for(size_t j=0; j<opt::nSamp; j++)
            xMean[i]+=x[j][i];
        xMean[i] /= 1.0*opt::nSamp;
    }

    double f[65536];
    for(size_t i=0; i<65536; i++) f[i]=0;
    f[1]= -1.0*xMean[1]/(xMean[0]*(log(xMean[0])-opt::rBits*log(2)));
    for(size_t i=2; i<65536; i++) {
        double sum=0.0;
        for(size_t j=1; j<i; j++)
            sum+=j*xMean[i-j]*f[j];
        f[i]=-1.0*xMean[i]/(xMean[0]*(log(xMean[0])-opt::rBits*log(2)))-sum/(i*xMean[0]);
    }
    double F0= (opt::rBits*log(2)-log(xMean[0])) * 1.0* ((size_t)1<<(opt::sBits+opt::rBits));
    std::cout << "F1\t" << totalKmers << "\n";
    std::cout << "F0\t" << (long long)F0 << "\n";
    for(size_t i=1; i<= opt::covMax; i++) {
        std::cout << "f" << i << "\t" << abs((long long)(f[i]*F0)) << "\n";
    }
    std::cerr << "Runtime(sec): " <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    return 0;
}

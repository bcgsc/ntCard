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

#define PROGRAM "ntCard"

static const char VERSION_MESSAGE[] =
    PROGRAM " 1.1.0 \n"
    "Written by Hamid Mohamadi.\n"
    "Copyright 2018 Hamid Mohamadi, Licensed under MIT License\n";

static const char USAGE_MESSAGE[] =
    "Usage: " PROGRAM " [OPTION]... FILE(S)...\n"
    "Estimates k-mer coverage histogram in FILE(S).\n\n"
    "Acceptable file formats: fastq, fasta, sam, bam and in compressed formats gz, bz, zip, xz.\n"
    "A list of files containing file names in each row can be passed with @ prefix.\n"

    "\n"
    " Options:\n"
    "\n"
    "  -t, --threads=N	use N parallel threads [1] (N>=2 should be used when input files are >=2)\n"
    "  -k, --kmer=N	the length of kmer \n"
    "  -c, --cov=N	the maximum coverage of kmer in output [1000]\n"
    "  -p, --pref=STRING    the prefix for output file name(s)\n"
    "  -o, --output=STRING	the name for output file name (used when output should be a single file)\n"
    "      --help	display this help and exit\n"
    "      --version	output version information and exit\n"
    "\n"
    "Report bugs to https://github.com/bcgsc/ntCard/issues\n";


using namespace std;

namespace opt {
unsigned nThrd=1;
unsigned kmLen=64;
size_t rBuck;
unsigned rBits=27;
unsigned sBits=11;
unsigned sMask;
unsigned covMax=1000;
size_t nSamp=2;
size_t nK=0;
string prefix;
string output;
bool samH=true;
}

static const char shortopts[] = "t:s:r:k:c:l:p:f:o:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads",	required_argument, NULL, 't' },
    { "kmer",	required_argument, NULL, 'k' },
    { "cov",	required_argument, NULL, 'c' },
    { "rbit",	required_argument, NULL, 'r' },
    { "sbit",	required_argument, NULL, 's' },
    { "output",	required_argument, NULL, 'o' },
    { "pref",	required_argument, NULL, 'p' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

size_t getInf(const char* inFile) {
    std::ifstream in(inFile, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

bool isNumber(const string &seq) {
    string::const_iterator itr = seq.begin();
    while(itr!= seq.end() && isdigit(*itr)) ++itr;
    return (itr==seq.end() && !seq.empty());
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
    std::istringstream alnSec(hseq);
    std::string s1,s2,s3,s4,s5,s6,s7,s8,s9,s10, s11;
    alnSec>>s1>>s2>>s3>>s4>>s5>>s6>>s7>>s8>>s9>>s10>>s11;
    if(isNumber(s2) && isNumber(s5)) {
        opt::samH=false;
        samSeq=hseq;
        return 2;
    }
    return 3;
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

inline void ntRead(const string &seq, const std::vector<unsigned> &kList, uint16_t *t_Counter, size_t totKmer[]) {
    for(unsigned k=0; k<kList.size(); k++) {
        ntHashIterator itr(seq, kList[k]);
        while (itr != itr.end()) {
            ntComp((*itr),t_Counter+k*opt::nSamp*opt::rBuck);
            ++itr;
            ++totKmer[k];
        }
    }
}

void getEfq(std::ifstream &in, const std::vector<unsigned> &kList, uint16_t *t_Counter, size_t totKmer[]) {
    bool good = true;
    for(string seq, hseq; good;) {
        good = static_cast<bool>(getline(in, seq));
        good = static_cast<bool>(getline(in, hseq));
        good = static_cast<bool>(getline(in, hseq));
        if(good)
            ntRead(seq, kList, t_Counter, totKmer);
        good = static_cast<bool>(getline(in, hseq));
    }
}

void getEfa(std::ifstream &in, const std::vector<unsigned> &kList, uint16_t *t_Counter, size_t totKmer[]) {
    bool good = true;
    for(string seq, hseq; good;) {
        string line;
        good = static_cast<bool>(getline(in, seq));
        while(good&&seq[0]!='>') {
            line+=seq;
            good = static_cast<bool>(getline(in, seq));
        }
        ntRead(line, kList, t_Counter, totKmer);
    }
}

void getEsm(std::ifstream &in, const std::vector<unsigned> &kList, const std::string &samSeq, uint16_t *t_Counter, size_t totKmer[]) {
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
        ntRead(seq, kList, t_Counter, totKmer);
    } while(getline(in,samLine));
}

void compEst(const uint16_t *t_Counter, double &F0Mean, double fMean[]) {
    unsigned p[opt::nSamp][65536];
    for(size_t i=0; i<opt::nSamp; i++)
        for(size_t j=0; j<65536; j++)
            p[i][j]=0;

    for(size_t i=0; i<opt::nSamp; i++)
        for(size_t j=0; j<opt::rBuck; j++)
            ++p[i][t_Counter[i*opt::rBuck+j]];

    double pMean[65536];
    for(size_t i=0; i<65536; i++) pMean[i]=0.0;
    for(size_t i=0; i<65536; i++) {
        for(size_t j=0; j<opt::nSamp; j++)
            pMean[i]+=p[j][i];
        pMean[i] /= 1.0*opt::nSamp;
    }

    F0Mean= (ssize_t)((opt::rBits*log(2)-log(pMean[0])) * 1.0* ((size_t)1<<(opt::sBits+opt::rBits)));
    for(size_t i=0; i<65536; i++) fMean[i]=0;
    fMean[1]= -1.0*pMean[1]/(pMean[0]*(log(pMean[0])-opt::rBits*log(2)));
    for(size_t i=2; i<65536; i++) {
        double sum=0.0;
        for(size_t j=1; j<i; j++)
            sum+=j*pMean[i-j]*fMean[j];
        fMean[i]=-1.0*pMean[i]/(pMean[0]*(log(pMean[0])-opt::rBits*log(2)))-sum/(i*pMean[0]);
    }
    for(size_t i=1; i<65536; i++)
        fMean[i]=abs((ssize_t)(fMean[i]*F0Mean));
}

void outDefault(const std::vector<unsigned> &kList, const size_t totalKmers[], const uint16_t *t_Counter) {
    std::ofstream histFiles[opt::nK];
    for (unsigned k=0; k<opt::nK; k++) {
        std::stringstream hstm;
        hstm << opt::prefix << "_k" << kList[k] << ".hist";
        histFiles[k].open((hstm.str()).c_str());
    }
    #pragma omp parallel for schedule(dynamic)
    for(unsigned k=0; k<kList.size(); k++) {
        double F0Mean=0.0;
        double fMean[65536];
        compEst(t_Counter+k*opt::nSamp*opt::rBuck, F0Mean, fMean);
        histFiles[k] << "F1\t" << totalKmers[k] << "\n";
        histFiles[k] << "F0\t" << (size_t)F0Mean << "\n";
        for(size_t i=1; i<= opt::covMax; i++)
            histFiles[k] << i << "\t" << (size_t)fMean[i] << "\n";
    }
    for (unsigned k=0; k<opt::nK; k++)
        histFiles[k].close();
}

void outCompact(const std::vector<unsigned> &kList, const size_t totalKmers[], const uint16_t *t_Counter) {
    ofstream histFile(opt::output.c_str());
    histFile << "k\tf\tn\n";
    for(unsigned k=0; k<kList.size(); k++) {
        double F0Mean=0.0;
        double fMean[65536];
        compEst(t_Counter+k*opt::nSamp*opt::rBuck, F0Mean, fMean);
        std::cerr << "k=" << kList[k] << "\tF1\t" << totalKmers[k] << "\n";
        std::cerr << "k=" << kList[k] << "\tF0\t" << (size_t)F0Mean << "\n";
        for(size_t i=1; i<= opt::covMax; i++)
            histFile << kList[k] << "\t" << i << "\t" << (size_t)fMean[i] << "\n";
    }
    histFile.close();
}

int main(int argc, char** argv) {

    double sTime = omp_get_wtime();

    vector<unsigned> kList;
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
            if(opt::covMax>65535)
                opt::covMax = 65535;
            break;
        case 'p':
            arg >> opt::prefix;
            break;
        case 'o':
            arg >> opt::output;
            break;
        case 'k':
        {
            std::string token;
            while(getline(arg, token, ',')) {
                unsigned myK;
                std::stringstream ss(token);
                ss >> myK;
                kList.push_back(myK);
                ++opt::nK;
            }
            break;
        }
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

    if (opt::nK == 0) {
        std::cerr << PROGRAM ": missing argument -k ... \n";
        die = true;
    }

    if (opt::prefix.empty() && opt::output.empty()) {
        std::cerr << PROGRAM ": missing argument -p/-o ... \n";
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

    size_t totalKmers[kList.size()];
    for(unsigned k=0; k<kList.size(); k++) totalKmers[k]=0;

    opt::rBuck = ((size_t)1) << opt::rBits;
    opt::sMask = (((size_t)1) << (opt::sBits-1))-1;
    uint16_t *t_Counter = new uint16_t [opt::nK*opt::nSamp*opt::rBuck]();

#ifdef _OPENMP
    omp_set_num_threads(opt::nThrd);
#endif

    #pragma omp parallel for schedule(dynamic)
    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
        size_t totKmer[kList.size()];
        for(unsigned k=0; k<kList.size(); k++) totKmer[k]=0;
        std::ifstream in(inFiles[file_i].c_str());
        std::string samSeq;
        unsigned ftype = getftype(in,samSeq);
        if(ftype==0)
            getEfq(in, kList, t_Counter, totKmer);
        else if (ftype==1)
            getEfa(in, kList, t_Counter, totKmer);
        else if (ftype==2)
            getEsm(in, kList, samSeq, t_Counter, totKmer);
        else {
            std::cerr << "Error in reading file: " << inFiles[file_i] << std::endl;
            exit(EXIT_FAILURE);
        }
        in.close();
        for(unsigned k=0; k<kList.size(); k++)
            #pragma omp atomic
            totalKmers[k]+=totKmer[k];
    }

    if(opt::output.empty())
        outDefault(kList, totalKmers, t_Counter);
    else
        outCompact(kList, totalKmers, t_Counter);

    delete [] t_Counter;

    std::cerr << "Runtime(sec): " <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    return 0;
}

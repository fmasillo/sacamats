
#ifndef SACAMATS_H
#define SACAMATS_H

#include <cstdint>
#include <vector>
#include <string>
#include <omp.h>

using data_type = uint8_t;
using filelength_type = uint64_t;

// struct containing command line parameters and other globals
struct Args {
  std::string filename = "";
  std::string outname = "";
  long long unsigned int prefixLength = UINT64_MAX; // prefix length for input collection file
  long long unsigned int nThreads = 0; // number of threads to use
  bool check = 0; // check the correctness of the output
};


const std::string lzInitialize(char *refFileName, char *collFileName, uint64_t prefixLength);
//void lzInitialize(data_type *ax, unsigned int an, bool isMismatchingSymbolNeeded, char *refFileName, char *collFileName);
//void lzInitialize(data_type *ax, unsigned int *sa, unsigned int an, bool isMismatchingSymbolNeeded, char *refFileName, char *collFileName);

//returns the length of the RLZ parsing of sx relative to ax
int lzFactorize(const std::string _sx);
//int lzFactorize(char *fileToParse, int seqno, char* outputfilename, bool verbose);

void computeGSA(char *refFileName, char *filename, uint64_t prefixLength, uint32_t nThreads, std::vector<std::pair<uint32_t, int32_t>> &MSGSA);

#endif

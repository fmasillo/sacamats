
#ifndef SACAMATS_H
#define SACAMATS_H

#include <cstdint>

using data_type = uint8_t;
using filelength_type = uint64_t;

void lzInitialize(data_type *ax, unsigned int an, bool isMismatchingSymbolNeeded, char *refFileName, char *collFileName);
//void lzInitialize(data_type *ax, unsigned int *sa, unsigned int an, bool isMismatchingSymbolNeeded, char *refFileName, char *collFileName);

//returns the length of the RLZ parsing of sx relative to ax
int lzFactorize(char *fileToParse, int seqno, char* outputfilename, bool verbose);

#endif

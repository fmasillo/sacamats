#include <stdlib.h>
#include <algorithm>
#include <chrono>
#include <ctime>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <vector>
#include <string.h>
#include <cstdlib>
#include <fstream>
#include <omp.h>
#ifndef _OPENMP
   #define _OPENMP ; 
#endif
#include "libsais/src/libsais.h"
#include "sacamats.h"
#include "rmq_tree.h"
#include "utils.h"
#include "match.h"
#include "predecessor.h"

#include "induce_parallel.cpp"

#undef max
#define sequenceSeparator 2

void computeLZFactorAt(const std::string &_sx, filelength_type i, filelength_type *pos, filelength_type *len, int64_t &leftB, int64_t &rightB, int64_t &match, bool &isSmaller, unsigned char &mismatchingSymbol);
inline int64_t binarySearchLB(int64_t lo, int64_t hi, filelength_type offset, data_type c);
inline int64_t binarySearchRB(int64_t lo, int64_t hi, filelength_type offset, data_type c);

static std::string outputFile;

//static data_type *_x;
static std::string _x;
static int *_SA;
static uint32_t *_ISA;
static int32_t *_LCP;
static int32_t *_PLCP;
static rmq_tree *_rmq;
static uint32_t _n;
//static data_type *_sx;
static filelength_type _sn;
static bool _isMismatchingSymbolNeeded;
static std::vector<uint64_t> docBoundaries;
static std::vector<uint32_t> headBoundaries;


uint64_t tiesCounter = 0;

uint64_t _maxLCP = 0;
uint64_t _numberOfShortFactors = 0;
uint64_t _numberOfSingleMatches = 0;

uint64_t maxFactorLength = 0;
int lenZeroFactors = 0;

bool verbose;


uint16_t sizeChars = 256;
uint64_t D = 0;
uint64_t *charBkts = new uint64_t[sizeChars]();

//std::vector<std::pair<uint32_t, int32_t>> MSGSA;

//std::vector<std::pair<uint32_t,uint32_t> > phrases;
std::vector<Match> phrases;
std::vector<MatchSA> headsSA;

void constructISA(int32_t *sa, uint32_t *isa, uint32_t n){
   fprintf(stderr,"\tComputing ISA...\n");
   for(uint32_t i=0; i<n; i++){
      isa[sa[i]] = i;
   }
}

std::pair<int,int> adjustInterval(int lo, int hi, int offset) {
   int psv = _rmq->psv(lo,offset);
   if(psv == -1){
      psv = 0;
   }
   int nsv = _rmq->nsv(hi+1,offset);
   if(nsv == -1){
      nsv = _n-1;
   }else{
      nsv--;
   }
   return std::make_pair(psv,nsv);
}

std::pair<int,int> contractLeft(int lo, int hi, int offset){
   uint32_t suflo = _SA[lo];
   uint32_t sufhi = _SA[hi];
   if(suflo == _n-1 || sufhi == _n-1){ //if true we must be at depth 1
      return std::make_pair(0,_n-1); //root
   }
   uint32_t tmplo = _ISA[suflo+1]; 
   uint32_t tmphi = _ISA[sufhi+1];
   return adjustInterval(tmplo, tmphi, offset);
}

void radixSortPos(const std::vector<SufSStar>::iterator a, const uint32_t start, const uint32_t end, const uint8_t radixIter){
   const int RADIX = 0x100;
   std::vector<SufSStar> result;
   result.resize(end-start);  
   uint32_t *buckets = new uint32_t[RADIX](); 
   uint32_t *startIndex = new uint32_t[RADIX]();
   data_type key = 0;
   for(uint32_t flag = 0; flag < radixIter; flag++){
      //std::cerr << "flag: " << flag << "\n";
      for(std::vector<SufSStar>::iterator it = a+start; it < a+end; ++it) {
         //key = (it->pos & (MASK << flag)) >> flag;
         key = it->diffLen.bytearray[flag];
         ++buckets[key];
      }
      startIndex[0] = 0;
      for(uint32_t j = 1; j < RADIX; ++j) startIndex[j] = startIndex[j - 1] + buckets[j - 1];
      for(std::vector<SufSStar>::iterator it = a+end-1; it >= a+start; --it){
         //key = (it->pos & (MASK << flag)) >> flag;
         key = it->diffLen.bytearray[flag];
         result[startIndex[key] + --buckets[key]].assign(it->idx, it->doc, it->head, it->diffLen.val);
      }
      std::copy(result.begin(), result.end(), a+start);
      //flag += MASK_BIT_LENGTH;
   }
   delete[] buckets;
   delete[] startIndex;
}


bool compareMatchSA(const MatchSA &a, const MatchSA &b){
   if((a.smaller == 0) & (b.smaller == 0)){
      return (a.next < b.next)*(a.len == b.len) + (a.len < b.len)*(a.len != b.len);
   }
   else if((a.smaller == 1) & (b.smaller == 1)){
      return (a.next < b.next)*(a.len == b.len) + (a.len > b.len)*(a.len != b.len);
   }
   return a.smaller < b.smaller;
}

bool compareSuf(const SufSStar &a, const SufSStar &b){
   //finalSuffCounter++;
   MatchSA headA = headsSA[a.head];
   MatchSA headB = headsSA[b.head];
   //if(headA.len - a.diffLen.val != headB.len - b.diffLen.val){
   if(headA.len - a.idx != headB.len - b.idx){
      return headA.smaller*((headA.len - a.idx) < (headB.len - b.idx)) + 
            !headB.smaller*((headA.len - a.idx) > (headB.len - b.idx));
   }
   else{
      return (headA.pos != headB.pos)*(headA.pos < headB.pos) + (headA.pos == headB.pos)*(headA.start < headB.start);
   }
}

uint64_t checkHeadsSA(std::vector<MatchSA> &GSA, uint64_t n, std::string &_sx){
   uint64_t err = 0;
   for(size_t i = 0; i < n; i++){
       if(verbose) std::cerr << "i=" << i << ": " << phrases[GSA[i].start].start << " " << GSA[i].len << " " << "\n";//MSGSA[i].head <<"\n";
      //  if(GSA[i].len == 0 || GSA[i-1].len == 0){
      //     std::cerr << "There was an empty entry\n";
      //     continue;
      //  }
       const std::string _slice_sx = _sx.substr(phrases[GSA[i].start].start + docBoundaries[GSA[i].len - 1]);
       std::string _slice_prev;
       uint32_t maxIdx;
       if(i > 0){
         _slice_prev = _sx.substr(phrases[GSA[i-1].start].start + docBoundaries[GSA[i-1].len - 1]);
         maxIdx = std::min(docBoundaries[GSA[i].len] - (phrases[GSA[i].start].start + docBoundaries[GSA[i].len - 1]), docBoundaries[GSA[i-1].len] - (phrases[GSA[i-1].start].start + docBoundaries[GSA[i-1].len - 1]));
       } 
       else{
          _slice_prev = "$";
          maxIdx = 1;
       }
       if(verbose) std::cerr << "suf_i-1 " << _slice_prev;
       if(verbose) std::cerr << "suf_i " << _slice_sx << "\n";
       
      if(memcmp(&_slice_sx[0], &_slice_prev[0], maxIdx) < 0){
         if(verbose) std::cerr << "PROBLEM with " << i-1 << " (" << phrases[GSA[i-1].start].start << "," << GSA[i-1].len << ") and " << i << " (" << phrases[GSA[i].start].start << "," << GSA[i].len << ")\n"; 
         err++;
      }
    }
    return err;
}

uint64_t checkHeadsSA_no_copy(std::vector<MatchSA> &GSA, uint64_t n, const std::string &_sx){
   uint64_t err = 0;
   for(size_t i = 1; i < n; i++){
      if(verbose) std::cerr << "i=" << i << ": " << phrases[GSA[i].start].start << " " << GSA[i].len << " " << "\n";//MSGSA[i].head <<"\n";
      //  if(GSA[i].len == 0 || GSA[i-1].len == 0){
      //     std::cerr << "There was an empty entry\n";
      //     continue;
      //  }
      uint32_t maxIdx = std::min(docBoundaries[GSA[i].len] - (phrases[GSA[i].start].start + docBoundaries[GSA[i].len - 1]), docBoundaries[GSA[i-1].len] - (phrases[GSA[i-1].start].start + docBoundaries[GSA[i-1].len - 1]));
      
      if(memcmp(&_sx[phrases[GSA[i].start].start + docBoundaries[GSA[i].len - 1]], &_sx[phrases[GSA[i-1].start].start + docBoundaries[GSA[i-1].len - 1]], maxIdx) < 0){
         if(verbose) std::cerr << "PROBLEM with " << i-1 << " (" << phrases[GSA[i-1].start].start << "," << GSA[i-1].len << ") and " << i << " (" << phrases[GSA[i].start].start << "," << GSA[i].len << ")\n"; 
         err++;
      }
   }
   return err;
}

uint64_t checkGSA(std::vector<SufSStar> GSA, uint64_t n, const std::string &_sx){
   uint64_t err = -1;
   //#pragma omp parallel for
   for(size_t i = 0; i < n; i++){
      if(verbose) std::cerr << "i=" << i << ": " << GSA[i].idx << " " << GSA[i].doc << " " << "\n";//MSGSA[i].head <<"\n";
      if(i > 0)
         if(GSA[i].doc == 0 || GSA[i-1].doc == 0){
            std::cerr << "There was an empty entry\n";
            //#pragma omp atomic
            err++;
            continue;
         }
      const std::string _slice_sx = _sx.substr(GSA[i].idx + docBoundaries[GSA[i].doc - 1]);
      std::string _slice_prev;
      uint32_t maxIdx;
      if(i > 0){
         _slice_prev = _sx.substr(GSA[i-1].idx + docBoundaries[GSA[i-1].doc - 1]);
         maxIdx = std::min(docBoundaries[GSA[i].doc] - (GSA[i].idx + docBoundaries[GSA[i].doc - 1]), docBoundaries[GSA[i-1].doc] - (GSA[i-1].idx + docBoundaries[GSA[i-1].doc - 1]));
      } 
      else{
         _slice_prev = sequenceSeparator;
         maxIdx = 1;
      }
      
      if(memcmp(&_slice_sx[0], &_slice_prev[0], maxIdx) < 0){
         std::cerr << "PROBLEM with " << i-1 << " (" << GSA[i-1].idx << "," << GSA[i-1].doc << ") and " << i << " (" << GSA[i].idx << "," << GSA[i].doc << ")\n"; 
         //#pragma omp atomic
         err++;
      }
    }
    return err;
}

uint64_t checkGSA_no_copy(std::vector<SufSStar> &GSA, uint64_t n, const std::string &_sx){
   uint64_t err = 0;
   //GSA[30] = GSA[n/2]; //just to check if errors can be found
   #pragma omp parallel for
   for(size_t i = 1; i < n; i++){
      if(verbose) std::cerr << "i=" << i << ": " << GSA[i].idx << " " << GSA[i].doc << " " << "\n";//MSGSA[i].head <<"\n";
      //  if(GSA[i].len == 0 || GSA[i-1].len == 0){
      //     std::cerr << "There was an empty entry\n";
      //     continue;
      //  }
      uint32_t maxIdx = std::min(docBoundaries[GSA[i].doc] - (GSA[i].idx + docBoundaries[GSA[i].doc - 1]), docBoundaries[GSA[i-1].doc] - (GSA[i-1].idx + docBoundaries[GSA[i-1].doc - 1]));
  
      if(memcmp(&_sx[GSA[i].idx + docBoundaries[GSA[i].doc - 1]], &_sx[GSA[i-1].idx + docBoundaries[GSA[i-1].doc - 1]], maxIdx) < 0){
         std::cerr << "PROBLEM with " << i-1 << " (" << GSA[i-1].idx << "," << GSA[i-1].doc << ") and " << i << " (" << GSA[i].idx << "," << GSA[i].doc << ")\n"; 
         err++;
      }
   }
   return err;
}

const std::string lzInitialize_omp(char *refFileName, char *collFileName, uint64_t prefixLength, uint32_t nThreads) {
//void lzInitialize(data_type *ax, unsigned int an, bool isMismatchingSymbolNeeded, char *refFileName, char *collFileName) {
   auto t1 = std::chrono::high_resolution_clock::now();
   omp_set_num_threads(nThreads);
   errno = 0;
   FILE *infileRef = fopen(refFileName, "r");
   if (!infileRef) {
      fprintf(stderr, "Error opening file of base sequence %s, errno=%d\n", refFileName,errno);
      exit(1);
   }
   fprintf(stderr, "About to read ref\n");

   unsigned int n = 0;
   //data_type *x = new data_type[n + 1];
   fseek(infileRef, 0, SEEK_END);
   n = ftell(infileRef) / sizeof(data_type);
   std::cerr << "n = " << n << '\n';
   fseek(infileRef, 0, SEEK_SET);
   if(n){
      char *firstChar = new char[1];
      int o = fread(firstChar, sizeof(char), 1, infileRef);
      //std::cerr << "firstChar: " << firstChar << '\n';
      if(firstChar[0] == '>'){
         _x.reserve(n);
         //std::cerr << "Yes FASTA\n";
         std::ifstream streamInfileRef(refFileName);
         std::string line, content;

         while(std::getline(streamInfileRef, line).good()){
            if( line.empty() || line[0] == '>' ){
               _x += content;
               std::string().swap(content);
            }
            else if (!line.empty()) {
               content += line;
            }
         }
         if(content.size()) _x += content;
         std::string().swap(content);
         _x.resize(_x.size());
         streamInfileRef.close();
      }
      else{
         fseek(infileRef, 0, SEEK_SET);
         _x.resize(n);
         if (n != fread(&_x[0], sizeof(data_type), n, infileRef)) {
            fprintf(stderr, "Error reading %u bytes from file %s\n", n, refFileName);
            exit(1);
         }
      }
   }
   else{
      std::cerr << "Reference file is empty!\n";
      exit(1);
   }
   fclose(infileRef);
   fprintf(stderr, "Reference (size = %lu):\n\t", _x.size());
   //std::cerr << "reference: " << _x << "\n";
   auto t01 = std::chrono::high_resolution_clock::now();
   //_x = std::string(reinterpret_cast<char const *>(x), n);
   if((_x[_x.size()-1] == '\n') | (_x[_x.size()-1] == '\r') | (_x[_x.size()-1] == 0)){
      _x.erase(_x.size()-1);
   }
   if(_x[_x.size()-1] == '$'){
      _x.erase(_x.size()-1);
   }

   FILE *infile = fopen(collFileName, "r");
   if(!infile){
      fprintf(stderr, "Error opening file of sequence (%s)\n", collFileName);
      exit(1);
   }
   filelength_type sn = 0;
   fseek(infile, 0L, SEEK_END);
   std::cerr << "ftello(infile): " << ftello(infile) << '\n';
   sn = std::min(ftello(infile) / sizeof(data_type), prefixLength);
   std::cerr << "sn: " << sn << '\n';
   fseek(infile, 0, SEEK_SET);
   std::string sx;
   docBoundaries.push_back(0);
   if(sn){
      char *firstChar = new char[1];
      int o = fread(firstChar, sizeof(char), 1, infile);
      fseek(infile, 0, SEEK_SET);
      //std::cerr << "firstChar: " << firstChar << '\n';
      if(firstChar[0] == '>'){
         sx.reserve(sn);
         std::cerr << "Yes FASTA\n";
         std::ifstream streamInfile(collFileName,std::ios::in);
         std::string line, content;
         
         while(std::getline(streamInfile, line).good()){
            if( line.empty() || line[0] == '>' ){
               sx += content;
               sx += sequenceSeparator;
               docBoundaries.push_back(sx.size());
               std::string().swap(content);
            }
            else if (!line.empty()) {
               content += line;
            }
            if(sx.size() > sn){ //if string is filled up to sn (useful for prefixLength)
               docBoundaries.pop_back();
               break;
            }
         }
         if(content.size() > 1){
            sx += content;
            sx += sequenceSeparator;
         }
         std::string().swap(content);
         
         sx.shrink_to_fit();
         sx.resize(std::min(sn,sx.size()));
         std::cerr << "sx[last]: " << sx[sx.size()-1] << '\n';
         if(sx[sx.size()-1] != sequenceSeparator) sx += sequenceSeparator;
         if(docBoundaries[docBoundaries.size()-1] != sx.size()) docBoundaries.push_back(sx.size());
         sn = sx.size();
         streamInfile.close();
         std::cerr << "File size after FASTA cleaning: " << sn << '\n';
      }
      else{
         std::cerr << "No FASTA\n";
         fseek(infile, 0, SEEK_SET);
         sx.resize(sn);
         if (sn != fread(&sx[0], sizeof(data_type), sn, infile)) {
            fprintf(stderr, "Error reading %u bytes from file %s\n", n, collFileName);
            exit(1);
         }
         if (sx[sn-1] != sequenceSeparator){
            sx += sequenceSeparator;
            sn++;
         }
         sx.shrink_to_fit();
         sn = sx.size();
      }
   }
   else{
      std::cerr << "Collection file is empty!\n";
      exit(1);
   }
   fclose(infile);

   std::cerr << "Collection read.\n";
   uint64_t *maxRunsReference = new uint64_t[sizeChars]();
   data_type c = '<';
   uint64_t runLen = 0;
   for(uint32_t i = 0; i < _x.size(); i++){
      if(_x[i] != c){
         if(runLen > maxRunsReference[c]){
            maxRunsReference[c] = runLen;
         }
         runLen = 1;
         c = _x[i];
         continue;
      }
      runLen++;
   }

   uint32_t *maxRunsCollection = new uint32_t[sizeChars]();
   #pragma omp parallel for
   for(uint32_t d = 1; d < docBoundaries.size(); d++){
      uint32_t *privateMaxRunsCollection = new uint32_t[sizeChars]();
      uint64_t *privateCharBkts = new uint64_t[sizeChars]();
      data_type c = (data_type)0;
      uint32_t runLen = 0;
      for(uint64_t i = docBoundaries[d-1]; i < docBoundaries[d]; i++){
         privateCharBkts[sx[i]]++;
         if(sx[i] != c){
            if(runLen > privateMaxRunsCollection[c]){
               privateMaxRunsCollection[c] = runLen;
            }
            runLen = 1;
            c = sx[i];
            continue;
         }
         runLen++;
      }

      #pragma omp critical
      {
         for(uint16_t j = 0; j < sizeChars; j++){
            charBkts[j] += privateCharBkts[j];
            maxRunsCollection[j] = maxRunsCollection[j] < privateMaxRunsCollection[j] ? privateMaxRunsCollection[j] : maxRunsCollection[j];
         }
         delete[] privateCharBkts;
         delete[] privateMaxRunsCollection;
      }
   }
   D = docBoundaries.size();

   for(uint16_t i = 0; i < sizeChars; i++){
      if((maxRunsReference[i] == 0) & (maxRunsCollection[i] != 0) & (i != sequenceSeparator)){
         _x.append(maxRunsCollection[i], (data_type)i);
      }
   }
   _x += (char)1;
   _x += (char)0;
   //docBoundaries.reserve(maxRunsCollection[sequenceSeparator]);
   headBoundaries.reserve(maxRunsCollection[sequenceSeparator]);
   //std::cerr << refAug << "\n";
   auto t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Augmenting reference done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t01 = std::chrono::high_resolution_clock::now();
   int *sa = new int[_x.size()];
   libsais_omp(reinterpret_cast<unsigned char*>(const_cast<char*>(_x.c_str())), sa, _x.size(), 0, NULL, nThreads);
   //sais(reinterpret_cast<unsigned char*>(const_cast<char*>(_x.c_str())), sa, _x.size());
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing SA done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   //std::cerr << sx << "\n";
   //_sx = reinterpret_cast<unsigned char*>(const_cast<char*>(sx.c_str()));
   const std::string _sx{sx};
   if(_sx[sn-1] != sequenceSeparator){
      _sn = sn - 1;
   }
   else{
      _sn = sn;
   }
   _n = _x.size();

   t01 = std::chrono::high_resolution_clock::now();
   _SA = sa;
   _ISA = new uint32_t[_n];
   _PLCP = new int32_t[_n];
   _LCP = new int32_t[_n];
   t01 = std::chrono::high_resolution_clock::now();
   constructISA(_SA,_ISA,_n);
   libsais_plcp_omp(reinterpret_cast<unsigned char*>(const_cast<char*>(_x.c_str())), _SA, _PLCP, _n, nThreads);
   libsais_lcp_omp(_PLCP, _SA, _LCP, _n, nThreads);
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing ISA done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   //constructLCP(_x,_n,_SA,_LCP,_PLCP);
   #pragma omp parallel for
   for(uint32_t i = 0; i < _n; i++){
      _PLCP[i] = std::max(_LCP[_ISA[i]],_LCP[_ISA[i]+1]);
   }
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing LCP and PLCP done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t01 = std::chrono::high_resolution_clock::now();
   fprintf(stderr,"\tComputing RMQ...\n");
   _rmq = new rmq_tree((int *)_LCP, _n, 7);
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing RMQ done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   //_isMismatchingSymbolNeeded = isMismatchingSymbolNeeded;
   std::cerr << "Finished pre-processing\n";

   auto t2 = std::chrono::high_resolution_clock::now();
   uint64_t preprocTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Preprocessing done in " << preprocTime << " ms\n";
   return _sx;
}

const std::string lzInitialize(char *refFileName, char *collFileName, uint64_t prefixLength) {
//void lzInitialize(data_type *ax, unsigned int an, bool isMismatchingSymbolNeeded, char *refFileName, char *collFileName) {
   auto t1 = std::chrono::high_resolution_clock::now();
   errno = 0;
   FILE *infileRef = fopen(refFileName, "r");
   if (!infileRef) {
      fprintf(stderr, "Error opening file of base sequence %s, errno=%d\n", refFileName,errno);
      exit(1);
   }
   fprintf(stderr, "About to read ref\n");

   unsigned int n = 0;
   fseek(infileRef, 0, SEEK_END);
   n = ftell(infileRef) / sizeof(data_type);
   std::cerr << "n = " << n << '\n';
   fseek(infileRef, 0, SEEK_SET);
   if(n){
      char *firstChar = new char[1];
      int o = fread(firstChar, sizeof(char), 1, infileRef);
      if(firstChar[0] == '>'){
         _x.reserve(n);
         std::ifstream streamInfileRef(refFileName);
         std::string line, content;
         while(std::getline(streamInfileRef, line).good()){
            if( line.empty() || line[0] == '>' ){
               _x += content;
               std::string().swap(content);
            }
            else if (!line.empty()) {
               content += line;
            }
         }
         if(content.size()) _x += content;
         std::string().swap(content);
         _x.resize(_x.size());
         streamInfileRef.close();
      }
      else{
         //std::cerr << "No FASTA\n";
         fseek(infileRef, 0, SEEK_SET);
         _x.resize(n);
         if (n != fread(&_x[0], sizeof(data_type), n, infileRef)) {
            fprintf(stderr, "Error reading %u bytes from file %s\n", n, refFileName);
            exit(1);
         }
      }
      //x[n] = 0;
   }
   else{
      std::cerr << "Reference file is empty!\n";
      exit(1);
   }
   fclose(infileRef);
   fprintf(stderr, "Reference (size = %lu):\n\t", _x.size());

   auto t01 = std::chrono::high_resolution_clock::now();
   //_x = std::string(reinterpret_cast<char const *>(x), n);
   if((_x[_x.size()-1] == '\n') | (_x[_x.size()-1] == '\r') | (_x[_x.size()-1] == 0)){
      _x.erase(_x.size()-1);
   }
   if(_x[_x.size()-1] == '$'){
      _x.erase(_x.size()-1);
   }
   FILE *infile = fopen(collFileName, "r");
   if(!infile){
      fprintf(stderr, "Error opening file of sequence (%s)\n", collFileName);
      exit(1);
   }
   filelength_type sn = 0;
   fseek(infile, 0L, SEEK_END);
   std::cerr << "ftello(infile): " << ftello(infile) << '\n';
   sn = std::min(ftello(infile) / sizeof(data_type), prefixLength);
   std::cerr << "sn: " << sn << '\n';
   fseek(infile, 0, SEEK_SET);
   std::string sx;
   if(sn){
      char *firstChar = new char[1];
      int o = fread(firstChar, sizeof(char), 1, infile);
      //std::cerr << "firstChar: " << firstChar << '\n';
      if(firstChar[0] == '>'){
         sx.reserve(sn);
         std::cerr << "Yes FASTA\n";
         std::ifstream streamInfile(collFileName);
         std::string line, content;
         while(std::getline(streamInfile, line).good()){
            if( line.empty() || line[0] == '>' ){
               sx += content;
               sx += sequenceSeparator;
               //std::cerr << _x << "size: " << _x.size() << "\n";
               std::string().swap(content);
            }
            else if (!line.empty()) {
               content += line;
            }
            if(sx.size() >= sn){ //if string is filled up to sn (useful for prefixLength)
               break;
            }
         }
         if(content.size()){
            sx += content;
            sx += sequenceSeparator;
         }
         //std::cerr << _x << "size: " << _x.size() << "\n";
         std::string().swap(content);
         sx.resize(sn);
         sx += sequenceSeparator;
         sn = sx.size();
         streamInfile.close();
         std::cerr << "File size after FASTA cleaning: " << sn << '\n';
      }
      else{
         std::cerr << "No FASTA\n";
         fseek(infile, 0, SEEK_SET);
         sx.resize(sn);
         if (sn != fread(&sx[0], sizeof(data_type), sn, infile)) {
            fprintf(stderr, "Error reading %u bytes from file %s\n", n, collFileName);
            exit(1);
         }
         if (sx[sn-1] != sequenceSeparator){
            sx += sequenceSeparator;
            sn++;
         }
         
      }
   }
   else{
      std::cerr << "Collection file is empty!\n";
      exit(1);
   }
   fclose(infile);

   std::cerr << "Collection read.\n";

   // uint8_t *_ssx = new uint8_t[sn+1];
   // for(uint64_t i = 0; i < sn; i++){
   //    _ssx[i] = sx[i];
   // }
   // _ssx[sn] = 0;
   // std::string().swap(sx);
   // uint64_t r = count_runs_in_bwt(_ssx, sn);
   // std::cerr << "Number of runs in BWT: " << r << "\n";
   // exit(0);

   uint64_t *maxRunsReference = new uint64_t[sizeChars]();
   data_type c = '<';
   uint64_t runLen = 0;
   for(uint32_t i = 0; i < _x.size(); i++){
      if(_x[i] != c){
         if(runLen > maxRunsReference[c]){
            maxRunsReference[c] = runLen;
         }
         runLen = 1;
         c = _x[i];
         continue;
      }
      runLen++;
   }

   c = '<';
   runLen = 0;
   docBoundaries.push_back(0);
   uint64_t *maxRunsCollection = new uint64_t[sizeChars]();
   for(uint64_t i = 0; i < sn; i++){
      charBkts[sx[i]]++;
      if(sx[i] == sequenceSeparator){
         D++;
         docBoundaries.push_back(i+1);
      }
      if(sx[i] != c){
         if(runLen > maxRunsCollection[c]){
            maxRunsCollection[c] = runLen;
         }
         runLen = 1;
         c = sx[i];
         continue;
      }
      runLen++;
   }
   docBoundaries.push_back(_sn);

   //std::string refAug(_x);
   for(uint16_t i = 0; i < sizeChars; i++){
      if(verbose) std::cerr << (char)i << " ref: " << maxRunsReference[i] << ", coll: " << maxRunsCollection[i] << "\n";
      if((maxRunsReference[i] == 0) & (maxRunsCollection[i] != 0) & (i != sequenceSeparator)){
         for(uint64_t x = 0; x < maxRunsCollection[i]; x++){
            _x += (char)i;
         }
      }
   }
   _x += (char)1;
   _x += (char)0;
   headBoundaries.reserve(maxRunsCollection[sequenceSeparator]);
   auto t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Augmenting reference done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t01 = std::chrono::high_resolution_clock::now();
   int32_t *sa = new int32_t[_x.size()];
   libsais(reinterpret_cast<unsigned char*>(const_cast<char*>(_x.c_str())), sa, _x.size(), 0, NULL);
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing SA done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   const std::string _sx{sx};
   if(_sx[sn-1] != sequenceSeparator){
      _sn = sn - 1;
   }
   else{
      _sn = sn;
   }
   std::cerr << "Last pos " << _sx[_sn - 1] << "\n";
   _n = _x.size();

   t01 = std::chrono::high_resolution_clock::now();
   _SA = sa;
   _ISA = new uint32_t[_n];
   _PLCP = new int32_t[_n];
   _LCP = new int32_t[_n];
   t01 = std::chrono::high_resolution_clock::now();
   constructISA(_SA,_ISA,_n);
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing ISA done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   libsais_plcp(reinterpret_cast<unsigned char*>(const_cast<char*>(_x.c_str())), _SA, _PLCP, _x.size());
   libsais_lcp(_PLCP, _SA, _LCP, _x.size());
   for(uint32_t i = 0; i < _n; i++){
      _PLCP[i] = std::max(_LCP[_ISA[i]],_LCP[_ISA[i]+1]);
   }
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing LCP and PLCP done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t01 = std::chrono::high_resolution_clock::now();
   fprintf(stderr,"\tComputing RMQ...\n");
   _rmq = new rmq_tree((int *)_LCP, _n, 7);
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing RMQ done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   std::cerr << "Finished pre-processing\n";

   auto t2 = std::chrono::high_resolution_clock::now();
   uint64_t preprocTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Preprocessing done in " << preprocTime << " ms\n";
   return _sx;
}

int lzFactorize_omp(const std::string &_sx, std::vector<std::pair<uint32_t, int32_t>> &MSGSA, uint32_t nThreads) {
//int lzFactorize(char *fileToParse, int seqno, char* outputfilename, const bool v) {
   verbose = 0;
   omp_set_num_threads(nThreads);
   std::cerr << "_x.len: " <<  _x.size() << '\n';
   std::cerr << "_sx.len: " <<  _sx.size() << '\n';
   auto t1 = std::chrono::high_resolution_clock::now();

   unsigned int numfactors = 0;

   unsigned int inc = 100000000;
   uint64_t mark = inc;

   std::cerr << "About to start main parsing loop...\n";

   uint32_t ndoc = docBoundaries.size();
   std::cerr << "ndoc: " << ndoc << '\n';
   uint64_t *bucketLengthStar = new uint64_t[_n]();
   uint64_t *bucketLengthStarChar = new uint64_t[sizeChars]();
   uint64_t *bucketLengthsHeads = new uint64_t[_n]();

   uint64_t nStar = 0;
   headBoundaries.push_back(0);
   uint64_t i = 0;
   std::vector<std::vector<Match>> phrasesBuffer;
   phrasesBuffer.resize(docBoundaries.size());
   //uint64_t sumLen = 0, maxLen = 0;
   //uint64_t heurYes = 0, heurNo = 0;
   #pragma omp parallel for schedule(dynamic)
   for(uint32_t d = 1; d < docBoundaries.size(); d++){
      uint64_t j = docBoundaries[d-1];
      uint64_t startDoc = docBoundaries[d-1];
      uint64_t lengthSequence = docBoundaries[d];
      phrasesBuffer[d].reserve((lengthSequence-j)/10);
  
      int64_t leftB = 0;
      int64_t rightB = _n-1;
      int64_t match;
      bool isSmallerThanMatch;
      unsigned char mismatchingSymbol;
      int64_t prevPos = -2;
      uint64_t pos = 0, len = 0;
      //uint32_t iCurrentDoc = 0;
      while(j < lengthSequence){
         if(_sx[j] == sequenceSeparator){ //new doc found
            //std::cerr << "NEWDOC\n";
            leftB = 0;
            rightB = _n-1;
            len = 0;
            pos = _n - 1;
            phrasesBuffer[d].push_back(Match(j-startDoc, pos, len, 0, sequenceSeparator));
            break;
         }
         else{
            computeLZFactorAt(_sx, j, &pos, &len, leftB, rightB, match, isSmallerThanMatch, mismatchingSymbol);
            if((int64_t)pos != prevPos+1){
               phrasesBuffer[d].push_back(Match(j-startDoc, pos, len, isSmallerThanMatch, mismatchingSymbol));
            }
            //j++;
            len--;
            
            if(leftB == rightB){
               while(len > _PLCP[pos+1]){
                  j++;
                  len--;
                  pos++;
               }
               leftB = rightB = _ISA[pos];
            }
            std::pair<int,int> interval = contractLeft(leftB,rightB,len);
            leftB = interval.first;
            rightB = interval.second;
         }
         j++;
         prevPos = pos;
      }
   }
   uint64_t totalSizeCMS = 0;
   uint64_t maxValue = 0;
   for(uint32_t d = 1; d < docBoundaries.size(); d++){
      totalSizeCMS += phrasesBuffer[d].size();
      headBoundaries.push_back(totalSizeCMS);
      maxValue = std::max(maxValue, (docBoundaries[d]-docBoundaries[d-1])-1);
   }
   std::cerr << "totalSizeCMS: " << totalSizeCMS << '\n';
   delete[] _SA;
   delete[] _LCP;
   delete[] _PLCP;
   delete _rmq;
   
   std::string().swap(_x);
   std::cerr << headBoundaries.size() << ", " << ndoc << "\n";
   // docBoundaries.push_back(_sn);
   if(verbose) std::cerr << "Printing docBoundaries" << "\n";
   if(verbose) for(size_t i = 0; i < docBoundaries.size(); i++){ std::cerr << docBoundaries[i] << ", letter: " << _sx[docBoundaries[i]] << "\n";}
   auto t2 = std::chrono::high_resolution_clock::now();
   uint64_t lzFactorizeTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Time to compute matching statistics: " << lzFactorizeTime << " milliseconds\n";

   t1 = std::chrono::high_resolution_clock::now();
   std::cerr << "Start Sorting procedure for MSGSA\n";

   std::vector<SufSStar> sStar;
   std::vector<SufSStar>::iterator begStar;
   uint32_t *prefSumBucketLengthsStar = new uint32_t[_n + 1]();
   prefSumBucketLengthsStar[0] = 0;
   
   if(_n < 256*256){
      std::cerr << "Bucketing S* with small reference\n";
      std::vector<std::vector<uint64_t>> partialBucketLengthStarChar(nThreads, std::vector<uint64_t>(sizeChars,0)); //= new uint64_t[n_threads][sizeChars]();
      std::vector<std::vector<uint64_t>> partialBucketLengthStar(nThreads+1, std::vector<uint64_t>(_n,0));
      #pragma omp parallel for schedule(static)
      for(uint32_t d = 1; d < docBoundaries.size(); d++){
         const uint32_t thread_id = omp_get_thread_num();
         uint64_t partialNStar = 0;

         partialBucketLengthStar[thread_id+1][_ISA[phrasesBuffer[d][phrasesBuffer[d].size()-1].pos]]++;
         partialBucketLengthStarChar[thread_id][sequenceSeparator]++;
         partialNStar++;
         bool currentType = 0;
         Match firstHead, secondHead;
         for(int32_t i = phrasesBuffer[d].size()-1; i > 0; i--){
            secondHead = phrasesBuffer[d][i];
            firstHead = phrasesBuffer[d][i-1];
            for(int32_t start = secondHead.start - firstHead.start - 1; start >= 0; start--){
               if(_sx[firstHead.start + docBoundaries[d - 1] + start] > _sx[firstHead.start + docBoundaries[d - 1] + start + 1]){
                  currentType = 1;
               }
               else if(_sx[firstHead.start + docBoundaries[d - 1] + start] < _sx[firstHead.start + docBoundaries[d - 1] + start + 1]){
                  currentType = 0;  
               }
               if(firstHead.start + start != 0){
                  if((currentType == 0) & (_sx[firstHead.start + docBoundaries[d - 1] + start - 1] > _sx[firstHead.start + docBoundaries[d - 1] + start])){
                     partialBucketLengthStar[thread_id+1][_ISA[firstHead.pos + start]]++;
                     partialBucketLengthStarChar[thread_id][_sx[firstHead.start + docBoundaries[d - 1] + start]]++;
                     partialNStar++;
                  }
               }
            }
         }
         #pragma omp critical
         {
            nStar += partialNStar;
         }
      }
      for(uint32_t thread_id = 0; thread_id < nThreads; thread_id++){
         #pragma omp parallel for schedule(static)
         for(uint32_t i = 0; i < sizeChars; i++){
            bucketLengthStarChar[i] += partialBucketLengthStarChar[thread_id][i];
         }
         #pragma omp parallel for schedule(static)
         for(uint32_t i = 0; i < _n; i++){
            bucketLengthStar[i] += partialBucketLengthStar[thread_id+1][i];
         }
      }
      std::vector<std::vector<uint64_t>>().swap(partialBucketLengthStarChar);
      #pragma omp parallel for schedule(static)
      for(uint32_t i = 1; i < _n; i++){
         for(uint32_t thread_id = 1; thread_id < nThreads; thread_id++){
            partialBucketLengthStar[thread_id][i] += partialBucketLengthStar[thread_id-1][i];
         }
      }
      
      std::cerr << "After parallel counting for buckets\n";
      std::cerr << "nStar: " << nStar << "\n";
      sStar.resize(nStar);
      begStar = sStar.begin();
      std::cerr << "After resize\n";
      if(verbose) std::cerr << prefSumBucketLengthsStar[0] << "\n";
      for(uint32_t i = 1; i < _n; i++){
         prefSumBucketLengthsStar[i] = prefSumBucketLengthsStar[i-1] + bucketLengthStar[i-1];
      }
      prefSumBucketLengthsStar[_n] = nStar;

      #pragma omp parallel for schedule(static)
      for(uint32_t d = 1; d < docBoundaries.size(); d++){
         const uint32_t thread_id = omp_get_thread_num();
         bool currentType = 0;
         Match firstHead, secondHead;

         sStar[d-1].assign(phrasesBuffer[d][phrasesBuffer[d].size()-1].start, d, phrasesBuffer[d].size()-1+headBoundaries[d-1], 0);
         uint32_t idx;
         for(uint64_t i = phrasesBuffer[d].size()-1; i > 0; i--){
            secondHead = phrasesBuffer[d][i];
            firstHead = phrasesBuffer[d][i-1];
            for(int32_t start = secondHead.start - firstHead.start - 1; start >= 0; start--){
               if(_sx[firstHead.start + docBoundaries[d - 1] + start] > _sx[firstHead.start + docBoundaries[d - 1] + start + 1]){
                  currentType =  1;
               }
               else if(_sx[firstHead.start + docBoundaries[d - 1] + start] < _sx[firstHead.start + docBoundaries[d - 1] + start + 1]){
                  currentType = 0;
               }
               if(firstHead.start + start != 0){
                  if((currentType == 0) & (_sx[firstHead.start + docBoundaries[d - 1] + start - 1] > _sx[firstHead.start + docBoundaries[d - 1] + start])){
                     idx = prefSumBucketLengthsStar[_ISA[firstHead.pos + start]] + partialBucketLengthStar[thread_id][_ISA[firstHead.pos + start]]++;

                     sStar[idx].assign(firstHead.start + start, d, i-1+headBoundaries[d-1], start);
                  }
               }
            }
         }
      }

      prefSumBucketLengthsStar[0] = 0;
      for(uint32_t i = 1; i < _n; i++){
         prefSumBucketLengthsStar[i] = prefSumBucketLengthsStar[i-1] + bucketLengthStar[i-1];
      }
      prefSumBucketLengthsStar[_n] = nStar;
      std::vector<std::vector<uint64_t>>().swap(partialBucketLengthStar);
   }
   else{
      //ndoc = 0;
      uint32_t max = _n;
      uint32_t shift = 0;
		for(uint8_t i = 0; i < 48; i++){
			if(max < nOfDigitsBig[i]){
				shift = i;
				break;
			}
		} 
      uint8_t radixIter = 2;
      if(shift < 8){
         radixIter--;
      }
      uint64_t *preBuckets = new uint64_t[65536]();

      #pragma omp parallel for schedule(static)
      for(uint32_t d = 1; d < docBoundaries.size(); d++){
         uint32_t *partialBucketLengthStarChar = new uint32_t[sizeChars]();
         uint32_t *partialBucketLengthStar = new uint32_t[65536]();
         uint64_t partialNStar = 0;
         
         partialBucketLengthStar[_ISA[phrasesBuffer[d][phrasesBuffer[d].size()-1].pos] >> shift]++;
         partialBucketLengthStarChar[sequenceSeparator]++;
         partialNStar++;
         bool currentType = 0;
         Match firstHead, secondHead;
         for(int32_t i = phrasesBuffer[d].size()-1; i > 0; i--){
            secondHead = phrasesBuffer[d][i];
            firstHead = phrasesBuffer[d][i-1];
            for(int32_t start = secondHead.start - firstHead.start - 1; start >= 0; start--){
               if(_sx[firstHead.start + docBoundaries[d - 1] + start] > _sx[firstHead.start + docBoundaries[d - 1] + start + 1]){
                  currentType = 1;
               }
               else if(_sx[firstHead.start + docBoundaries[d - 1] + start] < _sx[firstHead.start + docBoundaries[d - 1] + start + 1]){
                  currentType = 0;
               }
               if(firstHead.start + start != 0){
                  if((currentType == 0) & (_sx[firstHead.start + docBoundaries[d - 1] + start - 1] > _sx[firstHead.start + docBoundaries[d - 1] + start])){
                     partialBucketLengthStar[_ISA[firstHead.pos + start] >> shift]++;
                     partialBucketLengthStarChar[_sx[firstHead.start + docBoundaries[d - 1] + start]]++;
                     partialNStar++;
                  }
               }
            }
         }
         #pragma omp critical
         {
            for(uint32_t i = 0; i < sizeChars; i++){
               bucketLengthStarChar[i] += partialBucketLengthStarChar[i];
            }
            for(uint32_t i = 0; i < 65536; i++){
               preBuckets[i] += partialBucketLengthStar[i];
            }
            nStar += partialNStar;
            delete[] partialBucketLengthStarChar;
            delete[] partialBucketLengthStar;
         }
      }

      sStar.resize(nStar);
      begStar = sStar.begin();
      uint64_t *prefSumPreBuckets = new uint64_t[65536];
      prefSumPreBuckets[0] = 0;
      for(uint32_t i = 1; i < 65536; i++){
         prefSumPreBuckets[i] = prefSumPreBuckets[i-1] + preBuckets[i-1];
      }

      #pragma omp parallel for schedule(static)
      for(uint32_t d = 1; d < docBoundaries.size(); d++){
         bool currentType = 0;
         Match firstHead, secondHead;
         uint32_t idx;

         #pragma omp atomic capture
         idx = prefSumPreBuckets[_ISA[phrasesBuffer[d][phrasesBuffer[d].size()-1].pos] >> shift]++;
         
         sStar[idx].assign(phrasesBuffer[d][phrasesBuffer[d].size()-1].start, d, phrasesBuffer[d].size()-1+headBoundaries[d-1], 0);
         for(uint64_t i = phrasesBuffer[d].size()-1; i > 0; i--){
            secondHead = phrasesBuffer[d][i];
            firstHead = phrasesBuffer[d][i-1];
            for(int32_t start = secondHead.start - firstHead.start - 1; start >= 0; start--){
               if(_sx[firstHead.start + docBoundaries[d - 1] + start] > _sx[firstHead.start + docBoundaries[d - 1] + start + 1]){
                  currentType =  1;
               }
               else if(_sx[firstHead.start + docBoundaries[d - 1] + start] < _sx[firstHead.start + docBoundaries[d - 1] + start + 1]){
                  currentType = 0;
               }
               if(firstHead.start + start != 0){
                  if((currentType == 0) & (_sx[firstHead.start + docBoundaries[d - 1] + start - 1] > _sx[firstHead.start + docBoundaries[d - 1] + start])){
                     #pragma omp atomic capture
                     idx = prefSumPreBuckets[_ISA[firstHead.pos + start] >> shift]++;

                     sStar[idx].assign(firstHead.start + start, d, i-1+headBoundaries[d-1], _ISA[firstHead.pos + start]);;
                  }
               }
            }
         }
      }
      std::cerr << "Bucketed S*-suffixes\n";
      prefSumPreBuckets[0] = 0;
      for(uint32_t i = 1; i < 65536; i++){
         prefSumPreBuckets[i] = prefSumPreBuckets[i-1] + preBuckets[i-1];
      }
      std::cerr << "Radix iter: " << (uint32_t)radixIter << "\n";
      #pragma omp parallel for schedule(dynamic)
      for(uint32_t i = 1; i < 65536; i++){
         if(prefSumPreBuckets[i]-prefSumPreBuckets[i-1]){
            radixSortPos(begStar, prefSumPreBuckets[i-1], prefSumPreBuckets[i], radixIter);
         }
      }
      
      for(std::vector<SufSStar>::iterator it = sStar.begin(); it < sStar.end(); it++){
         bucketLengthStar[it->diffLen.val]++;
      }
      if(verbose) std::cerr << prefSumBucketLengthsStar[0] << "\n";
      for(uint32_t i = 1; i < _n; i++){
         prefSumBucketLengthsStar[i] = prefSumBucketLengthsStar[i-1] + bucketLengthStar[i-1];
      }
      prefSumBucketLengthsStar[_n] = nStar;
      delete[] preBuckets;
      delete[] prefSumPreBuckets;
   }

   delete[] bucketLengthStar;
   phrases.reserve(totalSizeCMS);
   for(uint32_t d = 1; d < docBoundaries.size(); d++){
      phrases.insert(phrases.end(), phrasesBuffer[d].begin(), phrasesBuffer[d].end());
   }
   std::vector<std::vector<Match>>().swap(phrasesBuffer);
   std::cerr << "phrases.size() = " << phrases.size() << "\n";
   if(verbose) for(uint64_t i = 0; i < phrases.size(); i++) std::cerr << phrases[i].start << "," << phrases[i].pos << "," << phrases[i].len << "," << phrases[i].smaller << '\n';
   predecessor2 *pHeads = new predecessor2(phrases, headBoundaries, ndoc, maxValue);
   for(uint64_t i = 0; i < phrases.size(); i++){
      bucketLengthsHeads[_ISA[phrases[i].pos]]++;
   }
   std::cerr << "Prefix sum done\n";
   std::cerr << "NStar: " << nStar << "\n";

   t2 = std::chrono::high_resolution_clock::now();
   uint64_t bucketingTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Time to bucket suffixes: " << bucketingTime << " milliseconds\n";

   t1 = std::chrono::high_resolution_clock::now();
   //Sort suffixes corrisponding to heads
   uint32_t *prefSumBucketLengthsHeads = new uint32_t[_n];
   uint32_t *prefSumBucketLengthsHeadsCopy = new uint32_t[_n+1];
   prefSumBucketLengthsHeads[0] = 0;
   prefSumBucketLengthsHeadsCopy[0] = 0; 
   for(size_t i = 1; i < _n; i++){
      prefSumBucketLengthsHeads[i] = prefSumBucketLengthsHeads[i-1] + bucketLengthsHeads[i-1];
      prefSumBucketLengthsHeadsCopy[i] = prefSumBucketLengthsHeads[i];
   }
   prefSumBucketLengthsHeadsCopy[_n] = phrases.size();
   
   headsSA.resize(phrases.size());
   i = 0, ndoc = 0;
   auto t01 = std::chrono::high_resolution_clock::now();
   for(std::vector<Match>::iterator it = phrases.begin(); it < phrases.end(); it++){
      if(it->start == 0){
         ndoc++;
      }
      headsSA[prefSumBucketLengthsHeads[_ISA[it->pos]]++] = MatchSA(i++, it->pos, it->len, !it->smaller, it->next);
   }
   auto t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Bucketing heads took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   if(verbose) std::cerr << "Outputting headsSA bucketed (size=" << headsSA.size() << ")\n";
   std::vector<MatchSA>::iterator begHeads = headsSA.begin();
   t01 = std::chrono::high_resolution_clock::now();
   #pragma omp parallel for schedule(dynamic)
   for(size_t i = 1; i < _n + 1; i++){
      std::sort(begHeads + prefSumBucketLengthsHeadsCopy[i-1], begHeads + prefSumBucketLengthsHeadsCopy[i], compareMatchSA);
   }
   delete[] prefSumBucketLengthsHeads;
   delete[] prefSumBucketLengthsHeadsCopy;
   delete[] bucketLengthsHeads;
   
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Head sort took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";
   
   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to re-rank heads\n";
   int32_t *stringHeads = new int32_t[headsSA.size()+1]; 
   uint64_t rank = 1;
   int64_t prevLen = -1, prevSmall = -1, prevNext = -1;
   int64_t prevPos = -1;
   i = 0;
   for(std::vector<MatchSA>::iterator it = headsSA.begin(); it < headsSA.end(); it++){
      if(i >= headBoundaries.size()-1){
         if(it->pos != prevPos || it->len != prevLen || it->smaller != prevSmall || it->next != prevNext){
            rank++;
            prevPos = it->pos;
            prevLen = it->len;
            prevSmall = it->smaller;
            prevNext = it->next;
         }
      }
      stringHeads[i++] = rank;
   }
   std::cerr << "Re-ranked heads\n";
   if(verbose) for(uint64_t i = 0; i < phrases.size(); i++){
      std::cerr << headsSA[i].pos << "," << headsSA[i].len << " --> " << stringHeads[i] << "\n";
   }
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Reranking took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to invert stringHeads\n";
   
   int64_t hLen = headsSA.size();
   int32_t *stringHeads2 = new int32_t[headsSA.size()+1]; 
   #pragma omp parallel for
   for(uint32_t i = 0; i < headsSA.size(); i++){
      stringHeads2[headsSA[i].start] = stringHeads[i];
   }
   /*
   int temp, next;
   for (int i = 0; i < hLen; i++) {
      next = i;

      // Check if it is already
      // considered in cycle
      while(headsSA[next].start < hLen) {

         // Swap the current element according
         // to the order in headsSA
         std::swap(stringHeads[i], stringHeads[headsSA[next].start]);
         temp = headsSA[next].start;

         // Subtract size of headsSA from an entry in headsSA
         // to make it negative which indicates
         // the corresponding move
         // has been performed
         headsSA[next].start -= hLen;
         next = temp;
      }
   }
   */
   #pragma omp parallel for
   for(uint32_t i = 1; i < headsSA.size()+1; i++){
      stringHeads[i] = stringHeads2[i];
   }
   delete[] stringHeads2;
   std::vector<MatchSA>().swap(headsSA);
   stringHeads[hLen] = 0;
   std::cerr << "Inverted stringHeads\n";
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Inverting stringHeads took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   uint32_t bufferForLibSais = 60000;
   int32_t *indicesSA = new int32_t[hLen+1+bufferForLibSais]();

   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to compute suffix array of MS-heads\n";
   //sacak_int(stringHeads, indicesSA, hLen+1, rank+1);
   libsais_int_omp(stringHeads, indicesSA, hLen+1, rank+1, bufferForLibSais, nThreads);
   std::cerr << "Computed suffix array of MS-heads\n";
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing MS-heads SA took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   if(verbose) for(uint64_t i = 0; i < phrases.size()+1; i++){
      std::cerr << stringHeads[i] << ", " << indicesSA[i] << "\n";
   }
   delete[] stringHeads;
   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to permute headsSA\n";

   const std::vector<Match>::iterator begPhrases = phrases.begin();
   //auto firstHead = std::begin(phrases);
   std::cerr << "New for loop\n";
   uint64_t currentSequence = 1;
   int32_t *indicesSA2 = new int32_t[hLen+1];
   #pragma omp parallel for
   for(uint32_t i = 0; i < hLen+1; i++){
      indicesSA2[indicesSA[i]] = i;
   }
   #pragma omp parallel for
   for(uint64_t seq = 1; seq < headBoundaries.size(); seq++){
      for(uint64_t i = headBoundaries[seq-1]; i < headBoundaries[seq]; i++){
         if(phrases[i].len == 0){
            continue;
         }
         uint32_t nextStart = phrases[i].start + phrases[i].len;
         const Match headNextStart = Match(nextStart, 0, seq);
         std::vector<Match>::iterator nextHead = pHeads->predQuery(headNextStart, begPhrases);
         phrases[i].pos = _ISA[nextHead->pos + (nextStart - nextHead->start)];
         indicesSA2[i] = indicesSA2[nextHead - begPhrases] - 1;
      }
   }

   std::cerr << "New for loop end\n";
   headsSA.resize(hLen);

   #pragma omp parallel for
   for(int64_t i = 1; i < hLen+1; i++) {
      headsSA[i-1] = MatchSA(indicesSA2[indicesSA[i]], phrases[indicesSA[i]].pos, phrases[indicesSA[i]].start + phrases[indicesSA[i]].len, phrases[indicesSA[i]].smaller, phrases[indicesSA[i]].next);
   }
   std::cerr << "headBoundaries.size(): " << headBoundaries.size() << "\n";
   std::cerr << "headsSA.size(): " << headsSA.size() << "\n";
   std::vector<uint32_t>().swap(headBoundaries);
   std::cerr << "Permuted headsSA\n";
   //delete[] indicesSA;
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Permuting MS-heads took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t2 = std::chrono::high_resolution_clock::now();
   uint64_t headSortTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Time to sort heads: " << headSortTime << " milliseconds\n";

   if(verbose){
      uint64_t errHeads = checkHeadsSA_no_copy(headsSA, phrases.size(), _sx);
      std::cerr << "N. errors on headsSA: " << errHeads << "\n";
   }
   t1 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing nextPos in headSA\n";

   #pragma omp parallel for
   for(uint32_t i = 0; i < hLen+1; i++){
      indicesSA2[indicesSA[i]] = i;
   }
   #pragma omp parallel for
   for(size_t i = 0; i < nStar; i++){
      sStar[i].head = indicesSA2[sStar[i].head] - 1;
   }
   std::vector<Match>().swap(phrases);
   delete pHeads;
   delete[] indicesSA2;
   delete[] indicesSA;
   // Now: 
   // headsSA[i] = (nextHeadRank, nextHeadISA[pos], length, smaller, mmchar);

   std::cerr << "Computing change len and start in headSA DONE\n";
   t2 = std::chrono::high_resolution_clock::now();
   uint64_t changeHeadsTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Change infomation in heads in: " << changeHeadsTime << " milliseconds\n";

   std::cerr << "Starting to sort S*-suffixes\n";
   t1 = std::chrono::high_resolution_clock::now();   

   #pragma omp parallel for schedule(dynamic)
   for(uint32_t i = 1; i < _n + 1; i++){
      std::sort(begStar + prefSumBucketLengthsStar[i-1], begStar + prefSumBucketLengthsStar[i], compareSuf);
   }
   std::vector<MatchSA>().swap(headsSA);
   delete[] prefSumBucketLengthsStar;

   t2 = std::chrono::high_resolution_clock::now();
   uint64_t sortTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "S*-suffixes sorted in: " << sortTime << " milliseconds\n";
   

   if(verbose){
      std::cerr << "Going to check errors in sStar\n";
      uint64_t errSStar = checkGSA_no_copy(sStar, nStar, _sx);
      std::cerr << "N. errors on sStar: " << errSStar << "\n";
   }

   t1 = std::chrono::high_resolution_clock::now();
   uint64_t *prefSumCharBkts = new uint64_t[sizeChars + 1];
   uint64_t *prefSumCharBktsEnds = new uint64_t[sizeChars + 1];
   //uint32_t t_sum = 0;
   prefSumCharBkts[0] = 0;
   prefSumCharBktsEnds[0] = 0;
   if(verbose) std::cerr << prefSumCharBkts[0] << "\n";
   for(size_t i = 1; i < sizeChars; i++){
      prefSumCharBkts[i] = prefSumCharBkts[i-1] + charBkts[i-1];
   }
   prefSumCharBkts[sizeChars] = _sn;
   for(size_t i = 0; i < sizeChars; i++){
      prefSumCharBktsEnds[i] = prefSumCharBkts[i+1];
   }

   std::cerr << "Going to convert types\n";
         
   std::vector<std::pair<uint32_t,int32_t>> sStarSuf;
   sStarSuf.resize(nStar);
   //MSGSA.resize(nStar);
   #pragma omp parallel for
   for(uint64_t i = 0; i < nStar; i++){
      sStarSuf[i] = std::pair<uint32_t,int32_t>{sStar[i].idx, sStar[i].doc};
   }
    
   std::vector<SufSStar>().swap(sStar);
   MSGSA.reserve(_sn);
   uint64_t diffLen;
   std::reverse(sStarSuf.begin(), sStarSuf.end());

   for(uint32_t x = 1; x < sizeChars; x++){
      diffLen = (prefSumCharBkts[x] - prefSumCharBkts[x-1]) - bucketLengthStarChar[x-1];
      for(uint64_t p = 0; p < diffLen; p++){
         MSGSA.push_back(std::pair<uint32_t,int32_t>(0,0));
      }
      for(uint32_t p = 0; p < bucketLengthStarChar[x-1]; p++){
         MSGSA.push_back(sStarSuf.back());
         sStarSuf.pop_back();
      }
      if(bucketLengthStarChar[x-1]) sStarSuf.shrink_to_fit();
   }

   std::vector<std::pair<uint32_t,int32_t>>().swap(sStarSuf);
   
   delete[] bucketLengthStarChar;
   
   std::vector<uint64_t> prefSumCharBktsV(prefSumCharBkts, prefSumCharBkts + sizeChars + 1);
   delete[] prefSumCharBkts;
   psais::induce_sort<uint64_t, uint32_t, std::pair<uint32_t, int32_t>>(_sx, prefSumCharBktsV, MSGSA, docBoundaries, nThreads);

   t2 = std::chrono::high_resolution_clock::now();
   uint64_t induceTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Induced in: " << induceTime << " milliseconds\n";

   int checkGSA = 1;
   if(checkGSA){
      //exit(0);
      uint64_t err = 0;
      #pragma omp parallel for
      for(uint64_t i = 1; i < MSGSA.size(); i++){
         exit(0);
         //std::cerr << "i=" << i << ": " << MSGSA[i].first << " " << MSGSA[i].second << " " << "\n";//MSGSA[i].head <<"\n";
         uint32_t maxIdx = std::min(docBoundaries[MSGSA[i].second] - (MSGSA[i].first + docBoundaries[MSGSA[i].second - 1]), docBoundaries[MSGSA[i-1].second] - (MSGSA[i-1].first + docBoundaries[MSGSA[i-1].second - 1]));
   
         if(memcmp(&_sx[MSGSA[i].first + docBoundaries[MSGSA[i].second - 1]], &_sx[MSGSA[i-1].first + docBoundaries[MSGSA[i-1].second - 1]], maxIdx) < 0){
            std::cerr << "PROBLEM with " << i-1 << " (" << MSGSA[i-1].first << "," << MSGSA[i-1].second << ") and " << i << " (" << MSGSA[i].first << "," << MSGSA[i].second << ")\n"; 
            err++;
         }
      }
      std::cerr << "n. errors " << err << "\n"; 
   } 

   return numfactors;
}

int lzFactorize(const std::string &_sx, std::vector<std::pair<uint32_t, int32_t>> &MSGSA) {
   verbose = 0;
   std::cerr << "File is in memory\n";
   auto t1 = std::chrono::high_resolution_clock::now();

   unsigned int numfactors = 0;

   unsigned int inc = 100000000;
   uint64_t mark = inc;

   std::cerr << "About to start main parsing loop...\n";

   int64_t leftB = 0;
   int64_t rightB = _n-1;
   int64_t match;
   bool isSmallerThanMatch;
   unsigned char mismatchingSymbol;
   int64_t prevPos = -2;
   uint64_t pos = 0, len = 0;
   uint32_t ndoc = 1;
   uint32_t iCurrentDoc = 0;
   uint64_t *bucketLengthsStar = new uint64_t[_n]();
   uint64_t *bucketLengthStarChar = new uint64_t[sizeChars]();
   uint32_t *bucketLengthsHeads = new uint32_t[_n]();

   uint64_t nStar = 0;
   headBoundaries.push_back(0);
   uint64_t i = 0;
   
   uint64_t maxValue = 0;
   phrases.reserve(_sn / 10);
   while(i < _sn){
      if(verbose) std::cerr << "i: " << i << "\n";
      if(i > mark){
         fprintf(stderr,"i = %lu; lpfRuns = %ld\n",i,phrases.size());
         mark = mark + inc;
      }
      if(_sx[i] == sequenceSeparator){ //new doc found
         leftB = 0;
         rightB = _n-1;
         len = 0;
         pos = _n - 1;
         phrases.push_back(Match(iCurrentDoc, pos, len, 0, sequenceSeparator));
         headBoundaries.push_back(phrases.size());
         if(maxValue < iCurrentDoc) maxValue = iCurrentDoc;
         iCurrentDoc = 0;
         ndoc++;
      }
      else{
         computeLZFactorAt(_sx, i, &pos, &len, leftB, rightB, match, isSmallerThanMatch, mismatchingSymbol);
         if((int64_t)pos != prevPos+1){
            phrases.push_back(Match(iCurrentDoc, pos, len, isSmallerThanMatch, mismatchingSymbol));
         }
         iCurrentDoc++;
         len--;
         
         if(leftB == rightB){
            while(len > _PLCP[pos+1]){
               i++;
               iCurrentDoc++;
               len--;
               pos++;
            }
            leftB = rightB = _ISA[pos];
         }
         std::pair<int,int> interval = contractLeft(leftB,rightB,len);
         leftB = interval.first;
         rightB = interval.second;
      }
      i++;
      prevPos = pos;
   }
   delete[] _SA;
   delete[] _LCP;
   delete[] _PLCP;
   delete _rmq;
   
   phrases.shrink_to_fit();
   std::string().swap(_x);
   std::cerr << headBoundaries.size() << ", " << ndoc << "\n";
   std::cerr << "phrases.size() = " << phrases.size() << "\n";
   if(verbose) for(uint64_t i = 0; i < phrases.size(); i++) std::cerr << phrases[i].start << "," << phrases[i].pos << "," << phrases[i].len << "," << phrases[i].smaller << '\n';
   predecessor2 *pHeads = new predecessor2(phrases, headBoundaries, ndoc, maxValue);
   for(uint64_t i = 0; i < phrases.size(); i++){
      bucketLengthsHeads[_ISA[phrases[i].pos]]++;
   }
   // docBoundaries.push_back(_sn);
   if(verbose) std::cerr << "Printing docBoundaries" << "\n";
   if(verbose) for(size_t i = 0; i < docBoundaries.size(); i++){ std::cerr << docBoundaries[i] << ", letter: " << _sx[docBoundaries[i]] << "\n";}
   auto t2 = std::chrono::high_resolution_clock::now();
   uint64_t lzFactorizeTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Time to compute matching statistics: " << lzFactorizeTime << " milliseconds\n";

   t1 = std::chrono::high_resolution_clock::now();
   std::cerr << "Start Sorting procedure for MSGSA\n";

   std::vector<SufSStar> sStar;
   std::vector<SufSStar>::iterator begStar;
   uint64_t *prefSumBucketLengthsStar = new uint64_t[_n + 1]();
   prefSumBucketLengthsStar[0] = 0;
   Match firstHead, secondHead;
   ndoc = 0;
   if(_n < 256*256){
      //std::cerr << "Bucketing S* with small reference\n";
      bucketLengthsStar[_ISA[phrases[phrases.size()-1].pos]]++;
      bucketLengthStarChar[sequenceSeparator]++;
      nStar++;
      bool currentType = 0;
      ndoc = D;
      for(uint64_t i = phrases.size()-1; i > 0; i--){
         //std::cerr << i << '\n';
         secondHead = phrases[i];
         firstHead = phrases[i-1];
         if(secondHead.start == 0){
            bucketLengthsStar[_ISA[firstHead.pos]]++;
            bucketLengthStarChar[sequenceSeparator]++;
            currentType = 0;
            nStar++;
            ndoc--; 
         }
         else{
            //std::cerr << "This head wasn't followed by $\n";
            for(int64_t start = secondHead.start - firstHead.start - 1; start >= 0; start--){
               if(_sx[firstHead.start + docBoundaries[ndoc - 1] + start] > _sx[firstHead.start + docBoundaries[ndoc - 1] + start + 1]){
                  currentType =  1;
               }
               else if(_sx[firstHead.start + docBoundaries[ndoc - 1] + start] < _sx[firstHead.start + docBoundaries[ndoc - 1] + start + 1]){
                  currentType = 0;
               }
               if(firstHead.start + start != 0){
                  if((currentType == 0) & (_sx[firstHead.start + docBoundaries[ndoc - 1] + start - 1] > _sx[firstHead.start + docBoundaries[ndoc - 1] + start])){
                     bucketLengthsStar[_ISA[firstHead.pos + start]]++;
                     bucketLengthStarChar[_sx[firstHead.start + docBoundaries[ndoc - 1] + start]]++;
                     nStar++;
                  }
               }
            }
         }
      }
      std::cerr << "Buckets filled\n";
      sStar.resize(nStar);
      begStar = sStar.begin();

      if(verbose) std::cerr << prefSumBucketLengthsStar[0] << "\n";
      for(uint32_t i = 1; i < _n; i++){
         prefSumBucketLengthsStar[i] = prefSumBucketLengthsStar[i-1] + bucketLengthsStar[i-1];
      }
      prefSumBucketLengthsStar[_n] = nStar;

      sStar[prefSumBucketLengthsStar[_ISA[phrases[phrases.size()-1].pos]]++].assign(phrases[phrases.size()-1].start, D, phrases.size() - 1, 0);
      currentType = 0;
      ndoc = D;
      for(uint64_t i = phrases.size()-1; i > 0; i--){
         secondHead = phrases[i];
         firstHead = phrases[i-1];
         if(secondHead.start == 0){
            ndoc--; 
            sStar[prefSumBucketLengthsStar[_ISA[firstHead.pos]]++].assign(firstHead.start, ndoc, i-1, 0);
            currentType = 0;
         }
         else{
            for(int32_t start = secondHead.start - firstHead.start - 1; start >= 0; start--){
               if(_sx[firstHead.start + docBoundaries[ndoc - 1] + start] > _sx[firstHead.start + docBoundaries[ndoc - 1] + start + 1]){
                  currentType =  1;
               }
               else if(_sx[firstHead.start + docBoundaries[ndoc - 1] + start] < _sx[firstHead.start + docBoundaries[ndoc - 1] + start + 1]){
                  currentType = 0;
               }
               if(firstHead.start + start != 0){
                  if((currentType == 0) & (_sx[firstHead.start + docBoundaries[ndoc - 1] + start - 1] > _sx[firstHead.start + docBoundaries[ndoc - 1] + start])){
                     sStar[prefSumBucketLengthsStar[_ISA[firstHead.pos + start]]++].assign(firstHead.start + start, ndoc, i-1, start);
                  }
               }
            }
         }
      }

      prefSumBucketLengthsStar[0] = 0;
      for(uint32_t i = 1; i < _n; i++){
         prefSumBucketLengthsStar[i] = prefSumBucketLengthsStar[i-1] + bucketLengthsStar[i-1];
      }
      prefSumBucketLengthsStar[_n] = nStar;
   }
   else{
      uint32_t max = _n;
      uint32_t shift = 0;
		for(uint8_t i = 0; i < 48; i++){
			if(max < nOfDigitsBig[i]){
				shift = i;
				break;
			}
		} 
      uint8_t radixIter = 2;
      if(shift < 8){
         radixIter--;
      }
      uint64_t *preBuckets = new uint64_t[65536]();

      preBuckets[_ISA[phrases[phrases.size()-1].pos] >> shift]++;
      bucketLengthStarChar[sequenceSeparator]++;
      nStar++;
      bool currentType = 0;
      ndoc = D;
      for(uint64_t i = phrases.size()-1; i > 0; i--){
         secondHead = phrases[i];
         firstHead = phrases[i-1];
         if(secondHead.start == 0){
            preBuckets[_ISA[firstHead.pos] >> shift]++;
            bucketLengthStarChar[sequenceSeparator]++;
            currentType = 0;
            nStar++;
            ndoc--; 
         }
         else{
            for(int32_t start = secondHead.start - firstHead.start - 1; start >= 0; start--){
               if(_sx[firstHead.start + docBoundaries[ndoc - 1] + start] > _sx[firstHead.start + docBoundaries[ndoc - 1] + start + 1]){
                  currentType =  1;
               }
               else if(_sx[firstHead.start + docBoundaries[ndoc - 1] + start] < _sx[firstHead.start + docBoundaries[ndoc - 1] + start + 1]){
                  currentType = 0;
               }
               if(firstHead.start + start != 0){
                  if((currentType == 0) & (_sx[firstHead.start + docBoundaries[ndoc - 1] + start - 1] > _sx[firstHead.start + docBoundaries[ndoc - 1] + start])){
                     preBuckets[_ISA[firstHead.pos + start] >> shift]++;
                     bucketLengthStarChar[_sx[firstHead.start + docBoundaries[ndoc - 1] + start]]++;
                     nStar++;
                  }
               }
            }
         }
      }
      sStar.resize(nStar);
      begStar = sStar.begin();
      uint64_t *prefSumPreBuckets = new uint64_t[65536];
      prefSumPreBuckets[0] = 0;
      for(uint32_t i = 1; i < 65536; i++){
         prefSumPreBuckets[i] = prefSumPreBuckets[i-1] + preBuckets[i-1];
      }

      sStar[prefSumPreBuckets[_ISA[phrases[phrases.size()-1].pos] >> shift]++].assign(phrases[phrases.size()-1].start, D, phrases.size() - 1, _ISA[phrases[phrases.size()-1].pos]);
      currentType = 0;
      ndoc = D;
      for(uint64_t i = phrases.size()-1; i > 0; i--){
         secondHead = phrases[i];
         firstHead = phrases[i-1];
         if(secondHead.start == 0){
            ndoc--; 
            sStar[prefSumPreBuckets[_ISA[firstHead.pos] >> shift]++].assign(firstHead.start, ndoc, i-1, _ISA[firstHead.pos]);
            currentType = 0;
         }
         else{
            for(int32_t start = secondHead.start - firstHead.start - 1; start >= 0; start--){
               if(_sx[firstHead.start + docBoundaries[ndoc - 1] + start] > _sx[firstHead.start + docBoundaries[ndoc - 1] + start + 1]){
                  currentType =  1;
               }
               else if(_sx[firstHead.start + docBoundaries[ndoc - 1] + start] < _sx[firstHead.start + docBoundaries[ndoc - 1] + start + 1]){
                  currentType = 0;
               }
               if(firstHead.start + start != 0){
                  if((currentType == 0) & (_sx[firstHead.start + docBoundaries[ndoc - 1] + start - 1] > _sx[firstHead.start + docBoundaries[ndoc - 1] + start])){
                     sStar[prefSumPreBuckets[_ISA[firstHead.pos + start] >> shift]++].assign(firstHead.start + start, ndoc, i-1, _ISA[firstHead.pos + start]);
                  }
               }
            }
         }
      }
      std::cerr << "Bucketed S*-suffixes\n";
      prefSumPreBuckets[0] = 0;
      for(uint32_t i = 1; i < 65536; i++){
         prefSumPreBuckets[i] = prefSumPreBuckets[i-1] + preBuckets[i-1];
      }
      std::cerr << "Radix iter: " << (uint32_t)radixIter << "\n";
      
      for(uint32_t i = 1; i < 65536; i++){
         if(prefSumPreBuckets[i]-prefSumPreBuckets[i-1])
            radixSortPos(begStar, prefSumPreBuckets[i-1], prefSumPreBuckets[i], radixIter);
      }
      
      for(std::vector<SufSStar>::iterator it = sStar.begin(); it < sStar.end(); it++){
         bucketLengthsStar[it->diffLen.val]++;
      }
      if(verbose) std::cerr << prefSumBucketLengthsStar[0] << "\n";
      for(uint32_t i = 1; i < _n; i++){
         prefSumBucketLengthsStar[i] = prefSumBucketLengthsStar[i-1] + bucketLengthsStar[i-1];
      }
      prefSumBucketLengthsStar[_n] = nStar;
      delete[] preBuckets;
      delete[] prefSumPreBuckets;
   }
   delete[] bucketLengthsStar;
   std::cerr << "Prefix sum done\n";
   std::cerr << "NStar: " << nStar << "\n";

   t2 = std::chrono::high_resolution_clock::now();
   uint64_t bucketingTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Time to bucket suffixes: " << bucketingTime << " milliseconds\n";

   t1 = std::chrono::high_resolution_clock::now();
   //Sort suffixes corrisponding to heads
   uint32_t *prefSumBucketLengthsHeads = new uint32_t[_n];
   uint32_t *prefSumBucketLengthsHeadsCopy = new uint32_t[_n+1];
   prefSumBucketLengthsHeads[0] = 0;
   prefSumBucketLengthsHeadsCopy[0] = 0; 
   for(size_t i = 1; i < _n; i++){
      prefSumBucketLengthsHeads[i] = prefSumBucketLengthsHeads[i-1] + bucketLengthsHeads[i-1];
      prefSumBucketLengthsHeadsCopy[i] = prefSumBucketLengthsHeads[i];
   }
   prefSumBucketLengthsHeadsCopy[_n] = phrases.size();
   
   headsSA.resize(phrases.size());
   i = 0, ndoc = 0;
   auto t01 = std::chrono::high_resolution_clock::now();
   for(std::vector<Match>::iterator it = phrases.begin(); it < phrases.end(); it++){
      if(it->start == 0){
         ndoc++;
      }
      headsSA[prefSumBucketLengthsHeads[_ISA[it->pos]]++] = MatchSA(i++, it->pos, it->len, !it->smaller, it->next);
   }
   auto t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Bucketing heads took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   if(verbose) std::cerr << "Outputting headsSA bucketed (size=" << headsSA.size() << ")\n";
   std::vector<MatchSA>::iterator begHeads = headsSA.begin();
   t01 = std::chrono::high_resolution_clock::now();
   
   for(size_t i = 1; i < _n + 1; i++){
      std::sort(begHeads + prefSumBucketLengthsHeadsCopy[i-1], begHeads + prefSumBucketLengthsHeadsCopy[i], compareMatchSA);
   }
   delete[] prefSumBucketLengthsHeads;
   delete[] prefSumBucketLengthsHeadsCopy;
   delete[] bucketLengthsHeads;
   
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Head sort took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";
   
   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to re-rank heads\n";
   int32_t *stringHeads = new int32_t[headsSA.size()+1]; 
   uint64_t rank = 1;
   int64_t prevLen = -1, prevSmall = -1, prevNext = -1;
   prevPos = -1;
   i = 0;
   for(std::vector<MatchSA>::iterator it = headsSA.begin(); it < headsSA.end(); it++){
      if(i >= headBoundaries.size()-1){
         if(it->pos != prevPos || it->len != prevLen || it->smaller != prevSmall || it->next != prevNext){
            rank++;
            prevPos = it->pos;
            prevLen = it->len;
            prevSmall = it->smaller;
            prevNext = it->next;
         }
      }
      stringHeads[i++] = rank;
   }
   std::cerr << "Re-ranked heads\n";
   if(verbose) for(uint64_t i = 0; i < phrases.size(); i++){
      std::cerr << headsSA[i].pos << "," << headsSA[i].len << " --> " << stringHeads[i] << "\n";
   }
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Reranking took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to invert stringHeads\n";
   int64_t hLen = headsSA.size();
   /*
   int temp, next;
   for (int i = 0; i < hLen; i++) {
      next = i;

      // Check if it is already
      // considered in cycle
      while(headsSA[next].start < hLen) {

         // Swap the current element according
         // to the order in headsSA
         std::swap(stringHeads[i], stringHeads[headsSA[next].start]);
         temp = headsSA[next].start;

         // Subtract size of headsSA from an entry in headsSA
         // to make it negative which indicates
         // the corresponding move
         // has been performed
         headsSA[next].start -= hLen;
         next = temp;
      }
   }
   */
   int32_t *stringHeads2 = new int32_t[headsSA.size()+1]; 
   for(uint32_t i = 0; i < headsSA.size(); i++){
      stringHeads2[headsSA[i].start] = stringHeads[i];
   }
   for(uint32_t i = 1; i < headsSA.size()+1; i++){
      stringHeads[i] = stringHeads2[i];
   }
   delete[] stringHeads2;
   std::vector<MatchSA>().swap(headsSA);
   stringHeads[hLen] = 0;
   std::cerr << "Inverted stringHeads\n";
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Inverting stringHeads took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   uint32_t bufferLibSais = 60000;
   int32_t *indicesSA = new int32_t[hLen+1+bufferLibSais]();

   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to compute suffix array of MS-heads\n";
   libsais_int(stringHeads, indicesSA, hLen+1, rank+1, bufferLibSais);
   //sacak_int(stringHeads, indicesSA, hLen+1, rank+1);
   std::cerr << "Computed suffix array of MS-heads\n";
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing MS-heads SA took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   if(verbose) for(uint64_t i = 0; i < phrases.size()+1; i++){
      std::cerr << stringHeads[i] << ", " << indicesSA[i] << "\n";
   }
   delete[] stringHeads;
   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to permute headsSA\n";
   headsSA.resize(hLen);

   const std::vector<Match>::iterator begPhrases = phrases.begin();
   std::cerr << "New for loop\n";
   uint64_t currentSequence = 1;
   int32_t *indicesSA2 = new int32_t[hLen+1];
   for(uint32_t i = 0; i < hLen+1; i++){
      indicesSA2[indicesSA[i]] = i;
   }
   for(size_t i = 0; i < headsSA.size(); i++){
      if(phrases[i].len == 0){
         currentSequence++;
         continue;
      }
      uint32_t nextStart = phrases[i].start + phrases[i].len;
      const Match headNextStart = Match(nextStart, 0, currentSequence);
      std::vector<Match>::iterator nextHead = pHeads->predQuery(headNextStart, begPhrases);
      phrases[i].pos = _ISA[nextHead->pos + (nextStart - nextHead->start)];
      indicesSA2[i] = indicesSA2[nextHead - begPhrases] - 1;
   }
   std::cerr << "New for loop end\n";

   for(int64_t i = 1; i < hLen+1; i++) {
      //headsSA[i-1] = MatchSA(indicesSA[i], phrases[indicesSA[i]].pos, std::upper_bound(std::begin(headBoundaries), std::end(headBoundaries), indicesSA[i]) - std::begin(headBoundaries), phrases[indicesSA[i]].smaller, phrases[indicesSA[i]].next);
      headsSA[i-1] = MatchSA(indicesSA2[indicesSA[i]], phrases[indicesSA[i]].pos, phrases[indicesSA[i]].start + phrases[indicesSA[i]].len, phrases[indicesSA[i]].smaller, phrases[indicesSA[i]].next);
      //if(verbose) std::cerr << headsSA[i-1].start << "," << _ISA[headsSA[i-1].pos] << "," << headsSA[i-1].len << "," << headsSA[i-1].smaller << "," << headsSA[i-1].next <<"\n";
   }
   std::cerr << "headBoundaries.size(): " << headBoundaries.size() << "\n";
   std::vector<uint32_t>().swap(headBoundaries);
   std::cerr << "Permuted headsSA\n";
   //delete[] indicesSA;
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Permuting MS-heads took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t2 = std::chrono::high_resolution_clock::now();
   uint64_t headSortTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Time to sort heads: " << headSortTime << " milliseconds\n";

   t1 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing nextPos in headSA\n";
   
   for(uint32_t i = 0; i < hLen+1; i++){
      //std::cerr << "indicesSA[" << i << "] = " << indicesSA[i] << "\n";
      indicesSA2[indicesSA[i]] = i;
   }
   for(size_t i = 0; i < nStar; i++){
      //sStar[i].diffLen.val = sStar[i].idx - phrases[sStar[i].head].start;
      //if( sStar[i].diffLen.val) std::cerr << "diffLen: " << sStar[i].diffLen.val << " idx: " << sStar[i].idx << " start: " << phrases[sStar[i].head].start <<"\n";
      sStar[i].head = indicesSA2[sStar[i].head] - 1;
   }
   std::vector<Match>().swap(phrases);
   delete pHeads;
   delete[] indicesSA;
   delete[] indicesSA2;
   // Now: 
   // headsSA[i] = (nextHeadRank, nextHeadISA[pos], length, smaller, mmchar);

   std::cerr << "Computing change len and start in headSA DONE\n";
   t2 = std::chrono::high_resolution_clock::now();
   uint64_t changeHeadsTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Change infomation in heads in: " << changeHeadsTime << " milliseconds\n";

   std::cerr << "Starting to sort S*-suffixes\n";
   t1 = std::chrono::high_resolution_clock::now();   

   for(uint32_t i = 1; i < _n + 1; i++){
      std::sort(begStar + prefSumBucketLengthsStar[i-1], begStar + prefSumBucketLengthsStar[i], compareSuf);
   }
   std::vector<MatchSA>().swap(headsSA);
   delete[] prefSumBucketLengthsStar;

   t2 = std::chrono::high_resolution_clock::now();
   uint64_t sortTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "S*-suffixes sorted in: " << sortTime << " milliseconds\n";
   

   if(verbose){
      uint64_t errSStar = checkGSA_no_copy(sStar, nStar, _sx);
      std::cerr << "N. errors on sStar: " << errSStar << "\n";
   }
   //exit(0);

   t1 = std::chrono::high_resolution_clock::now();
   uint64_t *prefSumCharBkts = new uint64_t[sizeChars + 1];
   uint64_t *prefSumCharBktsEnds = new uint64_t[sizeChars + 1];
   //uint32_t t_sum = 0;
   prefSumCharBkts[0] = 0;
   prefSumCharBktsEnds[0] = 0;
   if(verbose) std::cerr << prefSumCharBkts[0] << "\n";
   for(size_t i = 1; i < sizeChars; i++){
      //t_sum += bucketLengths[i-1];
      prefSumCharBkts[i] = prefSumCharBkts[i-1] + charBkts[i-1];
      //if(verbose) std::cerr << prefSumBucketLengths[i] << "\n";
   }
   prefSumCharBkts[sizeChars] = _sn;
   for(size_t i = 0; i < sizeChars; i++){
      prefSumCharBktsEnds[i] = prefSumCharBkts[i+1];
   }

   std::cerr << "Going to convert types\n";
   
   std::vector<std::pair<uint32_t,int32_t>> sStarSuf;
   sStarSuf.resize(nStar);
   //MSGSA.resize(nStar);
   for(uint64_t i = 0; i < nStar; i++){
      sStarSuf[i] = {sStar.back().idx, sStar.back().doc};
      sStar.pop_back();
      //sStarSuf.push_back(std::pair<uint32_t,int32_t>{sStar[i].idx, sStar[i].doc});
      //MSGSA[i] = std::pair<uint32_t,int32_t>{sStar[i].idx, sStar[i].doc};
   }
    
   std::vector<SufSStar>().swap(sStar);
   MSGSA.reserve(_sn);
   uint64_t diffLen;
   
   for(uint32_t x = 1; x < sizeChars; x++){
      diffLen = (prefSumCharBkts[x] - prefSumCharBkts[x-1]) - bucketLengthStarChar[x-1];
      MSGSA.insert(MSGSA.end(), diffLen, std::pair<uint32_t,int32_t>(0,0));
      for(uint64_t p = 0; p < bucketLengthStarChar[x-1]; p++){
         MSGSA.push_back(sStarSuf.back());
         sStarSuf.pop_back();
      }
      if(bucketLengthStarChar[x-1]) sStarSuf.shrink_to_fit();
   }

   std::vector<std::pair<uint32_t,int32_t>>().swap(sStarSuf);
   
   delete[] bucketLengthStarChar;
   
   std::cerr << "Start inducing L-types\n";
   data_type c1 = _sx[_sn-1];
   data_type c0;
   std::pair<uint32_t, int32_t> j;
   std::vector<std::pair<uint32_t, int32_t>>::iterator begGSA = MSGSA.begin();
   std::vector<std::pair<uint32_t, int32_t>>::iterator b = MSGSA.begin() + prefSumCharBkts[c1];
   *b++ = ((0 < _sn-1) && (_sx[_sn - 2] < c1)) ? std::pair<uint32_t, int32_t>(MSGSA[0].first, ~MSGSA[0].second) : MSGSA[0];
   for(uint64_t i = 0; i < _sn; ++i){
      if(MSGSA[i].first > 0){
         j = MSGSA[i]; 
         MSGSA[i].second = ~MSGSA[i].second;

         if(0 < j.second){
            if((c0 = _sx[j.first + docBoundaries[j.second - 1] - 1]) != c1) { prefSumCharBkts[c1] = b - begGSA; b = begGSA + prefSumCharBkts[c1 = c0]; }
            *b++ = ((0 < j.first - 1) && (_sx[j.first + docBoundaries[j.second - 1] - 2] < c1)) ? std::pair<uint32_t, int32_t>(j.first - 1, ~j.second) : std::pair<uint32_t, int32_t>(j.first - 1, j.second);
         }
      }
      //__builtin_prefetch(&_sx[MSGSA[i+(i<_sn)].first + docBoundaries[MSGSA[i+(i<_sn)].second - 1] - 1]);
   }
   std::cerr << "\nAfter inducing L-types\n";
   
   std::cerr << "Start inducing S-types\n";
   b = MSGSA.begin() + prefSumCharBktsEnds[c1 = 0];
   for(uint64_t i = _sn - 1; i < _sn; i--){
      if(MSGSA[i].first > 0){
         if(0 < MSGSA[i].second){
            j = MSGSA[i];
            if((c0 = _sx[j.first + docBoundaries[j.second - 1] - 1]) != c1) { prefSumCharBktsEnds[c1] = b - begGSA; b = begGSA + prefSumCharBktsEnds[c1 = c0]; }
            *--b = ((0 == j.first - 1) || (_sx[j.first + docBoundaries[j.second - 1] - 2] > c1)) ? std::pair<uint32_t, int32_t>(j.first - 1, ~j.second) : std::pair<uint32_t, int32_t>(j.first - 1, j.second);
         }
         else{
            MSGSA[i].second = ~MSGSA[i].second;
         }
      }
      else if(0 > MSGSA[i].second) MSGSA[i].second = ~MSGSA[i].second;; 
      //__builtin_prefetch(&_sx[MSGSA[i-(i>0)].first + docBoundaries[MSGSA[i-(i>0)].second - 1] - 1]);
   }
   std::cerr << "\nAfter inducing S-types\n";
   t2 = std::chrono::high_resolution_clock::now();
   uint64_t induceTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Induced in: " << induceTime << " milliseconds\n";

   int checkGSA = 0;
   if(checkGSA){
      uint64_t err = 0;
      //#pragma omp parallel for
      for(uint64_t i = 1; i < MSGSA.size(); i++){
         uint32_t maxIdx = std::min(docBoundaries[MSGSA[i].second] - (MSGSA[i].first + docBoundaries[MSGSA[i].second - 1]), docBoundaries[MSGSA[i-1].second] - (MSGSA[i-1].first + docBoundaries[MSGSA[i-1].second - 1]));
   
         if(memcmp(&_sx[MSGSA[i].first + docBoundaries[MSGSA[i].second - 1]], &_sx[MSGSA[i-1].first + docBoundaries[MSGSA[i-1].second - 1]], maxIdx) < 0){
            std::cerr << "PROBLEM with " << i-1 << " (" << MSGSA[i-1].first << "," << MSGSA[i-1].second << ") and " << i << " (" << MSGSA[i].first << "," << MSGSA[i].second << ")\n"; 
            err++;
         }
      }
      std::cerr << "n. errors " << err << "\n"; 
   } 

   return numfactors;
}

void computeGSA(char* refFileName, char* collFileName, uint64_t prefixLength, uint32_t nThreads, std::vector<std::pair<uint32_t, int32_t>> &MSGSA){
   if(nThreads > 1){
      const std::string _sx = lzInitialize_omp(refFileName, collFileName, prefixLength, nThreads);
      lzFactorize_omp(_sx, MSGSA, nThreads);
   }
   else{
      const std::string _sx = lzInitialize(refFileName, collFileName, prefixLength);
      lzFactorize(_sx, MSGSA);
   }
}

void computeLZFactorAt(const std::string &_sx, filelength_type i, filelength_type *pos, filelength_type *len, int64_t & leftB, int64_t & rightB, int64_t & maxMatch, bool & isSmallerThanMaxMatch, unsigned char &mismatchingSymbol) {
   filelength_type offset = *len;
   filelength_type j = i + offset;

   //int64_t nlb = 0, nrb = _n - 1;
   int64_t nlb = leftB, nrb = rightB;
   unsigned int match = _SA[nlb];
   while (j < _sn) { //scans the string from j onwards until a maximal prefix is found between _x and _sx
      if (nlb == nrb) { 
         //std::cout << "Search finished with nlb == nrb " << nlb << " " << nrb << "\n";
         if (_x[_SA[nlb] + offset] != _sx[j]) {
            //std::cout << "mismatch " << _x[_SA[nlb] + offset] << " and " << _sx[j] << "\n";
            //fprintf(stderr,"Breaking from 1\n");
            isSmallerThanMaxMatch = (_x[_SA[nlb] + offset] > _sx[j]);
            mismatchingSymbol = _sx[j];
            break;
         }
         leftB = nlb;
         rightB = nrb;
         maxMatch = nlb;
      } else { //refining the bucket in which the match is found, from left and then from right
         //renormalizations++;
         nlb = binarySearchLB(nlb, nrb, offset, _sx[j]);
         if (nlb < 0) {
            //no match, the game is up
            //fprintf(stderr,"Breaking from 2; offset = %lu; _sx[%lu] = %u\n",offset,j,_sx[j]);
            maxMatch = -(nlb)-1;
            isSmallerThanMaxMatch = true;
            mismatchingSymbol = _sx[j];
            if(maxMatch == nrb+1){
               maxMatch--;
               isSmallerThanMaxMatch = false;
            }
            match = _SA[maxMatch];
            break;
         }
         nrb = binarySearchRB(nlb, nrb, offset, _sx[j]);
         leftB = nlb;
         rightB = nrb;
      }
      //std::cerr << "After if nlb: " << nlb << "\n";
      match = _SA[nlb];
      j++;
      offset++;
   }
   *pos = match;
   *len = offset;
   //std::cout << "Match " << match << "\n";
   //std::cout << "Len " << offset << "\n";
}

//Returns the leftmost occurrence of the element if it is present or (if not present)
//then returns -(x+1) where x is the index at which the key would be inserted into the 
//array: i.e., the index of the first element greater than the key, or hi+1 if all elements 
//in the array are less than the specified key.
inline int64_t binarySearchLB(int64_t lo, int64_t hi, filelength_type offset, data_type c) {
   int64_t low = lo, high = hi;
   while (low <= high) {
      int64_t mid = (low + high) >> 1;
      data_type midVal = _x[_SA[mid] + offset];
      if (midVal < c)
         low = mid + 1;
      else if (midVal > c)
         high = mid - 1;
      else { //midVal == c
         if (mid == lo)
            return mid; // leftmost occ of key found
         data_type midValLeft = _x[_SA[mid - 1] + offset];
         if (midValLeft == midVal) {
            high = mid - 1; //discard mid and the ones to the right of mid
         } else { //midValLeft must be less than midVal == c
            return mid; //leftmost occ of key found
         }
      }
   }
   return -(low + 1);  // key not found.
}

inline int64_t binarySearchRB(int64_t lo, int64_t hi, filelength_type offset, data_type c) {
   int64_t low = lo, high = hi;
   while (low <= high) {
      int64_t mid = (low + high) >> 1;
      data_type midVal = _x[_SA[mid] + offset];
      if (midVal < c)
         low = mid + 1;
      else if (midVal > c)
         high = mid - 1;
      else { //midVal == c
         if (mid == hi)
            return mid; // rightmost occ of key found
         data_type midValRight = _x[_SA[mid + 1] + offset];
         if (midValRight == midVal) {
            low = mid + 1; //discard mid and the ones to the left of mid
         } else { //midValRight must be greater than midVal == c
            return mid; //rightmost occ of key found
         }
      }
   }
   return -(low + 1);  // key not found.
}


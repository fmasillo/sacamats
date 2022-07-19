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

#include "sacamats.h"
#include "rmq_tree.h"
#include "utils.h"
#include "match.h"
#include "predecessor.h"
#include "gsa-is/gsacak.h"
#include "sais-lite-2.4.1/sais.h"

#undef max
//#include "inplace-radixxx/inplace_radixxx.h"

//#define likely(x)       __builtin_expect((x),1)

void computeLZFactorAt(const std::string &_sx, filelength_type i, filelength_type *pos, filelength_type *len, int64_t &leftB, int64_t &rightB, int64_t &match, bool &isSmaller, unsigned char &mismatchingSymbol);
inline int64_t binarySearchLB(int64_t lo, int64_t hi, filelength_type offset, data_type c);
inline int64_t binarySearchRB(int64_t lo, int64_t hi, filelength_type offset, data_type c);

static std::string outputFile;

//static data_type *_x;
static std::string _x;
static int *_SA;
static uint32_t *_ISA;
static uint32_t *_LCP;
static uint32_t *_PLCP;
static rmq_tree *_rmq;
static uint32_t _n;
//static data_type *_sx;
static filelength_type _sn;
static bool _isMismatchingSymbolNeeded;
static std::vector<uint32_t> docBoundaries;
static std::vector<uint32_t> headBoundaries;
static predecessor2 pHeads;

uint64_t tiesCounter = 0;

uint64_t _maxLCP = 0;
uint64_t _numberOfShortFactors = 0;
uint64_t _numberOfSingleMatches = 0;

int renormalizations = 0;

uint64_t maxFactorLength = 0;
int lenZeroFactors = 0;

bool verbose;
uint64_t maxCounter = 0;
uint64_t denCounter = 0;
uint64_t sumCounter = 0;
uint64_t diffLenCounter = 0;
uint64_t finalSuffCounter = 0;

uint64_t diffSufPos;
uint64_t diffSufLen;
uint64_t diffSufHeads;

uint16_t sizeChars = 256;
uint64_t D = 0;

//std::vector<std::pair<uint32_t, int32_t>> MSGSA;

//std::vector<std::pair<uint32_t,uint32_t> > phrases;
std::vector<Match> phrases;
std::vector<MatchSA> headsSA;

//LCP array construction method of J. Kärkkäinen, G. Manzini, and S. J. Puglisi, CPM 2009
void constructLCP(std::string t, int32_t n, int32_t *sa, uint32_t *lcp, uint32_t *temp) {
//void constructLCP(unsigned char *t, int32_t n, uint32_t *sa, uint32_t *lcp, uint32_t *temp) {
   fprintf(stderr,"\tComputing LCP...\n");
   int32_t *phi = (int32_t *)lcp, *plcp = (int32_t *)temp, l = 0;
   for (int32_t i = 1; i < n; ++i)
     phi[sa[i]] = sa[i-1];
   phi[sa[0]] = -1;
   for (int32_t i = 0; i < n; ++i) {
     int32_t j = phi[i];
     if (j == -1) { plcp[i] = 0; continue; }
     else {
       while (i + l < n && j + l < n && t[i + l] == t[j + l]) ++l;
       plcp[i] = l;
       l = std::max(l - 1, 0);
     }
   }
   for (int32_t i = 0; i < n; ++i)
     lcp[i] = plcp[sa[i]];
}

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
  }else{
    //psv++;
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
   if(headA.len - a.diffLen.val != headB.len - b.diffLen.val){
      //diffSufLen++;
      return headA.smaller*((headA.len - a.diffLen.val) < (headB.len - b.diffLen.val)) + 
            !headB.smaller*((headA.len - a.diffLen.val) > (headB.len - b.diffLen.val));
   }
   else{
      return (headA.pos != headB.pos)*(headA.pos < headB.pos) + (headA.pos == headB.pos)*(headA.start<headB.start);
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

uint64_t checkGSA(std::vector<SufSStar> GSA, uint64_t n, const std::string &_sx){
   uint64_t err = -1;
   #pragma omp parallel for
   for(size_t i = 0; i < n; i++){
       if(verbose) std::cerr << "i=" << i << ": " << GSA[i].idx << " " << GSA[i].doc << " " << "\n";//MSGSA[i].head <<"\n";
       if(GSA[i].doc == 0 || GSA[i-1].doc == 0){
          std::cerr << "There was an empty entry\n";
          #pragma omp atomic
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
          _slice_prev = "$";
          maxIdx = 1;
       }
       if(verbose) std::cerr << "suf_i-1 " << _slice_prev;
       if(verbose) std::cerr << "suf_i " << _slice_sx << "\n";
       
      if(memcmp(&_slice_sx[0], &_slice_prev[0], maxIdx) < 0){
         std::cerr << "PROBLEM with " << i-1 << " (" << GSA[i-1].idx << "," << GSA[i-1].doc << ") and " << i << " (" << GSA[i].idx << "," << GSA[i].doc << ")\n"; 
         #pragma omp atomic
         err++;
      }
    }
    return err;
}

const std::string lzInitialize(char *refFileName, char *collFileName) {
//void lzInitialize(data_type *ax, unsigned int an, bool isMismatchingSymbolNeeded, char *refFileName, char *collFileName) {
   auto t1 = std::chrono::high_resolution_clock::now();
   errno = 0;
   FILE *infileRef = fopen(refFileName, "r");
   //FILE *infile = fopen("~/Desktop/Simon/data/themisto_data/genomes64.concat.da.dict.16K.1K", "r");
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
               //std::cerr << _x << "size: " << _x.size() << "\n";
               std::string().swap(content);
            }
            else if (!line.empty()) {
               content += line;
            }
         }
         if(content.size()) _x += content;
         //std::cerr << _x << "size: " << _x.size() << "\n";
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
   sn = ftello(infile) / sizeof(data_type);
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
               sx += "%";
               //std::cerr << _x << "size: " << _x.size() << "\n";
               std::string().swap(content);
            }
            else if (!line.empty()) {
               content += line;
            }
         }
         if(content.size()){
            sx += content;
            sx += "%";
         }
         //std::cerr << _x << "size: " << _x.size() << "\n";
         std::string().swap(content);
         sx.resize(sx.size());
         sn = sx.size();
         streamInfile.close();
      }
      else{
         std::cerr << "No FASTA\n";
         fseek(infile, 0, SEEK_SET);
         sx.resize(sn);
         if (sn != fread(&sx[0], sizeof(data_type), sn, infile)) {
            fprintf(stderr, "Error reading %u bytes from file %s\n", n, collFileName);
            exit(1);
         }
      }
   }
   else{
      std::cerr << "Collection file is empty!\n";
      exit(1);
   }

   // fseek(infile, 0L, SEEK_SET);
   // data_type *sx = new data_type[sn + 1];
   // if(sn != fread(sx, sizeof(data_type), sn, infile)){
   //    fprintf(stderr, "Error reading %lu bytes from file %s\n", sn, collFileName);
   //    exit(1);
   // }
   // sx[sn] = 1; //I don't think there is any reason to do this
   fclose(infile);

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
      if(sx[i] == '%'){
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
      if((maxRunsReference[i] == 0) & (maxRunsCollection[i] != 0) & (i != '%')){
         for(uint64_t x = 0; x < maxRunsCollection[i]; x++){
            _x += (char)i;
         }
      }
   }
   _x += '$';
   _x += (char)0;
   //docBoundaries.reserve(maxRunsCollection['%']);
   headBoundaries.reserve(maxRunsCollection['%']);
   //std::cerr << refAug << "\n";
   auto t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Augmenting reference done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t01 = std::chrono::high_resolution_clock::now();
   int *sa = new int[_x.size()];
   sais(reinterpret_cast<unsigned char*>(const_cast<char*>(_x.c_str())), sa, _x.size());
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing SA done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   //std::cerr << sx << "\n";
   //_sx = reinterpret_cast<unsigned char*>(const_cast<char*>(sx.c_str()));
   const std::string _sx{sx};
   if(_sx[sn-1] != '%'){
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
   _PLCP = new uint32_t[_n];
   _LCP = new uint32_t[_n];
   t01 = std::chrono::high_resolution_clock::now();
   constructISA(_SA,_ISA,_n);
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing ISA done in " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   constructLCP(_x,_n,_SA,_LCP,_PLCP);
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
   //sprintf(collFileName, "%s.gsa", collFileName);
   //outputFile = std::string(collFileName);
   return _sx;
}

int lzFactorize(const std::string _sx, std::vector<std::pair<uint32_t, int32_t>> &MSGSA) {
//int lzFactorize(char *fileToParse, int seqno, char* outputfilename, const bool v) {
   verbose = 0;
   //omp_set_num_threads(4);
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
   //uint32_t *bucketLengths = new uint32_t[_n]();
   uint32_t *bucketLengthsStar = new uint32_t[_n]();
   //uint32_t *bucketLengthsStar = new uint32_t[sizeChars]();
   uint32_t *bucketLengthStarChar = new uint32_t[sizeChars]();
   uint32_t *bucketLengthsHeads = new uint32_t[_n]();
   //uint32_t *bucketLengthsHeads = new uint32_t[sizeChars]();

   uint64_t nStar = 0;
   headBoundaries.push_back(0);
   uint64_t i = 0;
   uint64_t *charBkts = new uint64_t[sizeChars]();
   uint64_t maxValue = 0;
   phrases.reserve(_sn / 10);
   uint64_t sumLen = 0, maxLen = 0;
   //uint64_t heurYes = 0, heurNo = 0;
   while(i < _sn){
      //std::cerr << "i: " << i << "\n";
      //std::cerr << _sx[i] << "\n";
      if(verbose) std::cerr << "i: " << i << "\n";
      if(i > mark){
         fprintf(stderr,"i = %lu; lpfRuns = %ld\n",i,phrases.size());
         mark = mark + inc;
      }
      sumLen += len;
      if(len > maxLen) maxLen = len;
      charBkts[_sx[i]]++;
      if(_sx[i] == '%'){ //new doc found
         //std::cerr << "NEWDOC\n";
         leftB = 0;
         rightB = _n-1;
         len = 0;
         pos = _n - 1;
         phrases.push_back(Match(iCurrentDoc, pos, len, 0));
         headBoundaries.push_back(phrases.size());
         if(maxValue < iCurrentDoc) maxValue = iCurrentDoc;
         iCurrentDoc = 0;
         ndoc++;
      }
      else{
         computeLZFactorAt(_sx, i, &pos, &len, leftB, rightB, match, isSmallerThanMatch, mismatchingSymbol);
         if((int64_t)pos != prevPos+1){
            phrases.push_back(Match(iCurrentDoc, pos, len, isSmallerThanMatch));
         }
         iCurrentDoc++;
         len--;
         
         if(leftB == rightB){
            while(len > _PLCP[pos+1]){
               //sumLen += len;
               //if(len > maxLen) maxLen = len;
               i++;
               charBkts[_sx[i]]++;
               iCurrentDoc++;
               len--;
               pos++;
               //heurYes++;
            }
            leftB = rightB = _ISA[pos];
         }
         std::pair<int,int> interval = contractLeft(leftB,rightB,len);
         leftB = interval.first;
         rightB = interval.second;
         //heurNo++;
      }
      i++;
      prevPos = pos;
   }
   delete[] _SA;
   delete[] _LCP;
   delete[] _PLCP;
   
   //std::cerr << "heurYes " << heurYes << ", heurNo " << heurNo << "\n";
   //std::cerr << "maxLen " << maxLen << " meanLen " << (float)sumLen/_sn << "\n";
   phrases.shrink_to_fit();
   std::string().swap(_x);
   std::cerr << headBoundaries.size() << ", " << ndoc << "\n";
   std::cerr << "phrases.size() = " << phrases.size() << "\n";
   if(verbose) for(uint64_t i = 0; i < phrases.size(); i++) std::cerr << phrases[i].start << "," << phrases[i].pos << "," << phrases[i].len << "," << phrases[i].smaller << '\n';
   pHeads = predecessor2(phrases, headBoundaries, ndoc, maxValue);
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
   uint32_t *prefSumBucketLengthsStar = new uint32_t[_n + 1]();
   prefSumBucketLengthsStar[0] = 0;
   Match firstHead, secondHead;
   ndoc = 0;
   if(_n < 256*256){
      bucketLengthsStar[_ISA[phrases[phrases.size()-1].pos]]++;
      bucketLengthStarChar['%']++;
      nStar++;
      bool currentType = 0;
      ndoc = D;
      for(uint64_t i = phrases.size()-2; i > 0; i--){
         secondHead = phrases[i];
         firstHead = phrases[i-1];
         if(secondHead.start == 0){
            bucketLengthsStar[_ISA[firstHead.pos]]++;
            bucketLengthStarChar['%']++;
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
                     bucketLengthsStar[_ISA[firstHead.pos + start]]++;
                     bucketLengthStarChar[_sx[firstHead.start + docBoundaries[ndoc - 1] + start]]++;
                     nStar++;
                  }
               }
            }
         }
      }
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
      for(uint64_t i = phrases.size()-2; i > 0; i--){
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
      bucketLengthStarChar['%']++;
      nStar++;
      bool currentType = 0;
      ndoc = D;
      for(uint64_t i = phrases.size()-2; i > 0; i--){
         secondHead = phrases[i];
         firstHead = phrases[i-1];
         if(secondHead.start == 0){
            preBuckets[_ISA[firstHead.pos] >> shift]++;
            bucketLengthStarChar['%']++;
            currentType = 0;
            nStar++;
            ndoc--; 
            //std::cerr << "doc: " << ndoc << "\n";
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
      for(uint64_t i = phrases.size()-2; i > 0; i--){
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
      //#pragma omp parallel for schedule(dynamic)
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

   //put $ instead of X, otherwise the X characters does not lead to a correct comparison (because they are greater)
   //for(size_t i = 0; i < _sn; i++) {if(_sx[i] == '%') _sx[i] = "$";}

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
   //headsSA.reserve(phrases.size());
   i = 0, ndoc = 0;
   auto t01 = std::chrono::high_resolution_clock::now();
   for(std::vector<Match>::iterator it = phrases.begin(); it < phrases.end(); it++){
      if(it->start == 0){
         ndoc++;
      }
      headsSA[prefSumBucketLengthsHeads[_ISA[it->pos]]++] = MatchSA(i++, it->pos, it->len, !it->smaller, _sx[it->start + docBoundaries[ndoc-1] + it->len]);
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
   uint32_t *stringHeads = new uint32_t[headsSA.size()+1]; 
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
   std::vector<MatchSA>().swap(headsSA);
   stringHeads[hLen] = 0;
   std::cerr << "Inverted stringHeads\n";
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Inverting stringHeads took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   uint32_t *indicesSA = new uint32_t[hLen+1]();
   // for(uint64_t i = 0; i < phrases.size(); i++){
   //    headsSufArr[i] = phrases[i].pos << 31 + phrases[i].len;
   // }
   t01 = std::chrono::high_resolution_clock::now();
   std::cerr << "Going to compute suffix array of MS-heads\n";
   sacak_int(stringHeads, indicesSA, hLen+1, rank+1);
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
   for(int64_t i = 1; i < hLen+1; i++) {
      headsSA[i-1] = MatchSA(indicesSA[i], phrases[indicesSA[i]].pos, std::upper_bound(std::begin(headBoundaries), std::end(headBoundaries), indicesSA[i]) - std::begin(headBoundaries), phrases[indicesSA[i]].smaller);
      //if(verbose) std::cerr << headsSA[i-1].start << "," << _ISA[headsSA[i-1].pos] << "," << headsSA[i-1].len << "," << headsSA[i-1].smaller << "," << headsSA[i-1].next <<"\n";
   }
   std::cerr << "headBoundaries.size(): " << headBoundaries.size() << "\n";
   std::vector<uint32_t>().swap(headBoundaries);
   std::cerr << "Permuted headsSA\n";
   delete[] indicesSA;
   t02 = std::chrono::high_resolution_clock::now();
   std::cerr << "Permuting MS-heads took: " << std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01).count() << " ms\n";

   t2 = std::chrono::high_resolution_clock::now();
   uint64_t headSortTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Time to sort heads: " << headSortTime << " milliseconds\n";

   //if(verbose){
   //   uint64_t errHeads = checkHeadsSA(headsSA, phrases.size());
   //   std::cerr << "N. errors on headsSA: " << errHeads << "\n";
   //}
   //exit(0);
   t1 = std::chrono::high_resolution_clock::now();
   std::cerr << "Computing nextPos in headSA\n";
   std::vector<Match>::iterator begPhrases = phrases.begin();
   for(size_t i = 0; i < headsSA.size(); i++){
      //_ISA[headsSA[headA.pos].pos + (nextStartA - headA.start)]
      uint32_t nextStart = phrases[headsSA[i].start].start + phrases[headsSA[i].start].len;
      Match headNextStart = Match(nextStart, 0, headsSA[i].len);
      Match nextHead = *(pHeads.predQuery(headNextStart, begPhrases));
      headsSA[i].pos = _ISA[nextHead.pos + (nextStart - nextHead.start)];
      __builtin_prefetch(&phrases[headsSA[i+1].start]);
   }
   std::cerr << "Computing nextPos in headSA DONE\n";
   std::cerr << "Computing changeP in phrases\n";
   for(size_t i = 0; i < headsSA.size(); i++){
      phrases[headsSA[i].start].changeP(i);
   }
   std::cerr << "Computing changeP in phrases DONE\n";
   std::cerr << "Computing change heads in sStar\n";
   for(size_t i = 0; i < nStar; i++){
      sStar[i].diffLen.val = sStar[i].idx - phrases[sStar[i].head].start;
      //if( sStar[i].diffLen.val) std::cerr << "diffLen: " << sStar[i].diffLen.val << " idx: " << sStar[i].idx << " start: " << phrases[sStar[i].head].start <<"\n";
      sStar[i].head = phrases[sStar[i].head].pos;
   }
   std::cerr << "Computing change heads in sStar DONE\n";
   std::cerr << "Computing change len and start in headSA\n";
   for(size_t i = 0; i < headsSA.size(); i++){
      uint32_t nextStart = phrases[headsSA[i].start].start + phrases[headsSA[i].start].len;
      Match headNextStart = Match(nextStart, 0, headsSA[i].len);
      Match nextHead = *(pHeads.predQuery(headNextStart, begPhrases));
      headsSA[i].len = phrases[headsSA[i].start].len;
      headsSA[i].start = nextHead.pos;
      __builtin_prefetch(&phrases[headsSA[i+1].start]);
   }
   std::vector<Match>().swap(phrases);
   // Now: 
   // headsSA[i] = (nextHeadRank, nextHeadISA[pos], length, smaller, mmchar);

   std::cerr << "Computing change len and start in headSA DONE\n";
   std::cerr << "Starting to sort S*-suffixes\n";
   
   //#pragma omp parallel for schedule(dynamic)
   for(uint32_t i = 1; i < _n + 1; i++){
      std::sort(begStar + prefSumBucketLengthsStar[i-1], begStar + prefSumBucketLengthsStar[i], compareSuf);
   }
   std::vector<MatchSA>().swap(headsSA);
   delete[] prefSumBucketLengthsStar;

   t2 = std::chrono::high_resolution_clock::now();
   uint64_t sortTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "S*-suffixes sorted in: " << sortTime << " milliseconds\n";
   //std::cerr << "number of diffPos " << diffSufPos << "\n";
   //std::cerr << "number of diffLen " << diffSufLen << "\n";
   //std::cerr << "number of diffHeads " << diffSufHeads << "\n";

   if(verbose){
      uint64_t errSStar = checkGSA(sStar, nStar, _sx);
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
      //t_sum += bucketLengths[i-1];
      prefSumCharBkts[i] = prefSumCharBkts[i-1] + charBkts[i-1];
      //if(verbose) std::cerr << prefSumBucketLengths[i] << "\n";
   }
   prefSumCharBkts[sizeChars] = _sn;
   for(size_t i = 0; i < sizeChars; i++){
      prefSumCharBktsEnds[i] = prefSumCharBkts[i+1];
   }

   std::cerr << "Going to convert types\n";
   // //change back original positions in heads
   // for(size_t i = 0; i < headsSA.size(); i++){
   //    phrases[i].changeP(headsSA[phrases[i].pos].pos);
   // }
   MSGSA.reserve(_sn);
   uint64_t diffLen;
   std::reverse(sStar.begin(), sStar.end());
   
   for(uint32_t x = 1; x < sizeChars; x++){
      //diffLen = (prefSumCharBkts[x] - prefSumCharBkts[x-1]) - (prefSumBucketLengthsStar[x] - prefSumBucketLengthsStar[x-1]);
      diffLen = (prefSumCharBkts[x] - prefSumCharBkts[x-1]) - bucketLengthStarChar[x-1];
      for(uint64_t p = 0; p < diffLen; p++){
         //MSGSA.push_back(Suf());
         MSGSA.push_back(std::make_pair(0,0));
      }
      //for(uint32_t p = 0; p < (prefSumBucketLengthsStar[x] - prefSumBucketLengthsStar[x-1]); p++){
      for(uint32_t p = 0; p < bucketLengthStarChar[x-1]; p++){
         //MSGSA.push_back(Suf(sStar.back()));
         MSGSA.push_back(std::make_pair(sStar.back().idx, sStar.back().doc));
         sStar.pop_back();
      }
      if(bucketLengthStarChar[x-1]) sStar.shrink_to_fit();
   }

   std::vector<SufSStar>().swap(sStar);

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
      //__builtin_prefetch(&_sx[MSGSA[i+1].idx + docBoundaries[MSGSA[i+1].doc - 1] - 1]);
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
      //__builtin_prefetch(&_sx[MSGSA[i-1].idx + docBoundaries[MSGSA[i-1].doc - 1] - 1]);
   }
   std::cerr << "\nAfter inducing S-types\n";
   t2 = std::chrono::high_resolution_clock::now();
   uint64_t induceTime = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
   std::cerr << "Induced in: " << induceTime << " milliseconds\n";

   int checkGSA = 0;
   if(checkGSA){
      std::cerr << "Checking GSA\n"; 
      uint32_t err = 0;
      #pragma omp parallel for
      for(size_t i = 0; i < _sn; i++){
         if(verbose) std::cerr << "i=" << i << ": " << MSGSA[i].first << " " << MSGSA[i].second << " " << "\n";//MSGSA[i].head <<"\n";
         
         //data_type *_slice_sx = &_sx[MSGSA[i].first + docBoundaries[MSGSA[i].second - 1]];
         //data_type *_slice_prev;
         uint32_t maxIdx;
         if(i > 0){
            if(MSGSA[i].second == 0 || MSGSA[i-1].second == 0){
               std::cerr << "There was an empty entry\n";
               err++;
               continue;
            }
            //_slice_prev = &_sx[MSGSA[i-1].first + docBoundaries[MSGSA[i-1].second - 1]];
            maxIdx = std::min(docBoundaries[MSGSA[i].second] - (MSGSA[i].first + docBoundaries[MSGSA[i].second - 1]), docBoundaries[MSGSA[i-1].second] - (MSGSA[i-1].first + docBoundaries[MSGSA[i-1].second - 1]));
         } 
         else{
            //_slice_prev = (data_type *)"$";
            maxIdx = 1;
            continue;
         }
         //if(verbose) std::cerr << "suf_i-1 " << _slice_prev;
         //if(verbose) std::cerr << "suf_i " << _slice_sx << "\n";
         
         if(memcmp(&_sx[MSGSA[i].first + docBoundaries[MSGSA[i].second - 1]], &_sx[MSGSA[i-1].first + docBoundaries[MSGSA[i-1].second - 1]], maxIdx) < 0){
            if(verbose) std::cerr << "PROBLEM with " << i-1 << " (" << MSGSA[i-1].first << "," << MSGSA[i-1].second << ") and " << i << " (" << MSGSA[i].first << "," << MSGSA[i].second << ")\n"; 
            err++;
            //if(err) break;
         }
      }
      std::cerr << "n. errors " << err << "\n"; 
   } 
   // std::cerr << "maxCounter " << maxCounter << "\n";
   // std::cerr << "meanCounter " << sumCounter/(denCounter+0.00000001) << "\n";
   // std::cerr << "times it had to compare only one char " << diffLenCounter << "\n";
   // std::cerr << "times it had to compare more than one char " << denCounter << "\n";
   // std::cerr << "times it compared two suff in const time " << finalSuffCounter << "\n";
   // delete[] _LCP;
   // delete[] _ISA;
   // delete[] _SA;
   //std::vector<std::pair<uint32_t, int32_t>>().swap(MSGSA);

   return numfactors;
}

void computeGSA(char* refFileName, char* collFileName, std::vector<std::pair<uint32_t, int32_t>> &MSGSA){
   const std::string _sx = lzInitialize(refFileName, collFileName);
   lzFactorize(_sx, MSGSA);
   //return MSGSA;
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
            renormalizations++;
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


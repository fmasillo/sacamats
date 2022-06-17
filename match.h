#ifndef __MATCH_H
#define __MATCH_H
#include <stdint.h>

union Key
{
   uint32_t val;
   char bytearray[sizeof(uint32_t)];
   uint16_t doublebytearray[sizeof(uint32_t)];
};

struct Match{
   Match() : start(0),pos(0),len(0),smaller(false) { }
   //Match(uint32_t p, uint32_t l, unsigned char nxt){
   //   pos = p; len = l, next = nxt;
   //}
   Match(uint32_t s, uint32_t p, uint32_t l) : start(s), pos(p), len(l), smaller(false) { }
   Match(uint32_t s, uint32_t p, uint32_t l, bool sm) : start(s), pos(p), len(l), smaller(sm) { }
   void changeS(uint32_t s){
      start = s;
   }
   void changeP(uint32_t p){
      pos = p;
   }
   void changeD(uint32_t d){
      len = d;
   }

   //uint64_t suf; //position of suffix in collection
   uint32_t start; //position in the collection 
   uint32_t pos; //position of match in _reference
   uint32_t len; //length of match
   //data_type next; //symbol in the collection after the match
   bool smaller;
};

struct MatchSA
{
   MatchSA() : start(0),pos(0),len(0),smaller(false),next(0) { }
   //Match(uint32_t p, uint32_t l, unsigned char nxt){
   //   pos = p; len = l, next = nxt;
   //}
   MatchSA(uint32_t s, uint32_t p, uint32_t l, bool sm) : start(s), pos(p), len(l), smaller(sm), next(0) { }
   MatchSA(uint32_t s, uint32_t p, uint32_t l, bool sm, unsigned char n) : start(s), pos(p), len(l), smaller(sm), next(n) { }
   uint32_t start; //position in the collection 
   uint32_t pos; //position of match in _reference
   uint32_t len; //length of match
   bool smaller;
   unsigned char next; //symbol in the collection after the match
};


struct SufSStar{
   SufSStar() : head(0),  idx(0), doc(0) { diffLen.val = 0; }
   SufSStar(uint32_t i, uint32_t d) : head(0), idx(i), doc(d) { diffLen.val = 0; }
   SufSStar(uint32_t i, uint32_t d, uint32_t h) : head(h), idx(i), doc(d) { diffLen.val = 0; }
   SufSStar(uint32_t i, uint32_t d, uint32_t h, uint32_t dL) : head(h), idx(i), doc(d) { diffLen.val = dL; }

   bool operator<(const SufSStar &other){
      return diffLen.val < other.diffLen.val;
   }

   void assign(uint32_t i, uint32_t d, uint32_t h, uint32_t dL){ 
      head = h, diffLen.val = dL, idx = i, doc = d;
   }
   uint32_t head;
   Key diffLen;
   uint32_t idx;
   uint32_t doc;
   
};

struct Suf{
   Suf() : idx(0), doc(0) { }
   Suf(uint32_t i, uint32_t d) : idx(i), doc(d) { }
   Suf(SufSStar &s) : idx(s.idx), doc(s.doc) { }

   void changeDocSign(){
      doc = ~doc;
   }

   uint32_t idx;
   int32_t doc;
};

struct PsvNsv{
   PsvNsv() : psv(-1), nsv(-1) { }
   PsvNsv(uint32_t p, uint32_t n) : psv(p), nsv(n) { }
   uint32_t psv, nsv;
};

// inline std::vector<Match>::iterator findHead(const std::vector<Match>::iterator start, const uint32_t len, const uint32_t m){
//    std::vector<Match>::iterator base = start;
//    uint32_t l = len, half;
//    while(l > 1){
//       half = l / 2;
//       //__builtin_prefetch(&base[(len - half) / 2]);
//       //__builtin_prefetch(&base[half + (len - half) / 2]);
//       //base += half*(base[half].start <= m); 
//       base = ((base + half)->start > m ? base : base + half);
//       l -= half;
//    }
//    return base + ((base)->start <= m);

//    // int32_t l = 0, r = len - 1;
//    // while (l < r) {
//    //    int32_t mid = (l + r) / 2;
//    //    //__builtin_prefetch(&start[(r-mid)/2]);
//    //    //__builtin_prefetch(&start[(mid + (r - mid))/2]);
//    //    if (start[mid].start >= m)
//    //       r = mid;
//    //    else
//    //       l = mid + 1;
//    // }
//    // return start + l + (start[l].start <= m);
   
//    // size_t m_len = len;
//    // size_t low = -1;
//    // size_t high = m_len;
//    // //assert(m_len < std::numeric_limits<size_t>::max() / 2);
//    // while(high - low > 1){
//    //    uint32_t probe = (low + high) / 2;
//    //    if(start[probe].start > m)
//    //       high = probe;
//    //    else
//    //       low = probe;
//    // }
//    // return start + high;
//    // if (high == m_len)
//    //    return start + m_len;
//    // else
//    //    return start + high;

//    // uint32_t index = 0;
//    // uint32_t size = len;
//    // while (size > 1) {
//    //       size /= 2;
//    //       uint32_t probe = index + size;
//    //       if (start[probe].start <= m)
//    //          index = probe + 1;
//    // }
//    // return start + index + ((start)->start <= m);
// }



#endif
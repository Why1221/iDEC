#ifndef _IDEC_UTILS_H_
#define _IDEC_UTILS_H_
#include <cstdint>

// counting 1's inside `x`
inline int16_t countingOnes(uint64_t x) {
#if defined(__GNUC__) || defined(__GNUG__)
  return __builtin_popcountll(x);
#elif
  int16_t setBits = 0;

  while (x > 0) {
    setBits += x & 1;
    x >>= 1;
  }

  return setBits;
#endif
}

// calculate hamming distance of two uint64_t integers
inline int16_t hammingDistance(uint64_t n1,  // number one
                               uint64_t n2   // number two
) {
  return countingOnes(n1 ^ n2);
}

// calculate hamming distance of two uint64_t vectors with a length of `enc_dim`
template <typename ResType = int16_t>
inline ResType hammingDistance(const uint64_t *n1,  // vector one
                               const uint64_t *n2,  // vector two
                               unsigned enc_dim,    // vector length
                               ResType init = 0     // initial distance
) {
  for (unsigned i = 0; i < enc_dim; ++i) init += hammingDistance(n1[i], n2[i]);
  return init;
}

// calculate ToW sum
// argument order matters, the first must be the point (or coordinates of a
// points)
inline int16_t calcToWSum(uint64_t a,       // point (or partial of a point)
                          uint64_t b,       // random vector (or partial of a random vector)
                          int16_t init = 0  // initial value
) {
  return init + countingOnes(a & b) - countingOnes(a & (~b));
}

#endif  // _IDEC_UTILS_H_
#include "linear_scan_l1.h"
#include <cstdlib>
#include <algorithm>
#include <cassert>
#include <limits>

static float MAXFLOAT = std::numeric_limits<float>::max();

// // Hamming distance for two uint64 numbers
// inline int _hammingU64(uint64_t n1,  // the first number
//                        uint64_t n2   // the second
// ) {
// #if defined(__GNUC__) || defined(__GNUG__)
//   // use popcnt to speedup
//   return __builtin_popcountll(n1 ^ n2);
// #else
//   uint64_t x = n1 ^ n2;
//   int setBits = 0;
//   while (x > 0) {
//     setBits += x & 1;
//     x >>= 1;
//   }
//   return setBits;
// #endif
// }

// L1 distance for two high-dimension points 
// format.
inline float _L1Distance(const float *_p1,  // point one
                              const float *_p2,  // point two
                              int _dim          // original dimension 
) {
  float dist = 0.0f;
  for (int i = 0; i < _dim; ++i) {
    dist += std::abs(_p1[i]- _p2[i]);
  }
  return dist;
}

// Linear scan for 1nn
static std::pair<size_t, float> linear_scan_1nn(
    const float *dataset,       // ponits in database
    size_t n,                      // database size (# of points)
    const float *query,         // query point
    unsigned _dim,              // dimension 
    float current_best = MAXFLOAT  // current best result
) {
  auto distance = current_best;
  size_t index = -1;
  for (size_t i = 0; i < n; ++i) {
    float score = _L1Distance(dataset + i * _dim, query, _dim);
    if (score < distance) {
      index = i;
      distance = score;
    }
  }

  return {index, distance};
}

namespace {
using Pair = std::pair<size_t, float>;
class PairCmp {
 public:
  bool operator()(const Pair &a, const Pair &b) const {
    return a.second < b.second;
  }
};
}  // namespace

// Linear Scan for knn
void linear_scan(
    const float *dataset,  // points in database
    size_t n,                 // database size (# of points)
    const float *query,    // query point
    unsigned _dim,         // dimension
    unsigned k,               // # of nearest neighbors
    std::vector<Pair> &ans,   // resutls
    bool ct  // continue to search (not to clear existing results)
) {
  assert(k < n);
  if (k == 1) {
    float current_best = MAXFLOAT;
    if (ans.empty()) {
      ans.resize(1);
    } else {
      if (ct) current_best = ans.front().second;
    }

    ans[0] = linear_scan_1nn(dataset, n, query, _dim, current_best);
    return;
  }

  if (!ct) {
    ans.clear();
    ans.reserve(k);
  }

  for (size_t i = 0; i < n; ++i) {
    float score = _L1Distance(dataset + i * _dim, query, _dim);

    if (ans.size() < k - 1) {
      ans.emplace_back(i, score);
    } else if (ans.size() == k - 1) {
      ans.emplace_back(i, score);
      std::make_heap(ans.begin(), ans.end(), PairCmp{});
    } else {
      if (score < ans.front().second) {
        std::pop_heap(ans.begin(), ans.end(), PairCmp{});
        ans.pop_back();
        ans.emplace_back(i, score);
        std::push_heap(ans.begin(), ans.end(), PairCmp{});
      }
    }
  }

  if (!ct) std::sort_heap(ans.begin(), ans.end(), PairCmp{});
}

void linear_scan(
    const std::vector<float> &dataset,  // points in database
    const float *query,                 // query point
    unsigned _dim,         // dimension
    unsigned k,                            // # of nearest neighbors
    std::vector<Pair> &ans,                // resutls
    bool ct  // continue to search (not to clear existing results)
) {
  linear_scan(&dataset[0], dataset.size()/_dim, query, _dim, k,
                      ans, ct);
}

void linear_scan(
    const std::vector<float> &dataset,  // points in database
    std::vector<float> &query,          // query point
    unsigned k,                            // # of nearest neighbors
    std::vector<Pair> &ans,                // resutls
    bool ct  // continue to search (not to clear existing results)
) {
  unsigned dim = query.size();
  auto n = dataset.size() / dim;
  linear_scan(&dataset[0], n, &query[0], dim, k, ans, ct);
}
#include "linear_scan_compact.h"
#include <algorithm>
#include <cassert>
#include <limits>

static float MAXFLOAT = std::numeric_limits<float>::max();

// Hamming distance for two uint64 numbers
inline int _hammingU64(uint64_t n1,  // the first number
                       uint64_t n2   // the second
) {
#if defined(__GNUC__) || defined(__GNUG__)
  // use popcnt to speedup
  return __builtin_popcountll(n1 ^ n2);
#else
  uint64_t x = n1 ^ n2;
  int setBits = 0;
  while (x > 0) {
    setBits += x & 1;
    x >>= 1;
  }
  return setBits;
#endif
}

// Hamming distance for two high-dimension points encoded in cpmpact
// format.
inline float _hammingDistance(const uint64_t *_p1,  // point one
                              const uint64_t *_p2,  // point two
                              int _enc_dim          // dimension after encoding
) {
  float dist = 0.0f;
  for (int i = 0; i < _enc_dim; ++i) {
    dist += _hammingU64(_p1[i], _p2[i]);
  }
  return dist;
}

// Linear scan for 1nn
static std::pair<size_t, float> linear_scan_compact_1nn(
    const uint64_t *dataset,       // ponits in database
    size_t n,                      // database size (# of points)
    const uint64_t *query,         // query point
    unsigned enc_dim,              // dimension after encoding
    float current_best = MAXFLOAT  // current best result
) {
  auto distance = current_best;
  size_t index = -1;
  for (size_t i = 0; i < n; ++i) {
    float score = _hammingDistance(dataset + i * enc_dim, query, enc_dim);
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
void linear_scan_compact(
    const uint64_t *dataset,  // points in database
    size_t n,                 // database size (# of points)
    const uint64_t *query,    // query point
    unsigned enc_dim,         // dimension after compact encoding
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

    ans[0] = linear_scan_compact_1nn(dataset, n, query, enc_dim, current_best);
    return;
  }

  if (!ct) {
    ans.clear();
    ans.reserve(k);
  }

  for (size_t i = 0; i < n; ++i) {
    float score = _hammingDistance(dataset + i * enc_dim, query, enc_dim);

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

void linear_scan_compact(
    const std::vector<uint64_t> &dataset,  // points in database
    const uint64_t *query,                 // query point
    unsigned enc_dim,                      // dimension after compact encoding
    unsigned k,                            // # of nearest neighbors
    std::vector<Pair> &ans,                // resutls
    bool ct  // continue to search (not to clear existing results)
) {
  linear_scan_compact(&dataset[0], dataset.size() / enc_dim, query, enc_dim, k,
                      ans, ct);
}

void linear_scan_compact(
    const std::vector<uint64_t> &dataset,  // points in database
    std::vector<uint64_t> &query,          // query point
    unsigned k,                            // # of nearest neighbors
    std::vector<Pair> &ans,                // resutls
    bool ct  // continue to search (not to clear existing results)
) {
  unsigned enc_dim = query.size();
  auto n = dataset.size() / enc_dim;
  linear_scan_compact(&dataset[0], n, &query[0], enc_dim, k, ans, ct);
}
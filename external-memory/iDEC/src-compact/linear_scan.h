#include <cassert>
#include <limits>
#include <queue>
#include <vector>

namespace ls {
using Pair = std::pair<size_t, float>;
class PairCmp {
 public:
  bool operator()(const Pair& a, const Pair& b) const {
    return a.second < b.second;
  }
};
}  // namespace ls

//using num_t = int16_t;

template <typename num_t>
std::pair<size_t, float> linear_scan_1nn(const num_t* dataset, size_t n, const num_t* query,
                                         unsigned dim) {
  auto distance = std::numeric_limits<float>::max();
  size_t index = -1;
  for (size_t i = 0; i < n; ++i) {
    float score = 0.0f;
    for (unsigned j = 0; j < dim; ++j) {
      float tmp = dataset[i * dim + j] - query[j];
      score += tmp * tmp;
    }
    if (score < distance) {
      index = i;
      distance = score;
    }
  }

  return {index, distance};
}

template <typename num_t>
void linear_scan(const num_t* dataset, size_t n, const num_t* query,
                 unsigned dim, unsigned k,
                 std::vector<std::pair<size_t, float>>& ans) {
  assert(k < n);
  if (k == 1) {
    ans.resize(1);
    ans[0] = linear_scan_1nn(dataset, n, query, dim);
    return;
  }

  ans.reserve(k);

  for (size_t i = 0; i < n; ++i) {
    float score = 0.0f;
    for (unsigned j = 0; j < dim; ++j) {
      float tmp = dataset[i * dim + j] - query[j];
      score += tmp * tmp;
    }

    if (ans.size() < k - 1) {
      ans.emplace_back(i, score);
    } else if (ans.size() == k - 1) {
      ans.emplace_back(i, score);
      std::make_heap(ans.begin(), ans.end(), ls::PairCmp{});
    } else {
      if (score < ans.front().second) {
        std::pop_heap(ans.begin(), ans.end(), ls::PairCmp{});
        ans.pop_back();
        ans.emplace_back(i, score);
        std::push_heap(ans.begin(), ans.end(), ls::PairCmp{});
      }
    }
  }

  std::sort_heap(ans.begin(), ans.end(), ls::PairCmp{});
}

template <typename num_t>
void linear_scan(const std::vector<num_t>& dataset,
                 const std::vector<num_t> query, unsigned k,
                 std::vector<std::pair<size_t, float>>& ans) {
  assert(dataset.size() % query.size() == 0);
  linear_scan(&dataset[0], dataset.size() / query.size(), &query[0],
              query.size(), k, ans);
}


void linear_scan_compact(const uint64_t* dataset, size_t n, const uint64_t* query,
                 unsigned dim, unsigned k,
                 std::vector<std::pair<size_t, float>>& ans) {
  assert(k < n);
  if (k == 1) {
    ans.resize(1);
    ans[0] = linear_scan_1nn(dataset, n, query, dim);
    return;
  }

  ans.reserve(k);

  for (size_t i = 0; i < n; ++i) {
    float score = 0.0f;
    for (unsigned j = 0; j < dim; ++j) {
      float tmp = dataset[i * dim + j] - query[j];
      score += tmp * tmp;
    }

    if (ans.size() < k - 1) {
      ans.emplace_back(i, score);
    } else if (ans.size() == k - 1) {
      ans.emplace_back(i, score);
      std::make_heap(ans.begin(), ans.end(), ls::PairCmp{});
    } else {
      if (score < ans.front().second) {
        std::pop_heap(ans.begin(), ans.end(), ls::PairCmp{});
        ans.pop_back();
        ans.emplace_back(i, score);
        std::push_heap(ans.begin(), ans.end(), ls::PairCmp{});
      }
    }
  }

  std::sort_heap(ans.begin(), ans.end(), ls::PairCmp{});
}

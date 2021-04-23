#include <chrono>
#include <iostream>
#include <random>
#include "linear_scan.h"
#include "ptR-tree.h"

typedef ptRTree *RtreePtr;

const int SAMPLES_DIM = 10;

inline void generateRandomPoint(int16_t *point, size_t dim,
                                const int16_t min_range = -128,
                                const int16_t max_range = 128) {
  std::default_random_engine eng{
      (unsigned)std::chrono::system_clock::now().time_since_epoch().count()};

  std::uniform_int_distribution<> distr(min_range, max_range);

  for (size_t d = 0; d < dim; d++) point[d] = distr(eng);
}

void generateRandomPointCloud(std::vector<int16_t> &samples, const size_t N,
                              const size_t dim, const int16_t min_range = -128,
                              const int16_t max_range = 128) {
  std::cout << "Generating " << N << " random points...";
  samples.resize(N * dim);
  for (size_t i = 0; i < N; i++)
    generateRandomPoint(&samples[i * dim], dim, min_range, max_range);
  std::cout << "done\n";
}

void rtree_demo(const size_t nSamples, const int dim) {
  std::vector<int16_t> samples;
  // Generate points:
  generateRandomPointCloud(samples, nSamples, dim);

  // Query point:
  std::vector<int16_t> query_pt(dim);
  generateRandomPoint(&query_pt[0], dim);

  RtreePtr rt = new ptRTree();

  // Build index
  std::cout << "Building r-tree index ...";
  const int B = 128;
  rt->init((char *)"./rt-example.rt" /* file to store the rtree */,
           2 * B /* block size */, NULL /* cache */, dim /* dimension */);

  {
    for (size_t j = 0; j < nSamples; ++j) {
      Entry *e = NULL;

      e = new ptEntry();
      e->init(rt);
      e->level = 0;
      e->son = j; /* id */

      for (int i = 0; i < dim; ++i) {
        e->bounces[2 * i] = e->bounces[2 * i + 1] = samples[j * dim + i];
      }
      rt->insert(e);
    }
  }
  std::cout << "done\n";

  // do a knn search
  const size_t num_results = 3;

  std::vector<std::pair<size_t, float_t> > indices_dists;
  rt->kNN(&query_pt[0], num_results, indices_dists);

  std::cout << "knnSearch(nn=" << num_results << "): \n";

  for (size_t i = 0; i < num_results; i++)
    std::cout << "ret_index[" << i << "]=" << indices_dists[i].first
              << " out_dist_sqr=" << indices_dists[i].second << std::endl;

  // ground truth
  std::vector<std::pair<size_t, float_t> > gnd;
  linear_scan<int16_t>(samples, query_pt, num_results, gnd);

  std::cout << "\n\nknnSearch(nn=" << num_results << ") ground truth: \n";
  for (size_t i = 0; i < num_results; i++)
    std::cout << "ret_index[" << i << "]=" << gnd[i].first
              << " out_dist_sqr=" << gnd[i].second << std::endl;
}

int main() { rtree_demo(10000 /* samples */, SAMPLES_DIM /* dim */); }

#define USE_SEARCH_HISTORY 1
#include "KDTreeFlatVectorAdaptor.h"

#include <iostream>
#include <random>
#include <chrono>
#include <memory>

const int SAMPLES_DIM = 15;
unsigned G_SEED = 0;
std::unique_ptr<std::default_random_engine> G_GEN;

inline void generateRandomPoint(float *point, size_t dim, const float max_range = 10.0){

    std::default_random_engine eng{(unsigned)std::chrono::system_clock::now().time_since_epoch().count() };
    std::uniform_real_distribution<float> distr(0, max_range);

    for (size_t d = 0; d < dim; d++)
        point[d] = distr(G_SEED == 0 ? eng : *G_GEN);
}

inline void generateRandomPoint(std::vector<float> &point, size_t dim, const float max_range = 10.0){
    if (point.size() != dim) point.resize(dim);
    generateRandomPoint(&point[0], dim, max_range);
}

void generateRandomPointCloud(std::vector<float> &samples,
        const size_t N, const size_t dim,
        const float max_range = 10.0)
{
    std::cout << "Generating "<< N << " random points...";
    samples.resize(N * dim);
    for (size_t i = 0; i < N; i++)  generateRandomPoint(&samples[i * dim], dim, max_range);
    std::cout << "done\n";
}

void kdtree_demo(const size_t nSamples, const size_t dim)
{
    std::vector<float>  samples;

    const double max_range = 20;

    // Generate points:
    generateRandomPointCloud(samples, nSamples,dim, max_range);

    // Query point:
    std::vector<float> query_pt(dim);
    generateRandomPoint(query_pt, dim, max_range);

    // construct a kd-tree index:
    // Dimensionality set at run-time (default: L2)
    // ------------------------------------------------------------
    typedef KDTreeFlatVectorAdaptor< std::vector<float>, float >  my_kd_tree_t;

    my_kd_tree_t   mat_index(dim /*dim*/, samples, 10 /* max leaf */ );
    mat_index.index->buildIndex();

    // do a knn search
    const size_t num_results = 3;
    std::vector<size_t>   ret_indexes(num_results);
    std::vector<float> out_dists_sqr(num_results);

    nanoflann::KNNResultSet<float> resultSet(num_results);

    resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
    mat_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

    std::cout << "knnSearch(nn="<<num_results<<"): \n";
    for (size_t i = 0; i < num_results; i++)
        std::cout << "ret_index["<<i<<"]=" << ret_indexes[i] << " out_dist_sqr=" << out_dists_sqr[i] << std::endl;

    auto hist = mat_index.index->getHistory();

    for (size_t i = 0;i < hist.size();++ i) {
        std::cout << "index = " << hist[i].second << ", dist = " << hist[i].first << "\n"; 
    }
}

int main(int argc, char** argv)
{
    if (argc == 2) {
        G_SEED = atoi(argv[1]);
        G_GEN = std::unique_ptr<std::default_random_engine>(new std::default_random_engine(G_SEED));
    }
    kdtree_demo(1000 /* samples */, SAMPLES_DIM /* dim */);
}
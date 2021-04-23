#include <cstdint>
#include <utility>
#include <vector>


void linear_scan(
    const float *dataset,  // points in database
    size_t n,                 // database size (# of points)
    const float *query,    // query point
    unsigned _dim,         // dimension
    unsigned k,               // # of nearest neighbors
    std::vector<std::pair<size_t, float>> &ans,  // resutls
    bool ct = false  // continue to search (not to clear existing results)
);

void linear_scan(
    const std::vector<float> &dataset,  // points in database
    const float *query,                 // query point
    unsigned _dim,                      // dimension 
    unsigned k,                            // # of nearest neighbors
    std::vector<std::pair<size_t, float>> &ans,  // resutls
    bool ct = false  // continue to search (not to clear existing results)
);

void linear_scan(
    const std::vector<float> &dataset,        // points in database
    std::vector<float> &query,                // query point
    unsigned k,                                  // # of nearest neighbors
    std::vector<std::pair<size_t, float>> &ans,  // resutls
    bool ct = false  // continue to search (not to clear existing results)
);

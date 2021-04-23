#include <cstdint>
#include <utility>
#include <vector>


void linear_scan_compact(
    const uint64_t *dataset,  // points in database
    size_t n,                 // database size (# of points)
    const uint64_t *query,    // query point
    unsigned enc_dim,         // dimension after compact encoding
    unsigned k,               // # of nearest neighbors
    std::vector<std::pair<size_t, float>> &ans,  // resutls
    bool ct = false  // continue to search (not to clear existing results)
);

void linear_scan_compact(
    const std::vector<uint64_t> &dataset,  // points in database
    const uint64_t *query,                 // query point
    unsigned enc_dim,                      // dimension after compact encoding
    unsigned k,                            // # of nearest neighbors
    std::vector<std::pair<size_t, float>> &ans,  // resutls
    bool ct = false  // continue to search (not to clear existing results)
);

void linear_scan_compact(
    const std::vector<uint64_t> &dataset,        // points in database
    std::vector<uint64_t> &query,                // query point
    unsigned k,                                  // # of nearest neighbors
    std::vector<std::pair<size_t, float>> &ans,  // resutls
    bool ct = false  // continue to search (not to clear existing results)
);


#ifndef _IDEC_H_
#define _IDEC_H_

#include <cstdint>  // for uint64_t
#include <vector>

#include "KDTreeFlatVectorAdaptor.h"

using num_t = int16_t;
// using iDECKDTree = KDTreeFlatVectorAdaptor<std::vector<num_t>, int16_t>;
using iDECKDTree = KDTreeFlatVectorAdaptor<std::vector<int16_t>, int16_t, -1, float>;

class iDEC {
 public:
  // restore/construct an iDEC object from previous build
  // throws on error
  // to make this constructor work, previous run must call saveIndex with
  // argument `indexFolder`
  explicit iDEC(
      const std::vector<uint64_t>& raw_data,  // raw data in compact flat vector
      const std::string& indexFolder          // folder for index
  );
  // construct an iDEC object fron scratch
  // NOTE: this constructor only initilizes paramters (i.e., it does not build
  // the index)
  //       to build the index, please call build().
  explicit iDEC(const std::vector<uint64_t>& raw_data,  // raw data
                unsigned d,                             // original dimension
                unsigned m,        // dimension of the projection
                unsigned nt = 1,   // # of trees
                unsigned seed = 0  // random seed (0 -- using default seed)
  );

  // get hash value for a point
  inline int16_t getHashValue(const uint64_t* point,  // point
                              unsigned mid,           // project coordinate
                              unsigned tid = 0u       // tree id
                              ) const;
  // build the index
  void build();
  // save index to the folder `indexFolder`
  void saveIndex(const std::string& indexFolder) const;
  // perform knn search for query point `query`
  void knn(const uint64_t* query,  // query point
           unsigned k,             // # of nearest neighbors
           float t,  // maximum number (ratio) of points to be returned fron the
                     // projection space (MUST be from 0 to 1.0)
           std::vector<std::pair<float /* distance */, size_t /* index */>>&
               ans  // knn answers
           ) const;
  // index size (in bytes)
  size_t size() const { return _tot_mem; }

 private:
  // one-nearest neighbor query
  void _one_nn(const uint64_t* query, float t, float& dist,
               size_t& index) const;
  // generate random walk vectors
  // throws on error
  void genRandomWalkVectors();
  // load parameters and index from the folder `indexFolder`
  void loadIndex(const std::string& indexFolder);
  // sace all parameters to file
  void saveParam(const std::string& indexFolder) const;
  // load parameters from file
  void loadParam(const std::string& indexFolder);
  // get filename for tree with id `tid`
  inline std::string _getTreeFile(const std::string& indexFolder,
                                  unsigned tid) const;
  inline std::string _getProjDataFile(const std::string &indexFolder) const;
  // get filename for parameters
  inline std::string _getParamFile(const std::string& indexFolder) const;

  unsigned _seed;     // random seed
  unsigned _d;        // dimensionality
  unsigned _enc_dim;  // encoded dimenstion
  unsigned _n;        // # of points
  unsigned _m;        // projected dimension
  unsigned _nt;       // # of trees

  size_t _tot_mem;  // total memory

  std::vector<uint64_t> _a_array;              // random variables
  std::vector<std::vector<num_t>> _proj_data;  // data points in project space
  const std::vector<uint64_t>& _raw_data;      // reference to raw_data
  std::vector<iDECKDTree> _trees;
};
#endif  // _IDEC_H_

#include "iDEC.h"

#include <algorithm>
#include <bitset>
#include <fstream>
#include <random>

#include "IdecException.h"
#include "IdecUtils.h"

/* Local Constants */
static const size_t ONE_KB = 1024;                // 1 KB
static const size_t ONE_MB = ONE_KB * 1024;       // 1 MB
static const size_t ONE_GB = ONE_MB * 1024;       // 1 GB
static const size_t MAX_MEM = ONE_GB * 100;       // 100 GB
static const unsigned MAX_LEAF = 10;              // Maximum number of leaf
static const unsigned DEFAULT_SEED = 111650782u;  // default random seed

// The following was produced by script/chi2inv_calc.m
// m,delta,x
// 4,0.90,7.779440
// 6,0.90,10.644641
// 8,0.90,13.361566
// 10,0.90,15.987179
// 12,0.90,18.549348
// 14,0.90,21.064144
// 16,0.90,23.541829
// m,delta,x
// 4,0.95,9.487729
// 6,0.95,12.591587
// 8,0.95,15.507313
// 10,0.95,18.307038
// 12,0.95,21.026070
// 14,0.95,23.684791
// 16,0.95,26.296228
// m,delta,x
// 4,0.99,13.276704
// 6,0.99,16.811894
// 8,0.99,20.090235
// 10,0.99,23.209251
// 12,0.99,26.216967
// 14,0.99,29.141238
// 16,0.99,31.999927

static const float _CHI2INV_ = 10.644641f;

iDEC::iDEC(const std::vector<uint64_t> &raw_data,
           const std::string &indexFolder)
    : _seed(0),
      _d(0),
      _enc_dim(0),
      _n(0),
      _m(0),
      _nt(0),
      _tot_mem(0),
      _a_array(),
      _proj_data(),
      _raw_data(raw_data),
      _trees() {
  loadIndex(indexFolder);
}

iDEC::iDEC(const std::vector<uint64_t> &raw_data, unsigned d, unsigned m,
           unsigned nt, unsigned seed)
    : _seed(seed),
      _d(d),
      _enc_dim(d / 64),
      _n(raw_data.size() / _enc_dim),
      _m(m),
      _nt(nt),
      _tot_mem(0),
      _a_array(),
      _proj_data(),
      _raw_data(raw_data),
      _trees() {
  IDEC_REQUIRED_MSG(_d == _enc_dim * 64,
                    "Only support dimension of multiples of 64");
  IDEC_REQUIRED(raw_data.size() == _n * _enc_dim);
  if (_seed == 0) _seed = DEFAULT_SEED;
}

void iDEC::genRandomWalkVectors() {
  size_t len = (size_t)_enc_dim * _nt * _m;
  _a_array.reserve(len);
  std::bitset<64> bits;

  std::mt19937_64 gen(_seed);
  std::uniform_int_distribution<> dist(0, 1);

  for (size_t i = 0, j = 0; i < len;) {
    bits[j] = dist(gen);
    if (j == 63) {
      _a_array.push_back(bits.to_ullong());
      j = 0;
      ++i;
    } else {
      ++j;
    }
  }

  IDEC_REQUIRED(_a_array.size() == len);
}

int16_t iDEC::getHashValue(const uint64_t *point,  // point
                           unsigned mid,           // project coordinate
                           unsigned tid            // tree id
                           ) const {
  int16_t res = 0;
  size_t start_ind = tid * (_m * _enc_dim) + mid * _enc_dim;
  IDEC_REQUIRED(start_ind + _enc_dim <= _a_array.size());
  const uint64_t *random_vec = &_a_array[start_ind];
  for (unsigned i = 0u; i < _enc_dim; i++) {
    res = calcToWSum(point[i], random_vec[i], res);
  }
  return res;
}

void iDEC::build() {
  if (_a_array.empty()) genRandomWalkVectors();

  _proj_data.resize(_nt);
  const size_t sz_each = (uint64_t) _n * _m; // Unsigned type is not supported for 1B dataset.
  for (unsigned tid = 0; tid < _nt; ++tid) {
#ifndef DISABLE_VERBOSE
    printf("\tBuild tree %u (out of %u)\n", tid + 1, _nt);
    printf("\t\tProjecting data\n");
    printf("\t\tAllocate memory for project data\n");
#endif
    _proj_data[tid].reserve(sz_each);
    for (size_t pid = 0; pid < _n; ++pid) {
      for (unsigned mid = 0; mid < _m; ++mid)
        _proj_data[tid].push_back(
            getHashValue(&_raw_data[pid * _enc_dim], mid, tid));
    }
    _tot_mem += sizeof(int16_t) * sz_each;
#ifndef DISABLE_VERBOSE
    printf("\t\tConstruct a new tree\n");
#endif
    _trees.emplace_back(_m /* dimension */, _proj_data[tid] /* data */,
                        MAX_LEAF);
    //_trees.back().index->buildIndex();
    _tot_mem += _trees.back().index->usedMemory(*(_trees.back().index));
  }
}

std::string iDEC::_getTreeFile(const std::string &indexFolder,
                               unsigned tid) const {
  return indexFolder + "/tree_" + std::to_string(tid) + ".bin";
  ;
}

std::string iDEC::_getProjDataFile(const std::string &indexFolder) const {
  return indexFolder + "/idec_proj_data.flat";
}

std::string iDEC::_getParamFile(const std::string &indexFolder) const {
  return indexFolder + "/idec_param.txt";
}

void iDEC::saveParam(const std::string &indexFolder) const {
  const std::string param_fn = _getParamFile(indexFolder);
  FILE *fp = fopen(param_fn.c_str(), "w");
  IDEC_REQUIRED(fp != NULL);
  fprintf(fp,
          "n = %u\nm = %u\nnumber of trees = %u\nseed = %u\nd = %u\nenc_dim = "
          "%u\n",
          _n, _m, _nt, _seed, _d, _enc_dim);

  fprintf(fp, "%lu", _a_array.front());

  for (size_t i = 1; i < _a_array.size(); ++i) fprintf(fp, " %lu", _a_array[i]);

  fprintf(fp, "\n");
  IDEC_REQUIRED(fclose(fp) == 0);
}

void iDEC::loadParam(const std::string &indexFolder) {
  const std::string param_fn = _getParamFile(indexFolder);
  FILE *fp = fopen(param_fn.c_str(), "r");
  IDEC_REQUIRED(fp != NULL);

  IDEC_REQUIRED(
      fscanf(
          fp,
          "n = %u\nm = %u\nnumber of trees = %u\nseed = %u\nd = %u\nenc_dim = "
          "%u\n",
          &_n, &_m, &_nt, &_seed, &_d, &_enc_dim) == 6);

  uint64_t temp = 0;
  IDEC_REQUIRED(fscanf(fp, "%lu", &temp) == 1);
  const size_t sz = _nt * _m * _enc_dim;
  _a_array.reserve(sz);
  _a_array.push_back(temp);

  while (!feof(fp) && _a_array.size() < sz) {
    IDEC_REQUIRED(fscanf(fp, " %lu", &temp) == 1);
    _a_array.push_back(temp);
  }

  IDEC_REQUIRED(fscanf(fp, "\n") == 0);
  IDEC_REQUIRED(feof(fp) && _a_array.size() == sz);
  IDEC_REQUIRED(fclose(fp) == 0);
}

void iDEC::saveIndex(const std::string &indexFolder) const {
  saveParam(indexFolder);
  {  // save projected
    const std::string pd_fn = _getProjDataFile(indexFolder);
    FILE *fp = fopen(pd_fn.c_str(), "wb");
    IDEC_REQUIRED(fp != NULL);

    auto sz = fwrite(&_proj_data[0][0], sizeof(int16_t), _proj_data[0].size(), fp);
    IDEC_REQUIRED(sz == _proj_data[0].size());
    IDEC_REQUIRED(fclose(fp) == 0);
  }
  for (unsigned tid = 0; tid < _nt; ++tid) {
    std::string tree_fn = _getTreeFile(indexFolder, tid);
    FILE *fp = fopen(tree_fn.c_str(), "wb");
    IDEC_REQUIRED(fp != NULL);
    _trees[tid].index->saveIndex(fp);
    if (fp != NULL) IDEC_REQUIRED(fclose(fp) == 0);
  }
}

void iDEC::loadIndex(const std::string &indexFolder) {
  loadParam(indexFolder);

  if (_a_array.empty()) genRandomWalkVectors();

  _proj_data.resize(_nt);
  const size_t sz_each = _n * _m;

  for (unsigned tid = 0; tid < _nt; ++tid) {
    std::string tree_fn = _getTreeFile(indexFolder, tid);
    FILE *fp = fopen(tree_fn.c_str(), "rb");
    IDEC_REQUIRED(fp != NULL);
#ifndef DISABLE_VERBOSE
    printf("\tLoad tree %u (out of %u)\n", tid + 1, _nt);
    printf("\t\tProjecting data\n");
    printf("\t\tAllocate memory for project data\n");
#endif
    _proj_data[tid].reserve(sz_each);
    for (size_t pid = 0; pid < _n; ++pid) {
      for (unsigned mid = 0; mid < _m; ++mid)
        _proj_data[tid].push_back(
            getHashValue(&_raw_data[pid * _enc_dim], mid, tid));
    }
    _tot_mem += sizeof(int16_t) * sz_each;
#ifndef DISABLE_VERBOSE
    printf("\t\tLoad a new tree from \"%s\"\n", tree_fn.c_str());
#endif
    _trees.emplace_back(_m, _proj_data[tid], MAX_LEAF);
    //_trees.back().index->loadIndex(fp);
    _tot_mem += _trees.back().index->usedMemory(*(_trees.back().index));
    if (fp != NULL) IDEC_REQUIRED(fclose(fp) == 0);
#ifndef DISABLE_VERBOSE
    printf("\t\tLoading succeeded\n");
#endif
  }
}

void iDEC::_one_nn(const uint64_t *query, float t, float &dist,
                   size_t &index) const {
  dist = _d;
  std::vector<int16_t> q_proj;
  q_proj.resize(_m);

  unsigned nn_kd = std::round(t * _n);

#ifdef SAVE_PROJ_DATA
  // only for debug purpose
  FILE *fp = fopen("query_proj.dat", "ab+");
  IDEC_REQUIRED(fp != NULL);
#endif
  for (unsigned tid = 0; tid < _nt; ++tid) {
    // calculate projection for query
    for (unsigned mid = 0; mid < _m; ++mid) {
      q_proj[mid] = getHashValue(query, mid, tid);
    }
#ifdef SAVE_PROJ_DATA
    IDEC_REQUIRED(fwrite(&q_proj[0], sizeof(int16_t), q_proj.size(), fp) ==
                  q_proj.size());
    IDEC_REQUIRED(fclose(fp) == 0);
#endif
    std::vector<size_t> ret_indexes(nn_kd);
    // std::vector<int16_t> out_dists_sqr(nn_kd);
    std::vector<float> out_dists_sqr(nn_kd);
    _trees[tid].query(&q_proj[0], nn_kd, &ret_indexes[0], &out_dists_sqr[0]);
    for (unsigned kid = 0; kid < nn_kd; ++kid) {
      float dist_tmp = hammingDistance(
          query, &_raw_data[ret_indexes[kid] * _enc_dim], _enc_dim);
      if (dist_tmp < dist) {
        dist = dist_tmp;
        index = ret_indexes[kid];
      }
    }
#ifdef USE_SEARCH_HISTORY
    // process search history
    for (const auto &p : _trees[tid].index->getHistory()) {
      if (p.first < _CHI2INV_ * dist) {
        // with a prob of
        float dist_tmp =
            hammingDistance(query, &_raw_data[p.second * _enc_dim], _enc_dim);
        if (dist_tmp < dist) {
          dist = dist_tmp;
          index = p.second;
        }
      }
    }
#endif  // USE_SEARCH_HISTORY
  }
}

void iDEC::knn(const uint64_t *query, unsigned k, float t,
               std::vector<std::pair<float, size_t>> &ans) const {
  if (k == 1) {
    ans.resize(1);
    _one_nn(query, t, ans[0].first, ans[0].second);
    return;
  }
  ans.clear();
  ans.reserve(k);
  std::vector<int16_t> q_proj;
  q_proj.resize(_m);

  unsigned nn_kd = std::round(t * _n);
  IDEC_REQUIRED(nn_kd > 0 && nn_kd < _n);

  for (unsigned tid = 0; tid < _nt; ++tid) {
    // calculate projection for query
    for (unsigned mid = 0; mid < _m; ++mid) {
      q_proj[mid] = getHashValue(query, mid, tid);
    }
    std::vector<size_t> ret_indexes(nn_kd);
    // std::vector<int16_t> out_dists_sqr(nn_kd);
    std::vector<float> out_dists_sqr(nn_kd);
    _trees[tid].query(&q_proj[0], nn_kd, &ret_indexes[0], &out_dists_sqr[0]);

    for (unsigned j = 0; j < nn_kd; ++j) {
      float dist_tmp = hammingDistance(
          query, &_raw_data[ret_indexes[j] * _enc_dim], _enc_dim);

      if (ans.size() < k - 1) {
        ans.emplace_back(dist_tmp, ret_indexes[j]);
      } else if (ans.size() == k - 1) {
        ans.emplace_back(dist_tmp, ret_indexes[j]);
        std::make_heap(ans.begin(), ans.end());
      } else {
        if (dist_tmp < ans.front().first) {
          std::pop_heap(ans.begin(), ans.end());
          ans.pop_back();
          ans.emplace_back(dist_tmp, ret_indexes[j]);
          std::push_heap(ans.begin(), ans.end());
        }
      }
    }

#ifdef USE_SEARCH_HISTORY
    // process search history
    // TODO: using threshold
    auto _cur_worst = ans.front().first;
    for (const auto &p : _trees[tid].index->getHistory()) {
      if (p.first < _CHI2INV_ * _cur_worst) {
        // with a prob of
        float dist_tmp =
            hammingDistance(query, &_raw_data[p.second * _enc_dim], _enc_dim);
        if (dist_tmp < _cur_worst) {
          std::pop_heap(ans.begin(), ans.end());
          ans.pop_back();
          ans.emplace_back(dist_tmp, p.second);
          std::push_heap(ans.begin(), ans.end());
          _cur_worst = ans.front().first;
        }
      }
    }
#endif  // USE_SEARCH_HISTORY
  }

  std::sort_heap(ans.begin(), ans.end());
}

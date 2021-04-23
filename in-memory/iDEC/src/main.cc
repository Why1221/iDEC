#include "iDEC.h"

#include "IdecException.h"

#include <AnnResultWriter.hpp>
#include <Timer.hpp>

#include<MemUtils.h>

#include <cstring>
#include <iostream>

const unsigned DEFAULT_NUM_TREES = 1u;

// print parameters to stdout
void show_params(const char* fmt, ...) {
  va_list args;
  va_start(args, fmt);

  while (*fmt != '\0') {
    char* name = va_arg(args, char*);
    if (*fmt == 'i') {
      int val = va_arg(args, int);
      printf("%s: %d\n", name, val);
    } else if (*fmt == 'c') {
      int val = va_arg(args, int);
      printf("%s: \'%c\'\n", name, val);
    } else if (*fmt == 'f') {
      double val = va_arg(args, double);
      printf("%s: %f\n", name, val);
    } else if (*fmt == 's') {
      char* val = va_arg(args, char*);
      printf("%s: \"%s\"\n", name, val);
    } else {
      fprintf(stderr, "Unsupported format");
      return;
    }
    ++fmt;
  }

  va_end(args);
}

void usage() {
  printf("iDEC (v1.0)\n");
  printf("Options\n");
  printf("-d {value}     \trequired \tdimensionality\n");
  printf("-ds {string}   \trequired \tdataset file\n");
  printf("-n {value}     \trequired \tnumber of points in database\n");
  printf("-m {value}     \trequired \tdimension of projection space\n");
  printf("-qn {value}     \trequired (only for query) \tnumber of queries\n");
  printf(
      "-L {value}     \toptional (for indexing) \tnumber of trees (default: "
      "1)\n");
  printf(
      "-t {value}     \trequired (only for query) \tparameter t (maximum "
      "checked ratio))\n");
  printf("-if {string}   \trequired \tfile path for iDEC index folder\n");
  printf(
      "-k {value}     \toptional (only for query) \tnumber of neighbors "
      "(default: 1)\n");
  printf(
      "-gt {string}   \trequired (only for query) \tfile of exact results\n");
  printf("-rf {string}   \trequired (only for query) \tresult file\n");
  printf("-qs {string}   \trequired (only for query) \tfile of query set\n");

  printf("\n");
  printf("Build an iDEC index\n");
  printf("-d -n -ds -if -m [-L]\n");
  printf("Use iDEC to answer a knn workload\n");
  printf("-d -n -ds -qn -qs -if -rf -gt -t [-k]\n");

  printf("\n");
}

// load data from binary file
void load_data(const char* filename,         // filename for data file
               std::vector<uint64_t>& data,  // data (output)
               unsigned num,                 // number of points
               unsigned enc_dim              // dimension after encoding
               )  {
  FILE* fp = fopen(filename, "rb");
  IDEC_REQUIRED(fp != NULL);
  auto tot_n = (size_t)num * enc_dim;
  data.resize(tot_n);
  IDEC_REQUIRED(fread(&data[0], sizeof(uint64_t), tot_n, fp) == tot_n);
  IDEC_REQUIRED(fclose(fp) == 0);
}

void indexing(const char* ds_filename,  // database filename
              const char* i_filename,   // index folder
              unsigned num,             // # of points
              unsigned dim,             // dimension
              unsigned m,
              unsigned num_trees  // # of trees
) {
  std::vector<uint64_t> train;
  auto enc_dim = dim / 64;
  load_data(ds_filename, train, num, enc_dim);

  iDEC idec(train, dim, m, num_trees);

  std::string perf_filename = std::string(i_filename) + "/idec-indexing.txt";

  AnnResultWriter writer(perf_filename);
  writer.writeRow(
      "s", "dsname,#n,#dim,#trees,m,index_size(bytes),construction_time(us),Peak_mem(bytes)");
  const char* fmt = "siiiiifi";

  HighResolutionTimer timer;
  timer.restart();
  idec.build();
  auto e = timer.elapsed();
  auto peak_mem = getPeakRSS();
  
  idec.saveIndex(i_filename);

  auto isz = idec.size();
  
  writer.writeRow(fmt, ds_filename, num, dim, num_trees, m, isz, e,peak_mem);
}

void knn(const char* ds_filename,  // database filename
         const char* q_filename,   // query filename
         const char* i_filename,   // index filename
         const char* gt_filename,  // ground truth filename
         const char* r_filename,   // result filename
         unsigned num,             // # of points in database
         unsigned dim,             // dimensionality
         unsigned qn,              // # of queries
         unsigned K,               // # of NNs
         float t                   // parameter t
) {
  std::vector<uint64_t> train, test;
  auto enc_dim = dim / 64;
  load_data(ds_filename, train, num, enc_dim);
  load_data(q_filename, test, qn, enc_dim);

  iDEC idec(train, i_filename);

  unsigned r_qn, r_maxk;
  FILE* fp = fopen(gt_filename, "r");
  IDEC_REQUIRED(fp != NULL);
  IDEC_REQUIRED(fscanf(fp, "%d %d\n", &r_qn, &r_maxk) >= 0);
  IDEC_REQUIRED(r_qn >= qn && r_maxk >= K);

  std::vector<float> gt(qn * r_maxk, -1.0f);

  for (unsigned i = 0; i < qn; ++i) {
    unsigned j;
    IDEC_REQUIRED(fscanf(fp, "%d", &j) >= 0);
    IDEC_REQUIRED(j == i);
    for (j = 0; j < r_maxk; ++j) {
      IDEC_REQUIRED(fscanf(fp, " %f", &gt[i * r_maxk + j]) >= 0);
#ifdef DEBUG
      printf("(%d,%d): %.6f\n", i, j, gt[i * r_maxk + j]);
#endif
    }
    IDEC_REQUIRED(fscanf(fp, "\n") >= 0);
  }
  IDEC_REQUIRED(fclose(fp) == 0);

  HighResolutionTimer timer;
  AnnResultWriter writer(r_filename);

  writer.writeRow("s", AnnResults::_DEFAULT_HEADER_I_);

  #ifdef DEBUG
  printf("write header finished\n");
  printf("qn: %u, # of queries: %lu, K: %d, enc_dim: %u\n", qn, test.size(), K, enc_dim);
  #endif

  for (unsigned i = 0; i < qn; i++) {
    std::vector<std::pair<float, size_t>> res;
    if(std::round(t * num) == 0) {
      fprintf(stderr, "t is too small\n");
      continue;
    }
    res.resize(K);
    timer.restart();
    idec.knn(&test[i * enc_dim], K, t, res);
    auto query_time = timer.elapsed();

    for (unsigned j = 0; j < K; ++j) {
      float dist = res[j].first;
      float gdist = std::round(gt[i * r_maxk + j] * gt[i * r_maxk + j]);
      writer.writeRow(AnnResults::_DEFAULT_FMT_I_, i, j, res[j].second, (int)dist,
                      (int)gdist, (dist / gdist), query_time);
    }
  }
}

/*
 * Get the index of next unblank char from a string.
 */
int GetNextChar(char* str) {
  int rtn = 0;

  // Jump over all blanks
  while (str[rtn] == ' ') {
    rtn++;
  }

  return rtn;
}

/*
 * Get next word from a string.
 */
void GetNextWord(char* str, char* word) {
  // Jump over all blanks
  while (*str == ' ') {
    str++;
  }

  while (*str != ' ' && *str != '\0') {
    *word = *str;
    str++;
    word++;
  }

  *word = '\0';
}

int main(int argc, char** argv) {
  // These two are global variables
  unsigned nPoints = 0;          // the number of points
  unsigned pointsDimension = 0;  // the dimensionality of points

  int qn = -1;  // the number of queries
  int k = 1;    // the k of k-NN

  char ds[200] = "";   // the file path of dataset
  char qs[200] = "";   // the file path of query set
  char gt[200] = "";   // the file path of ground truth
  char rf[200] = "";   // the folder path of results
  char inf[200] = "";  // the folder path of index

  unsigned num_trees = DEFAULT_NUM_TREES;
  unsigned m = 0;
  float t = -1;

  int cnt = 1;
  bool failed = false;
  char* arg;
  int i;
  char para[10];

  std::string err_msg;

  while (cnt < argc && !failed) {
    arg = argv[cnt++];
    if (cnt == argc) {
      failed = true;
      break;
    }

    i = GetNextChar(arg);
    if (arg[i] != '-') {
      failed = true;
      err_msg = "Wrong format!";
      break;
    }

    GetNextWord(arg + i + 1, para);

    arg = argv[cnt++];

    if (strcmp(para, "n") == 0) {
      nPoints = atoi(arg);
      if (nPoints <= 0) {
        failed = true;
        err_msg = "n should a positive integer!";
        break;
      }
    } else if (strcmp(para, "d") == 0) {
      pointsDimension = atoi(arg);
      if (pointsDimension <= 0) {
        failed = true;
        err_msg = "d should a positive integer!";
        break;
      }
    } else if (strcmp(para, "qn") == 0) {
      qn = atoi(arg);
      if (qn <= 0) {
        failed = true;
        err_msg = "qn should a positive integer!";
        break;
      }
    } else if (strcmp(para, "k") == 0) {
      k = atoi(arg);
      if (k <= 0) {
        failed = true;
        err_msg = "k should a positive integer!";
        break;
      }
    } else if (strcmp(para, "L") == 0) {
      num_trees = atoi(arg);
      if (num_trees <= 0) {
        failed = true;
        err_msg = "L should a positive integer!";
        break;
      }
    } else if (strcmp(para, "m") == 0) {
      m = atoi(arg);
      if (m <= 0) {
        failed = true;
        err_msg = "R should a positive integer!";
        break;
      }
    } else if (strcmp(para, "t") == 0) {
      t = atof(arg);
      if (!(t > 0.0f && t < 1.0f)) {
        failed = true;
        err_msg = "t should a positive number within (0,1)!";
        break;
      }
    } else if (strcmp(para, "ds") == 0) {
      GetNextWord(arg, ds);

    } else if (strcmp(para, "qs") == 0) {
      GetNextWord(arg, qs);

    } else if (strcmp(para, "gt") == 0) {
      GetNextWord(arg, gt);

    } else if (strcmp(para, "if") == 0) {
      GetNextWord(arg, inf);

    } else if (strcmp(para, "rf") == 0) {
      GetNextWord(arg, rf);

    } else {
      failed = true;
      fprintf(stderr, "Unknown option -%s!\n\n", para);
    }
  }

  if (failed) {
    fprintf(stderr, "%s:%d: %s\n\n", __FILE__, __LINE__, err_msg.c_str());
    usage();
    return EXIT_FAILURE;
  }

  int nargs = (cnt - 1) / 2;

  if (!(nargs == 5 || nargs == 6 || nargs == 9 || nargs == 10)) {
    fprintf(stderr, "%s:%d: %s\n\n", __FILE__, __LINE__,
            "Wrong number of arguements!");
    usage();
    return EXIT_FAILURE;
  }

#ifndef DISABLE_VERBOSE
  printf("=====================================================\n");
  show_params("iiiiiifsssss", "# of points", nPoints, "dimension",
              pointsDimension, "# of queries", qn, "number of trees", num_trees,
              "projection dimension", m, "# of nns (i.e., k)", k,
              "maximum checked ratio", t, "dataset filename", ds,
              "index folder", inf, "result filename", rf,
              "ground truth filename", gt, "query filename", qs);
  printf("=====================================================\n");
#endif

  try {
    if (nargs == 5 || nargs == 6) {
      IDEC_REQUIRED(strlen(ds) != 0);
      IDEC_REQUIRED(strlen(inf) != 0);
      IDEC_REQUIRED(num_trees > 0 && m > 0 && nPoints > 0 &&
                    pointsDimension > 0);
      indexing(ds, inf, nPoints, pointsDimension, m, num_trees);
    } else {
      IDEC_REQUIRED(strlen(ds) != 0);
      IDEC_REQUIRED(strlen(qs) != 0);
      IDEC_REQUIRED(strlen(gt) != 0);
      IDEC_REQUIRED(strlen(inf) != 0);
      IDEC_REQUIRED(strlen(rf) != 0);
      IDEC_REQUIRED(nPoints > 0 && pointsDimension > 0 && qn > 0 && k > 0 &&
                    t > 0);
      knn(ds, qs, inf, gt, rf, nPoints, pointsDimension, qn, k, t);
    }
  } catch (const IdecException& e) {
    std::cerr << e.what() << std::endl;
  } catch (const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }

  return EXIT_SUCCESS;
}
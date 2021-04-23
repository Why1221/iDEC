
#include <errno.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include <Exception.h>
#include <AnnResultWriter.hpp>
#include <Timer.hpp>

#include "Hdf5File.h"
#include "linear_scan_l1.h"

const char *_HEADER_E_ = "#qid,#kid,qtime(us),#io";
const char *_FMT_E_ = "iifi";
const char *_HEADER_I_ = "#qid,#kid,qtime(us)";
const char *_FMT_I_ = "iif";

/*
 * Prints the usage of the linear_scan_main.
 */
void usage(const char *progname) {
  printf("Linear Scan Usgae: %s options\n", progname);
  printf("Option\n");
  printf("-n  {integer}\t required\t the number of points\n");
  printf("-qn {integer}\t required\t the number of queries\n");
  printf("-d  {integer}\t required\t the dimensionality of points\n");
  printf("-ds {string}\t required\t the file path of dataset\n");
  printf("-k  {integer}\t required\t the k of k-NN\n");
  printf(
      "-b  {integer}\t required (only for external-memory)\t the page size\n");

  printf("-rf {string}\t required\t the file path for storing distance results\n");
  printf("-ri {string}\t required\t the file path for storing index results\n");
  printf("-pf {string}\t required\t the file path for storing performance\n");
  printf("\n");
}

/*
 * Get the index of next unblank char from a string.
 */
int GetNextChar(char *str) {
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
void GetNextWord(char *str, char *word) {
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

void linear_scan_ram(
    const std::string &ds_filename,    // dataset file (MUST be a hdf5 file)
    const std::string &res_filename,   // result file
    const std::string &res_index_filename,   // result file for index 
    const std::string &perf_filename,  // performance file
    int n,                             // dataset size
    int dim,      // dataset dimension 
    int qn,       // number of queries
    int k         // number of nearest neighbors
) {
  std::vector<float> train, query;

  //   try
  {
    Hdf5File h5(ds_filename, Hdf5File::Mode::ReadOnly);

    h5.read<float>(train, "train");
    h5.read<float>(query, "test");

    const int n_ = train.size() / dim;
    const int qn_ = query.size() / dim;

    NPP_ASSERT(n_ == n);
    NPP_ASSERT(qn_ == qn);
  }

  std::vector<std::pair<size_t, float>> results;
  HighResolutionTimer timer;
  AnnResultWriter writer(perf_filename);
  writer.writeRow("s",_HEADER_I_);

  FILE *fp = fopen(res_filename.c_str(), "w");
  if (!fp)
    throw npp::Exception("Failed to open file " + res_filename +
                             ". Error: " + std::string(strerror(errno)),
                         __FILE__, __LINE__);

  fprintf(fp, "%d %d\n", qn, k);

  FILE *fp_2 = fopen(res_index_filename.c_str(), "w");
  if (!fp_2)
    throw npp::Exception("Failed to open file " + res_index_filename +
                             ". Error: " + std::string(strerror(errno)),
                         __FILE__, __LINE__);
  fprintf(fp_2, "%d %d\n", qn, k); 

  for (int i = 0; i < qn; ++i) {
    results.clear();
    timer.restart();
    linear_scan(train, &query[i * dim], dim, k, results);
    auto exec_time = timer.elapsed();
    fprintf(fp, "%d", i);
    fprintf(fp_2, "%d", i);
    for (int j = 0; j < k; ++j) {
      writer.writeRow(_FMT_I_, (i + 1), (j + 1), (j == k - 1 ? exec_time : -1));
      fprintf(
          fp, " %.2f",
          results[j].second);  // Here, we stored L_1 distance
      fprintf(
          fp_2, " %lu",
          results[j].first);  // Here, we stored raw index
    }
    fprintf(fp, "\n");
    fprintf(fp_2, "\n");
  }
  NPP_ASSERT(fclose(fp) == 0);
  NPP_ASSERT(fclose(fp_2) == 0);
}

namespace {
using Pair = std::pair<size_t, float>;
class PairCmp {
 public:
  bool operator()(const Pair &a, const Pair &b) const {
    return a.second < b.second;
  }
};
}

int main(int argc, char **argv) {
  char *progname;
  char *p;
  progname = ((p = strrchr(argv[0], '/')) ? ++p : argv[0]);

  int n = -1;
  int dim = -1;
  int qn = -1;
  int k = -1;
  int b = -1;
  char ds[200] = "";  // the file path of dataset
  char pf[200] = "";  // the folder path of
  char rf[200] = "";  // the folder path of results
  char ri[200] = "";  // the folder path of index results

  int cnt = 1;
  bool failed = false;
  char *arg;
  int i;
  char para[40];

  while (cnt < argc && !failed) {
    arg = argv[cnt++];
    if (cnt == argc) {
      failed = true;
      break;
    }

    i = GetNextChar(arg);
    if (arg[i] != '-') {
      failed = true;
      break;
    }

    GetNextWord(arg + i + 1, para);

    arg = argv[cnt++];

    if (strcmp(para, "n") == 0) {
      n = atoi(arg);
      if (n <= 0) {
        failed = true;
        break;
      }
    } else if (strcmp(para, "d") == 0) {
      dim = atoi(arg);
      if (dim <= 0) {
        failed = true;
        break;
      }
    } else if (strcmp(para, "qn") == 0) {
      qn = atoi(arg);
      if (qn <= 0) {
        failed = true;
        break;
      }
    } else if (strcmp(para, "k") == 0) {
      k = atoi(arg);
      if (k <= 0) {
        failed = true;
        break;
      }
    } else if (strcmp(para, "b") == 0) {
      b = atoi(arg);
      if (b <= 0) {
        failed = true;
        break;
      }
    } else if (strcmp(para, "ds") == 0) {
      GetNextWord(arg, ds);
    } else if (strcmp(para, "pf") == 0) {
      GetNextWord(arg, pf);

    }  else if (strcmp(para, "ri") == 0) {
      GetNextWord(arg, ri);

    } 
    else if (strcmp(para, "rf") == 0) {
      GetNextWord(arg, rf);

    } else {
      failed = true;
      printf("Unknown option -%s!\n\n", para);
    }
  }

  auto nargs = cnt / 2;
  if (failed || (!(nargs == 8 || nargs == 9))) {
    usage(progname);
    return -1;
  }

  try {
    linear_scan_ram(std::string(ds), std::string(rf),std::string(ri), std::string(pf), n,
                      dim, qn, k);
  } 
  catch (const npp::Exception &e) {
    std::cerr << e.what() << std::endl;
    return -1;
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return -1;
  }

  return 0;
}
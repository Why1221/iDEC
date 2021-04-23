/*
 * File:  main-hamming.cpp
 * By:    yifangs
 */

#include "gendef.h"
#include "ptR-tree.h"
#include "srs-compact.h"
#include <cassert>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <time.h>

using namespace std;

/*****************************************************************
builds SRS tree(s) by incremental insertions

-para-
dfname		dataset file
ffolder		forest folder
n			    cardinality
d			    dimensionality
B			    page size
L			    number of srs-trees
m         dimensionality of the projected space

-return-
0			success
1			failure
*****************************************************************/

int build(char *_dfname, char *_ffolder, int _n, int _d, int _B, int _L,
          int _m) {
  int ret = 0;
  SRS *srs = NULL;

  printf("Build SRS trees.\n\n");

  srs = new SRS();
  srs->init(_d, _n, _B, _L, _m);

  if (srs->buildFromFile(_dfname, _ffolder)) {
    ret = 1;
    goto recycle;
  }

recycle:
  if (srs) {
    delete srs;
    srs = NULL;
  }

  return ret;
}

/*****************************************************************
continue builds a srs tree/forest by incremental insertions
*****************************************************************/

int ctbuild(char *_dfname, char *_ffolder, int _n, int _d, int _B, int _L,
            int _m) {
  int ret = 0;
  SRS *srs = NULL;

  printf("Build SRS trees.\n\n");

  srs = new SRS();

  // srs->quiet = 5;
  srs->init(_d, _n, _B, _L, _m);

  if (srs->ctBuildFromFile(_dfname, _ffolder)) {
    ret = 1;
    goto recycle;
  }

recycle:
  if (srs) {
    delete srs;
    srs = NULL;
  }

  return ret;
}

/*****************************************************************
builds binary data file

-para-
dfname    dataset file
ffolder   forest folder
n         cardinality
d         dimensionality
B         page size

-return-
0     success
1     failure
*****************************************************************/

int buildBin(char *_dfname, char *_ffolder, int _n, int _d, int _B) {
  int ret = 0;
  SRS *srs = NULL;

  printf("Build binary file.\n\n");

  srs = new SRS();
  srs->init(_d, _n, _B, 1, 1);

  if (srs->buildBinFromFile(_dfname, _ffolder)) {
    ret = 1;
    goto recycle;
  }

recycle:
  if (srs) {
    delete srs;
    srs = NULL;
  }

  return ret;
}

/*****************************************************************
this function finds the exact k NN distances of a workload of
queries. these results are written to a file for computing the
approximation ratios of approximate methods.

-para-
dfname			dataset file
n				cardinality
d				dimensionality
qfname			query file
k				number of neighbors
ofname			output file

-return-
0		success
1		failure

-by-
yufei tao
*****************************************************************/

int seqscan(char *_dfname, int _n, int _d, char *_qfname, int _k,
            char *_ofname) {
  int ret = 0;

  FILE *dfp = NULL;
  FILE *qfp = NULL;
  FILE *ofp = NULL;
  float dist = -1;
  float *knndist = NULL;
  int c = -1;
  // int * ds = NULL;
  uint64_t *ds = NULL;
  int i = -1;
  int ii = -1;
  int j = -1;
  int rsltCnt = -1;
  // int * q = NULL;
  uint64_t *q = NULL;
  int qn = -1;
  unsigned enc_dim = _d / 64;

  dfp = fopen(_dfname, "rb");
  if (!dfp) {
    printf("Could not open the data file.\n");

    ret = 1;
    goto recycle;
  }

  qfp = fopen(_qfname, "r");
  if (!qfp) {
    printf("Could not open the query file.\n");

    ret = 1;
    goto recycle;
  }

  ofp = fopen(_ofname, "w");
  if (!ofp) {
    printf("Could not create the output file.\n");

    ret = 1;
    goto recycle;
  }

  printf("Finding the accurate results (nearest neighbor)\n");

  // ds = new int[_n * _d];
  ds = new uint64_t[_n * enc_dim];

  c = fread(ds, sizeof(uint64_t), _n * enc_dim, dfp) / enc_dim;

  for (int i = 0;i < 5;++ i) {
    for (int j = 0;j < enc_dim;++ j) fprintf(stdout, "%llu ", ds[i*enc_dim +j]);
    fprintf(stdout, "\n");
  }

  // while (!feof(dfp) && c < _n) {
  //   fscanf(dfp, "%d", &j);

  //   for (i = 0; i < _d; i++)
  //     fscanf(dfp, " %d", &(ds[c * _d + i]));

  //   fscanf(dfp, "\n");

  //   c++;
  // }

  // if (!feof(dfp) && c == _n) {
  uint64_t temp;
  if (c == _n && fread(&temp, sizeof(uint64_t), 1, dfp) != 0) {
    printf("Dataset larger than you said.\n");

    ret = 1;
    goto recycle;
  } else if (feof(dfp) && c < _n) {
    printf("Set the dataset size to %d and try again.\n", c);

    ret = 1;
    goto recycle;
  }

  // q = new int[_d];
  q = new uint64_t[enc_dim];

  knndist = new float[_k];

  qn = 0;

  while (!feof(qfp)) {
    fscanf(qfp, "%d", &c);

    // for (i = 0; i < _d; i++)
    // fscanf(qfp, " %d", &(q[i]));

    for (i = 0; i < enc_dim; i++)
      fscanf(qfp, " %lu", &(q[i]));
    fscanf(qfp, "\n");

    qn++;
  }

  fprintf(ofp, "%d %d\n", qn, _k);

  if (fseek(qfp, 0, SEEK_SET)) {
    printf("Could not rewind to the beginning of the query file.\n");

    ret = 1;
    goto recycle;
  }

  c = 0;
  while (!feof(qfp)) {
    fscanf(qfp, "%d", &c);

    // for (i = 0; i < _d; i++)
    //   fscanf(qfp, " %d", &(q[i]));

    for (i = 0; i < enc_dim; i++)
      fscanf(qfp, " %lu", &(q[i]));
    fscanf(qfp, "\n");

    for (i = 0; i < _k; i++)
      knndist[i] = (float)MAXREAL;

    for (i = 0; i < _n; i++) {
      // dist = l2_dist_int(&(ds[i * _d]), q, _d);
      dist = l2_dist_int(&(ds[i * enc_dim]), q, enc_dim);

      for (j = 0; j < _k; j++) {
        if (compfloats(dist, knndist[j]) == -1)
          break;
      }

      if (j < _k) {
        for (ii = _k - 1; ii >= j + 1; ii--)
          knndist[ii] = knndist[ii - 1];

        knndist[j] = dist;
      }
    }

    fprintf(ofp, "%d", c);

    for (i = 0; i < _k; i++) {
      fprintf(ofp, " %f", knndist[i]);
    }

    fprintf(ofp, "\n");

    ++c;
  }

recycle:
  if (dfp)
    fclose(dfp);

  if (qfp)
    fclose(qfp);

  if (ofp)
    fclose(ofp);

  if (ds)
    delete[] ds;

  if (q)
    delete[] q;

  if (knndist)
    delete[] knndist;

  return ret;
}

/*****************************************************************
answer a workload of knn queries

-para-
ffolder			folder of the lsb-forest
numTrees		number of lsb-trees used for querying
rfname			file of accurate results
qfname			query file
k				number of neighbors
ofname			output file

-return-
0		success
1		failure
*****************************************************************/

int srsknn(char *_ffolder, int _numTrees, char *_rfname, char *_qfname, int _k,
           char *_ofname, int _start_id, float _c, float _tau, float _t,
           bool early_stop) {
  int ret = 0;

  FILE *rfp = NULL;
  FILE *ofp = NULL;
  FILE *qfp = NULL;
  float dist = -1;
  float *knndist = NULL;
  float overallRatioSum = -1;
  float thisOverallRatio = -1;
  float *R = NULL;
  float p = -1;
  int c = -1;
  int costSum = -1;
  int *ds = NULL;
  int i = -1;
  int j = -1;
  int maxk = -1;
  // int * q = NULL;
  uint64_t *q = NULL;
  int qn = -1;
  int rsltCnt = -1;
  int thiscost = -1;
  unsigned enc_dim = 0;
  SRS *srs = NULL;
  SRS_Hentry *rslt = NULL;

  srs = new SRS();

  if (srs->restore(_ffolder, _start_id, _numTrees)) {
    ret = 1;
    goto recycle;
  }

  rfp = fopen(_rfname, "r");
  if (!rfp) {
    printf("Could not open the file of accurate results.\n");
    ret = 1;
    goto recycle;
  }

  qfp = fopen(_qfname, "r");
  if (!qfp) {
    printf("Could not open the query file.\n");
    ret = 1;
    goto recycle;
  }

  ofp = fopen(_ofname, "w");
  if (!ofp) {
    printf("Could not create the output file.\n");
    ret = 1;
    goto recycle;
  }

  fscanf(rfp, "%d %d\n", &qn, &maxk);

  if (_k > maxk) {
    printf("Result file good for k <= %d only.\n", maxk);

    ret = 1;
    goto recycle;
  }

  R = new float[qn * maxk];

  for (i = 0; i < qn; i++) {
    fscanf(rfp, "%d", &j);

    for (j = 0; j < maxk; j++)
      fscanf(rfp, "%f", &(R[i * maxk + j]));
  }

  rslt = new SRS_Hentry[_k];

  for (i = 0; i < _k; i++) {
    rslt[i].d = srs->d;
  }

  overallRatioSum = 0;

  c = 0;

  costSum = 0;

  assert(srs->d % 64 == 0 && "Only support mutiples of 64");

  enc_dim = srs->d / 64;
  q = new uint64_t[srs->d];

  fprintf(ofp, "cost\toverall Ratio\n");

  while (!feof(qfp)) {
    if (c == qn) {
      printf("Query file longer than expected. I thought there were %d queries "
             "only.\n",
             qn);

      ret = 1;
      goto recycle;
    }

    fscanf(qfp, "%d", &i);

    // for (i = 0; i < srs->d; i++)
    //   fscanf(qfp, "%d", &(q[i]));

    for (i = 0; i < enc_dim; i++)
      fscanf(qfp, "%llu", &(q[i]));

    fscanf(qfp, "\n");

    thiscost =
        srs->knn(q, _k, rslt, _start_id, _numTrees, _c, _tau, _t, early_stop);

    thisOverallRatio = 0;

    for (i = 0; i < _k; i++) {
      p = rslt[i].dist / R[c * maxk + i];
      thisOverallRatio += p;
    }
    thisOverallRatio /= _k;

    fprintf(ofp, "%d\t%f\n", thiscost, thisOverallRatio);
    // printf("%d\t%f\n", thiscost, thisOverallRatio);
    overallRatioSum += thisOverallRatio;
    costSum += thiscost;

    c++;
  }

  fprintf(ofp, "------------------------------------------\n");
  fprintf(ofp, "%d\t%f\n", costSum / c, overallRatioSum / c);

recycle:
  if (rfp)
    fclose(rfp);

  if (qfp)
    fclose(qfp);

  if (ofp)
    fclose(ofp);

  if (srs)
    delete srs;

  if (R)
    delete[] R;

  if (q)
    delete[] q;

  if (rslt) {
    delete[] rslt;
  }

  return ret;
}

void usage() {
  printf("SRS (v1.1)\n");
  printf("Options\n");
  printf("-b {value}\tpage size in bytes (need to be a multiple of 4)\n");
  printf("-d {value}\tdimensionality\n");
  printf("-ds {string}\tdataset file\n");
  // printf("-l {value}\tnumber of srs-trees in the index (default 1)\n");
  printf("-n {value}\tcardinality\n");
  printf("-m {value}\tdimensionality of the projected space\n");
  printf("-f {string}\tfolder for srs index\n");
  printf("-k {value}\tnumber of neighbors wanted\n");
  printf("-r {string}\tfile of exact results\n");
  printf("-q {string}\tfile of query set\n");
  printf("-o {string}\toutput file\n");
  printf("-c {value}\tapproximate ratio (>= 1.0, default 4.0)\n");
  printf("-tau {value}\tthreshold of early termination condition (default "
         "0.1809)\n");
  printf("-t {value}\tmaximum percentage of examed points (default 0.00242)\n");
  printf(
      "-es [0|1]\t1: enable early stop; 0: disable early stop (default 1)\n");
  // printf("-si {value}\tstart index in query processing (default 0)\n");
  printf("-ct [0|1]\t1: insert into old index, otherwise build new index\n");
  printf("\n");
  printf("Build a srs index\n");
  printf("-b -d -ds -n -m -f\n");
  printf("Continue build a index\n");
  printf("-b -d -ds -n -m -f -ct 1\n");
  printf("Generate the binary data file\n");
  printf("-b -d -ds -n -f\n");
  printf("Exact knn results for a query workload\n");
  printf("-d -ds -k -n -r -q\n");
  printf("Use SRS to answer a knn workload\n");
  printf("-f -k -o -r -q [-c -tau -t -es]\n");
}

int main(int argc, char **argv) {

  bool failed = false;
  char c = -1;
  char *cp = NULL;
  char dfname[100] = "";
  char ffolder[100] = "";
  char ofname[100] = "";
  char qfname[100] = "";
  char rfname[100] = "";
  char para[100] = "";
  clock_t stTime = -1;
  clock_t edTime = -1;
  int cnt = -1;
  int Bb = -1;
  int d = -1;
  int k = -1;
  int L = 1;
  int n = -1;
  int nt = 1;
  int m = -1;
  int startID = 0;
  float tau = 0.1809;
  float ratio = 4.0;
  float t = 0.00242;
  int es = 1;
  int ct = 0;

  cnt = 1;

  while (cnt < argc && !failed) {
    cp = argv[cnt];

    c = getnextChar(&cp);
    if (c != '-') {
      failed = true;
      break;
    }

    getnextWord(&cp, para);

    cnt++;
    if (cnt == argc) {
      failed = true;
      break;
    }

    cp = argv[cnt];

    if (strcmp(para, "ds") == 0) {
      getnextWord(&cp, dfname);
    } else if (strcmp(para, "f") == 0) {
      getnextWord(&cp, ffolder);
    } else if (strcmp(para, "q") == 0) {
      getnextWord(&cp, qfname);
    } else if (strcmp(para, "o") == 0) {
      getnextWord(&cp, ofname);
    } else if (strcmp(para, "r") == 0) {
      getnextWord(&cp, rfname);
    } else if (strcmp(para, "n") == 0) {
      n = atoi(cp);
      if (n <= 0) {
        failed = true;
        break;
      }
    } else if (strcmp(para, "d") == 0) {
      d = atoi(cp);
      if (d <= 0) {
        failed = true;
        break;
      }
    } else if (strcmp(para, "b") == 0) {
      Bb = atoi(cp);
      if (Bb <= 0 || (Bb / 4 * 4 != Bb)) {
        printf("Illegal page size.\n\n");
        failed = true;
        break;
      }
    } else if (strcmp(para, "t") == 0) {
      t = atof(cp);
      if (t <= 0) {
        failed = true;
        break;
      }
    } else if (strcmp(para, "k") == 0) {
      k = atoi(cp);
      if (k <= 0) {
        failed = true;
        break;
      }
    } else if (strcmp(para, "l") == 0) {
      L = atoi(cp);
      if (L <= 0) {
        failed = true;
        break;
      }
    } else if (strcmp(para, "m") == 0) {
      m = atoi(cp);
      if (m <= 0) {
        failed = true;
        break;
      }
    } else if (strcmp(para, "tau") == 0) {
      tau = atof(cp);
      if (tau <= 0.0 || tau >= 1.0) {
        failed = true;
        break;
      }
    } else if (strcmp(para, "c") == 0) {
      ratio = atof(cp);
      if (ratio < 1.0) {
        failed = true;
        break;
      }
    } else if (strcmp(para, "ct") == 0) {
      ct = atoi(cp);
    } else if (strcmp(para, "es") == 0) {
      es = atoi(cp);
      if (es != 0 && es != 1) {
        failed = true;
        break;
      }
    } else {
      printf("Unknown option '-%s'.\n\n", para);
      failed = true;
      break;
    }
    cnt++;
  }

  stTime = clock();

  if (!failed) {
    if (dfname[0] != '\0' && ffolder[0] != '\0' && n != -1 && d != -1 &&
        Bb != -1 && m != -1 && L != -1 && ct == 0) {
      build(dfname, ffolder, n, d, Bb / 4, L, m);
      goto end;
    } else if (dfname[0] != '\0' && ffolder[0] != '\0' && n != -1 && d != -1 &&
               Bb != -1 && m != -1 && L != -1 && ct == 1) {
      ctbuild(dfname, ffolder, n, d, Bb / 4, L, m);
      goto end;
    }
    if (dfname[0] != '\0' && ffolder[0] != '\0' && n != -1 && d != -1 &&
        Bb != -1) {
      buildBin(dfname, ffolder, n, d, Bb / 4);
      goto end;
    } else if (dfname[0] != '\0' && qfname[0] != '\0' && rfname[0] != '\0' &&
               n != -1 && d != -1 && k != -1) {
      seqscan(dfname, n, d, qfname, k, rfname);
      goto end;
    } else if (ffolder[0] != '\0' && rfname[0] != '\0' && qfname[0] != '\0' &&
               ofname[0] != '\0' && k != -1) {
      bool early_stop = true;
      if (es == 0) {
        early_stop = false;
      }
      srsknn(ffolder, nt, rfname, qfname, k, ofname, startID, ratio, tau, t,
             early_stop);
      goto end;
    } else {
      failed = true;
    }
  }

end:
  if (failed) {
    usage();
  } else {
    edTime = clock();

    printf("Running time %.1f seconds\n",
           ((float)edTime - stTime) / CLOCKS_PER_SEC);
  }

  return 0;
}

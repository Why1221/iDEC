#include "idec.h"

// added by long @02132020 for measuring IO time
#ifdef MEASURE_IO_TIME
#include <Timer.hpp>
// too avoid changes in existing APIs/classes, we used gflags
//#include <gflags/gflags.h>
DEFINE_double(G_IO_TIME, 0.0 /* initial value */, "IO time");
DEFINE_int32(G_RIO_CNT, 0 /* initial value */, "IO for raw data");
#endif

#include <string.h>
#include <cmath>
#include "ChiSquared.h"
#include "binheap.h"
#include "gendef.h"
#include "rand.h"

//#include <boost/math/distributions/chi_squared.hpp>



// to handle large files
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#define O_LARGEFILE 0

IDEC::IDEC() {
  d = -1;
  n = -1;
  L = -1;
  B = -1;
  m = -1;
  numUsedTrees = -1;
  trees = NULL;
  a_array = NULL;
}

IDEC::~IDEC() {
  if (a_array) {
    delete[] a_array;
    a_array = NULL;
  }

  int l = L;
  if (numUsedTrees != -1) {
    l = numUsedTrees;
  }

  if (trees) {
    for (int i = 0; i < l; i++) {
      trees[i]->close();
      delete trees[i];
      trees[i] = NULL;
    }

    delete[] trees;
    trees = NULL;
  }
}

void IDEC::init(int _d, int _n, int _B, int _L, int _m) {
  d = _d;
  n = _n;
  B = _B;
  L = _L;
  m = _m;

  gen_vectors();
}

void IDEC::gen_vectors() {
  int i = -1;

  // a_array = new float[m * L * d];
  a_array = new int[m * L * d];
  for (i = 0; i < m * L * d; i++)
    a_array[i] = random_walk();  // generate random variables which equals +1 or
                                 // -1 (with equal prob)
}

/*****************************************************************
builds IDEC tree(s) from a dataset. data format:
id coordinate_1 _2 ... _dim

-para-
dsPath			the dataset path
indexPath	  the folder containing the trees

-return-
0				success
1				failure

-prior-
all tree parameters properly set.
*****************************************************************/

int IDEC::buildFromFile(char *_dsPath, char *_indexPath) {
  int ret = 0;

  char fname[100];
  int cnt = -1;
  int i = -1;
  //   int * key = NULL;
  uint64_t *key = NULL;
  int son = -1;
  FILE *fp = NULL;

  // fp = fopen64(_dsPath, "r");
  fp = fopen64(_dsPath, "rb");
  BlockFile *d_blockFile = NULL;
  char *blk = NULL;
  int blk_pos = -1;
  unsigned enc_dim = d / 64;

  if (!fp) {
    printf("Could not open the source file.\n");
    ret = 1;

    goto recycle;
  }

  fclose(fp);

  strcpy(dsPath, _dsPath);
  getFNameFromPath(dsPath, dsName);

  strcpy(indexPath, _indexPath);
  if (indexPath[strlen(indexPath) - 1] != '/') strcat(indexPath, "/");

  strcpy(dName, indexPath);
  strcat(dName, "data");

  strcpy(fname, indexPath);
  strcat(fname, "para");

  if (writeParaFile(fname)) {
    ret = 1;
    goto recycle;
  }

  key = new uint64_t[enc_dim];
  trees = new RtreePtr[L];

  for (i = 0; i < L; i++) {
    printf("Tree %d (out of %d)\n", i + 1, L);

    trees[i] = new ptRTree();

    getTreeFname(i, fname);
    trees[i]->init(fname, 4 * B, NULL, m);

    cnt = 0;
    // if (isBinary)
    //   fp = fopen64(dsPath, "rb");
    // else
    // fp = fopen64(dsPath, "r");
    fp = fopen64(dsPath, "rb");

    while (!feof(fp) && cnt < n) {
      // freadNextEntry(fp, &son, key);
      son = cnt + 1;  // index starting from 1: not sure why
      memset(key, sizeof(uint64_t), enc_dim);
      if (fread(key, sizeof(uint64_t), enc_dim, fp) != enc_dim) {
        perror("IDEC::buildFromFile -- fread() failed");
      }
      //   if (son <= 0) {
      //     printf("Sorry no negative id please.\n");

      //     ret = 1;
      //     goto recycle;
      //   }
      insert(i, son, key);
      cnt++;

      // if (cnt % 1000 == 0) {
      //   printf("\tInserted %d (%d%%)\n", cnt, cnt * 100 / n);
      // }
    }

    // if (cnt == n && !feof(fp)) {
    uint64_t temp;
    if (cnt == n && (fread(&temp, sizeof(uint64_t), 1, fp) != 0)) {
      printf("The dataset is larger than you said. Giving up...\n", n);
      fclose(fp);
      break;
    }
    fclose(fp);
  }

  for (i = 0; i < L; i++) {
    trees[i]->close();
    delete trees[i];
    trees[i] = NULL;
  }

  delete[] trees;
  trees = NULL;

recycle:
  delete d_blockFile;
  delete[] blk;
  delete[] key;

  return ret;
}

/*****************************************************************
Build a binary file for a dataset. data format:
id coordinate_1 _2 ... _dim

-para-
dsPath      the dataset path
indexPath  the folder containing the binary file (same as the one containing the
trees)

-return-
0       success
1       failure
*****************************************************************/

int IDEC::buildBinFromFile(char *_dsPath, char *_indexPath) {
  int ret = 0;
  char fname[100];
  int cnt = -1;
  int i = -1;
  //   int * key = NULL;
  uint64_t *key = NULL;
  int son = -1;
  FILE *fp = NULL;
  unsigned enc_dim = d / 64;

  //   fp = fopen64(_dsPath, "r");
  fp = fopen64(_dsPath, "rb");
  BlockFile *d_blockFile = NULL;
  char *blk = NULL;
  int blk_pos = -1;

  if (!fp) {
    printf("Could not open the source file.\n");
    ret = 1;

    goto recycle;
  }

  fclose(fp);

  strcpy(dsPath, _dsPath);
  getFNameFromPath(dsPath, dsName);

  strcpy(indexPath, _indexPath);
  if (indexPath[strlen(indexPath) - 1] != '/') strcat(indexPath, "/");

  strcpy(dName, indexPath);
  strcat(dName, "data");

  d_blockFile = new BlockFile(dName, 4 * B);
  blk = new char[d_blockFile->blocklength];
  blk_pos = 0;

  strcpy(fname, indexPath);
  strcat(fname, "para");

  //   key = new int[d];
  key = new uint64_t[d];
  cnt = 0;
  //   fp = fopen64(dsPath, "r");
  fp = fopen64(dsPath, "rb");

  while (!feof(fp) && cnt < n) {
    // freadNextEntry(fp, &son, key);
    // son = cnt;
    memset(key, sizeof(uint64_t), enc_dim);
    if (fread(key, sizeof(uint64_t), enc_dim, fp) != enc_dim) {
      perror("IDEC::buildBinFromFile -- fread() failed");
    }
    // for (int j = 0; j < d; ++j) {
    for (int j = 0; j < enc_dim; ++j) {
      //   memcpy(&blk[blk_pos], &key[j], sizeof (int));
      memcpy(&blk[blk_pos], &key[j], sizeof(uint64_t));
      blk_pos += sizeof(uint64_t);  // sizeof (int);
      if (blk_pos == d_blockFile->blocklength) {
        blk_pos = 0;
        d_blockFile->append_block(blk);
      }
    }

    // if (son <= 0) {
    //   printf("Sorry no negative id please.\n");

    //   ret = 1;
    //   goto recycle;
    // }
    cnt++;
  }
  if (blk_pos != 0) {
    d_blockFile->append_block(blk);  // Thus the last block is not clear!!!
  }

  // if (cnt == n && !feof(fp)) {
  uint64_t temp;
  if (cnt == n && (fread(&temp, sizeof(uint64_t), 1, fp) != 0)) {
    printf("The dataset is larger than you said. Giving up...\n", n);

    ret = 1;
    goto recycle;
  }

  fclose(fp);

recycle:
  delete d_blockFile;
  delete[] blk;

  // it seems this function has mem leak.
  return ret;
}

/*****************************************************************
continue builds IDEC tree(s) from a dataset. data format:
id coordinate_1 _2 ... _dim

-para-
dsPath      the dataset path
indexPath   the folder containing the trees

-return-
0       success
1       failure

-prior-
all tree parameters properly set
*****************************************************************/

int IDEC::ctBuildFromFile(char *_dsPath, char *_indexPath) {
  int ret = 0;

  char fname[100];
  int cnt = -1;
  int i = -1;
  //   int * key = NULL;
  uint64_t *key = NULL;
  // float key = -0.0;
  int son = -1;
  FILE *fp = NULL;

  //   fp = fopen64(_dsPath, "r");
  fp = fopen64(_dsPath, "rb");
  BlockFile *d_blockFile = NULL;
  char *blk = NULL;

  int sz_each_point = d;  // sizeof(uint64_t) * (d / 64)
  int enc_dim = d / 64;
  //   uint64_t* buf;

  if (!fp) {
    printf("Could not open the source file.\n");

    ret = 1;

    goto recycle;
  }

  fclose(fp);

  strcpy(dsPath, _dsPath);
  getFNameFromPath(dsPath, dsName);

  strcpy(indexPath, _indexPath);
  if (indexPath[strlen(indexPath) - 1] != '/') strcat(indexPath, "/");

  strcpy(dName, indexPath);
  strcat(dName, "data");

  strcpy(fname, indexPath);
  strcat(fname, "para");

  key = new uint64_t[d];
  trees = new RtreePtr[L];

  for (i = 0; i < L; i++) {
    printf("Tree %d (out of %d)\n", i + 1, L);

    trees[i] = new ptRTree();

    getTreeFname(i, fname);
    trees[i]->init(fname, NULL);

    cnt = 0;
    fp = fopen64(dsPath, "r");

    while (!feof(fp) && cnt < n) {
      // freadNextEntry(fp, &son, key);
      son = cnt;
      memset(key, sizeof(64), enc_dim);
      if (fread(key, sizeof(uint64_t), enc_dim, fp) != enc_dim) {
        perror("fread() failed");
      }

      //   if (son <= 0) {
      //     printf("Sorry no negative id please.\n");

      //     ret = 1;
      //     goto recycle;
      //   }
      insert(i, son, key);
      cnt++;

      if (cnt % 1000 == 0) {
        printf("\tInserted %d (%d%%)\n", cnt, cnt * 100 / n);
      }
    }

    // if (cnt == n && !feof(fp)) {
    uint64_t temp;
    if (cnt == n && (fread(&temp, sizeof(uint64_t), 1, fp) != 0)) {
      printf("The dataset is larger than you said. Giving up...\n", n);

      ret = 1;
      goto recycle;
    }

    fclose(fp);
  }

  for (i = 0; i < L; i++) {
    trees[i]->close();
    delete trees[i];
    trees[i] = NULL;
  }

  delete[] trees;
  trees = NULL;

recycle:
  delete d_blockFile;
  delete[] blk;

  return ret;
}

/*****************************************************************
write the para file to the disk

-para-
fname		path of the file

-return-
0			success
1			failure
*****************************************************************/

int IDEC::writeParaFile(char *_fname) {
  int ret = 0;

  int i = -1;
  int u = -1;
  FILE *fp = NULL;
  int *aVector = NULL;

  fp = fopen64(_fname, "r");

  if (fp) {
    printf("index exists. ");

    ret = 1;

    goto recycle;
  }

  fp = fopen64(_fname, "w");
  if (!fp) {
    printf("I could not create %s.\n", _fname);
    printf("Perhaps no such folder %s?\n", indexPath);

    ret = 1;

    goto recycle;
  }

  fprintf(fp, "%s\n", dsPath);
  fprintf(fp, "B = %d\n", B);
  fprintf(fp, "n = %d\n", n);
  fprintf(fp, "d = %d\n", d);
  fprintf(fp, "l = %d\n", L);
  fprintf(fp, "m = %d\n", m);

  for (i = 0; i < m * L; i++) {
    getHashPara(i, &aVector);
    // fprintf(fp, "%f", aVector[0]);
    fprintf(fp, "%d", aVector[0]);
    for (u = 1; u < d; u++) {
      // fprintf(fp, " %f", aVector[u]);
      fprintf(fp, " %d", aVector[u]);
    }
    fprintf(fp, "\n");
  }

recycle:
  if (fp) fclose(fp);

  return ret;
}

/*****************************************************************
get the file name of a IDEC-tree

-para-
i		  which tree it is, from 0 to L- 1
fname	(out) the file name
*****************************************************************/

void IDEC::getTreeFname(int _i, char *_fname) {
  char c[100];

  strcpy(_fname, indexPath);
  sprintf(c, "%d", _i);  // modified by Yifang
  strcat(_fname, c);
}

/*****************************************************************
insert a point into the IDEC-trees.

-para-
treeID	into which tree are we inserting (from 0 to L-1)
son		  the object's id
key		  its coordinates

 -return-
0		success
1		failure
*****************************************************************/

int IDEC::insert(int _treeID, int _son, uint64_t *_key) {
  int ret = 0;
  Entry *e = NULL;

  e = new ptEntry();
  e->init(trees[_treeID]);
  e->level = 0;
  e->son = _son;

  for (int i = 0; i < m; ++i) {
    e->bounces[2 * i] = e->bounces[2 * i + 1] =
        getHashValue(_treeID * m + i, _key);
  }
  trees[_treeID]->insert(e);

  return ret;
}

/*****************************************************************
returns the u-th random projection vector.

-para-
u			    see above
a_vector	(out) starting address of the a-vector, which is an array of
size d.
*****************************************************************/

void IDEC::getHashPara(int _u, int **_a_vector) {
  (*_a_vector) = &(a_array[_u * d]);
}

/*****************************************************************
gets the inner product of key and the u-th random vector.

-para-
u		  see above
key		raw coordinates

-return-
see above
*****************************************************************/

// float IDEC::getHashValue(int _u, int *_key) {
//   float ret = 0;
//   float * a_vector = NULL;

//   getHashPara(_u, &a_vector);

//   for (int i = 0; i < d; i++) {
//     ret += a_vector[i] * _key[i];
//   }
//   return ret;
// }
// inline bool isKthBitSet(uint64_t n, unsigned k) {
//     return ((n >> (k - 1)) & 1llu) != 0llu;
// }

int16_t IDEC::getHashValue(int _u, uint64_t *_key) {
  int16_t ret = 0;
  int *a_vector = NULL;

  getHashPara(_u, &a_vector);

  for (int i = 0, j = 0, c = 0; i < d; i++) {
    //
    // Note that the absolute order itself does not matter, as long as
    // all points use the same order.
    // However for debug purpose (i.e., maintaining the same order as un-compact
    // verion), we should pay attention to order issue
    if (isKthBitSet(_key[c], 64 - j)) ret += a_vector[i];
    if (j == 63) {
      ++c;
      j = 0;
    } else {
      ++j;
    }
  }
  return ret;
}

/*****************************************************************
load an existing IDEC index.

-para-
paraPath		    The full path of the parameter file
start_tree_ID   the first tree that we are going to use (default: 0)
numTrees        the number of trees that we are going to use (default: 1)

-return-
0		success
1		failure
*****************************************************************/

int IDEC::restore(char *_paraPath, int _start_tree_ID, int _numTrees) {
  int ret = 0;

  char fname[100];
  int i = -1;
  int len = -1;
  numUsedTrees = _numTrees;

  strcpy(indexPath, _paraPath);

  len = strlen(indexPath);
  if (indexPath[len - 1] != '/' && indexPath[len - 1] != '\\')
    strcat(indexPath, "/");

  strcpy(dName, indexPath);
  strcat(dName, "data");

  strcpy(fname, indexPath);
  strcat(fname, "para");

  if (readParaFile(fname)) {
    ret = 1;
    goto recycle;
  }

  if (_start_tree_ID + _numTrees > L) {
    ret = 1;
    goto recycle;
  }

  trees = new RtreePtr[_numTrees];

  for (i = 0; i < _numTrees; i++) {
    getTreeFname(_start_tree_ID + i, fname);

    trees[i] = new ptRTree();
    trees[i]->init(fname, NULL);
  }

recycle:
  return ret;
}

/*****************************************************************
read the parameter file.

-para-
fname		full path of the para file

-return-
0		success
1		failure
*****************************************************************/

int IDEC::readParaFile(char *_fname) {
  int ret = 0;

  FILE *fp = NULL;
  int cnt = 0;
  int i = -1;
  int j = -1;
  int k = -1;

  fp = fopen64(_fname, "r");

  if (!fp) {
    printf("Could not open %s.\n", _fname);

    ret = 1;
    goto recycle;
  }

  fscanf(fp, "%s\n", dsPath);
  getFNameFromPath(dsPath, dsName);

  fscanf(fp, "B = %d\n", &B);
  fscanf(fp, "n = %d\n", &n);
  fscanf(fp, "d = %d\n", &d);
  // fscanf(fp, "t = %d\n", &t);
  fscanf(fp, "l = %d\n", &L);
  fscanf(fp, "m = %d\n", &m);

  a_array = new int[m * L * d];

  cnt = 0;

  for (i = 0; i < L * m; i++) {
    for (k = 0; k < d; k++) {
      // fscanf(fp, "%f", &a_array[cnt]);
      fscanf(fp, "%d", &a_array[cnt]);
      cnt++;
    }
    fscanf(fp, "\n");
  }

recycle:

  if (fp) fclose(fp);

  return ret;
}

/*****************************************************************
reads the next pt from the original data file.

-para-
fp		handle of the file
son		(out) pt id
key		pt coordinates
*****************************************************************/

void IDEC::freadNextEntry(FILE *_fp, int *_son, int *_key) {
  int i;
  fscanf(_fp, "%d", _son);
  for (i = 0; i < d; i++) {
    fscanf(_fp, " %d", &(_key[i]));
  }

  fscanf(_fp, "\n");
}

/*****************************************************************
reads the pt with given id from the binary data file.

-para-
pt    (out) pt coordinates
bf   block file
id   pt id

-return-
number of I/O that cost in this step
*****************************************************************/

// int IDEC::readData(int *_pt, BlockFile *_bf, int _id) {
// _id starts from 1?
int IDEC::readData(uint64_t *_pt, BlockFile *_bf, int _id) {
  int ret = 0;
  unsigned enc_dim = d / 64;
  //   int blk_pos = ((long long)d * sizeof(int) * (_id - 1)) %
  //   _bf->blocklength; int blk_id = (long long)d * sizeof(int) * (_id - 1) /
  //   _bf->blocklength;

  int blk_pos =
      ((long long)enc_dim * sizeof(uint64_t) * (_id - 1)) % _bf->blocklength;
  int blk_id =
      (long long)enc_dim * sizeof(uint64_t) * (_id - 1) / _bf->blocklength;
  char *blk = new char[_bf->blocklength];

  // printf("%d %d\n",blk_pos,blk_id);

  _bf->read_block(blk, blk_id);
  ret++;  // io count
          //   for (int i = 0; i < d; ++i) {
  for (int i = 0; i < enc_dim; ++i) {
    if (blk_pos == _bf->blocklength) {
      blk_id++;
      _bf->read_block(blk, blk_id);
      ret++;
      blk_pos = 0;
    }
    // memcpy(&_pt[i], &blk[blk_pos], sizeof(int));
    // blk_pos += sizeof(int);
    memcpy(&_pt[i], &blk[blk_pos], sizeof(uint64_t));
    blk_pos += sizeof(uint64_t);
  }

  delete[] blk;
  return ret;
}

/*----------------------------------------------------------------
  auxiliary function called by IDEC:knn.
  ----------------------------------------------------------------*/

int IDEC_hcomp(const void *_e1, const void *_e2) {
  IDEC_Hentry *e1 = (IDEC_Hentry *)_e1;
  IDEC_Hentry *e2 = (IDEC_Hentry *)_e2;

  int ret = 0;

  if (e1->gap < e2->gap)
    ret = -1;
  else if (e1->gap > e2->gap)
    ret = 1;

  return ret;
}

/*----------------------------------------------------------------
  auxiliary function called by IDEC:knn.
  ----------------------------------------------------------------*/

void IDEC_hdestroy(const void *_e) {
  IDEC_Hentry *e = (IDEC_Hentry *)_e;

  delete e;
  e = NULL;
}

/*----------------------------------------------------------------
  auxiliary function called by IDEC:knn.
  ----------------------------------------------------------------*/

float IDEC::updateknn(IDEC_Hentry *_rslt, IDEC_Hentry *_he, int _k) {
  float ret = -1;

  int i = -1;
  int pos = -1;
  bool alreadyIn = false;

  for (i = 0; i < _k; i++) {
    if (_he->id == _rslt[i].id) {
      alreadyIn = true;
      break;
    } else if (compfloats(_he->dist, _rslt[i].dist) == -1) {
      break;
    }
  }

  pos = i;

  if (!alreadyIn && pos < _k) {
    for (i = _k - 1; i > pos; i--) _rslt[i].setto(&(_rslt[i - 1]));

    _rslt[pos].setto(_he);
  }
  ret = _rslt[_k - 1].dist;

  return ret;
}

/*****************************************************************
finds the top-k c-approximate nearest neighbors.

-para-
q             query
k             k of top-k
rslt          (out) top-k cANN
numTrees      number of trees that used (default 1) (only report the case of
"numTrees = 1" in the paper) start_tree_id id of the first tree (default 0) c
approximation ratio that will be used in the early termination condition tau
threshold in Algorithm 2 in the paper t             maximum proportion of points
that will be verified (i.e., T/n) early_stop    whether the early termination
condition is applied or not

-return-
number of I/O
*****************************************************************/
// int IDEC::knn(int *_q, int _k, IDEC_Hentry *_rslt, int _start_tree_ID,
//              int _numTrees, float _c, float _tau, float _t, bool early_stop)
//              {
int IDEC::knn(uint64_t *_q, int _k, IDEC_Hentry *_rslt, int _start_tree_ID,
              int _numTrees, float _c, float _tau, float _t, bool early_stop) {
  int ret = 0;  // total I/O
  int cnt = 0;
  int treeId = 0;
  int son = 0;
  bool again = true;
  // float *qk = NULL;
  int16_t *qk = NULL;
  BinHeap *hp = NULL;
  BinHeapEntry *bhe = NULL;
  BinHeapEntry *bhe_top = NULL;
  IDEC_Hentry *he = NULL;
  BlockFile *d_blockFile = NULL;
  float knnDist = -1;
  // boost::math::chi_squared chi(m);
  float *d_array = NULL;
  int limit = -1;
  const unsigned enc_dim = d / 64;

  IDEC_Hentry *inner_rslt = new IDEC_Hentry[_k];
  for (int i = 0; i < _k; i++) {
    inner_rslt[i].d = d;
  }

  // inner_rslt is the top-k NN found through R-tree; rslt also consider the
  // points that in the same page when reading a point from the binary file
  for (int i = 0; i < _k; i++) {
    _rslt[i].id = -1;
    _rslt[i].dist = (float)MAXREAL;
    inner_rslt[i].id = -1;
    inner_rslt[i].dist = (float)MAXREAL;
  }

  hp = new BinHeap();
  hp->compare_func = &IDEC_hcomp;
  hp->destroy_data = &IDEC_hdestroy;

  d_blockFile = new BlockFile(dName, 4 * B);
  knnDist = (float)MAXREAL;

  limit = (int)ceil(_t * n) + _k - 1;
  if (limit > n) limit = n;

  qk = new int16_t[m * _numTrees];

  // d_array records all the read points, so if the r-tree returns a point that
  // is already read, we can save a I/O
  //
  d_array = new float[n + 1];
  for (int i = 0; i <= n; ++i) {
    d_array[i] = -1.0;
  }

  for (int i = 0; i < _numTrees; ++i) {
    he = new IDEC_Hentry();
    he->id = trees[i]->root;
    he->gap = 0;
    he->level = 1;
    he->d = d;
    he->treeId = i;
    he->dist = -1;

    bhe = new BinHeapEntry();
    bhe->data = he;
    hp->insert(bhe);

    for (int j = 0; j < m; ++j) {
      qk[i * m + j] = getHashValue((_start_tree_ID + i) * m + j, _q);
#ifdef DEBUG
      printf("%.6f ", qk[i * m + j]);
#endif
    }
#ifdef DEBUG
    printf("\n ");
#endif
  }

#ifdef MEASURE_IO_TIME
  HighResolutionTimer io_timer;
  FLAGS_G_IO_TIME = 0.0;
  FLAGS_G_RIO_CNT = 0;
#endif
  while (again) {
    bhe_top = hp->remove();
    if (!bhe_top) {
      again = false;
    } else {
      he = (IDEC_Hentry *)bhe_top->data;
      treeId = he->treeId;
      son = he->id;

      if (he->level == 0) {
        cnt++;

        if (d_array[son] < 0) {
          //   int *pt = new int[d];
          uint64_t *pt = new uint64_t[enc_dim];
          int ret_temp = 0;
          // locate the page that contains point son
          //   int blk_pos = ((long long)d * sizeof(int) * (son - 1)) %
          //                 d_blockFile->blocklength;
          //   int blk_id =
          //       (long long)d * sizeof(int) * (son - 1) /
          //       d_blockFile->blocklength;
          ///// NOTE THAT son starts from 1
          int blk_pos = ((long long)enc_dim * sizeof(uint64_t) * (son - 1)) %
                        d_blockFile->blocklength;
          int blk_id = (long long)enc_dim * sizeof(uint64_t) * (son - 1) /
                       d_blockFile->blocklength;
          char *blk = new char[d_blockFile->blocklength];

#ifdef MEASURE_IO_TIME
          io_timer.restart();
#endif
          d_blockFile->read_block(blk, blk_id);
#ifdef MEASURE_IO_TIME
          FLAGS_G_IO_TIME += io_timer.elapsed();
          ++FLAGS_G_RIO_CNT;
#endif
          ret++;

          // read points in the page before son
          for (int i = 1;; i++) {
            // int tmp_pos = blk_pos - d * sizeof(int) * i;
            int tmp_pos = blk_pos - enc_dim * sizeof(uint64_t) * i;
            if (tmp_pos < 0 || son - i < 1) {
              break;
            }

            // for (int j = 0; j < d; ++j) {
            //   memcpy(&pt[j], &blk[tmp_pos], sizeof(int));
            //   tmp_pos += sizeof(int);
            // }

            for (int j = 0; j < enc_dim; ++j) {
              memcpy(&pt[j], &blk[tmp_pos], sizeof(uint64_t));
              tmp_pos += sizeof(uint64_t);
            }

            d_array[son - i] =
                l2_dist_int(pt, _q, enc_dim);  // l2_dist_int(pt, _q, d);

#ifndef DEBUG
            /* We comment out this block for debugging purpose. As the size of
               each point is different in compact mode and un-compact mode, the
               results for them under the same t is different (more data points
               would be checked under compact mode==>compact mode would be
               better). To make sure they will produce the same result when
               using the same parameter, we comment out this block. */
            IDEC_Hentry *tmp_e = new IDEC_Hentry();
            tmp_e->id = son - i;
            tmp_e->dist = d_array[son - i];
            updateknn(_rslt, tmp_e, _k);
            delete tmp_e;
#endif
          }

          // read point son
          //   for (int i = 0; i < d; ++i) {
          for (int i = 0; i < enc_dim; ++i) {
            if (blk_pos == d_blockFile->blocklength) {
              blk_id++;
#ifdef MEASURE_IO_TIME
              io_timer.restart();
#endif
              d_blockFile->read_block(blk, blk_id);
#ifdef MEASURE_IO_TIME
              FLAGS_G_IO_TIME += io_timer.elapsed();
              ++FLAGS_G_RIO_CNT;
#endif
              ret++;
              blk_pos = 0;
            }
            // memcpy(&pt[i], &blk[blk_pos], sizeof(int));
            // blk_pos += sizeof(int);
            memcpy(&pt[i], &blk[blk_pos], sizeof(uint64_t));
            blk_pos += sizeof(uint64_t);
          }
          d_array[son] = l2_dist_int(pt, _q, enc_dim);  // l2_dist_int(pt, _q,
                                                        // d);

          // read point in the page after son
          for (int i = 1;; ++i) {
            // int tmp_pos = blk_pos + d * sizeof(int) * (i - 1);
            int tmp_pos = blk_pos + enc_dim * sizeof(uint64_t) * (i - 1);
            // if ((tmp_pos + d * sizeof(int)) > d_blockFile->blocklength ||
            //     son + i > n) {
            if ((tmp_pos + enc_dim * sizeof(uint64_t)) >
                    d_blockFile->blocklength ||
                son + i > n) {
              break;
            }

            // for (int j = 0; j < d; ++j) {
            for (int j = 0; j < enc_dim; ++j) {
              //   memcpy(&pt[j], &blk[tmp_pos], sizeof(int));
              //   tmp_pos += sizeof(int);
              memcpy(&pt[j], &blk[tmp_pos], sizeof(uint64_t));
              tmp_pos += sizeof(uint64_t);
            }
            d_array[son + i] =
                l2_dist_int(pt, _q, enc_dim);  // l2_dist_int(pt, _q, d);
            if (son + i > n) {
              printf("%d %d\n", son, i);
            }

#ifndef DEBUG
            /* We comment out this block for debugging purpose. As the size of
               each point is different in compact mode and un-compact mode, the
               results for them under the same t is different (more data points
               would be checked under compact mode==>compact mode would be
               better). To make sure they will produce the same result when
               using the same parameter, we comment out this block. */
            IDEC_Hentry *tmp_e = new IDEC_Hentry();
            tmp_e->id = son + i;
            tmp_e->dist = d_array[son + i];
            updateknn(_rslt, tmp_e, _k);
            delete tmp_e;
#endif
          }
          delete[] pt;
          delete[] blk;
        }

        he->dist = d_array[son];
        updateknn(_rslt, he, _k);
        knnDist = updateknn(inner_rslt, he, _k);

        // early termination condition
        //        if (early_stop && (compfloats(knnDist, 0.0) == 0 || 1 - pow(1
        //        - boost::math::cdf(chi, _c * _c * he->gap / knnDist /
        //        knnDist), (double) _numTrees) >= _tau)) {
        if (early_stop &&
            (compfloats(knnDist, 0.0) == 0 ||
             (knnDist * knnDist) > (he->gap * chi2inv(_tau, m) / _c / _c))) {
          again = false;
        }

      } else {
        ret++;
        RTNode *child = new ptRTNode();

#ifdef MEASURE_IO_TIME
        io_timer.restart();
#endif
        child->init(trees[treeId], son);
#ifdef MEASURE_IO_TIME
        FLAGS_G_IO_TIME += io_timer.elapsed();
#endif
        for (int i = 0; i < child->num_entries; i++) {
          he = new IDEC_Hentry();
          he->id = child->entries[i]->son;
          he->level = child->level;
          he->gap = MINDIST(qk, child->entries[i]->bounces, m, treeId);
          he->treeId = treeId;
          he->d = d;

          bhe = new BinHeapEntry();
          bhe->data = he;
          hp->insert(bhe);
        }
        delete child;
      }
      // normal termination condition (i.e., verified T points)
      if (cnt >= limit) {
        again = false;
      }
    }
  }

recycle:
  if (d_blockFile) {
    delete d_blockFile;
  }

  if (qk) {
    delete[] qk;
  }

  if (d_array) {
    delete[] d_array;
  }

  delete[] inner_rslt;

  if (hp->root) hp->root->recursive_data_wipeout(hp->destroy_data);
  delete hp;

  return ret;
}

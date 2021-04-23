#include "gendef.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//////////////////////////////////////////////////////////////////////////////
// globals
//////////////////////////////////////////////////////////////////////////////

#ifdef UNIX
void strupr(char *_msg) {
  int dist = 'A' - 'a';
  char *c_ptr = _msg;

  while (*c_ptr) {
    if (*c_ptr >= 'a') *c_ptr += dist;
    c_ptr++;
  }
}
#endif

void error(char *t, bool ex) {
  fprintf(stderr, t);
  if (ex) exit(0);
}

float area(int dimension, float *mbr) {
  int i;
  float sum;

  sum = 1.0;
  for (i = 0; i < dimension; i++) sum *= mbr[2 * i + 1] - mbr[2 * i];

  return sum;
}

float margin(int dimension, float *mbr) {
  float *ml, *mu, *m_last, sum;

  sum = 0.0;
  m_last = mbr + 2 * dimension;
  ml = mbr;
  mu = ml + 1;
  while (mu < m_last) {
    sum += *mu - *ml;
    ml += 2;
    mu += 2;
  }

  return sum;
}

bool inside(float &p, float &lb, float &ub) { return (p >= lb && p <= ub); }

bool check_point_in_rec(float *v, float *mbr, int dimension) {
  int i;

  for (i = 0; i < dimension; i++)
    if (!inside(v[i], mbr[2 * i], mbr[2 * i + 1])) return FALSE;

  return TRUE;
}

// calculates the overlapping rectangle between r1 and r2
// if they do not overlap, returns null
float *overlapRect(int dimension, float *r1, float *r2) {
  float *overlap = new float[2 * dimension];
  for (int i = 0; i < dimension; i++) {
    if ((r1[i * 2] > r2[i * 2 + 1]) ||
        (r1[i * 2 + 1] < r2[i * 2]))  // non overlapping
    {
      delete[] overlap;
      return NULL;
    }
    overlap[2 * i] = max(r1[i * 2], r2[i * 2]);
    overlap[2 * i + 1] = min(r1[i * 2 + 1], r2[i * 2 + 1]);
  }

  return overlap;
}

float overlap(int dimension, float *r1, float *r2)
// calculates the overlapping area of r1 and r2
{
  float sum;
  float *r1pos, *r2pos, *r1last, r1_lb, r1_ub, r2_lb, r2_ub;

  sum = 1.0;
  r1pos = r1;
  r2pos = r2;
  r1last = r1 + 2 * dimension;

  while (r1pos < r1last) {
    r1_lb = *(r1pos++);
    r1_ub = *(r1pos++);
    r2_lb = *(r2pos++);
    r2_ub = *(r2pos++);

    // calculate overlap in this dimension

    if (inside(r1_ub, r2_lb, r2_ub))
    // upper bound of r1 is inside r2
    {
      if (inside(r1_lb, r2_lb, r2_ub))
        // and lower bound of r1 is inside
        sum *= (r1_ub - r1_lb);
      else
        sum *= (r1_ub - r2_lb);
    } else {
      if (inside(r1_lb, r2_lb, r2_ub))
        // and lower bound of r1 is inside
        sum *= (r2_ub - r1_lb);
      else {
        if (inside(r2_lb, r1_lb, r1_ub) && inside(r2_ub, r1_lb, r1_ub))
          // r1 contains r2
          sum *= (r2_ub - r2_lb);
        else
          // r1 and r2 do not overlap
          sum = 0.0;
      }
    }
  }

  return sum;
}

void enlarge(int dimension, float **mbr, float *r1, float *r2)
// enlarge r to tightly contain s
{
  int i;

  *mbr = new float[2 * dimension];
  for (i = 0; i < 2 * dimension; i += 2) {
    (*mbr)[i] = min(r1[i], r2[i]);

    (*mbr)[i + 1] = max(r1[i + 1], r2[i + 1]);
  }
}

bool section(int dimension, float *mbr1, float *mbr2) {
  int i;

  for (i = 0; i < dimension; i++) {
    if (mbr1[2 * i] > mbr2[2 * i + 1] || mbr1[2 * i + 1] < mbr2[2 * i])
      return FALSE;
  }
  return TRUE;
}

int sort_lower_mbr(const void *d1, const void *d2) {
  SortMbr *s1, *s2;
  float erg;
  int dimension;

  s1 = (SortMbr *)d1;
  s2 = (SortMbr *)d2;
  dimension = s1->dimension;
  erg = s1->mbr[2 * dimension] - s2->mbr[2 * dimension];
  if (erg < 0.0)
    return -1;
  else if (erg == 0.0)
    return 0;
  else
    return 1;
}

int sort_upper_mbr(const void *d1, const void *d2) {
  SortMbr *s1, *s2;
  float erg;
  int dimension;

  s1 = (SortMbr *)d1;
  s2 = (SortMbr *)d2;
  dimension = s1->dimension;
  erg = s1->mbr[2 * dimension + 1] - s2->mbr[2 * dimension + 1];
  if (erg < 0.0)
    return -1;
  else if (erg == 0.0)
    return 0;
  else
    return 1;
}

int sort_center_mbr(const void *d1, const void *d2) {
  SortMbr *s1, *s2;
  int i, dimension;
  float d, e1, e2;

  s1 = (SortMbr *)d1;
  s2 = (SortMbr *)d2;
  dimension = s1->dimension;

  e1 = e2 = 0.0;
  for (i = 0; i < dimension; i++) {
    d = ((s1->mbr[2 * i] + s1->mbr[2 * i + 1]) / 2.0) - s1->center[i];
    e1 += d * d;
    d = ((s2->mbr[2 * i] + s2->mbr[2 * i + 1]) / 2.0) - s2->center[i];
    e2 += d * d;
  }

  if (e1 < e2)
    return -1;
  else if (e1 == e2)
    return 0;
  else
    return 1;
}

float ptDIST(float *p1, float *p2, int _dim) {
  float summe = 0;
  int i;

  for (i = 0; i < _dim; i++) summe += pow(p1[i] - p2[i], 2);

  return (sqrt(summe));
  // return(summe);
}

/*****************************************************************
this function returns the maxdist of 2 mbrs
para:
m1: the bounces of the 1st mbr
m2: the bounces of the 2nd mbr
dim: dimensionality
*****************************************************************/

float MbrMAXDIST(float *_m1, float *_m2, int _dim) {
  float dist = 0;
  for (int i = 0; i < _dim; i++) {
    float d1 = fabs(_m1[2 * i] - _m2[2 * i + 1]);
    float d2 = fabs(_m1[2 * i + 1] - _m2[2 * i]);
    float d = max(d1, d2);
    dist += pow(d, 2);
  }
  return dist;
}

/*****************************************************************
this function returns the mindist of 2 mbrs
para:
m1: the bounces of the 1st mbr
m2: the bounces of the 2nd mbr
dim: dimensionality
*****************************************************************/

float MbrMINDIST(float *_m1, float *_m2, int _dim) {
  float dist = 0;
  for (int i = 0; i < _dim; i++) {
    if (_m1[2 * i] > _m2[2 * i + 1])
      dist += pow(_m1[2 * i] - _m2[2 * i + 1], 2);
    else if (_m1[2 * i + 1] < _m2[2 * i])
      dist += pow(_m1[2 * i + 1] - _m2[2 * i], 2);
  }
  return dist;
}

float MINDIST(float *p, float *bounces, int _dim, int _start_id) {
  float summe = 0.0;
  float r;
  int i;

  for (i = 0; i < _dim; i++) {
    if (p[i + _start_id * _dim] < bounces[2 * i])
      r = bounces[2 * i];
    else {
      if (p[i + _start_id * _dim] > bounces[2 * i + 1])
        r = bounces[2 * i + 1];
      else
        r = p[i + _start_id * _dim];
    }

    summe += pow(p[i + _start_id * _dim] - r, 2);
  }
  return (summe);
}

float MINDIST(float *p, float *bounces, int _dim) {
  //
  // Berechne die kuerzeste Entfernung zwischen einem Punkt Point
  // und einem MBR bounces (Lotrecht!)
  //

  float summe = 0.0;
  float r;
  int i;

  for (i = 0; i < _dim; i++) {
    if (p[i] < bounces[2 * i])
      r = bounces[2 * i];
    else {
      if (p[i] > bounces[2 * i + 1])
        r = bounces[2 * i + 1];
      else
        r = p[i];
    }

    summe += pow(p[i] - r, 2);
  }
  return (summe);
}

float MAXDIST(float *p, float *bounces, int dim) {
  float summe = 0.0;
  float r;
  int i;

  for (i = 0; i < dim; i++) {
    if (p[i] < bounces[2 * i])
      r = bounces[2 * i + 1];
    else {
      if (p[i] > bounces[2 * i + 1])
        r = bounces[2 * i];
      else if (p[i] - bounces[2 * i] > bounces[2 * i + 1] - p[i])
        r = bounces[2 * i];
      else
        r = bounces[2 * i + 1];
    }

    summe += pow(p[i] - r, 2);
  }

  return (summe);
}

float MINMAXDIST(float *p, float *bounces, int _dim) {
  // Berechne den kleinsten maximalen Abstand von einem Punkt Point
  // zu einem MBR bounces.
  // Wird benutzt zur Abschaetzung von Abstaenden bei NearestNarborQuery.
  // Kann als Supremum fuer die aktuell kuerzeste Distanz:
  // Alle Punkte mit einem Abstand > MINMAXDIST sind keine Kandidaten mehr
  // fuer den NearestNarbor
  // vgl. Literatur:
  // Nearest Narbor Query v. Roussopoulos, Kelley und Vincent,
  // University of Maryland

  float summe = 0;
  float minimum = 1.0e20;
  float S = 0;

  float rmk, rMi;
  int k, i;

  for (i = 0; i < _dim; i++) {
    rMi = (p[i] >= (bounces[2 * i] + bounces[2 * i + 1]) / 2)
              ? bounces[2 * i]
              : bounces[2 * i + 1];
    S += pow(p[i] - rMi, 2);
  }

  for (k = 0; k < _dim; k++) {
    rmk = (p[k] <= (bounces[2 * k] + bounces[2 * k + 1]) / 2)
              ? bounces[2 * k]
              : bounces[2 * k + 1];

    summe = pow(p[k] - rmk, 2);

    rMi = (p[k] >= (bounces[2 * k] + bounces[2 * k + 1]) / 2)
              ? bounces[2 * k]
              : bounces[2 * k + 1];

    summe += S - pow(p[k] - rMi, 2);

    minimum = min(minimum, summe);
  }

  return (minimum);
}

/*****************************************************************
compares two float numbers

para:
- v1:
- v2:

return:
- -1: v1 smaller
  0: tie
  1: v2 smaller

Coded by Yufei Tao, 31 july 08
 *****************************************************************/

int compfloats(float _v1, float _v2) {
  int ret = 0;

  if (_v1 - _v2 < -FLOATZERO)
    ret = -1;
  else if (_v1 - _v2 > FLOATZERO)
    ret = 1;
  else
    ret = 0;

  return ret;
}

/*****************************************************************
checks if v is a power of 2

para:
- v

return:
- true or false

Coded by Yufei Tao, 31 july 08
 *****************************************************************/

bool is_pow_of_2(int _v) {
  int x = (int)(log((float)_v) / log(2.0));
  int y = (int)pow(2, x);

  return (_v == y);
}

/*****************************************************************
get the part of a path after the last blackslash.

e.g., given path = "./ex/ex/ex.h", return "ex.h"

para
- path:			the given path.
- (out) fname:	the returned part (usually a file name)

Coded by Yufei Tao, 4 aug 08
 *****************************************************************/

void getFNameFromPath(char *_path, char *_fname) {
  int i;
  int len = strlen(_path);
  int pos = -1;

  for (i = len - 1; i >= 0; i--) {
    if (_path[i] == '/') {
      pos = i;
      break;
    }
  }

  pos++;

  for (i = pos; i <= len; i++) {
    _fname[i - pos] = _path[i];
  }

  _fname[i - pos] = '\0';
}

/*****************************************************************
this function gets the part of the given path up to the last folder.
e.g, given ./ex/ex/1.zip, the function returns ./ex/ex/

para:
- path
- (out) folder:

Coded by Yufei Tao, 7 aug 08
 *****************************************************************/

void get_leading_folder(char *_path, char *_folder) {
  int len = strlen(_path);
  int pos = -1;

  for (int i = len - 1; i >= 0; i--) {
    if (_path[i] == '/') {
      pos = i;
      break;
    }
  }

  int i = 0;             // modified by Yifang
  for (; i <= pos; i++)  // modified by Yifang
  {
    _folder[i] = _path[i];
  }

  _folder[i] = '\0';
}

/*****************************************************************
just gets n tabs on the screen.

para:
- n

Coded by Yufei Tao, 7 aug 08
 *****************************************************************/

void blank_print(int _n) {
  for (int i = 0; i < _n; i++) printf("   ");
}

/*****************************************************************
Coded by Yufei Tao, 7 aug 08
 *****************************************************************/

float l2_dist_int(int *_p1, int *_p2, int _dim) {
  //    long long p = 0;
  //    long long q = 0;
  //    long long pq = 0;
  //    for (int i = 0; i < _dim; ++i){
  //        p += (long)_p1[i] * (long)_p1[i];
  //    }
  //    for (int i = 0; i < _dim; ++i){
  //        q += (long)_p2[i] * (long)_p2[i];
  //    }
  //    for (int i = 0; i < _dim; ++i){
  //        pq += (long)_p2[i] * (long)_p1[i];
  //    }
  //    if ((double)pq / sqrt((double)p) / sqrt((double)q) <= 0.95)
  //        printf("%f\t", (double)pq / sqrt((double)p) / sqrt((double)q));
  float ret = 0;
  for (int i = 0; i < _dim; i++) {
    float dif = (float)(_p1[i] - _p2[i]);
    ret += dif * dif;
    // tt += (long)(_p1[i] - _p2[i]) * (long)(_p1[i] - _p2[i]);
  }
  //    printf("%lld\t%f\n", t, ret);
  /// x += tt;
  // ret = sqrt(ret);
  // printf("%f\t", (double)x / (2*(double)t));
  return sqrt(ret); // ret
}

/*----------------------------------------------------------------
DESCRIPTION:
Get the next word from a string.

PARA:

_s (in/out):		the string. At finish, _s will point to the
                    char right after the word returned.
_w (out):			the word returned.


RETURN:
the character.

KNOWN ISSUES:
None.

AUTHOR:
Yufei Tao

LAST MODIFIED:
5 Mar. 2009.
----------------------------------------------------------------*/

void getnextWord(char **_s, char *_w) {
  while (**_s == ' ') {
    (*_s)++;
  }

  while (**_s != ' ' && **_s != '\0') {
    *_w = **_s;
    (*_s)++;
    _w++;
  }

  *_w = '\0';
}

/*----------------------------------------------------------------
DESCRIPTION:
Get the next non-white character from a string.

PARA:
_s (in/out):		the string. when the function finishes, _s
                                        will point to the character after the
one returned.

RETURN:
the character.

KNOWN ISSUES:
None.

AUTHOR:
Yufei Tao

LAST MODIFIED:
5 Mar. 2009.
----------------------------------------------------------------*/

char getnextChar(char **_s) {
  char ret;

  while (**_s == ' ') {
    (*_s)++;
  }

  ret = **_s;
  (*_s)++;

  return ret;
}

/*****************************************************************
print the last m bits of an integer n

-para-
n
m

-by-
yf tao

-last touch-
19 aug 09
 *****************************************************************/

void printINT_in_BIN(int _n, int _m) {
  int i = -1;
  int mask = -1;
  int n = -1;

  mask = 1 << _m;
  mask--;

  n = _n & mask;

  mask = 1 << (_m - 1);

  for (i = 0; i < _m; i++) {
    if (mask & n)
      printf("1");
    else
      printf("0");

    mask >>= 1;
  }
  printf("\n");
}

bool isKthBitSet(uint64_t n, unsigned k) {
  return ((n >> (k - 1)) & 1llu) != 0llu;
}

// Function to calculate hamming distance
inline int hammingDistance(uint64_t n1, uint64_t n2) {
  uint64_t x = n1 ^ n2;
  int setBits = 0;

  while (x > 0) {
    setBits += x & 1;
    x >>= 1;
  }

  return setBits;
}

float l2_dist_int(uint64_t *_p1, uint64_t *_p2, int _dim) {
  float dist = 0.0f;
  for (int i = 0; i < _dim; ++i) {
#if defined(__GNUC__) || defined(__GNUG__)
    // use popcnt to speedup
    dist += __builtin_popcountll((_p1[i] ^ _p2[i]));
#else
    dist += hammingDistance(_p1[i], _p2[i]);
#endif
  }
  return sqrt(dist); // dist
}

float hammingDistance(uint64_t *_p1, uint64_t *_p2, int _dim) {
  float dist = 0.0f;
  for (int i = 0; i < _dim; ++i) {
    dist += hammingDistance(_p1[i], _p2[i]);
  }
  return dist;
}

float area(int dimension, int16_t *mbr) {
  int i;
  float sum;

  sum = 1.0;
  for (i = 0; i < dimension; i++) sum *= mbr[2 * i + 1] - mbr[2 * i];

  return sum;
}

float margin(int dimension, int16_t *mbr) {
  // float *ml, *mu, *m_last, sum;

  int16_t *ml, *mu, *m_last;
  float sum;
  sum = 0.0;
  m_last = mbr + 2 * dimension;
  ml = mbr;
  mu = ml + 1;
  while (mu < m_last) {
    sum += *mu - *ml;
    ml += 2;
    mu += 2;
  }

  return sum;
}

bool inside(int16_t &p, int16_t &lb, int16_t &ub) {
  return (p >= lb && p <= ub);
}

bool check_point_in_rec(int16_t *v, int16_t *mbr, int dimension) {
  int i;

  for (i = 0; i < dimension; i++)
    if (!inside(v[i], mbr[2 * i], mbr[2 * i + 1])) return FALSE;

  return TRUE;
}

int16_t *overlapRect(int dimension, int16_t *r1, int16_t *r2) {
  // float *overlap = new float[2 * dimension];
  int16_t *overlap = new int16_t[2 * dimension];
  for (int i = 0; i < dimension; i++) {
    if ((r1[i * 2] > r2[i * 2 + 1]) ||
        (r1[i * 2 + 1] < r2[i * 2]))  // non overlapping
    {
      delete[] overlap;
      return NULL;
    }
    overlap[2 * i] = max(r1[i * 2], r2[i * 2]);
    overlap[2 * i + 1] = min(r1[i * 2 + 1], r2[i * 2 + 1]);
  }

  return overlap;
}

float overlap(int dimension, int16_t *r1, int16_t *r2)
// calculates the overlapping area of r1 and r2
{
  float sum;
  int16_t *r1pos, *r2pos, *r1last, r1_lb, r1_ub, r2_lb, r2_ub;

  sum = 1.0;
  r1pos = r1;
  r2pos = r2;
  r1last = r1 + 2 * dimension;

  while (r1pos < r1last) {
    r1_lb = *(r1pos++);
    r1_ub = *(r1pos++);
    r2_lb = *(r2pos++);
    r2_ub = *(r2pos++);

    // calculate overlap in this dimension

    if (inside(r1_ub, r2_lb, r2_ub))
    // upper bound of r1 is inside r2
    {
      if (inside(r1_lb, r2_lb, r2_ub))
        // and lower bound of r1 is inside
        sum *= (r1_ub - r1_lb);
      else
        sum *= (r1_ub - r2_lb);
    } else {
      if (inside(r1_lb, r2_lb, r2_ub))
        // and lower bound of r1 is inside
        sum *= (r2_ub - r1_lb);
      else {
        if (inside(r2_lb, r1_lb, r1_ub) && inside(r2_ub, r1_lb, r1_ub))
          // r1 contains r2
          sum *= (r2_ub - r2_lb);
        else
          // r1 and r2 do not overlap
          sum = 0.0;
      }
    }
  }

  return sum;
}

void enlarge(int dimension, int16_t **mbr, int16_t *r1, int16_t *r2)
// enlarge r to tightly contain s
{
  int i;

  // *mbr = new float[2 * dimension];
  *mbr = new int16_t[2 * dimension];
  for (i = 0; i < 2 * dimension; i += 2) {
    (*mbr)[i] = min(r1[i], r2[i]);

    (*mbr)[i + 1] = max(r1[i + 1], r2[i + 1]);
  }
}

bool section(int dimension, int16_t *mbr1, int16_t *mbr2) {
  int i;

  for (i = 0; i < dimension; i++) {
    if (mbr1[2 * i] > mbr2[2 * i + 1] || mbr1[2 * i + 1] < mbr2[2 * i])
      return FALSE;
  }
  return TRUE;
}

float ptDIST(int16_t *p1, int16_t *p2, int _dim) {
  float summe = 0;
  int i;

  for (i = 0; i < _dim; i++) summe += pow(p1[i] - p2[i], 2);

  return (sqrt(summe));
  // return (summe);
}

float MbrMAXDIST(int16_t *_m1, int16_t *_m2, int _dim) {
  float dist = 0;
  for (int i = 0; i < _dim; i++) {
    float d1 = fabs(_m1[2 * i] - _m2[2 * i + 1]);
    float d2 = fabs(_m1[2 * i + 1] - _m2[2 * i]);
    float d = max(d1, d2);
    dist += pow(d, 2);
  }
  return dist;
}

float MbrMINDIST(int16_t *_m1, int16_t *_m2, int _dim) {
  float dist = 0;
  for (int i = 0; i < _dim; i++) {
    if (_m1[2 * i] > _m2[2 * i + 1])
      dist += pow(_m1[2 * i] - _m2[2 * i + 1], 2);
    else if (_m1[2 * i + 1] < _m2[2 * i])
      dist += pow(_m1[2 * i + 1] - _m2[2 * i], 2);
  }
  return dist;
}

float MINDIST(int16_t *p, int16_t *bounces, int _dim, int _start_id) {
  float summe = 0.0;
  float r;
  int i;

  for (i = 0; i < _dim; i++) {
    if (p[i + _start_id * _dim] < bounces[2 * i])
      r = bounces[2 * i];
    else {
      if (p[i + _start_id * _dim] > bounces[2 * i + 1])
        r = bounces[2 * i + 1];
      else
        r = p[i + _start_id * _dim];
    }

    summe += pow(p[i + _start_id * _dim] - r, 2);
  }
  return (summe);
}

float MINDIST(int16_t *p, int16_t *bounces, int _dim) {
  //
  // Berechne die kuerzeste Entfernung zwischen einem Punkt Point
  // und einem MBR bounces (Lotrecht!)
  //

  float summe = 0.0;
  float r;
  int i;

  for (i = 0; i < _dim; i++) {
    if (p[i] < bounces[2 * i])
      r = bounces[2 * i];
    else {
      if (p[i] > bounces[2 * i + 1])
        r = bounces[2 * i + 1];
      else
        r = p[i];
    }

    summe += pow(p[i] - r, 2);
  }
  return (summe);
}

float MAXDIST(int16_t *p, int16_t *bounces, int dim) {
  float summe = 0.0;
  float r;
  int i;

  for (i = 0; i < dim; i++) {
    if (p[i] < bounces[2 * i])
      r = bounces[2 * i + 1];
    else {
      if (p[i] > bounces[2 * i + 1])
        r = bounces[2 * i];
      else if (p[i] - bounces[2 * i] > bounces[2 * i + 1] - p[i])
        r = bounces[2 * i];
      else
        r = bounces[2 * i + 1];
    }

    summe += pow(p[i] - r, 2);
  }

  return (summe);
}

float MINMAXDIST(int16_t *p, int16_t *bounces, int _dim) {
  // Berechne den kleinsten maximalen Abstand von einem Punkt Point
  // zu einem MBR bounces.
  // Wird benutzt zur Abschaetzung von Abstaenden bei NearestNarborQuery.
  // Kann als Supremum fuer die aktuell kuerzeste Distanz:
  // Alle Punkte mit einem Abstand > MINMAXDIST sind keine Kandidaten mehr
  // fuer den NearestNarbor
  // vgl. Literatur:
  // Nearest Narbor Query v. Roussopoulos, Kelley und Vincent,
  // University of Maryland

  float summe = 0;
  float minimum = 1.0e20;
  float S = 0;

  float rmk, rMi;
  int k, i;

  for (i = 0; i < _dim; i++) {
    rMi = (p[i] >= (bounces[2 * i] + bounces[2 * i + 1]) / 2)
              ? bounces[2 * i]
              : bounces[2 * i + 1];
    S += pow(p[i] - rMi, 2);
  }

  for (k = 0; k < _dim; k++) {
    rmk = (p[k] <= (bounces[2 * k] + bounces[2 * k + 1]) / 2)
              ? bounces[2 * k]
              : bounces[2 * k + 1];

    summe = pow(p[k] - rmk, 2);

    rMi = (p[k] >= (bounces[2 * k] + bounces[2 * k + 1]) / 2)
              ? bounces[2 * k]
              : bounces[2 * k + 1];

    summe += S - pow(p[k] - rMi, 2);

    minimum = min(minimum, summe);
  }

  return (minimum);
}

int compfloats(int16_t _v1, int16_t _v2) {
  int ret = 0;

  if (_v1 - _v2 < 0)
    ret = -1;
  else if (_v1 - _v2 > 0)
    ret = 1;
  else
    ret = 0;

  return ret;
}

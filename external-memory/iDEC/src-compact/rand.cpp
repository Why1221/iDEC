//--------------------------------------------
// this file contains necessary random number
// generators
// 
// collected by Yufei Tao
//--------------------------------------------
#include <math.h>
#include <stdlib.h>
#include "rand.h"

/************************************************************
 ***  Given a mean and a standard deviation, gaussian       **
 **   generates a normally distributed random number        **
 **   Algorithm:  Polar Method, p.  104, Knuth, vol. 2      **
 ************************************************************/

float gaussian(float mean, float sigma) {
    float v1, v2;
    float s;
    float x;

    do {
        v1 = 2 * uniform(0, 1) - 1;
        v2 = 2 * uniform(0, 1) - 1;
        s = v1 * v1 + v2*v2;
    } while (s >= 1.);

    x = v1 * (float) sqrt(-2. * log(s) / s);

    /*  x is normally distributed with mean 0 and sigma 1.  */
    x = x * sigma + mean;

    return (x);
}

/************************************************************
 ** Generates a random number between _min and _max         **
 ** uniformly                                               **
   By Yufei Tao
 ************************************************************/

float uniform(float _min, float _max) {
    int int_r = rand();
    long base = RAND_MAX - 1;
    float f_r = ((float) int_r) / base;
    return (_max - _min) * f_r + _min;
}

/*************************************************************/

/*  zipf generates a random number that follows Zipf         **
 **  distribution and lies between x1 and x2.                 **
 **  original code by Christos Faloutsos, 1995

 **  The original node outputs discrete data only. The current**
 **  function remedies the problem.			                 **
 **  Modified by Yufei Tao (08/Dec/02)                         **
 **************************************************************/
float zipf(float x1, float x2, double p) {

    float x;
    float i;
    double r, HsubV, sum;
    int V = 100;

    //float uniform();

    /* calculate the V-th harmonic number HsubV. WARNING: V>1 */
    HsubV = 0.0;
    for (i = 1; i <= V; i++)
        HsubV += 1.0 / pow((double) i, p);

    r = uniform(0., 1.) * HsubV;
    sum = 1.0;
    i = uniform(0, 1);
    while (sum < r) {
        //i++;  //commented by Yufei Tao
        i += uniform(1, 2);
        sum += 1.0 / pow((double) i, p);
    }

    /* i follows Zipf distribution and lies between 1 and V */

    /* x lies between 0. and 1. and then between x1 and x2 */
    x = ((float) i - (float) 1.) / ((float) V - (float) 1.);
    x = (x2 - x1) * x + x1;

    return (x);
}

float new_uniform(int _d_num) {
    float base = 1;
    float sum = 0;
    for (int i = 0; i < _d_num; i++) {
        int digit = (int) uniform(0, 10);
        if (digit == 10) digit = 9;
        sum += base*digit;
        base *= 10;
    }
    return sum;
}

float new_uniform(float _min, float _max) {
    float ran_base = (float) 9999999999;
    float ran = new_uniform(10);
    return ran / ran_base * (_max - _min) + _min;
}

/*****************************************************************
calculates the pdf of the guassian distribution

para:
- x:
- u: mean
- sigma: standard devication

return:
- pdf(x)

Coded by Yufei Tao, 2 aug 08
 *****************************************************************/

double normal_pdf(float _x, float _u, float _sigma) {
    double PI = 3.14159265;

    double ret = exp(-(_x - _u) * (_x - _u) / (2 * _sigma * _sigma));
    ret /= _sigma * sqrt(2 * PI);

    return ret;
}

/*****************************************************************
calculates the cdf of the normal distribution (mean 0, var 1)

para:
- x:
- step: controls the precision of the result. recommended value:
  0.001.

return:
- cdf(x)

Coded by Yufei Tao, 2 aug 08
 *****************************************************************/

double normal_cdf(float _x, float _step) {
    double ret = 0;

    for (double i = -10; i < _x; i += _step) {
        ret += _step * normal_pdf((float) i, 0.0f, 1.0f);
    }

    return ret;
}


int random_walk(float p) {
    if (uniform(0.0f, 1.0f) > p) return 1;
    return -1;
}
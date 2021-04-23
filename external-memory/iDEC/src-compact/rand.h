//--------------------------------------------
// this function contains necessary random number
// generators
// 
// collected by Yufei Tao
//--------------------------------------------

#ifndef RAND_H
#define RAND_H
//--------------------------------------------
float gaussian(float mean, float sigma);
float new_uniform(float _min, float _max);
double normal_pdf(float _x, float _u, float _sigma);
double normal_cdf(float _x, float _step);
float uniform(float _min, float _max);
float zipf(float x1, float x2, double p);
int random_walk(float p=0.5);
//--------------------------------------------
#endif


#ifndef __GENERAL_DEFINITION
#define __GENERAL_DEFINITION

#include <stdio.h>
#include <ctype.h>
#include <cstdint>

#define MAXREAL         1e20
#define MAXINT16        32767
#define MININT16        -32768
#define FLOATZERO       1e-7
#define MAX_DIMENSION   256

#define TRUE 1
#define FALSE 0

// NOTE: The following two macros will cause confliction 
#define min(a, b) (((a) < (b))? (a) : (b)  )
#define max(a, b) (((a) > (b))? (a) : (b)  )

//==========================================================
typedef float *floatptr;
typedef int16_t *int16ptr;

// struct SortMbr
// {
//     int dimension;
//     float *mbr;
//     float *center;
//     int index;
// };

struct SortMbr
{
    int dimension;
    int16_t *mbr;
    float *center;
    int index;
};

//General functions--------------------------------------
void error(char *_errmsg, bool _terminate);				//exit the program with an error message

//1d functions ------------------------------------------
bool inside(float &p, float &lb, float &ub);			//checks if a value falls in an interval
bool inside(int16_t &p, int16_t &lb, int16_t &ub);			//checks if a value falls in an interval

//Spatial functions--------------------------------------
float area(int dimension, float *mbr);					//calculates the area of a rectangle
float area(int dimension, int16_t *mbr);					//calculates the area of a rectangle

bool check_point_in_rec(float *v, float *mbr, int dimension);	//checks if a point falls in a rectangle
bool check_point_in_rec(int16_t *v, int16_t *mbr, int dimension);	//checks if a point falls in a rectangle

float margin(int dimension, float *mbr);				//the perimeter of a rectangle
float margin(int dimension, int16_t *mbr);				//the perimeter of a rectangle

float overlap(int dimension, float *r1, float *r2);		//the overlapping area of two rectangles
float overlap(int dimension, int16_t *r1, int16_t *r2);		//the overlapping area of two rectangles

float* overlapRect(int dimension, float *r1, float *r2);	//returns a rectangle which is the overlapping area of two rectangles
int16_t* overlapRect(int dimension, int16_t *r1, int16_t *r2);	//returns a rectangle which is the overlapping area of two rectangles

float ptDIST(float *p1, float *p2, int _dim);			//returns the L2 distance of two points
float ptDIST(int16_t *p1, int16_t *p2, int _dim);			//returns the L2 distance of two points

float MINMAXDIST(float *Point, float *bounces, int _dim);	//the minmaxdist between a point and a rectangle
float MINMAXDIST(int16_t *Point, int16_t *bounces, int _dim);	//the minmaxdist between a point and a rectangle

float MINDIST(float *Point, float *bounces, int Pdim);	//the mindist between a point and a rectangle
float MINDIST(int16_t *Point, int16_t *bounces, int Pdim);	//the mindist between a point and a rectangle

float MAXDIST(float *p, float *bounces, int dim);		//the maxdist between a point and a rectangle
float MAXDIST(int16_t *p, int16_t *bounces, int dim);		//the maxdist between a point and a rectangle

float MINDIST(float *Point, float *bounces, int Pdim, int _start_id); 
float MINDIST(int16_t *p, int16_t *bounces, int _dim, int _start_id);

float MbrMINDIST(float *_m1, float *_m2, int _dim);		//the mindist between two rectangles
float MbrMINDIST(int16_t *_m1, int16_t *_m2, int _dim);		//the mindist between two rectangles

float MbrMAXDIST(float *_m1, float *_m2, int _dim);		//the maxdist between two rectangles
float MbrMAXDIST(int16_t *_m1, int16_t *_m2, int _dim);		//the maxdist between two rectangles

bool section(int dimension, float *mbr1, float *mbr2);	//reports the topological relations between two rectangles
bool section(int dimension, int16_t *mbr1, int16_t *mbr2);	//reports the topological relations between two rectangles

void enlarge(int dimension, float **mbr, float *r1, float *r2);	//increases r1 to tightly enclose r2
void enlarge(int dimension, int16_t **mbr, int16_t *r1, int16_t *r2);	//increases r1 to tightly enclose r2

//functions used by the R*-split algorithms -------------
int sort_lower_mbr(const void *d1, const void *d2);
int sort_upper_mbr(const void *d1, const void *d2);
int sort_center_mbr(const void *d1, const void *d2);
int sortmindist(const void *element1, const void *element2);

void blank_print(int _n);
int compfloats(float _v1, float _v2);
//void error(char *_msg, bool _exit);
void get_leading_folder(char *_path, char *_folder);
void getFNameFromPath(char *_path, char *_fname);
char getnextChar(char **_s);
void getnextWord(char **_s, char *_w);
bool is_pow_of_2(int _v);
float l2_dist_int(int *_p1, int *_p2, int _dim);
void printINT_in_BIN(int _n, int _m);
float l2_dist_int(uint64_t *_p1, uint64_t *_p2, int _dim);
// only for test purpose
float hammingDistance(uint64_t *_p1, uint64_t *_p2, int _dim);

// copied from https://www.geeksforgeeks.org/check-whether-k-th-bit-set-not/
bool isKthBitSet(uint64_t n, unsigned k) ;
#ifdef UNIX
void strupr(char *_msg);
#endif
//-----------------------------------------------------------
#endif
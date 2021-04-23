/* 
 * File:   srs-hamming-compact.h
 * By:     yifangs
 */

#ifndef SRS_H
#define	SRS_H

#include "ptR-tree.h"
#include <cstdint>

typedef ptRTree * RtreePtr;

struct SRS_Hentry {
    int d; /* dimensionality */
    int id; /* id of the point */
    float gap; /* projected distance to query */
    int treeId;
    float dist; /* actual distance to query */
    int level;

    void setto(SRS_Hentry * _e) {
        d = _e->d;
        id = _e->id;
        gap = _e->gap;
        treeId = _e->treeId;
        dist = _e->dist;
        level = _e->level;
    }
};

class SRS {
public:
    //--=== on disk ===--
    int d; /* dimensionality */
    int n; /* cardinality */
    int B; /* page size in words */
    int m; /* m in paper */

    //--=== others ===--
    int numUsedTrees;
    char dsPath[100]; /* folder containing the dataset file */
    char dsName[100]; /* dataset file name */
    char indexPath[100]; /* folder containing the index */
    char dName[100]; /* data file name (binary file) */

    int L; /* number of srs-trees */

    float * a_array; /* each b-tree requires m hash functions, and each hash function requires a d-dimensional vector a, thus L*d in total*/
    RtreePtr *trees; /* rtrees */

    //--=== Functions ===--
    SRS();
    ~SRS();

    //--=== internal ===--
    virtual void freadNextEntry(FILE *_fp, int * _son, int * _key);
    virtual void gen_vectors();
    virtual void getTreeFname(int _i, char *_fname);
    virtual float getHashValue(int _tableID, uint64_t *_key);
    virtual void getHashPara(int _u, float **_a_vector);
    virtual int insert(int _treeID, int _son, uint64_t * _key);
    virtual int readParaFile(char *_fname);
    virtual float updateknn(SRS_Hentry * _rslt, SRS_Hentry *_he, int _k);
    virtual int writeParaFile(char *_fname);

//    //--=== external ===--
    virtual int buildFromFile(char *_dsPath, char *_indexPath);
    virtual int buildBinFromFile(char *_dsPath, char *_indexPath);
    virtual int ctBuildFromFile(char *_dsPath, char *_indexPath);
    virtual void init(int _d, int _n, int _B, int _L, int _m);
    // virtual int knn(int* _q, int _k, SRS_Hentry* _rslt, int _start_ID, int _numTrees, float _c, float _tau, float _t, bool useErwa);
    virtual int knn(uint64_t* _q, int _k, SRS_Hentry* _rslt, int _start_ID, int _numTrees, float _c, float _tau, float _t, bool useErwa);
    virtual int restore(char *_paraPath, int _start_ID, int _numTrees);
    // virtual int readData(int *_d, BlockFile *_bf, int _id);
    virtual int readData(uint64_t *_d, BlockFile *_bf, int _id);
};


#endif	/* SRS_H */


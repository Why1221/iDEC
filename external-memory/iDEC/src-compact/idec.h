#ifndef IDEC_H
#define	IDEC_H

// added by long @02132020 for measuring IO time
#ifdef MEASURE_IO_TIME
#include <gflags/gflags.h>
DECLARE_double(G_IO_TIME);
DECLARE_int32(G_RIO_CNT);
#endif 

#include "ptR-tree.h"
#include <cstdint>



typedef ptRTree * RtreePtr;

struct IDEC_Hentry {
    int d; /* dimensionality */
    int id; /* id of the point */
    float gap; /* projected distance to query */
    int treeId;
    float dist; /* actual distance to query */
    int level;

    void setto(IDEC_Hentry * _e) {
        d = _e->d;
        id = _e->id;
        gap = _e->gap;
        treeId = _e->treeId;
        dist = _e->dist;
        level = _e->level;
    }
};

class IDEC {
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

    int L; /* number of IDEC-trees */

    int * a_array; /* each b-tree requires m hash functions, and each hash function requires a d-dimensional vector a, thus L*d in total*/
    RtreePtr *trees; /* rtrees */

    //--=== Functions ===--
    IDEC();
    ~IDEC();

    //--=== internal ===--
    virtual void freadNextEntry(FILE *_fp, int * _son, int * _key);
    virtual void gen_vectors();
    virtual void getTreeFname(int _i, char *_fname);
    virtual int16_t getHashValue(int _tableID, uint64_t *_key);
    virtual void getHashPara(int _u, int **_a_vector);
    virtual int insert(int _treeID, int _son, uint64_t * _key);
    virtual int readParaFile(char *_fname);
    virtual float updateknn(IDEC_Hentry * _rslt, IDEC_Hentry *_he, int _k);
    virtual int writeParaFile(char *_fname);

//    //--=== external ===--
    virtual int buildFromFile(char *_dsPath, char *_indexPath);
    virtual int buildBinFromFile(char *_dsPath, char *_indexPath);
    virtual int ctBuildFromFile(char *_dsPath, char *_indexPath);
    virtual void init(int _d, int _n, int _B, int _L, int _m);
    virtual int knn(uint64_t* _q, int _k, IDEC_Hentry* _rslt, int _start_ID, int _numTrees, float _c, float _tau, float _t, bool useErwa);
    virtual int restore(char *_paraPath, int _start_ID, int _numTrees);
    virtual int readData(uint64_t *_d, BlockFile *_bf, int _id);
};


#endif	/* IDEC_H */


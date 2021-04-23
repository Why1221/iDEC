#ifndef __BINHEAP_H
#define __BINHEAP_H

class BinHeapEntry {
public:
    int leftcnt; //how many nodes in total in the left subtree.
    int rightcnt; //how many nodes in the right subtree.
    //When there is only 1 node in the heap (root), then both its leftcnt and rightcnt are 0;   
    BinHeapEntry *leftchild;
    BinHeapEntry *rightchild;
    BinHeapEntry *parent;
    void *data;

    //-----functions-----
    BinHeapEntry();
    ~BinHeapEntry();
    BinHeapEntry * crude_insert(BinHeapEntry * _new_e);
    BinHeapEntry * crude_delete();
    void recursive_data_wipeout(void (*_destroy_data)(const void *_e));
};

typedef BinHeapEntry * BinHeapEntryptr;

class BinHeap {
public:

    BinHeapEntry *root;
    int (*compare_func)(const void *_e1, const void *_e2);
    void (*destroy_data)(const void *_e);

    //-----functions-----
    BinHeap();
    ~BinHeap();
    void adjust_upward(BinHeapEntry *_he);
    void adjust_downward();
    virtual int compare(const void *_e1, const void *_e2);
    void insert(BinHeapEntry *_he);
    BinHeapEntry * remove();
    void swap_data(BinHeapEntry *_e1, BinHeapEntry *_e2);
};

#endif

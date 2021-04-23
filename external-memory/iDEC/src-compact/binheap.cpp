#include "binheap.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "gendef.h"

/*****************************************************************
Coded by Yufei Tao, 4 aug 08
 *****************************************************************/

BinHeapEntry::BinHeapEntry() {
    leftcnt = 0;
    rightcnt = 0;
    leftchild = NULL;
    rightchild = NULL;
    parent = NULL;

    data = NULL;
}

/*****************************************************************
Coded by Yufei Tao, 4 aug 08
 *****************************************************************/

BinHeapEntry::~BinHeapEntry() {
    if (leftchild) {
        delete leftchild;
        leftchild = NULL;
    }

    if (rightchild) {
        delete rightchild;
        rightchild = NULL;
    }

    if (data) {
        printf("warning: heap entry's data field not empty. possible memory leak. \n");
    }
}

/*****************************************************************
this function puts the new entry at the next position of the lowest 
level. first, if currently the binary tree is a perfect tree, put the
entry as the left child of the leftmost leaf. otherwise, put the 
entry as the left (or right if left exists) child of the leftmost 
node at the last-but-one-level node that doesn't have two children.

para:
- new_e: 

return:
- the final parent of new_e

Coded by Yufei Tao, 4 aug 08
 *****************************************************************/

BinHeapEntry * BinHeapEntry::crude_insert(BinHeapEntry * _new_e) {
    BinHeapEntry *ret = NULL;

    if (leftcnt == 0) {
        _new_e->parent = this;
        leftchild = _new_e;
        leftcnt++;
        ret = leftchild;
    } else if (rightcnt == 0) {
        _new_e->parent = this;
        rightchild = _new_e;
        rightcnt++;
        ret = rightchild;
    } else if (leftcnt == rightcnt) {
        leftcnt++;
        ret = leftchild->crude_insert(_new_e);
    } else {
        if (!is_pow_of_2(leftcnt + 1)) {
            leftcnt++;
            ret = leftchild->crude_insert(_new_e);
        } else {
            rightcnt++;
            ret = rightchild->crude_insert(_new_e);
        }
    }

    return ret;
}

/*****************************************************************
this function returns the rightmost leaf.

para:

return:
- the leaf

Coded by Yufei Tao, 4 aug 08
 *****************************************************************/

BinHeapEntry * BinHeapEntry::crude_delete() {
    BinHeapEntry *ret = NULL;

    if (leftcnt == 0) {
        ret = this;
    } else if (rightcnt == 0) {
        ret = leftchild;
        leftcnt--;
    } else if (leftcnt == rightcnt) {
        rightcnt--;
        ret = rightchild->crude_delete();
    } else {
        if (!is_pow_of_2(rightcnt + 1)) {
            rightcnt--;
            ret = rightchild->crude_delete();
        } else {
            leftcnt--;
            ret = leftchild->crude_delete();
        }
    }

    return ret;
}

/*****************************************************************
Coded by Yufei Tao, 4 aug 08
 *****************************************************************/

void BinHeapEntry::recursive_data_wipeout(void (*_destroy_data)(const void *_e)) {
    if (!_destroy_data)
        error("binheapentry::recurseive_data_wipeout - destroy_data is not given\n", true);

    if (leftchild)
        leftchild->recursive_data_wipeout(_destroy_data);
    if (rightchild)
        rightchild->recursive_data_wipeout(_destroy_data);

    _destroy_data(data);
    data = NULL;
}

/*****************************************************************
Coded by Yufei Tao, 4 aug 08
 *****************************************************************/

BinHeap::BinHeap() {
    root = NULL;
    compare_func = NULL;
    destroy_data = NULL;
}

/*****************************************************************
Coded by Yufei Tao, 4 aug 08
 *****************************************************************/

BinHeap::~BinHeap() {
    delete root;
    root = NULL;
}

/*****************************************************************
Coded by Yufei Tao, 4 aug 08
 *****************************************************************/

void BinHeap::insert(BinHeapEntry *_he) {
    _he->leftchild = _he->rightchild = NULL;
    _he->leftcnt = _he->rightcnt = 0;

    if (!_he->data) {
        printf("the data of the heap entry being inserted is null.\n");
        exit(1);
    }

    if (!root) {
        root = _he;
    } else {
        BinHeapEntry *leafhe = root->crude_insert(_he);
        adjust_upward(leafhe);
    }
}

/*****************************************************************
Coded by Yufei Tao, 4 aug 08
 *****************************************************************/

BinHeapEntry * BinHeap::remove() {
    BinHeapEntry *ret;

    if (!root) {
        ret = NULL;
    } else if (root->leftcnt == 0) {
        ret = root;
        root = NULL;
    } else {
        BinHeapEntry *leafhe = root->crude_delete();
        swap_data(root, leafhe);
        if (leafhe == leafhe->parent->leftchild)
            leafhe->parent->leftchild = NULL;
        else
            leafhe->parent->rightchild = NULL;
        leafhe->parent = NULL;
        ret = leafhe;
        adjust_downward();
    }

    return ret;
}

/*****************************************************************
this function fixes the data field to satisfy the heap property. 
start from a leaf and go upward

Coded by Yufei Tao, 4 aug 08
 *****************************************************************/

void BinHeap::adjust_upward(BinHeapEntry *_he) {
    BinHeapEntry *he = _he;

    while (he->parent) {
        int rslt = compare(he->parent->data, he->data);

        if (rslt == 1) {
            swap_data(he, he->parent);
            he = he->parent;
        } else
            break;
    }
}

/*****************************************************************
swap the data fields of two entries

para:
- e1:
- e2
  
Coded by Yufei Tao, 4 aug 08
 *****************************************************************/

void BinHeap::swap_data(BinHeapEntry *_e1, BinHeapEntry *_e2) {
    void *tmp;
    tmp = _e1->data;
    _e1->data = _e2->data;
    _e2->data = tmp;
}

/*****************************************************************
this function fixes the data field to satisfy the heap property. 
start from the root and go downward

Coded by Yufei Tao, 4 aug 08
 *****************************************************************/

void BinHeap::adjust_downward() {
    BinHeapEntry *he = root;

    while (he->leftchild) {
        if (!he->rightchild) {
            int rslt = compare(he->data, he->leftchild->data);
            if (rslt == 1) {
                swap_data(he, he->leftchild);
                he = he->leftchild;
            } else
                break;
        } else {
            int rslt1 = compare(he->data, he->leftchild->data);
            int rslt2 = compare(he->data, he->rightchild->data);

            if (rslt1 != 1 && rslt2 != 1)
                break;
            else {
                int rslt3 = compare(he->leftchild->data, he->rightchild->data);
                if (rslt3 == -1) {
                    swap_data(he, he->leftchild);
                    he = he->leftchild;
                } else {
                    swap_data(he, he->rightchild);
                    he = he->rightchild;
                }
            }
        }
    }
}

/*****************************************************************
this function fixes the data field to satisfy the heap property. 
start from the root and go downward

Coded by Yufei Tao, 4 aug 08
 *****************************************************************/

int BinHeap::compare(const void *_e1, const void *_e2) {
    if (!compare_func) {
        printf("compare func not set yet.\n");
        exit(1);
    }

    //you may wonder why I need to introduce such an intermediate ::compare member func, instead
    //of just use compare_func directly at all places of ::compare. there is actually an important
    //point to this. in a nutshell, it allows much better flexibility when compare_func is not
    //exactly the comparison functions we need. for example, this is the case for external sort.
    //see gadget.h.

    return compare_func(_e1, _e2);
}

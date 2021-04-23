#ifndef __RTNODE
#define __RTNODE
//------------------------------------------------------------
#include "rtree.h"
//------------------------------------------------------------
class SortedLinList;
class Entry;
class RTree;
class Heap;
//------------------------------------------------------------
typedef Entry *EntryPtr; //this is necessary to achieve polymophism

class RTNode
{
public:
//--===on disk===--
	char level; 
	int block;
	int num_entries;
	EntryPtr *entries;
//--===others===--
	bool dirty;
	int capacity;
    int dimension;
	RTree *my_tree;  
//--===functions===--
	RTNode();
    
	virtual ~RTNode();
	virtual int choose_subtree(Entry* _e);	//selects the subtree in which the given entry should be inserted
	virtual void enter(Entry *_de);			//registers a new entry into the code
	virtual R_DELETE delete_entry(Entry *e);	//deletes an object from the subtree
	virtual bool FindLeaf(Entry *_e);		//checks if an object exists in the tree
	virtual float *get_mbr();				//gets the MBR of the node
	virtual void init(int _level, RTree *_rt);	//creates a new node on the disk at the level specified
	virtual void init(RTree *_rt, int _block);	//reads an existing node from the disk
	virtual R_OVERFLOW insert(Entry *_d, RTNode **_sn);	//inserts an entry into the subtree of this node
	virtual Entry* new_one_entry();			//initiates an entry. Override this function to create an object
											//of a class that inherits Entry
	virtual RTNode* new_one_node();			//initiates a node. Override this function to create an object
											//of a class that inherits RTNode
	virtual int rangeQuery(float *_mbr, SortedLinList *_res, bool _ids_wanted);	//perform a range query
	virtual void read_from_buffer(char *_buffer);	//loads the content of a node from a disk buffer page
	virtual int mbr_split(float **_mbr, int **_distribution);	//this function implements the R*-split strategy 
																//(see the head of function for detailed explanations)
	virtual void split(RTNode *_sn);			//splits the node into itself and _sn
	virtual void write_to_buffer(char *_buffer);	//writes the content of a node to a disk buffer page

	//int get_num_of_data();
	//bool is_data_node() { return (level==0) ? TRUE : FALSE ;};
	//void NNSearch(float *QueryPoint, SortedLinList *res, float *nearest_distanz);
	//void print();
	//void rank_qry_inquiry(float *_weight, float _qscore, int *_rslt);
};

#endif
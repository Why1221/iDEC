#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "rtnode.h"
#include "rtree.h"
#include "entry.h"
#include "blk_file.h"
#include "cache.h"
#include "linlist.h"
#include "gendef.h"
//------------------------------------------------------------
RTNode::RTNode()
{
  entries = NULL; 
  my_tree = NULL;
}
//------------------------------------------------------------
//this function creates a new node on the disk

void RTNode::init(int _level, RTree *_rt)
{
  char *b;
  int header_size;
  Entry * d;
  int i;

  my_tree = _rt;
  dimension = _rt->dimension;
  num_entries = 0;
  dirty = TRUE;
  level = _level;

  d = new_one_entry();
  d->init(_rt);
  d->level = level;
  header_size = sizeof(char) + sizeof(int);  // level + num_entries
  capacity = (_rt -> file -> get_blocklength() - header_size) / d -> get_size();
  delete d;

  entries = new EntryPtr[capacity];
  for (i = 0; i < capacity; i ++)
  {
    entries[i] = new_one_entry();
    entries[i]->init(my_tree);
    entries[i]->level = level;
  }

  //assign a new block on the disk
  b = new char[_rt -> file -> get_blocklength()];
  block = _rt -> file -> append_block(b);
  delete [] b;
}
//------------------------------------------------------------
//reads an existing node from the disk

void RTNode::init(RTree *_rt, int _block)
{
  char *b;
  int header_size;
  Entry * d;
  int i;

  my_tree = _rt;
  dimension = _rt->dimension;
  num_entries = 0;
  dirty = FALSE;

  block = _block;
  b = new char[_rt -> file -> get_blocklength()];
  if (_rt -> cache == NULL) // no cache
    _rt -> file -> read_block(b, block);
  else
    _rt -> cache -> read_block(b, block, _rt);

  //here it is very important to get the level of the node first before initalizing the entries
  memcpy(&level, b, sizeof(char));

  d = new_one_entry();
  d->init(_rt);
  d->level = level;
  header_size = sizeof(char) + sizeof(int);
  capacity = (_rt -> file -> get_blocklength() - header_size) / d -> get_size();
  delete d;

  entries = new EntryPtr[capacity];
  for (i = 0; i < capacity; i ++)
  {
    entries[i] = new_one_entry();
    entries[i]->init(my_tree);
    entries[i]->level = level;
  }
    
  read_from_buffer(b);
  delete [] b;
}
//------------------------------------------------------------
RTNode::~RTNode()
{
  char *b;

  if (dirty)
  {
    b = new char[my_tree->file->get_blocklength()];
    write_to_buffer(b);

    if (my_tree->cache == NULL) // no cache
      my_tree->file->write_block(b, block);
    else
      my_tree->cache->write_block(b, block, my_tree);

    delete [] b;
    dirty = false;
  }


  for (int i = 0; i < capacity; i ++)
  {
    delete entries[i]; entries[i] = NULL;
  }


  /*
  //	Modified by Xiaokui Xiao, Jun 8, 2006
  for (int i = 0; i < capacity; ++i)
  {
  if (entries[i] != NULL)
  {
  delete entries[i];
  entries[i] = NULL;
  }
  }
  */
  delete [] entries;
}
//------------------------------------------------------------
//selects the subtree in which the given entry should be inserted

int RTNode::choose_subtree(Entry* _e)
{
  float *mbr = new float[2 * dimension];
  memcpy(mbr, _e->bounces, 2 * dimension * sizeof(float));

  int i, j, follow, minindex, *inside, inside_count, *over;
  float *bmbr, old_o, o, omin, a, amin, f, fmin;

  inside_count = 0;
  minindex = 0;
  inside = new int[num_entries];
  over = new int[num_entries];
  for (i = 0; i < num_entries; i++)
  {
    switch (entries[i]->section(mbr))
    {
      case INSIDE:
        inside[inside_count++] = i;
        break;
    }
  }

  if (inside_count == 1)
    // Case 1: There is exactly one dir_mbr that contains mbr
    follow = inside[0];
  else if (inside_count > 1)
    // Case 2: There are many dir_mbrs that contain mbr
    // choose the one with the minimum area
  {
    //    	fmin = MAXREAL;
    fmin=area(dimension, entries[inside[0]]->bounces); minindex=0;
    //changed for hi dimensionality

    for (i = 1; i < inside_count; i++)
    {
      f = area(dimension, entries[inside[i]]->bounces);
      if (f < fmin)
      {
        minindex = i;
        fmin = f;
      }
    }
    follow = inside[minindex];
  }
  else
    // Case 3: There are no dir_mbrs that contain mbr
    // choose the one for which insertion causes the minimun overlap if son_is_data
    // else choose the one for which insertion causes the minimun area enlargement
  {
    if (level == 1) // son_is_data
    {
      omin = MAXREAL;
      fmin = MAXREAL;
      amin = MAXREAL;
      for (i = 0; i < num_entries; i++)
      {
        enlarge(dimension, &bmbr, mbr, entries[i]->bounces);

        // calculate area and area enlargement
        a = area(dimension, entries[i]->bounces);
        f = area(dimension, bmbr) - a;

        // calculate overlap before enlarging entry_i
        old_o = o = 0.0;
        for (j = 0; j < num_entries; j++)
        {
          if (j != i)
          {
            old_o += overlap(dimension,
                             entries[i]->bounces,
                             entries[j]->bounces);
            o += overlap(dimension,
                         bmbr,
                         entries[j]->bounces);
          }
        }
        o -= old_o;

        // is this entry better than the former optimum ?
        if ((o < omin) ||
            (o == omin && f < fmin) ||
            (o == omin && f == fmin && a < amin))
        {
          minindex = i;
          omin = o;
          fmin = f;
          amin = a;
        }
        delete [] bmbr;
      }
    }
    else // son is not a data node
    {
      fmin = MAXREAL;
      amin = MAXREAL;
      for (i = 0; i < num_entries; i++)
      {
        enlarge(dimension, &bmbr, mbr, entries[i]->bounces);

        // calculate area and area enlargement
        a = area(dimension, entries[i]->bounces);
        f = area(dimension, bmbr) - a;

        // is this entry better than the former optimum ?
        if ((f < fmin) || (f == fmin && a < amin))
        {
          minindex = i;
          fmin = f;
          amin = a;
        }
        delete [] bmbr;
      }
    }

    follow = minindex;

    dirty = TRUE;
  }

  delete [] inside;
  delete [] over;
  delete [] mbr; //added by Jianbin

  return follow;
}
//------------------------------------------------------------
//deletes an object from the subtree

R_DELETE RTNode::delete_entry(Entry *_e)
{
  RTNode *succ;
  float *tmp;
  if (level > 0)
  {
    if (this == my_tree->root_ptr)
      //i.e. this is the root
    {
      for (int i = 0; i < num_entries; i++)
      {
        tmp = overlapRect(dimension, entries[i]->bounces, _e->bounces);
        if (tmp != NULL)
        {
          delete [] tmp;
          succ = entries[i]->get_son();
          R_DELETE del_ret;
          del_ret = succ -> delete_entry(_e);
          if (del_ret != NOTFOUND)
          {
            switch (del_ret)
            {
              case NORMAL:

                float *mbr;
                mbr = succ -> get_mbr();
                memcpy(entries[i]->bounces, mbr, sizeof(float) * 2 * dimension);
                dirty = true;
                delete [] mbr;

                delete entries[i]->son_ptr;
                entries[i]->son_ptr = NULL;

                return NORMAL;
                break;

              case ERASED:
                delete entries[i]->son_ptr;
                entries[i]->son_ptr = NULL;

                int j;
                for (j = i; j < num_entries - 1; j++)
                  entries[j]->set_equal_to(entries[j+1]);
                for (j = num_entries - 1; j < capacity; j++)
                  entries[j]->son_ptr = NULL;

                num_entries--;

                dirty = true;
                return NORMAL;
                break;
            }
          }
        }
      }
      //Not found;
      return NOTFOUND;
    }
    else//is not root and not leaf
    {
      for (int i = 0; i < num_entries; i++)
      {
        tmp = overlapRect(dimension, entries[i]->bounces, _e->bounces);
        if (tmp != NULL)
        {
          delete [] tmp;
          succ = entries[i]->get_son();
          R_DELETE del_ret;
          del_ret = succ->delete_entry(_e);
          if (del_ret != NOTFOUND)
          {
            switch (del_ret)
            {
              case NORMAL:

                float *mbr;
                mbr = succ -> get_mbr();
                memcpy(entries[i]->bounces, mbr, sizeof(float) * 2 * dimension);
                dirty = true;
                delete [] mbr;

                entries[i]->del_son();

                return NORMAL;
                break;

              case ERASED:

                entries[i]->del_son();

                int j;
                for (j = i; j < num_entries - 1; j++)
                  entries[j]->set_equal_to(entries[j + 1]);
                for (j = num_entries - 1; j < capacity; j++)
                  entries[j]->son_ptr = NULL;
							
                num_entries--;
                dirty = true;

                if (num_entries < (int)ceil(0.4 * capacity))
                {
                  for (int j = 0; j < num_entries; j++)
                  {
                    Linkable *e;
                    e = entries[j]->gen_Linkable();
                    e->level=level;
                    my_tree->deletelist->insert(e);
                  }

                  my_tree->num_of_inodes --;
                  return ERASED;
                }
                else
                  return NORMAL;
                break;
            }
          }
        }
      }
    }
  }
  else//it's a leaf
  {
    for (int i = 0; i < num_entries; i++)
    {
      if (entries[i]->check_equal(_e))
      {
        for (int j = i; j < num_entries-1; j++)
          entries[j]->set_equal_to(entries[j+1]);
				
        num_entries--;
        dirty = true;

        if (this != my_tree -> root_ptr && num_entries < (int)ceil(0.4 * capacity))
        {
          for (int k = 0; k < num_entries; k++)
          {
            Linkable *en;
            en = entries[k]->gen_Linkable();
            en->level = 0;
            my_tree->deletelist->insert(en);
          }

          my_tree -> num_of_dnodes --;
          return ERASED;
        }
        else
          return NORMAL;
      }
    }
    return NOTFOUND;
  }
}
//------------------------------------------------------------
//registers a new entry into the code

void RTNode::enter(Entry *_de)
{
  if (num_entries > (capacity-1))
    error("RTNode::enter: called, but node is full", TRUE);

  entries[num_entries]->set_equal_to(_de);
  num_entries++;
  dirty = true;
  _de->son_ptr = NULL;
  delete _de;
}
//------------------------------------------------------------
//checks if an object exists in the tree

bool RTNode::FindLeaf(Entry *_e)
{
  RTNode *succ;
  if (level > 0)
  {
    for (int i = 0; i < num_entries; i++)
    {
      float *f;
      f = overlapRect(my_tree -> dimension,
                      entries[i]->bounces, _e -> bounces);
      if (f != NULL)
      {
        delete [] f;
        succ = entries[i]->get_son();
        bool find;
        find = succ->FindLeaf(_e);
        entries[i]->del_son();
        if (find)
          return true;
      }
    }
    return false;
  }
  else
  {
    for (int i = 0; i < num_entries; i++)
    {
      if (entries[i]->check_equal(_e))
        return true;
    }
    return false;
  }
  return false;
}
//------------------------------------------------------------
//gets the MBR of the node

float* RTNode::get_mbr()
{
  int i, j;
  float *mbr;

  mbr = new float[2*dimension];
  for (i = 0; i < 2*dimension; i ++ )
    mbr[i] = entries[0]->bounces[i];

  for (j = 1; j < num_entries; j++)
  {
    for (i = 0; i < 2*dimension; i += 2)
    {
      mbr[i]   = min(mbr[i],   entries[j]->bounces[i]);
      mbr[i+1] = max(mbr[i+1], entries[j]->bounces[i+1]);
    }
  }

  return mbr;
}
//------------------------------------------------------------
//int RTNode::get_num_of_data()
//{
//    int i, sum;
//    RTNode* succ;
//
//    if (level == 0)
//        return num_entries;
//
//    sum = 0;
//    for (i = 0; i < num_entries ; i++)
//    {
//        succ = entries[i]->get_son();
//        sum += succ->get_num_of_data();
//		entries[i].del_son();
//    }
//
//    return sum;
//}
//------------------------------------------------------------
//inserts an entry into the subtree of this node

R_OVERFLOW RTNode::insert(Entry *_d, RTNode **_sn)
{
  int follow;
  RTNode *succ, *new_succ;
  RTNode *brother;
  Entry *de;
  R_OVERFLOW ret;
  float *mbr,*nmbr;

  int i, last_cand;
  float *center;
  SortMbr *sm;
  EntryPtr *new_entries;

  if (level > 0) // direcrtory node
  {
    if (level > _d -> level)
    {
      follow = choose_subtree(_d);
      succ = entries[follow]->get_son();
      ret = succ -> insert(_d, &new_succ);
    
      mbr = succ -> get_mbr();
      memcpy(entries[follow]->bounces, mbr, sizeof(float) * 2 * dimension);
      delete [] mbr;

      entries[follow]->del_son();

      if (ret == SPLIT)
        // node has split into itself and *new_succ
      {
        if (num_entries == capacity)
          error("RTNode::insert: maximum capacity violation", TRUE);

        de = new_one_entry();
        de->init(my_tree);
        de->level = level;
        nmbr = new_succ -> get_mbr();
        memcpy(de -> bounces, nmbr, 2 * dimension * sizeof(float));
        delete [] nmbr;
        de -> son = new_succ -> block;
        delete new_succ;
        de -> son_ptr = NULL;
        enter(de);

        if (num_entries == (capacity - 1))
        {
          brother = new_one_node();
          brother->init(level, my_tree);
				
          my_tree -> num_of_inodes++;
          brother -> level = level;
          split(brother);
          *_sn = brother;
          ret = SPLIT;
        }
        else
          ret = NONE;
      }
      dirty = TRUE;

      return ret;
    }
    else //level==d->level
    {
      enter(_d);    //note that d will be deleted on return
		    
      if (num_entries == (capacity - 1))
        // maximun no of entries --> Split
        // this happens already if the node is nearly filled
        // for the algorithms are more easy then
      {
        brother = new_one_node(); 
        brother->init(level, my_tree);
        my_tree -> num_of_inodes++;
        brother -> level = level;
        split(brother);
        *_sn = brother;
        ret = SPLIT;
      }
      else
        ret = NONE;

      dirty=true;
      return ret;
    }	
  }
  else // data (leaf) node
  {
    if (num_entries == capacity)
      error("RTDataNode::insert: maximum capacity violation", TRUE);

    enter(_d);

    dirty = TRUE;

    if (num_entries == (capacity - 1))
      // maximum # of entries --> Split
      // this happens already if the node is nearly filled
      // for the algorithms are more easy then
    {
      if (my_tree->re_level[0] == FALSE && my_tree -> root_ptr -> level != level)
        // there was no reinsert on level 0 during this insertion
        //-------------------------------------------------------
        //Here I changed the condition as if it is already root, no need
        //to reinsert.  Split directly in this case
        //-----------By TAO Yufei
      {
        // calculate center of page
        mbr = get_mbr();
        center = new float[dimension];
        for (i = 0; i < dimension; i++)
          center[i] = (mbr[2*i] + mbr[2*i+1]) / 2.0;

        new_entries = new EntryPtr[capacity];

        for (i = 0; i < capacity; i ++)
        {
          new_entries[i] = new_one_entry();
          new_entries[i]->init(my_tree);
          new_entries[i]->level = level;
        }

        sm = new SortMbr[num_entries];
        for (i = 0; i < num_entries; i++)
        {
          sm[i].index = i;
          sm[i].dimension = dimension;
          sm[i].mbr = entries[i]->bounces;
          sm[i].center = center;
        }

        qsort(sm, num_entries, sizeof(SortMbr), sort_center_mbr);

        last_cand = (int) ((float)num_entries * 0.30);

        // copy the nearest candidates to new array
        for (i = 0; i < num_entries - last_cand; i++)
          new_entries[i]->set_equal_to(entries[sm[i].index]);

        // insert candidates into reinsertion list
        for ( ; i < num_entries; i++)
        {
          Linkable *nd = entries[sm[i].index]->gen_Linkable();
          nd->level = level;
          my_tree -> re_data_cands -> insert(nd);
        }

        // free and copy data array
        delete [] entries;
        entries = new_entries;
				
        delete sm;
        delete [] mbr;
        delete [] center;
        my_tree -> re_level[0] = TRUE;

        // correct # of entries
        num_entries -= last_cand;

        // must write page
        dirty = TRUE;

        return REINSERT;
      }
      else  //there has been reinsertion on this level
      {
        *_sn = new_one_node();
        (*_sn)->init(level, my_tree);
        (*_sn)->level = level;
        my_tree->num_of_dnodes++;
        split((RTNode *) *_sn);
      }
      return SPLIT;
    }
    else
      return NONE;
  }
}
//------------------------------------------------------------
//void RTNode::NNSearch(float *QueryPoint, 
//					  SortedLinList *res,
//				      float *nearest_distanz)
//{
//	if (level > 0)
//	{
//		float *minmax_distanz;		// Array fuer MINMAXDIST aller Eintr"age
//		int *indexliste;		// Liste (for Sorting and Prunching)
//		int i,j,k,last,n;
//		float akt_min_dist;		// minimal distanz computed upto now 
//		float minmaxdist,mindist;
//
//		BranchList *activebranchList;
//    
//		n = num_entries;
//    
//		//k = res->get_num(); 	// wird haben eine k-nearest-Narbor-Query
//    
//		//*nearest_distanz = res->get(k-1)->distanz;  // der aktuell letzte 
//													// n"achste Nachbar wird
//													// versucht zu ersetzen.
//		if (res -> get_num() > 0)
//		{
//			if (*nearest_distanz != res -> get_first() -> distanz)
//			{
//				printf("testing...\n");
//				*nearest_distanz = res -> get_first() -> distanz;
//			}
//		}
//
//		activebranchList = new BranchList [n]; // Array erzeugen mit n Elementen
// 
//		for( i = 0; i < n; i++)
//		{
//			activebranchList[i].entry_number = i;
//			activebranchList[i].minmaxdist = MINMAXDIST(QueryPoint,entries[i]->bounces);
//			activebranchList[i].mindist = MINDIST(QueryPoint,entries[i]->bounces, dimension);	
//		}	
//
//		// sortbranchList
//		qsort(activebranchList,n,sizeof(BranchList),sortmindist); 
// 
//		// pruneBrunchList
//		last = pruneBrunchList(nearest_distanz,activebranchList,n);
//
//		for( i = 0; i < last; i++)
//		{
//			entries[activebranchList[i].entry_number]->get_son()->NNSearch(QueryPoint, res, nearest_distanz);
//			entries[i].del_son();
// 			
//			last = pruneBrunchList(nearest_distanz,activebranchList,last);
//		}
//
//		delete [] activebranchList;
//	}
//	else //level == 0
//	{
//		int i,j;
//		float nearest_dist,distanz;
//		bool t;
//		Linkable *element;
//
//		for (i = 0; i < num_entries; i++)
//		{
//			//distanz = objectDIST(QueryPoint,entries[i]->bounces);
//			//i change this line to cater for both points and rectangles
//			distanz = MINDIST(QueryPoint,entries[i]->bounces, dimension);
//			
//			if (distanz <= *nearest_distanz)
//			{
//				// l"osche letztes Elemente der res-Liste		
//				//res->get(k-1);
//				//t = res->erase();
//				if (res -> get_num() > 0)
//				{
//					Linkable *lin = res -> get_first();
//					res -> erase();
//				}
//				
//				// Erzeuge neuen Punkt in der res-Liste
//				element = entries[i].gen_Linkable();
//				element->distanz = distanz;
//				res->insert(element);
//				
//				//*nearest_distanz = res->get(k-1)->distanz;			
//				*nearest_distanz = distanz;
//			}
//		}
//	}
//}
//------------------------------------------------------------
//void RTNode::print()
//{
//    int i;
//
//	printf("level %d  Block: %d\n", level, block);
//	
//    for (i = 0; i < num_entries ; i++)
//    {
//        printf("(%4.1lf, %4.1lf, %4.1lf, %4.1lf)\n",
//	       entries[i]->bounces[0],
//	       entries[i]->bounces[1],
//	       entries[i]->bounces[2],
//	       entries[i]->bounces[3]);
//    }
//}
//------------------------------------------------------------
//perform a range query

int RTNode::rangeQuery(float *_mbr, SortedLinList *_res, bool _ids_wanted)
{
  int i, n;
  int ret = 0;
  SECTION s;
  RTNode *succ;

  n = num_entries;
  for (i = 0; i < n; i++)
  {
    s = entries[i]->section(_mbr);
    if (s == INSIDE || s == OVERLAP)
    {
      if (level == 0)
      {
        if (_ids_wanted)
        {
          Linkable *copy;
          copy = entries[i]->gen_Linkable();
          _res -> insert(copy);
        }

        ret ++;
      }
      else
      {
        succ = entries[i]->get_son();
        ret += succ->rangeQuery(_mbr, _res, _ids_wanted);
        entries[i]->del_son();
      }
    }
  }

  return ret;
}
//------------------------------------------------------------
//writes the content of an entry to a disk buffer page

void RTNode::read_from_buffer(char *_buffer)
{
  int i, j, s;

  // Level
  memcpy(&level, _buffer, sizeof(char));
  j = sizeof(char);

  // num_entries
  memcpy(&num_entries, &(_buffer[j]), sizeof(int));
  j += sizeof(int);

  s = entries[0]->get_size();
  for (i = 0; i < num_entries; i++)
  {
    entries[i]->read_from_buffer(&_buffer[j]);
    j += s;
  }
}
//------------------------------------------------------------
//this function implements the R*-split strategy. 

//"mbr" is an array of MBRs. When the function finishes, "array" distribution will be filled with the indexes of the 
//entries in the node, sorted according to the decided split manner. Furthermore, the function returns a value v
//which indicates that entries corresponding to the first v elements in "distribution" will be put in the first
//node, and the others in the second node. 

//For example, assume that the node being split has 4 entries, whose indexes in the "entries" array are 0, ..., 4
//respectively. Assume that this function yields the following "distribution" array:

//{4, 0, 2, 3, 1}

//and returns a value 3. It means that the entries with indexes 4, 0, 2 should be put in the first node, and 
//the others in the second.

int RTNode::mbr_split(float **_mbr, int **_distribution)
{
  bool lu = FALSE;
  int i, j, k, l, s, n, m1, dist, split_axis = 0;
  SortMbr *sml, *smu;
  float minmarg, marg, minover, mindead, dead, over, *rxmbr, *rymbr;

  n = num_entries;

  m1 = (int) ceil((float)n * 0.40);
  //added by yifang
  dist = m1;
  //end added by yifang

  sml = new SortMbr[n];
  smu = new SortMbr[n];
  rxmbr = new float[2*dimension];
  rymbr = new float[2*dimension];

  // choose split axis
  minmarg = MAXREAL;
  for (i = 0; i < dimension; i++)
    // for each axis
  {
    for (j = 0; j < n; j++)
    {
      sml[j].index = smu[j].index = j;
      sml[j].dimension = smu[j].dimension = i;
      sml[j].mbr = smu[j].mbr = _mbr[j];
    }

    // Sort by lower and upper value perpendicular axis_i
    qsort(sml, n, sizeof(SortMbr), sort_lower_mbr);
    qsort(smu, n, sizeof(SortMbr), sort_upper_mbr);

    marg = 0.0;
    // for all possible distributions of sml
    for (k = 0; k < n - 2 * m1 + 1; k++)
    {
      for (s = 0; s < 2 * dimension; s += 2)
      {
        rxmbr[s] =    MAXREAL;
        rxmbr[s+1] = -MAXREAL;
      }
      for (l = 0; l < m1 + k; l++)
      {
        for (s = 0; s < 2*dimension; s += 2)
        {
          rxmbr[s] =   min(rxmbr[s],   sml[l].mbr[s]);
          rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
        }
      }
      marg += margin(dimension, rxmbr);

      for (s = 0; s < 2 * dimension; s += 2)
      {
        rxmbr[s] =    MAXREAL;
        rxmbr[s+1] = -MAXREAL;
      }
      for ( ; l < n; l++)
      {
        for (s = 0; s < 2 * dimension; s += 2)
        {
          rxmbr[s] =   min(rxmbr[s],   sml[l].mbr[s]);
          rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
        }
      }
      marg += margin(dimension, rxmbr);
    }

    // for all possible distributions of smu
    for (k = 0; k < n - 2 * m1 + 1; k++)
    {
      // now calculate margin of R1
      // initialize mbr of R1
      for (s = 0; s < 2 * dimension; s += 2)
      {
        rxmbr[s] =    MAXREAL;
        rxmbr[s+1] = -MAXREAL;
      }
      for (l = 0; l < m1+k; l++)
      {
        // calculate mbr of R1
        for (s = 0; s < 2 * dimension; s += 2)
        {
          rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
          rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
        }
      }
      marg += margin(dimension, rxmbr);

      // now calculate margin of R2
      // initialize mbr of R2
      for (s = 0; s < 2 * dimension; s += 2)
      {
        rxmbr[s] =    MAXREAL;
        rxmbr[s+1] = -MAXREAL;
      }
      for ( ; l < n; l++)
      {
        // calculate mbr of R1
        for (s = 0; s < 2 * dimension; s += 2)
        {
          rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
          rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
        }
      }
      marg += margin(dimension, rxmbr);
    }

    if (marg < minmarg)
    {
      split_axis = i;
      minmarg = marg;
    }
  }

  // choose best distribution for split axis
  for (j = 0; j < n; j++)
  {
    sml[j].index = smu[j].index = j;
    sml[j].dimension = smu[j].dimension = split_axis;
    sml[j].mbr = smu[j].mbr = _mbr[j];
  }

  // Sort by lower and upper value perpendicular split axis
  qsort(sml, n, sizeof(SortMbr), sort_lower_mbr);
  qsort(smu, n, sizeof(SortMbr), sort_upper_mbr);

  minover = MAXREAL;
  mindead = MAXREAL;
  // for all possible distributions of sml and snu
  for (k = 0; k < n - 2 * m1 + 1; k++)
  {
    dead = 0.0;
    for (s = 0; s < 2 * dimension; s += 2)
    {
      rxmbr[s] =    MAXREAL;
      rxmbr[s+1] = -MAXREAL;
    }
    for (l = 0; l < m1 + k; l++)
    {
      for (s = 0; s < 2*dimension; s += 2)
      {
        rxmbr[s] =   min(rxmbr[s],   sml[l].mbr[s]);
        rxmbr[s+1] = max(rxmbr[s+1], sml[l].mbr[s+1]);
      }
      dead -= area(dimension, sml[l].mbr);
    }
    dead += area(dimension, rxmbr);
    //**************note**************
    //this does not compute the dead space for all the cases.  some overlapping
    //area may be subtrated twice.
    //********************************

    for (s = 0; s < 2*dimension; s += 2)
    {
      rymbr[s] =    MAXREAL;
      rymbr[s+1] = -MAXREAL;
    }
    for ( ; l < n; l++)
    {
      for (s = 0; s < 2*dimension; s += 2)
      {
        rymbr[s] =   min(rymbr[s],   sml[l].mbr[s]);
        rymbr[s+1] = max(rymbr[s+1], sml[l].mbr[s+1]);
      }
      dead -= area(dimension, sml[l].mbr);
    }
    dead += area(dimension, rymbr);

    over = overlap(dimension, rxmbr, rymbr);

    if ((over < minover) || (over == minover) && dead < mindead)
    {
      minover = over;
      mindead = dead;
      dist = m1+k;
      lu = TRUE;
    }

    //Now we do the same thing for smu
    dead = 0.0;
    for (s = 0; s < 2*dimension; s += 2)
    {
      rxmbr[s] =    MAXREAL;
      rxmbr[s+1] = -MAXREAL;
    }
    for (l = 0; l < m1+k; l++)
    {
      for (s = 0; s < 2*dimension; s += 2)
      {
        rxmbr[s] =   min(rxmbr[s],   smu[l].mbr[s]);
        rxmbr[s+1] = max(rxmbr[s+1], smu[l].mbr[s+1]);
      }
      dead -= area(dimension, smu[l].mbr);
    }
    dead += area(dimension, rxmbr);

    for (s = 0; s < 2*dimension; s += 2)
    {
      rymbr[s] =    MAXREAL;
      rymbr[s+1] = -MAXREAL;
    }
    for ( ; l < n; l++)
    {
      for (s = 0; s < 2*dimension; s += 2)
      {
        rymbr[s] =   min(rymbr[s],   smu[l].mbr[s]);
        rymbr[s+1] = max(rymbr[s+1], smu[l].mbr[s+1]);
      }
      dead -= area(dimension, smu[l].mbr);
    }
    //correcting errors
    dead += area(dimension, rymbr);

    over = overlap(dimension, rxmbr, rymbr);

    if ((over < minover) || (over == minover) && dead < mindead)
    {
      minover = over;
      mindead = dead;
      dist = m1+k;
      lu = FALSE;
    }
  }

  // calculate best distribution
  // the array distribution is deleted in split(rtnode *sn);
  *_distribution = new int[n];
  for (i = 0; i < n; i++)
  {
    if (lu)
      (*_distribution)[i] = sml[i].index;
    else
      (*_distribution)[i] = smu[i].index;
  }

  delete [] sml;
  delete [] smu;
  delete [] rxmbr;
  delete [] rymbr;

  return dist;
}
//------------------------------------------------------------
//splits the node into itself and _sn

void RTNode::split(RTNode *_sn)
{
  int i, *distribution, dist, n;
  float **mbr_array;
  EntryPtr *new_entries1, *new_entries2;

  n = num_entries;

  mbr_array = new floatptr[n];
  for (i = 0; i < n; i++)
    mbr_array[i] = entries[i]->bounces;

  dist = mbr_split(mbr_array, &distribution);

  new_entries1 = new EntryPtr[capacity];
  new_entries2 = new EntryPtr[capacity];

  for (i = 0; i < capacity; i ++)
  {
    new_entries1[i] = new_one_entry(); new_entries1[i]->init(my_tree);
    new_entries1[i]->level = level;
    new_entries2[i] = new_one_entry(); new_entries2[i]->init(my_tree);
    new_entries2[i]->level = level;
  }

  for (i = 0; i < dist; i++)
    new_entries1[i]->set_equal_to(entries[distribution[i]]);

  for (i = dist; i < n; i++)
  {
    //		cout << i << endl;
    new_entries2[i-dist]->set_equal_to(entries[distribution[i]]);
  }

  //destroy arrays "entries" and "_sn->entries" ---
  for (i = 0; i < n; i++)
  {
    entries[i]->son_ptr = NULL; delete entries[i];
    _sn->entries[i]->son_ptr = NULL; delete _sn->entries[i];
  }
  delete [] entries;
  delete [] _sn->entries;
  //---

  entries = new_entries1;
  _sn->entries = new_entries2;

  num_entries = dist;
  _sn->num_entries = n - dist;

  delete [] mbr_array;
  delete [] distribution;
}
//------------------------------------------------------------
//writes the content of a node to a disk buffer page

void RTNode::write_to_buffer(char *_buffer)
{
  int i, j, s;

  // Level
  memcpy(_buffer, &level, sizeof(char));
  j = sizeof(char);

  // num_entries
  memcpy(&_buffer[j], &num_entries, sizeof(int));
  j += sizeof(int);

  s = entries[0]->get_size();
  for (i = 0; i < num_entries; i++)
  {
    entries[i]->write_to_buffer(&_buffer[j]);
    j += s;
  }
}

/*****************************************************************
this function performs a rank-inquiry query with linear function using
the breadth-first search.
para:
weight: an array (with size equal to the dimensionality) storing
  the query vector
qscore: the query score
rslt: the link list storing the ids of the returned tuples (disabled
for the time being)

Coded by Yufei Tao 05/04/02
*****************************************************************/

//void RTNode::rank_qry_inquiry(float *_weight, float _qscore, int *_rslt)
//{
//	my_tree->na[level]++;
//	if (level==0)
//	{
//		return; //this is only a simulation
//	}
//	for (int i=0; i<num_entries; i++)
//	{
//		float range[2];
//		my_tree->get_score_range(&(entries[i]), _weight, dimension, range); //get the score range for this entry
//		if (range[0]<=_qscore && _qscore<=range[1])
//		{
//			RTNode *succ = entries[i]->get_son();
//			succ->rank_qry_inquiry(_weight, _qscore, _rslt);
//			entries[i].del_son();
//		}
//	}
//}

/*****************************************************************
initiates an entry. Override this function to create an object of a class that inherits Entry

Coded by Yufei Tao 05/31/06
*****************************************************************/

Entry* RTNode::new_one_entry()
{
  Entry *e = new Entry();
  return e;
}

/*****************************************************************
initiates a node. Override this function to create an object of a class that inherits RTNode

Coded by Yufei Tao 05/31/06
*****************************************************************/

RTNode* RTNode::new_one_node()
{
  RTNode *nd = new RTNode();
  return nd;
}


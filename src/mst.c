/*
 * =====================================================================================
 *
 *       Filename:  mst.c
 *
 *    Description: minimum spanning tree using union find 
 *
 *        Version:  1.0
 *        Created:  18/03/2020 15:54:41
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D.G.), dfguan9@gmail.com
 *   Organization:  Harbin Institute of Technology
 *
 * =====================================================================================
 */

// refer to https://www.geeksforgeeks.org/kruskals-minimum-spanning-tree-algorithm-greedy-algo-2/
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>


#include "mst.h"
#include "kvec.h"
#include "sdict.h"


int _find(subset_t *subsets, int i, int *_stk) 
{ 
    // find root and make root as parent of i  
    // (path compression) 
	/*while (subsets[i].p != i) */
		/*i = subsets[i].p;*/
	/*return i;*/
	int j = 0, k;
	while (subsets[i].p != i) _stk[j++] = i, i = subsets[i].p;
	for (k = 0; k < j; ++k) subsets[_stk[k]].p = i;
	/*if (subsets[i].p != i)*/
        /*subsets[i].p = _find(subsets, subsets[i].p);*/
    return i; //subsets[_stk[0]].p
} 

void _union(subset_t *subsets, int  x, int y, int *_stk)
{
	int xroot = _find(subsets, x, _stk); 
    int yroot = _find(subsets, y, _stk); 
  
    // Attach smaller rank tree under root of high  
    // rank tree (Union by Rank) 
    if (subsets[xroot].r < subsets[yroot].r) 
        subsets[xroot].p = yroot; 
    else if (subsets[xroot].r > subsets[yroot].r) 
        subsets[yroot].p = xroot; 
    else { 
		// If ranks are same, then make one as root and  
		// increment its rank by one 
        subsets[yroot].p = xroot; 
        subsets[xroot].r++; 
    } 	
	return ;
}

// nv: number of vertex
int cmp_met(const void *a, const void *b)
{
	mst_edge_t *p = (mst_edge_t *)a;
	mst_edge_t *q = (mst_edge_t *)b;
	if (p->wt < q->wt) return 1;
	else if (p->wt > q->wt) return -1;
	else return 0;	
}

int _kmst(mst_edge_t *met, int n, int nv)
{
	qsort(met, n, sizeof(mst_edge_t), cmp_met);	
  
    // Allocate memory for creating V ssubsets 
    subset_t *subsets = (subset_t*) malloc( nv * sizeof(subset_t)); 
	int *_stk = malloc(nv * sizeof(int));  
    // Create V subsets with single elements 
	int i;
    for (i = 0; i < nv; ++i) 
    { 
        subsets[i].p = i; 
        subsets[i].r = 0; 
    } 
  
    // Number of edges to be taken is equal to n 
	for (i = 0; i < n; ++i) {
		int x = _find(subsets, met[i].s, _stk);	
		int y = _find(subsets, met[i].d, _stk);	
        if (x != y) 
        { 
            met[i].in_mst = 1;
			_union(subsets, x, y, _stk); 
        } 
	}
	if (subsets) free(subsets);
	if (_stk) free(_stk);
	return 0;
}

#ifdef MST_MAIN
int main(int argc, char *argv[])
{
	sdict_t *sd = sd_init();
	FILE *fp = fopen(argv[1], "r");		
	char v[1024], w[1024];
	float wt; 
	char e;
	kvec_t(mst_edge_t) mets;
	kv_init(mets);
	mst_edge_t tmp;
	while (fscanf(fp, "%c\t%s\t%s\t%f\n", &e, v, w, &wt)!=EOF) {
		int vid = sd_put(sd, v);
		int wid = sd_put(sd, w);
		tmp = (mst_edge_t){vid, wid, 0, 0, wt};
		kv_push(mst_edge_t, mets, tmp);	
		/*fprintf(stderr, "%s\t%s\n", v, w);*/
	}	
	_kmst(mets.a, mets.n, sd->n_seq);
	int i;
	mst_edge_t *in_m = mets.a; 
	for (i = 0; i < mets.n; ++i) if (in_m[i].in_mst) fprintf(stderr, "%s\t%s\t%f\n", sd->seq[in_m[i].s].name, sd->seq[in_m[i].d].name, in_m[i].wt*10000);
	kv_destroy(mets);
	return 0;
}
#endif


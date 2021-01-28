/*
 * =====================================================================================
 *
 *       Filename:  mst.h
 *
 *    Description: header of minimum spanning tree 
 *
 *        Version:  1.0
 *        Created:  18/03/2020 18:10:22
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D.G.), dfguan9@gmail.com
 *   Organization:  Harbin Institute of Technology
 *
 * =====================================================================================
 */


#ifndef __MST_H
#define __MST_H
typedef struct {
	int32_t p, r;
} subset_t;

typedef struct {
	uint32_t s, d, aux:30, in_mst:1;
	float wt;
} mst_edge_t;

int _kmst(mst_edge_t *met, int n, int nv);

#endif

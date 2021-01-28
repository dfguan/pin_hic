/*
 * =====================================================================================
 *
 *       Filename:  build_graph.h
 *
 *    Description:  header for build graph 
 *
 *        Version:  1.0
 *        Created:  19/11/2018 19:53:35
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#ifndef BUILD_GRAPH_H
#define BUILD_GRAPH_H

#ifdef __cplusplus
extern "C" {
#endif
int main_bldg(int argc, char *argv[], int ishic);
int buildg(char *fn, char *edge_fn, int min_wt, int use_sat, int norm, float min_mdw, int mlc, char *out_fn);
int buildg_hic(char *fn, char *edge_fn, int min_wt, int use_sat, int norm, float min_mdw, int mlc, char *out_fn, int use_nw, int amode, int igm, int usep, int use_mst);
int buildg_hic_cont(int which, char *fn, char *edge_fn, int min_wt, int use_sat, int norm, float min_mdw, int mlc, char *out_fn, int use_nw, int amode, int igm, int usep, int use_mst);
#ifdef __cplusplus
}
#endif

#endif

/*
 * =====================================================================================
 *
 *       Filename:  col_hic_lnks.h
 *
 *    Description:  header for col_hic_lnks
 *
 *        Version:  1.0
 *        Created:  19/11/2018 20:52:33
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#ifndef _COL_HIC_LNKS_H
#define _COL_HIC_LNKS_H

#include <stdint.h>
#ifdef __cplusplus 
extern "C" {
#endif
int main_hic_lnks(int argc, char *argv[]);
int col_hic_lnks(char *sat_fn, char **bam_fn, int n_bam, int min_mq, uint32_t win_s, int use_min_dist, int ca, int igm, char *out_fn);
#ifdef __cplusplus
}
#endif


#endif


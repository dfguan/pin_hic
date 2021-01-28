/*
 * =====================================================================================
 *
 *       Filename:  make_brk.h
 *
 *    Description:  make_brk header
 *
 *        Version:  1.0
 *        Created:  14/11/2019 15:28:46
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D.G.), dfguan9@gmail.com
 *   Organization:  Harbin Institute of Technology
 *
 * =====================================================================================
 */

#ifndef MAKE_BRK_H
#define MAKE_BRK_H

#ifdef __cplusplus
extern "C" {
#endif
int main_brks(int argc, char *argv[]);
int main_brks_10x(int argc, char *argv[]);
int mk_brks(char *sat_fn, char *bam_fn[], int n_bams, int min_mq, float min_rat, char *out_dir, char *prefix);
#ifdef __cplusplus
}
#endif

#endif

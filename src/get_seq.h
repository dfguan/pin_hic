/*
 * =====================================================================================
 *
 *       Filename:  get_seq.h
 *
 *    Description:  header
 *
 *        Version:  1.0
 *        Created:  20/11/2018 13:21:06
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#ifndef _GET_SEQ_H
#define _GET_SEQ_H
#include <stdint.h>

#ifdef __cplusplus 
extern "C" {
#endif
int get_seq(char *sat_fn, char *seq_fn, uint32_t min_l, char *out_fn);
int main_get_seq(int argc, char *argv[]);
#ifdef __cplusplus
}
#endif


#endif

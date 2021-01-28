/*
 * =====================================================================================
 *
 *       Filename:  cdict.h
 *
 *    Description:  header for cdict.c
 *
 *        Version:  1.0
 *        Created:  21/10/2018 10:08:54
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#ifndef _CDICT_H
#define _CDICT_H

#include <stdio.h>
#include <stdint.h>

typedef struct {
	char	*name; 
	float cnt;
	uint32_t ucnt;
	//uint32_t cnt:31, is_added:1;	
	uint32_t snp_n:31, is_l:1;
} cd_cnt_t;

typedef struct {
	size_t		n_cnt, m_cnt;
	size_t		lim;
	cd_cnt_t	*cnts;
	void		*h;
}cdict_t;

typedef struct {
	char	*name; 
    float cnt[4]; // 00 01 10 11 --, -+, +-, ++
    float ncnt;//could be enough normalized cnt
} cd_cnt2_t;

typedef struct {
	size_t		n_cnt, m_cnt;
	size_t		lim;
	cd_cnt2_t	*cnts;
	void		*h;
}cdict2_t;

#ifdef __cplusplus
extern "C" {
#endif
	void cd_init(cdict_t *c);
	void cd_destroy(cdict_t *c);
	void cd_norm(cdict_t *c);

	void cd_add(cdict_t *c, const char *name, uint32_t is_l, uint32_t snp_n, float w);
void cd_add2(cdict_t *c, const char *name, uint32_t is_l, float cnt, uint32_t snp_n);
void cd_set_lim(cdict_t *c, uint32_t n, uint32_t min_wt, float min_rat, int max_cand, int norm);
	//void cd_add(cdict_t *c, char *name, int is_l);
	void cd_sort(cdict_t *c);
void cd_filt(cdict_t *c, uint32_t n, float min_rat);
//cdict2 operations	
	void cd2_init(cdict2_t *c);
	void cd2_destroy(cdict2_t *c);
    void cd2_add(cdict2_t *c, uint32_t is_l1, const char *name, uint32_t is_l2, float cnt);
	void cd2_sort(cdict2_t *c);
	void cd2_set_lim(cdict2_t *c, uint32_t n, int cann);
	int32_t cd2_get(cdict2_t *c, char *name);
#ifdef __cplusplus
}
#endif






#endif




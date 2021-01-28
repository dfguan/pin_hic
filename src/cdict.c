/*
 * =====================================================================================
 *
 *       Filename:  cdict.c
 *
 *    Description:  counter dictionary 
 *
 *        Version:  1.0
 *        Created:  21/10/2018 10:05:56
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include "khash.h"
#include "cdict.h"
#include "utls.h"


KHASH_MAP_INIT_STR(str, uint32_t)
typedef khash_t(str) shash_t;

void cd_init(cdict_t *c)
{
	c->h = kh_init(str);
}

void cd2_init(cdict2_t *c)
{
	c->h = kh_init(str);
}

void cd2_destroy(cdict2_t *c)
{
	if (c) {
		uint32_t i;
		if (c->h) kh_destroy(str, (shash_t *)c->h);
		for ( i = 0; i < c->n_cnt; i++)  {
			if (c->cnts[i].name) 
				free(c->cnts[i].name);	
		}
		if (c->cnts) free(c->cnts);
	}
}

void cd_destroy(cdict_t *c)
{
	if (c) {
		uint32_t i;
		if (c->h) kh_destroy(str, (shash_t *)c->h);
		for ( i = 0; i < c->n_cnt; i += 2)  {
			if (c->cnts[i].is_l && c->cnts[i].name) 
				free(c->cnts[i].name);	
		}
		if (c->cnts) free(c->cnts);
	}
}

void cd_filt(cdict_t *c, uint32_t n, float min_rat)
{
	cdict_t *t;
	uint32_t i, j;
	for (i = 0; i < n; ++i)  {
		t = c + i;
		//maximum weight
		/*t->lim = t->n_cnt ? 1 : 0;*/
		if (!t->n_cnt || !t->lim)  continue;

		uint32_t max_wt = t->cnts[0].cnt;
		/*fprintf(stderr, "%u\t", max);*/
		/*if (max <= 2) continue;*/
		/*max_wt >>= 2;*/
		uint32_t sum_wt = 0;
		for(j = 0; j < t->lim; ++j) sum_wt += t->cnts[j].cnt;
	  	if ((float) max_wt / sum_wt < min_rat) t->lim = 0;	
		fprintf(stderr, "%d\t%d\t%d\n", t->lim, max_wt, sum_wt);
	}
}


void cd2_set_lim(cdict2_t *c, uint32_t n, int cann)
{
	cdict2_t *t;

	uint32_t i, j;

	for (i = 0; i < n; ++i)  {
		t = c + i;
		//maximum weight
		/*t->lim = t->n_cnt ? 1 : 0;*/
		t->lim = cann < t->n_cnt? cann : t->n_cnt;

		/*if (!t->n_cnt)  continue;*/
	}
}


void cd_set_lim(cdict_t *c, uint32_t n, uint32_t min_wt, float min_rat, int max_cand, int norm)
{
	cdict_t *t;

	uint32_t i, j;

	for (i = 0; i < n; ++i)  {
		t = c + i;
		//maximum weight
		/*t->lim = t->n_cnt ? 1 : 0;*/
		t->lim = 0;

		if (!t->n_cnt)  continue;
		if (~max_cand) {
			t->lim = max_cand <  t->n_cnt ? max_cand : t->n_cnt;
		} else {
		
			/*fprintf(stderr, "%u\t", max);*/
			/*if (max <= 2) continue;*/
			/*max_wt >>= 2;*/
			
			for(j = 0; j < t->n_cnt; ++j) {
				if (t->cnts[j].cnt >= min_wt) 
					++t->lim;
				else
					break;
				/*if (t->cnts[j].cnt >= 0.1 * max_wt)*/
					/*++t->lim;*/
				/*else*/
					/*break;				*/
			}
		
		}
		//the first candiate is not significantly better than the second 
		if (!norm && t->n_cnt > 1 && norm_cdf(t->cnts[0].cnt, 0.5, t->cnts[1].cnt + t->cnts[0].cnt) <= min_rat && strcmp(t->cnts[0].name, t->cnts[1].name)) t->lim = 0; 
	}
}

int cmp_cnt2(const void *a, const void *b)
{
	cd_cnt2_t *f = (cd_cnt2_t *)a;
	cd_cnt2_t *h = (cd_cnt2_t *)b;
	if (f->ncnt > h->ncnt) return -1;
	else if (f->ncnt == h->ncnt) return 0; 
	else return 1;
}
void cd2_sort(cdict2_t *c)
{
    if (c) {
        qsort(c->cnts, c->n_cnt, sizeof(cd_cnt2_t), cmp_cnt2);
    }

}
int32_t cd2_get(cdict2_t *c, char *name)
{
	shash_t *h = (shash_t *)c->h;
	khint_t k;
	/*if (h) fprintf(stderr, "cd add");*/
	k = kh_get(str, h, name);
	return (k == kh_end(h)) ? -1 : kh_val(h, k);
}
void cd2_add(cdict2_t *c, uint32_t is_l1, const char *name, uint32_t is_l2, float cnt)
{
	shash_t *h = (shash_t *)c->h;
	khint_t k;
	int absent;
	/*if (h) fprintf(stderr, "cd add");*/
	k = kh_put(str, h, name, &absent);
	if (absent) {
		/*fprintf(stderr, "%u\n", c->n_cnt);*/
		if (c->n_cnt == c->m_cnt) {
			c->m_cnt = c->m_cnt ? c->m_cnt << 1 : 16;
			cd_cnt2_t *ncnts = calloc(c->m_cnt, sizeof(cd_cnt2_t));	
			if (c->cnts) memcpy(ncnts, c->cnts, sizeof(cd_cnt2_t) * c->n_cnt);
			if (c->cnts) free(c->cnts);
			c->cnts = ncnts;
		}
		kh_key(h, k) = c->cnts[c->n_cnt].name = strdup(name); //
		c->cnts[c->n_cnt].cnt[is_l1<<1 | is_l2] = cnt;
		kh_val(h, k) = c->n_cnt++;  
	} else {
		uint32_t ind = kh_val(h, k);
		/*fprintf(stderr, "%u\t %s exist\n", ind, name);*/
		c->cnts[ind].cnt[is_l1 << 1 | is_l2] = cnt;
	}
}
void cd_add2(cdict_t *c, const char *name, uint32_t is_l, float cnt, uint32_t snp_n)
{
	shash_t *h = (shash_t *)c->h;
	khint_t k;
	int absent;
	/*if (h) fprintf(stderr, "cd add");*/
	k = kh_put(str, h, name, &absent);
	if (absent) {
		/*fprintf(stderr, "%u\n", c->n_cnt);*/
		if (c->n_cnt == c->m_cnt) {
			c->m_cnt = c->m_cnt ? c->m_cnt << 1 : 16;
			cd_cnt_t *ncnts = calloc(c->m_cnt, sizeof(cd_cnt_t));	
			if (c->cnts) memcpy(ncnts, c->cnts, sizeof(cd_cnt_t) * c->n_cnt);
			if (c->cnts) free(c->cnts);
			c->cnts = ncnts;
		}
		kh_key(h, k) = c->cnts[c->n_cnt].name = strdup(name); //
		kh_val(h, k) = c->n_cnt;  
		c->cnts[c->n_cnt].is_l = 0;
		
		c->cnts[c->n_cnt+1].name = c->cnts[c->n_cnt].name; // init two 
		c->cnts[c->n_cnt+1].is_l = 1;	
		
		c->cnts[c->n_cnt | is_l].cnt = cnt;
		/*c->cnts[c->n_cnt | is_l].intcnt = cnt;*/
		c->cnts[c->n_cnt | is_l].snp_n = snp_n;
		c->n_cnt += 2;
	} else {
		uint32_t ind = kh_val(h, k);
		/*fprintf(stderr, "%u\t %s exist\n", ind, name);*/
		c->cnts[ind | is_l].cnt = cnt;
		/*c->cnts[ind | is_l].intcnt = cnt;*/
		c->cnts[ind | is_l].snp_n = snp_n;
	}
		
	/*fprintf(stderr, "%s\t%s\t%d\t%d\n", name, c->cnts[kh_val(h,k)].name, kh_val(h, k), absent);*/
	/*if (absent) {*/
		/*sd_seq_t *s;*/

		/*if (d->n_seq == d->m_seq) {*/
			/*d->m_seq = d->m_seq? d->m_seq<<1 : 16;*/
			/*d->seq = (sd_seq_t*)realloc(d->seq, d->m_seq * sizeof(sd_seq_t));*/
		/*}*/
		/*s = &d->seq[d->n_seq];*/
		/*kh_key(h, k) = s->name = strdup(name);*/
		/*kh_val(h, k) = d->n_seq++;*/
	/*} // TODO: test if len is the same;*/
	/*return kh_val(h, k);*/
}
/* 
void cd_add2(cdict_t *c, const char *name, uint32_t is_l, uint32_t cnt)
{
	shash_t *h = (shash_t *)c->h;
	khint_t k;
	int absent;
	k = kh_put(str, h, name, &absent);
	if (absent) {
		if (c->n_cnt == c->m_cnt) {
			c->m_cnt = c->m_cnt ? c->m_cnt << 1 : 16;
			cd_cnt_t *ncnts = calloc(c->m_cnt, sizeof(cd_cnt_t));	
			if (c->cnts) memcpy(ncnts, c->cnts, sizeof(cd_cnt_t) * c->n_cnt);
			if (c->cnts) free(c->cnts);
			c->cnts = ncnts;
		}
		kh_key(h, k) = c->cnts[c->n_cnt].name = strdup(name); //
		kh_val(h, k) = c->n_cnt;  
		c->cnts[c->n_cnt].is_l = 0;
		
		c->cnts[c->n_cnt+1].name = c->cnts[c->n_cnt].name; // init two 
		c->cnts[c->n_cnt+1].is_l = 1;	
		
		c->cnts[c->n_cnt | is_l].cnt = cnt;
		c->cnts[c->n_cnt | is_l].snp_n = 1;
		c->n_cnt += 2;
	} else {
		uint32_t ind = kh_val(h, k);
		c->cnts[ind | is_l].cnt = cnt;
		c->cnts[ind | is_l].snp_n = 1;
	}
} 
*/
void cd_add(cdict_t *c, const char *name, uint32_t is_l, uint32_t snp_n, float w)
{
	shash_t *h = (shash_t *)c->h;
	khint_t k;
	int absent;
	/*if (h) fprintf(stderr, "cd add");*/
	k = kh_put(str, h, name, &absent);
	if (absent) {
		/*fprintf(stderr, "%u\n", c->n_cnt);*/
		if (c->n_cnt == c->m_cnt) {
			c->m_cnt = c->m_cnt ? c->m_cnt << 1 : 16;
			cd_cnt_t *ncnts = calloc(c->m_cnt, sizeof(cd_cnt_t));	
			if (c->cnts) memcpy(ncnts, c->cnts, sizeof(cd_cnt_t) * c->n_cnt);
			if (c->cnts) free(c->cnts);
			c->cnts = ncnts;
		}
		kh_key(h, k) = c->cnts[c->n_cnt].name = strdup(name); //
		kh_val(h, k) = c->n_cnt; // necessary ? 
		c->cnts[c->n_cnt].is_l = 0;
		
		c->cnts[c->n_cnt+1].name = c->cnts[c->n_cnt].name; // init two 
		c->cnts[c->n_cnt+1].is_l = 1;	
		
		c->cnts[c->n_cnt | is_l].cnt = w;
		++c->cnts[c->n_cnt | is_l].ucnt;
		c->cnts[c->n_cnt | is_l].snp_n = snp_n;
		c->n_cnt += 2;
	} else {
		uint32_t ind = kh_val(h, k);
		c->cnts[ind | is_l].cnt += w;
		++c->cnts[ind | is_l].ucnt;
		c->cnts[ind | is_l].snp_n = snp_n;
	}
		
	/*if (absent) {*/
		/*sd_seq_t *s;*/

		/*if (d->n_seq == d->m_seq) {*/
			/*d->m_seq = d->m_seq? d->m_seq<<1 : 16;*/
			/*d->seq = (sd_seq_t*)realloc(d->seq, d->m_seq * sizeof(sd_seq_t));*/
		/*}*/
		/*s = &d->seq[d->n_seq];*/
		/*kh_key(h, k) = s->name = strdup(name);*/
		/*kh_val(h, k) = d->n_seq++;*/
	/*} // TODO: test if len is the same;*/
	/*return kh_val(h, k);*/
}

void cd_norm(cdict_t *c)
{
	if (c) {
		uint32_t i;
		cd_cnt_t *t = c->cnts;
		for ( i = 0; i < c->n_cnt; ++i) 
			if (t[i].snp_n)  t[i].cnt = t[i].cnt / t[i].snp_n;
			else t[i].cnt = 0;	
		
	}	
}

int cmp_cnt(const void *a, const void *b) 
{
	cd_cnt_t *f = (cd_cnt_t *)a;
	cd_cnt_t *h = (cd_cnt_t *)b;
	if (f->cnt > h->cnt) return -1;
	else if (f->cnt == h->cnt) {
		if (f->snp_n > h->snp_n) return 1;
		else if (f->snp_n == h->snp_n) return 0;
		else return -1;
	}
	else return 1;
}

void cd_sort(cdict_t *c)
{
	if (c) {
		qsort(c->cnts, c->n_cnt, sizeof(cd_cnt_t), cmp_cnt);
	}	
}




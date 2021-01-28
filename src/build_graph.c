/*
 * =====================================================================================
 *
 *       Filename:  build_graph.c
 *
 *    Description:  build graph with links information
 *
 *        Version:  1.0
 *        Created:  19/11/2018 19:39:30
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "bed.h"
#include "graph.h"
#include "sdict.h"
#include "cdict.h"
#include "utls.h"
#include "mst.h"
#include "kvec.h"



sdict_t *col_ctgs(char *fn)
{
	bed_file_t *bf = bed_open(fn);
	if (!bf) return 0;
	sdict_t *ctgs = sd_init();		
	bed_rec_t r;
	while (bed_read(bf, &r) >= 0) 
		sd_put2(ctgs, r.ctgn, r.len, r.le, r.rs, r.l_snp_n, r.r_snp_n);
	bed_close(bf);
	return ctgs;
}

sdict_t *col_ctgs_from_graph(graph_t *g)
{
	sdict_t *ctgs = sd_init();
	asm_t *ca = &g->as.asms[g->as.casm];
	path_t *pt = g->pt.paths;
	uint32_t n = ca->n;
	uint32_t i, j;
	for ( i = 0; i < n; ++i) {
		//get path length
		/*uint32_t path_len = 0;*/
		/*for (j = 0; j < pt[ca->pn[i]>>1].n; ++j) */
			/*path_len += ctgs->seq[sd_get(ctgs, pt[ca->pn[i]>>1].name)].len;	*/
		/*path_len += (pt[ca->pn[i]>>1].n - 1) * 200;	 //gap size = 200	*/
		sd_put4(ctgs, pt[ca->pn[i]>>1].name, pt[ca->pn[i]>>1].len, 0, 0, 0, 0, pt[ca->pn[i]>>1].is_circ);
	} 
	for (uint32_t i = 0; i < ctgs->n_seq; ++i) fprintf(stderr, "ctgs: %s %d\n", ctgs->seq[i].name, ctgs->seq[i].is_circ);
	return ctgs;
}

int get_links_hic(char *links_fn, cdict2_t *cds, sdict_t *ctgs)
{	
	bed_file_t *bf = bed_open(links_fn);
	if (!bf) return 0;
	lnk_rec_t r;
	uint32_t line_n = 0;
	while (lnk_read(bf, &r) >= 0) {
		uint32_t ind1;
		if (r.is_l) 
			ind1 = sd_put2(ctgs, r.ctgn, 0, 0, 0, r.llen, 0);		
		else
			ind1 = sd_put2(ctgs, r.ctgn, 0, 0, 0, 0, r.llen);		
		if (r.is_l2)
			sd_put2(ctgs, r.ctgn2, 0, 0, 0, r.rlen, 0);		
		else
			sd_put2(ctgs, r.ctgn2, 0, 0, 0, 0, r.rlen);		
		/*uint32_t ind2 = sd_put2(ctgs, r.ctgn, 0, 0, 0, r.llen, r.rlen);		*/
		line_n += 1;
		cd2_add(&cds[ind1], r.is_l, r.ctgn2, r.is_l2, r.fwt);	//this has been normalized	
	} 
	bed_close(bf);
	return 0;	
}
int get_links(char *links_fn, cdict_t *cds, sdict_t *ctgs)
{	
	bed_file_t *bf = bed_open(links_fn);
	if (!bf) return 0;
	lnk_rec_t r;
	uint32_t line_n = 0;
	while (lnk_read(bf, &r) >= 0) {
		uint32_t ind1;
		if (r.is_l) 
			ind1 = sd_put2(ctgs, r.ctgn, 0, 0, 0, r.rlen, 0); //wrong length add to contig
		else
			ind1 = sd_put2(ctgs, r.ctgn, 0, 0, 0, 0, r.rlen);		
		if (r.is_l2)
			sd_put2(ctgs, r.ctgn2, 0, 0, 0, r.llen, 0);		
		else
			sd_put2(ctgs, r.ctgn2, 0, 0, 0, 0, r.llen);		
		/*uint32_t ind2 = sd_put2(ctgs, r.ctgn, 0, 0, 0, r.llen, r.rlen);		*/
		line_n += 1;
		cd_add2(&cds[ind1<<1|r.is_l], r.ctgn2, r.is_l2, r.fwt, r.llen);	//this has been normalized	
	} 
	bed_close(bf);
	return 0;	
}
int anothernorm(cdict_t *cds, sdict_t *ctgs)
{
	uint32_t n_cds = ctgs->n_seq << 1;
	uint32_t i;
	cdict_t *c ;
	for ( i = 0; i < n_cds; ++i) {
		char *name1 = ctgs->seq[i>>1].name;
		uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;
		uint32_t j;
		c = cds + i;
		uint32_t icnt;
		for (j = 0; j < c->n_cnt; ++j) {
			char *name2 = c->cnts[j].name; 
			icnt = c->cnts[j].cnt ; 
			uint32_t snp2 = c->cnts[j].is_l ? ctgs->seq[sd_get(ctgs, name2)].l_snp_n:ctgs->seq[sd_get(ctgs,name2)].r_snp_n;
			fprintf(stderr, "%s\t%c\t%s\t%c\t%u\t%u\t%u\t%lf\n", name1, i&1?'+':'-', name2, c->cnts[j].is_l?'+':'-', icnt, snpn, snp2, 100000.0*(double)icnt/(snp2*snpn));
		}
	}
	return 0;
}

int print_cdict2(cdict2_t *cds, sdict_t *ctgs)
{
	uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c ;
	for ( i = 0; i < n_cds; ++i) {
		char *name1 = ctgs->seq[i].name;
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		uint32_t j;
		c = cds + i;
		uint32_t icnt;
        
		for (j = 0; j < c->n_cnt; ++j) {
            /*fprintf(stderr, "%s\n", c->cnts[j].name);*/
			char *name2 = c->cnts[j].name; 
            uint32_t ctg2_idx = sd_get(ctgs, name2);
			/*fprintf(stderr, "ctg_idx: %u\n", ctg2_idx);*/
			icnt = c->cnts[j].cnt[0] + c->cnts[j].cnt[1] + c->cnts[j].cnt[2] + c->cnts[j].cnt[3]; 
            /*c->cnts[j].ncnt = (float) icnt / ctgs->seq[ctg2_idx].len;*/
            /*c->cnts[j].ncnt = (float) icnt / (ctgs->seq[ctg2_idx].l_snp_n + ctgs->seq[ctg2_idx].r_snp_n);*/
            /*uint32_t z;*/
            /*for ( z = 0; z < 4; ++z) c->cnts[j].fcnt[z] = (float) c->cnts[j].cnt[z]/(z >> 1 ? ctgs->seq[i].l_snp_n : ctgs->seq[i].r_snp_n) / ( z & 0x1 ? ctgs->seq[ctg2_idx].l_snp_n : ctgs->seq[ctg2_idx].r_snp_n);  */
			/*uint32_t snp2 = c->cnts[j].is_l ? ctgs->seq[sd_get(ctgs, name2)].l_snp_n:ctgs->seq[sd_get(ctgs,name2)].r_snp_n;*/
			fprintf(stderr, "%s\t%s\t%f\t%f\t%f\t%f\t%f\t%u\n", name1, name2, c->cnts[j].ncnt, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3], ctgs->seq[ctg2_idx].len);
		}
	}
	return 0;
}

int print_cdict(cdict_t *cds, sdict_t *ctgs)
{
	uint32_t n_cds = ctgs->n_seq << 1;
	uint32_t i;
	cdict_t *c ;
	for ( i = 0; i < n_cds; ++i) {
		char *name1 = ctgs->seq[i>>1].name;
		uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;
		uint32_t j;
		c = cds + i;
		float icnt;
		uint32_t intcnt;
		for (j = 0; j < c->n_cnt; ++j) {
			char *name2 = c->cnts[j].name; 
			icnt = c->cnts[j].cnt ; 
			intcnt = 1; //c->cnts[j].intcnt;
			uint32_t snp2 = c->cnts[j].is_l ? ctgs->seq[sd_get(ctgs, name2)].l_snp_n:ctgs->seq[sd_get(ctgs,name2)].r_snp_n;
			fprintf(stderr, "%s\t%c\t%s\t%c\t%f\t%u\t%u\t%u\n", name1, i&1?'+':'-', name2, c->cnts[j].is_l?'+':'-', icnt, intcnt, snpn, snp2);
		}
	}
	return 0;
}

graph_t *build_graph(cdict_t *cds, sdict_t *ctgs)
{
	graph_t *g = graph_init();
	
	uint32_t n = ctgs->n_seq << 1;
	uint32_t i, j;
	//create nodes
	for ( i = 0; i < ctgs->n_seq; ++i) 
		add_node(g, ctgs->seq[i].name, 0, ctgs->seq[i].len, ctgs->seq[i].is_circ);
	//create edges	
	for ( i = 0; i < n; ++i) {
		char *name1 = ctgs->seq[i>>1].name;
		uint8_t is_l = i & 1;	
		cdict_t *c = cds + i;	
		if (!c->n_cnt) continue;	
		for (j = 0; j < c->lim; ++j) {
			char *name2 = c->cnts[j].name;
			if (strcmp(name1, name2) == 0) continue;
			//hsortand shaking
			/*fprintf(stderr, "try hand shaking\n");*/
			uint32_t ind = sd_get(ctgs, name2) << 1 | c->cnts[j].is_l;
			uint32_t k;
			uint8_t hand_shaking = 0;
			float ocnt = 0;
			for ( k = 0; k < cds[ind].lim; ++k) {
					if (strcmp(name1, cds[ind].cnts[k].name) == 0 && cds[ind].cnts[k].is_l == is_l) {
						hand_shaking = 1;
						ocnt = cds[ind].cnts[k].cnt;
						break;
					}
			}	
			/*if (hand_shaking) fprintf(stderr, "I hand shaking\n");*/
			ocnt = 1;
			if (hand_shaking) add_dedge(g, name1, is_l, name2, c->cnts[j].is_l, ocnt * c->cnts[j].cnt);	 //kinda residule cause index of name1 is the same as its index in ctgs but user doesn't know how the node is organized so better keep this.
		}		
	}	
	return g;
}

float dif(float *a, int n)
{
	int i;
	int maxi = 0; float maxv = a[0];
	for (i = 1; i < n; ++i) if (a[i] > maxv) maxv = a[i], maxi=i;
	float parv = a[n-1-maxi];//only works for 4 
	++maxv, ++parv; //in case parv is zero
	return 2.0;
	/*fprintf(stderr, "DIF\t%f\t%f\t%f\n", maxv, parv, maxv/parv);*/
	return maxv/parv;
}

int norm_links(cdict2_t *cds, sdict_t *ctgs, int norm, int igm, int usep, int min_wt)
{
	uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c ;
	for ( i = 0; i < n_cds; ++i) {
		/*char *name1 = ctgs->seq[i].name;*/
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		float len1 = ctgs->seq[i].len;
		uint32_t j;
		c = cds + i;
		float icnt;
        
		for (j = 0; j < c->n_cnt; ++j) {
            /*fprintf(stderr, "%s\n", c->cnts[j].name);*/
			char *name2 = c->cnts[j].name; 
            uint32_t ctg2_idx = sd_get(ctgs, name2);
			float len2 = ctgs->seq[ctg2_idx].len; 
			icnt = c->cnts[j].cnt[0] + c->cnts[j].cnt[1] + c->cnts[j].cnt[2] + c->cnts[j].cnt[3]; 
            /*c->cnts[j].ncnt = norm ? (float) icnt / ctgs->seq[ctg2_idx].len : (float) icnt;*/
            /*c->cnts[j].ncnt = norm ? icnt / (len1/2 + len2/2) : icnt;*/
			float mul = 2.0;
			if (igm) mul = 3;
			float div = len1 + len2; //float range 3.4E+38
			if (usep) mul *= mul, div = len1 * len2; 
			uint32_t k, ltmin = 0;
			for (k = 0; k < 4; ++k) 
				if (c->cnts[j].cnt[k] > min_wt) ltmin=1, c->cnts[j].cnt[k] = c->cnts[j].cnt[k] * mul / div;
				else c->cnts[j].cnt[k] = 0;

			c->cnts[j].ncnt = 0.0;
			if (ltmin) {
				for (k = 0; k < 4; ++k) {
					if (c->cnts[j].ncnt < c->cnts[j].cnt[k])
						c->cnts[j].ncnt = c->cnts[j].cnt[k];
				}
			}
			/*c->cnts[j].ncnt = 0.0;*/
			/*if (ltmin)  c->cnts[j].ncnt = icnt * mul / div; //require maximum larger than min_wt*/
			
			/*int sel = igm << 1 | usep;*/
			/*if (sel == 0) {*/
					/*c->cnts[j].ncnt = icnt * 2 / (len1 + len2);*/
			
			/*} else if (sel == 1) {*/
					/*c->cnts[j].ncnt = icnt * 4 / len1 / len2;*/
			/*} else if (sel == 2) {*/
					/*c->cnts[j].ncnt = icnt * 3 / len1 + len2;*/

			/*} else if (sel == 3) {*/
					/*c->cnts[j].ncnt = icnt * 9 / len1 / len2;*/
			
			
			/*}*/
			/*switch ((igm << 1)|usep) */
			/*{*/
				/*case 0:*/
					/*break;*/
				/*case 1:*/
					/*break;*/
				/*case 2:*/
					/*break;*/
				/*case 3:*/
					/*break;*/
			/*}*/
			/*else*/
				/*c->cnts[j].ncnt *= (dif(c->cnts[j].cnt, 4) - 1);*/

				/*c->cnts[j].ncnt = icnt / (len1/3 + len2/3);*/

            /*c->cnts[j].ncnt = (float) icnt / (ctgs->seq[ctg2_idx].l_snp_n + ctgs->seq[ctg2_idx].r_snp_n);*/
            /*uint32_t z;*/
            /*for ( z = 0; z < 4; ++z) c->cnts[j].fcnt[z] = (float) c->cnts[j].cnt[z]/(z >> 1 ? ctgs->seq[i].l_snp_n : ctgs->seq[i].r_snp_n) / ( z & 0x1 ? ctgs->seq[ctg2_idx].l_snp_n : ctgs->seq[ctg2_idx].r_snp_n);  */
			/*uint32_t snp2 = c->cnts[j].is_l ? ctgs->seq[sd_get(ctgs, name2)].l_snp_n:ctgs->seq[sd_get(ctgs,name2)].r_snp_n;*/
			/*fprintf(stderr, "%s\t%c\t%s\t%c\t%u\t%u\t%u\t%lf\n", name1, i&1?'+':'-', name2, c->cnts[j].is_l?'+':'-', icnt, snpn, snp2, 100000.0*(double)icnt/(snp2*snpn));*/
		}
	}
	return 0;
}

int norm_links2(cdict_t *cds, sdict_t *ctgs, int norm, int igm, int usep, int min_wt)
{
	uint32_t n_cds = ctgs->n_seq << 1;
	uint32_t i;
	cdict_t *c ;
	for ( i = 0; i < n_cds; ++i) {
		/*char *name1 = ctgs->seq[i].name;*/
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		float len1 = ctgs->seq[i>>1].len;
		uint32_t j;
		c = cds + i;
		float icnt;
        
		for (j = 0; j < c->n_cnt; ++j) {
            /*fprintf(stderr, "%s\n", c->cnts[j].name);*/
			char *name2 = c->cnts[j].name; 
            uint32_t ctg2_idx = sd_get(ctgs, name2);
			float len2 = ctgs->seq[ctg2_idx].len; 
			icnt = c->cnts[j].cnt;
			float mul = 2.0;
			if (igm) mul = 3;
			float div = len1 + len2; //float range 3.4E+38
			if (usep) mul *= mul, div = len1 * len2; 
			
			c->cnts[j].cnt = 0.0;
			if (icnt > min_wt)  c->cnts[j].cnt = icnt * mul / div; //require maximum larger than min_wt
		}
	}
	return 0;
}
int det_ori(float *ws, int isf, int min_wt, float min_mdw, int amode)
{
	float l, sl;
	l = ws[0];
	sl = 0;
	int i, maxi = 0;
	for (i = 1; i < 4; ++i ) 
		if (ws[i] >= l) 
				maxi = i, sl = l, l = ws[i];  
		else if (ws[i] > sl) sl = ws[i];
	/*if (l <= min_wt) return -1;*/
	if (amode) {
		if (!isf) return -1;
		else if (norm_cdf(l, 0.5, sl + l) <= min_mdw) return -1;
		else return maxi;	
	} else {
		if (sl != l)
			return maxi;
		else
			return -1;
	}
	/*if (norm_cdf(l, 0.5, sl + l) <= 0.95) */
			/*return -1; */
	/*else*/
		   /*return maxi;	*/
}


int cmp_v(const void *a, const void *b)
{
	mst_edge_t *p = (mst_edge_t *)a;
	mst_edge_t *q = (mst_edge_t *)b;
	if (p->s < q->s) return -1;
	else if (p->s > q->s) return 1;
	else return 0;	
}

graph_t *nns_mst(cdict2_t *cds, sdict_t *ctgs)
{
	graph_t *g = graph_init();
	
	uint32_t i, j;
	//create nodes
	for ( i = 0; i < ctgs->n_seq; ++i) {
		add_node(g, ctgs->seq[i].name, 0, ctgs->seq[i].len, ctgs->seq[i].is_circ);
	}
	//create edges
	kvec_t(mst_edge_t) mets;
	kv_init(mets);
	mst_edge_t tmp;
	for ( i = 0; i < ctgs->n_seq; ++i) {
		char *name1 = ctgs->seq[i].name;
		uint32_t len1 = ctgs->seq[i].len;
		cdict2_t *c = cds + i;	
		int isf = 1;
		for (j = 0; j < c->lim; ++j) {
			char *name2 = c->cnts[j].name;
			if (strcmp(name1, name2) == 0) continue;
			//hsortand shaking
			/*fprintf(stderr, "try hand shaking\n");*/
			uint32_t ind = sd_get(ctgs, name2);
			uint32_t k;
			uint8_t hand_shaking = 0;
			uint32_t is_l, is_l2;
			uint32_t len2 = ctgs->seq[ind].len;
			for ( k = 0; k < cds[ind].lim; ++k) {
					if (strcmp(name1, cds[ind].cnts[k].name) == 0) {
						hand_shaking = 1;
						break;
					}
			}	
			if (!hand_shaking) continue;
			if (i < ind) {
				tmp = (mst_edge_t){i, ind, j, 0, c->cnts[j].ncnt};
				kv_push(mst_edge_t, mets, tmp);
			}
			//vertex degrees

				/*fprintf(stderr, "E\t%s\t%s\t%f\n",name1, name2, c->cnts[j].ncnt);*/
			/*int idx;*/
			/*if (~(idx = det_ori(c->cnts[j].cnt, isf, min_wt, min_mdw, amode))) */
			
				// if build in accurate mode 
			/*if (~(idx = det_ori(c->cnts[j].cnt)))*/
			/*uint32_t hh = c->cnts[j].cnt[0];*/
			/*uint32_t ht = c->cnts[j].cnt[1];*/
			/*uint32_t th = c->cnts[j].cnt[2];*/
			/*uint32_t tt = c->cnts[j].cnt[3];*/
			/*uint32_t tl = c->cnts[j].cnt[2] + c->cnts[j].cnt[3];*/
			/*uint32_t hd = c->cnts[j].cnt[0] + c->cnts[j].cnt[1];	*/
			/*uint32_t hd2 = c->cnts[j].cnt[0] + c->cnts[j].cnt[2];*/
			/*uint32_t tl2 = c->cnts[j].cnt[1] + c->cnts[j].cnt[3];*/
			/*if (tl != hd) */
					/*is_l = tl > hd ? 0 : 1; */
			/*else */
					/*continue;*/
			/*if (hd2 != tl2) */
					/*is_l2 = tl2 > hd ? 0 : 1; 	*/
			/*else */
					/*continue;*/
			//when min_wt == -1; don't normalize weight 
			/*is_l = idx >> 1, is_l2 = idx & 1, add_dedge(g, name1, is_l, name2, is_l2, ~min_wt?c->cnts[j].cnt[idx] / (len1 / 2 + len2 / 2) : c->cnts[j].cnt[idx]);	 //kinda residule cause index of name1 is the same as its index in ctgs but user doesn't know how the node is organized so better keep this.*/
			/*is_l = idx >> 1, is_l2 = idx & 1, add_dedge(g, name1, is_l, name2, is_l2, use_nw?c->cnts[j].ncnt: c->cnts[j].cnt[idx] * 2/(len1 + len2));	 //kinda residule cause index of name1 is the same as its index in ctgs but user doesn't know how the node is organized so better keep this.*/
			/*isf = 0;*/
		}		
	}	
	mst_edge_t *in_mst = mets.a;
	_kmst(mets.a, mets.n, ctgs->n_seq);	
	for (i = 0; i < ctgs->n_seq; ++i) ctgs->seq[i].snp_n = 0;
	for (i = 0; i < mets.n; ++i) {
		if (in_mst[i].in_mst) {
			uint32_t v = in_mst[i].s;
			uint32_t w = in_mst[i].d;
			++ctgs->seq[v].snp_n, ++ctgs->seq[w].snp_n;
		}
	}
	// filter out the vertices with more than 2 degrees			
	uint64_t *idx = (uint64_t *)calloc(ctgs->n_seq, sizeof(uint64_t));
	for (i = 0, j = 0; i < mets.n; ++i) {
		if (in_mst[i].in_mst) {
			uint32_t v = in_mst[i].s;
			uint32_t w = in_mst[i].d;
			if (ctgs->seq[v].snp_n > 2 || ctgs->seq[w].snp_n > 2) 
				in_mst[i].in_mst = 0, ++idx[v], ++idx[w];
			else
				in_mst[j++] = in_mst[i]; 
		} 	
	}
	for (i = 0; i < ctgs->n_seq; ++i) ctgs->seq[i].snp_n -= idx[i];
	mets.n = j;
	for (i = 0; i < mets.n; ++i)  fprintf(stderr, "E\t%s\t%s\t%u\t%u\t%f\n", ctgs->seq[in_mst[i].s].name, ctgs->seq[in_mst[i].d].name, ctgs->seq[in_mst[i].s].len, ctgs->seq[in_mst[i].d].len, in_mst[i].wt*10000);
	// double the edges and index
	for (i = 0; i < j; ++i) {
		uint32_t v = mets.a[i].s;
		uint32_t w = mets.a[i].d;
		uint32_t aux = mets.a[i].aux;
		tmp = (mst_edge_t){w, v, aux, 0, 0};
		kv_push(mst_edge_t, mets, tmp);
	}	
	qsort(mets.a, mets.n, sizeof(mst_edge_t), cmp_v);	
	//index 
	memset(idx, 0, ctgs->n_seq * sizeof(uint64_t));
	in_mst = mets.a; // mets.a may have been changed
	for (i = 1, j = 0; i <= mets.n; ++i) {
		if (i==mets.n || in_mst[i].s != in_mst[j].s) {
			/*fprintf(stderr, "s: %u d: %u\n", in_mst[j].s, in_mst[j].d);*/
			idx[in_mst[j].s] = ((uint64_t)j << 32 | (i -j)) << 1 | 0;
			j = i;	
		}
	}
		
	// search for 1 degree node and use dynamic programming to find paths
	typedef struct {
		float	sc;
		int		bt;
	} vsc_t;
	vsc_t *vts = (vsc_t *)calloc(ctgs->n_seq<<1, sizeof(vsc_t)); //0 h 1 t
	kvec_t(uint32_t) endv;
	kv_init(endv);

	for (i = 0; i < ctgs->n_seq; ++i) {
		if (idx[i] & 0x1) continue; //visited before
		if (ctgs->seq[i].snp_n == 1) {
			//start node	
			// looking for next node
			uint32_t zprev = -1, z = i, pos;
			/*fprintf(stderr, "S: %s zprev:%d z:%d\n", ctgs->seq[z].name, zprev, z);*/
			while (1) {
				for (pos = idx[z]>>33; pos < (idx[z]>>33) + (uint32_t)(idx[z]>>1); ++pos) 
					if (in_mst[pos].d != zprev) break;
				uint32_t znext = in_mst[pos].d;
				/*fprintf(stderr, "M: %s zprev:%d z:%d pos: %d %d znext: %d\n", ctgs->seq[z].name, zprev, z,idx[z]>>1, idx[z+1]>>1, znext);*/
				uint32_t aux = in_mst[pos].aux;
				uint32_t lenz = ctgs->seq[z].len;
				uint32_t lenznext = ctgs->seq[znext].len;
				float hh, ht, th, tt, ort_cnt[4];
				if (in_mst[pos].in_mst) { // 
					float *_cnts = cds[z].cnts[aux].cnt;		
					/*fprintf(stderr, "%s\t%s\t%f\t%f\t%f\t%f\n", ctgs->seq[z].name, ctgs->seq[znext].name, _cnts[0], _cnts[1], _cnts[2], _cnts[3]);*/
					ort_cnt[0] = _cnts[0], ort_cnt[1] = _cnts[1], ort_cnt[2] = _cnts[2], ort_cnt[3] = _cnts[3];
					/*hh = _cnts[0]*2/(lenz + lenznext);*/
					/*ht = _cnts[1]*2/(lenz + lenznext);*/
					/*th = _cnts[2]*2/(lenz + lenznext);*/
					/*tt = _cnts[3]*2/(lenz + lenznext);*/
					/*hh = _cnts[0]*4/lenz / lenznext;*/
					/*ht = _cnts[1]*4/lenz / lenznext;*/
					/*th = _cnts[2]*4/lenz / lenznext;*/
					/*tt = _cnts[3]*4/lenz / lenznext;*/
				} else {
					float *_cnts = cds[znext].cnts[aux].cnt;		
					/*fprintf(stderr, "%s\t%s\t%f\t%f\t%f\t%f\n", ctgs->seq[z].name, ctgs->seq[znext].name, _cnts[0], _cnts[2], _cnts[1], _cnts[3]);*/
					ort_cnt[0] = _cnts[0], ort_cnt[1] = _cnts[2], ort_cnt[2] = _cnts[1], ort_cnt[3] = _cnts[3];
					/*hh = _cnts[0]*2/(lenz + lenznext);*/
					/*ht = _cnts[2]*2/(lenz + lenznext);*/
					/*th = _cnts[1]*2/(lenz + lenznext);*/
					/*tt = _cnts[3]*2/(lenz + lenznext);*/
					/*hh = _cnts[0]*2/(lenz * lenznext);*/
					/*ht = _cnts[2]*2/(lenz * lenznext);*/
					/*th = _cnts[1]*2/(lenz * lenznext);*/
					/*tt = _cnts[3]*2/(lenz * lenznext);*/
				}	
				
				if (vts[z<<1].sc + ort_cnt[2] > vts[z << 1 | 1].sc + ort_cnt[0]) {
					// -> ->  th                     <- -> hh
					vts[znext << 1].sc = vts[i<<1].sc + ort_cnt[2];
					vts[znext << 1].bt = ( z << 1 ) + 1;//plus 1 here to use 0 indicate an end
				} else {
					vts[znext << 1].sc = vts[i<<1|1].sc + ort_cnt[0];
					vts[znext << 1].bt = (z << 1 | 1 ) + 1;
				}

				if (vts[z<<1].sc + ort_cnt[3] > vts[z << 1 | 1].sc + ort_cnt[1]) {
					// -> <-  tt                <- <- ht
					vts[znext << 1 | 1].sc = vts[z<<1].sc + ort_cnt[3];
					vts[znext << 1 | 1].bt = (z << 1 ) + 1;
				} else {
					vts[znext << 1 | 1].sc = vts[z<<1|1].sc + ort_cnt[1];
					vts[znext << 1 | 1].bt = (z << 1 | 1) + 1;
				}
				zprev = z, z = znext;	
				idx[znext] |= 1; //set visited
				if (ctgs->seq[znext].snp_n == 1) {
				
					/*fprintf(stderr, "D: %s zprev:%d z:%d znext: %d\n", ctgs->seq[z].name, zprev, z, znext);*/
					if (vts[znext << 1 | 1].sc > vts[znext << 1].sc) 
						kv_push(uint32_t, endv, (znext<<1 | 1)+1);
					else
						kv_push(uint32_t, endv, (znext<<1) + 1);
					break;
				}	
			}
		} 
	}
	// traverse up the path 
	for (i = 0; i < endv.n; ++i) {
		uint32_t z;
		j = 0;
		for (z = endv.a[i]; z; z = vts[z - 1].bt) {
			/*fprintf(stderr, "z: %d\n", z >> 1);*/
			idx[j++] = z-1;
			/*fprintf(stderr, "%s%c\t", ctgs->seq[(z-1)>>1].name, (z-1) & 1 ? '+':'-');	*/
		}
		fprintf(stderr, "P\t");
		for (z = j - 1; z > 0; --z) {
		   add_udedge(g, ctgs->seq[idx[z]>>1].name, !(idx[z] & 1), ctgs->seq[idx[z-1]>>1].name, (idx[z-1] & 1), (vts[idx[z-1]].sc - vts[idx[z]].sc)*1000000);
			fprintf(stderr, "%s%c,",ctgs->seq[idx[z]>>1].name, idx[z] & 1 ? '+':'-');
		}
			fprintf(stderr, "%s%c\n",ctgs->seq[idx[z]>>1].name, idx[z] & 1 ? '+':'-');
			/*fprintf(stderr, "\n");	*/
			/*is_l = idx >> 1, is_l2 = idx & 1, add_dedge(g, name1, is_l, name2, is_l2, use_nw?c->cnts[j].ncnt: c->cnts[j].cnt[idx] * 2/(len1 + len2));	 //kinda residule cause index of name1 is the same as its index in ctgs but user doesn't know how the node is organized so better keep this.*/
	}	
	kv_destroy(endv);	
	free(vts);
	free(idx);
	kv_destroy(mets);
	return g;
}

graph_t *nns_straight(cdict2_t *cds, sdict_t *ctgs, int min_wt, float min_mdw, int use_nw, int amode)
{
	graph_t *g = graph_init();
	
	uint32_t i, j;
	//create nodes
	for ( i = 0; i < ctgs->n_seq; ++i) {
		/*fprintf(stderr, "add node %s %d\n", ctgs->seq[i].name, ctgs->seq[i].is_circ);*/
		add_node(g, ctgs->seq[i].name, 0, ctgs->seq[i].len, ctgs->seq[i].is_circ);
	}
	//create edges
		
	for ( i = 0; i < ctgs->n_seq; ++i) {
		char *name1 = ctgs->seq[i].name;
		uint32_t len1 = ctgs->seq[i].len;
		cdict2_t *c = cds + i;	
		int isf = 1;
		for (j = 0; j < c->lim; ++j) {
			char *name2 = c->cnts[j].name;
			if (strcmp(name1, name2) == 0) continue;
			//hsortand shaking
			/*fprintf(stderr, "try hand shaking\n");*/
			uint32_t ind = sd_get(ctgs, name2);
			uint32_t k;
			uint8_t hand_shaking = 0;
			uint32_t is_l, is_l2;
			uint32_t len2 = ctgs->seq[ind].len;
			for ( k = 0; k < cds[ind].lim; ++k) {
					if (strcmp(name1, cds[ind].cnts[k].name) == 0) {
						hand_shaking = 1;
						break;
					}
			}	
			/*if (hand_shaking) fprintf(stderr, "I hand shaking\n");*/
			//determine joint direction, here we use a very easy model
			if (!hand_shaking) continue;
			fprintf(stderr, "L\t%s\t%s\n", name1, name2);
			int idx;
			if (~(idx = det_ori(c->cnts[j].cnt, isf, min_wt, min_mdw, amode))) 
			// if build in accurate mode 
			/*if (~(idx = det_ori(c->cnts[j].cnt)))*/
			/*uint32_t hh = c->cnts[j].cnt[0];*/
			/*uint32_t ht = c->cnts[j].cnt[1];*/
			/*uint32_t th = c->cnts[j].cnt[2];*/
			/*uint32_t tt = c->cnts[j].cnt[3];*/
			/*uint32_t tl = c->cnts[j].cnt[2] + c->cnts[j].cnt[3];*/
			/*uint32_t hd = c->cnts[j].cnt[0] + c->cnts[j].cnt[1];	*/
			/*uint32_t hd2 = c->cnts[j].cnt[0] + c->cnts[j].cnt[2];*/
			/*uint32_t tl2 = c->cnts[j].cnt[1] + c->cnts[j].cnt[3];*/
			/*if (tl != hd) */
					/*is_l = tl > hd ? 0 : 1; */
			/*else */
					/*continue;*/
			/*if (hd2 != tl2) */
					/*is_l2 = tl2 > hd ? 0 : 1; 	*/
			/*else */
					/*continue;*/
			//when min_wt == -1; don't normalize weight 
			/*is_l = idx >> 1, is_l2 = idx & 1, add_dedge(g, name1, is_l, name2, is_l2, ~min_wt?c->cnts[j].cnt[idx] / (len1 / 2 + len2 / 2) : c->cnts[j].cnt[idx]);	 //kinda residule cause index of name1 is the same as its index in ctgs but user doesn't know how the node is organized so better keep this.*/
			is_l = idx >> 1, is_l2 = idx & 1, add_dedge(g, name1, is_l, name2, is_l2, c->cnts[j].cnt[idx]);	 //kinda residule cause index of name1 is the same as its index in ctgs but user doesn't know how the node is organized so better keep this.
			/*isf = 0;*/
		}		
	}	
	return g;
}

graph_t *nns_straight2(cdict_t *cds, sdict_t *ctgs, int min_wt, float min_mdw, int use_nw, int amode)
{
	graph_t *g = graph_init();
	
	uint32_t n = ctgs->n_seq << 1;
	uint32_t i, j;
	//create nodes
	for ( i = 0; i < ctgs->n_seq; ++i) 
		add_node(g, ctgs->seq[i].name, 0, ctgs->seq[i].len, ctgs->seq[i].is_circ);
	//create edges	
	for ( i = 0; i < n; ++i) {
		char *name1 = ctgs->seq[i>>1].name;
		uint8_t is_l = i & 1;	
		cdict_t *c = cds + i;	
		if (!c->n_cnt) continue;	
		for (j = 0; j < c->lim; ++j) {
			char *name2 = c->cnts[j].name;
			if (strcmp(name1, name2) == 0) continue;
			//hsortand shaking
			/*fprintf(stderr, "try hand shaking\n");*/
			uint32_t ind = sd_get(ctgs, name2) << 1 | c->cnts[j].is_l;
			uint32_t k;
			uint8_t hand_shaking = 0;
			float ocnt = 0;
			for ( k = 0; k < cds[ind].lim; ++k) {
					if (strcmp(name1, cds[ind].cnts[k].name) == 0 && cds[ind].cnts[k].is_l == is_l) {
						hand_shaking = 1;
						ocnt = cds[ind].cnts[k].cnt;
						break;
					}
			}	
			/*if (hand_shaking) fprintf(stderr, "I hand shaking\n");*/
			ocnt = 1;
			if (hand_shaking) add_dedge(g, name1, is_l, name2, c->cnts[j].is_l, ocnt * c->cnts[j].cnt);	 //kinda residule cause index of name1 is the same as its index in ctgs but user doesn't know how the node is organized so better keep this.
		}		
	}	
	return g;
}

graph_t *nns_mst2(cdict_t *cds, sdict_t *ctgs)
{
	graph_t *g = graph_init();
	
	uint32_t i, j;
	//create nodes
	for ( i = 0; i < ctgs->n_seq; ++i) {
		add_node(g, ctgs->seq[i].name, 0, ctgs->seq[i].len, ctgs->seq[i].is_circ);
	}
	//create edges
	kvec_t(mst_edge_t) mets;
	kv_init(mets);
	mst_edge_t tmp;
	float max_edge = 0;
	for ( i = 0; i < ctgs->n_seq << 1; ++i) {
		char *name1 = ctgs->seq[i>>1].name;
		uint8_t is_l = i & 1;	
		cdict_t *c = cds + i;	
		if (!c->n_cnt) continue;	
		for (j = 0; j < c->lim; ++j) {
			char *name2 = c->cnts[j].name;
			if (strcmp(name1, name2) == 0) continue;
			//hsortand shaking
			/*fprintf(stderr, "try hand shaking\n");*/
			uint32_t ind = sd_get(ctgs, name2) << 1 | c->cnts[j].is_l;
			uint32_t k;
			uint8_t hand_shaking = 0;
			float ocnt = 0;
			for ( k = 0; k < cds[ind].lim; ++k) {
					if (strcmp(name1, cds[ind].cnts[k].name) == 0 && cds[ind].cnts[k].is_l == is_l) {
						hand_shaking = 1;
						ocnt = cds[ind].cnts[k].cnt;
						break;
					}
			}	
			if (!hand_shaking) continue;
			if (i < ind) {
				if (c->cnts[j].cnt > max_edge) max_edge = c->cnts[j].cnt + 1;
				tmp = (mst_edge_t){i, ind, j, 0, c->cnts[j].cnt};
				kv_push(mst_edge_t, mets, tmp);
			}
		}		
	}	
		
	//add edges
	for (i = 0; i < ctgs->n_seq; ++i) {
		tmp = (mst_edge_t) {i << 1, (i << 1) | 1, 0, 0, max_edge};
		kv_push(mst_edge_t, mets, tmp);
	}
	mst_edge_t *in_mst = mets.a;
	_kmst(mets.a, mets.n, ctgs->n_seq << 1);	
	for (i = 0; i < ctgs->n_seq; ++i) ctgs->seq[i].l_snp_n = 0, ctgs->seq[i].r_snp_n = 0;
	for (i = 0; i < mets.n; ++i) {
		if (in_mst[i].in_mst) {
			uint32_t v = in_mst[i].s;
			uint32_t w = in_mst[i].d;
			if (v&1) ++ctgs->seq[v>>1].r_snp_n;
			else ++ctgs->seq[v>>1].l_snp_n; 
			if (w&1) ++ctgs->seq[w>>1].r_snp_n;
			else ++ctgs->seq[w>>1].l_snp_n;
		}
	}
	// filter out the vertices with more than 2 degrees, then break			0:left 1:right 
	uint64_t *idx = (uint64_t *)calloc(ctgs->n_seq << 1, sizeof(uint64_t));
	for (i = 0, j = 0; i < mets.n; ++i) {
		if (in_mst[i].in_mst) {
			uint32_t v = in_mst[i].s;
			uint32_t w = in_mst[i].d;
			uint32_t dv = v & 1 ? ctgs->seq[v>>1].r_snp_n: ctgs->seq[v>>1].l_snp_n;	
			uint32_t dw = w & 1 ? ctgs->seq[w>>1].r_snp_n: ctgs->seq[w>>1].l_snp_n;	
			
			if (dv > 2 || dw > 2) 
				in_mst[i].in_mst = 0, ++idx[v], ++idx[w];
			else
				in_mst[j++] = in_mst[i]; 
		} 	
	}
	for (i = 0; i < ctgs->n_seq << 1; i += 2) ctgs->seq[i>>1].l_snp_n -= idx[i];
	for (i = 1; i < ctgs->n_seq << 1; i += 2) ctgs->seq[i>>1].r_snp_n -= idx[i];

	mets.n = j;
	/*for (i = 0; i < mets.n; ++i)  fprintf(stderr, "E\t%s\t%s\t%u\t%u\t%f\n", ctgs->seq[in_mst[i].s>>1].name, ctgs->seq[in_mst[i].d>>1].name, ctgs->seq[in_mst[i].s>>1].len, ctgs->seq[in_mst[i].d>>1].len, in_mst[i].wt*10000);*/
	
	for (i = 0; i < mets.n; ++i) if ((in_mst[i].s | 1) != in_mst[i].d)add_udedge(g, ctgs->seq[in_mst[i].s>>1].name, in_mst[i].s & 1, ctgs->seq[in_mst[i].d>>1].name, in_mst[i].d & 1, in_mst[i].wt * 1000000); 
	// double the edges and index
	/*for (i = 0; i < j; ++i) {*/
		/*uint32_t v = mets.a[i].s;*/
		/*uint32_t w = mets.a[i].d;*/
		/*uint32_t aux = mets.a[i].aux;*/
		/*tmp = (mst_edge_t){w, v, aux, 0, 0};*/
		/*kv_push(mst_edge_t, mets, tmp);*/
	/*}	*/
	/*qsort(mets.a, mets.n, sizeof(mst_edge_t), cmp_v);	*/
	//index 
	/*memset(idx, 0, ctgs->n_seq * sizeof(uint64_t));*/
	/*in_mst = mets.a; // mets.a may have been changed*/
	/*for (i = 1, j = 0; i <= mets.n; ++i) {*/
		/*if (i==mets.n || in_mst[i].s != in_mst[j].s) {*/
			/*fprintf(stderr, "s: %u d: %u\n", in_mst[j].s, in_mst[j].d);*/
			/*idx[in_mst[j].s] = ((uint64_t)j << 32 | (i -j)) << 1 | 0;*/
			/*j = i;	*/
		/*}*/
	/*}*/
		
	free(idx);
	kv_destroy(mets);
	return g;
}
//single hic link
int buildg_hic2(char *fn, char *edge_fn, int min_wt, int use_sat, int norm, float min_mdw, int mlc, char *out_fn, int use_nw, int amode, int igm, int usep, int use_mst)
{
	graph_t *og; 
	sdict_t *ctgs = 0;
	if (use_sat) {
#ifdef VERBOSE
		fprintf(stderr, "[M::%s] collecting contigs from sat file\n", __func__);
#endif
		og = load_sat(fn);
		ctgs = col_ctgs_from_graph(og);
	} else {
		if (!fn) {
			fprintf(stderr, "[E::%s] please set reference index file with -c\n", __func__);
			return 1;
		}
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting contigs from faidx file\n", __func__);
#endif
		ctgs = col_ctgs(fn);	
		og = graph_init();	
	}
	if (!ctgs) return 1;
	uint32_t n_ctg = ctgs->n_seq;	
	/*fprintf(stderr, "%u\n", ctgs->n_seq);*/
	cdict_t* cds = (cdict_t*) calloc(n_ctg<<1, sizeof(cdict_t));
	uint32_t n_cds = n_ctg<<1;
	uint32_t i;
	for ( i = 0; i < n_cds; ++i) cd_init(cds+i);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting links\n", __func__);
#endif
	get_links(edge_fn, cds, ctgs);
	/*anothernorm(cds, ctgs);*/
	/*return 0;*/
	norm_links2(cds, ctgs, norm, igm, usep, min_wt);
	for ( i = 0; i < n_cds; ++i) cd_sort(cds+i);
	cd_set_lim(cds, n_cds, min_wt, min_mdw, mlc, norm);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] building graph\n", __func__);
#endif
	graph_t *g = 0;
	if (use_mst) 
		g = nns_mst2(cds, ctgs);
	else
		g = nns_straight2(cds, ctgs, min_wt, min_mdw, use_nw, amode);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing graph\n", __func__);
#endif
	process_graph(g, 0);
			
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] merging graph\n", __func__);
#endif
	merge_graph(og, g, 1);

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] output graph\n", __func__);
#endif

	dump_sat(og, out_fn);
		/*fprintf(stderr, "%s\n", ctgs->seq[i>>1].name);*/
	for (i = 0; i < n_cds; ++i) {
		/*fprintf(stderr, "%s\n", ctgs->seq[i>>1].name);*/
		cd_destroy(cds +i);	
	} 
	if (cds) free(cds);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] releasing memory\n", __func__);
#endif
	graph_destroy(og);
	return 0;
}

int buildg_hic(char *fn, char *edge_fn, int min_wt, int use_sat, int norm, float min_mdw, int mlc, char *out_fn, int use_nw, int amode, int igm, int usep, int use_mst)
{
	/*fprintf(stderr, "%s %s\n", fn, edge_fn);*/
	/*fprintf(stderr, "%d %d\n", min_wt, use_sat);*/
	/*fprintf(stderr, "%d %f\n", norm, min_mdw);*/
	/*fprintf(stderr, "%d\n", mlc);*/
	/*fprintf(stderr, "%p %s\n", out_fn, out_fn);*/
	graph_t *og; 
	sdict_t *ctgs = 0;
	if (use_sat) {
#ifdef VERBOSE
		fprintf(stderr, "[M::%s] collecting contigs from sat file\n", __func__);
#endif
		og = load_sat(fn);
		ctgs = col_ctgs_from_graph(og);
	} else {
		if (!fn) {
			fprintf(stderr, "[E::%s] please set reference index file with -c\n", __func__);
			return 1;
		}
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting contigs from faidx file\n", __func__);
#endif
		ctgs = col_ctgs(fn);	
		og = graph_init();	
	}
	if (!ctgs) return 1;
	uint32_t n_ctg = ctgs->n_seq;	
	/*fprintf(stderr, "%u\n", ctgs->n_seq);*/
	cdict2_t* cds = calloc(n_ctg, sizeof(cdict2_t)); 
	uint32_t i;
	for ( i = 0; i < n_ctg; ++i) cd2_init(cds+i); 
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting links\n", __func__);
#endif
	get_links_hic(edge_fn, cds, ctgs);
	/*anothernorm(cds, ctgs);*/
	/*return 0;*/
	/*if (norm) for (i = 0; i < n_cds; ++i) cd_norm(cds + i);*/
	/*print_cdict2(cds, ctgs);	*/
	norm_links(cds, ctgs, norm, igm, usep, min_wt);
	for ( i = 0; i < n_ctg; ++i) cd2_sort(cds+i); 
	/*print_cdict2(cds, ctgs);	*/
	cd2_set_lim(cds, n_ctg, mlc); 
	/*if (norm) */
	/*if (norm) cd_filt(cds, n_cds, min_rat); */
	/*if (norm) for (i = 0; i < n_cds; ++i) cd_norm(cds + i);*/
	/*if (norm) for (i = 0; i < n_cds; ++i) cd_norm(cds + i);*/
	/*for (i = 0; i < n_cds; ++i)	cd_norm(cds + i);*/
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] building graph\n", __func__);
#endif
	graph_t *g = 0;
	if (use_mst) 
		g = nns_mst(cds, ctgs);
	else
		g = nns_straight(cds, ctgs, min_wt, min_mdw, use_nw, amode);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing graph\n", __func__);
#endif
	process_graph(g, 0);
			
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] merging graph\n", __func__);
#endif
	merge_graph(og, g, 1);

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] output graph\n", __func__);
#endif

	dump_sat(og, out_fn);
	for (i = 0; i < n_ctg; ++i) {
		/*fprintf(stderr, "%s\n", ctgs->seq[i>>1].name);*/
		cd2_destroy(cds +i);	

	} 
	/*fprintf(stderr, "CheckPoint2\n");*/
	/*fprintf(stderr, "leave\n");*/
	if (cds) free(cds);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] releasing memory\n", __func__);
#endif
	graph_destroy(og);
	return 0;

}

int buildg(char *fn, char *edge_fn, int min_wt, int use_sat, int norm, float min_mdw, int mlc, char *out_fn)
{
	/*fprintf(stderr, "%s %s\n", fn, edge_fn);*/
	/*fprintf(stderr, "%d %d\n", min_wt, use_sat);*/
	/*fprintf(stderr, "%d %f\n", norm, min_mdw);*/
	/*fprintf(stderr, "%d\n", mlc);*/
	/*fprintf(stderr, "%p\n", out_fn);*/
	graph_t *og; 
	sdict_t *ctgs = 0;
	if (use_sat) {
#ifdef VERBOSE
		fprintf(stderr, "[M::%s] collecting contigs from sat file\n", __func__);
#endif
		og = load_sat(fn);
		ctgs = col_ctgs_from_graph(og);
	} else {
		if (!fn) {
			fprintf(stderr, "[E::%s] please set reference index file with -c\n", __func__);
			return 1;
		}
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting contigs from faidx file\n", __func__);
#endif
		ctgs = col_ctgs(fn);	
		og = graph_init();	
	}
	if (!ctgs) return 1;
	uint32_t n_ctg = ctgs->n_seq;	
	/*fprintf(stderr, "%u\n", ctgs->n_seq);*/
	cdict_t* cds = calloc(n_ctg<<1, sizeof(cdict_t)); 
	uint32_t n_cds = n_ctg<<1;
	uint32_t i;
	for ( i = 0; i < n_cds; ++i) cd_init(cds+i); 
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] collecting links\n", __func__);
#endif
	get_links(edge_fn, cds, ctgs);
	/*anothernorm(cds, ctgs);*/
	/*return 0;*/
	if (norm) for (i = 0; i < n_cds; ++i) cd_norm(cds + i);
	for ( i = 0; i < n_cds; ++i) cd_sort(cds+i); 
	print_cdict(cds, ctgs);	
	cd_set_lim(cds, n_cds, min_wt, min_mdw, mlc, norm); 
	/*if (norm) */
	/*if (norm) cd_filt(cds, n_cds, min_rat); */
	/*if (norm) for (i = 0; i < n_cds; ++i) cd_norm(cds + i);*/
	/*if (norm) for (i = 0; i < n_cds; ++i) cd_norm(cds + i);*/
	/*for (i = 0; i < n_cds; ++i)	cd_norm(cds + i);*/
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] building graph\n", __func__);
#endif
	graph_t *g = build_graph(cds, ctgs);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing graph\n", __func__);
#endif
	process_graph(g, 1);
			
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] merging graph\n", __func__);
#endif
	merge_graph(og, g, 1);

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] output graph\n", __func__);
#endif

	dump_sat(og, out_fn);

	for (i = 0; i < n_cds; ++i) {
		/*fprintf(stderr, "%s\n", ctgs->seq[i>>1].name);*/
		cd_destroy(cds +i);	

	} 
	/*fprintf(stderr, "leave\n");*/
	if (cds) free(cds);
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] releasing memory\n", __func__);
#endif
	graph_destroy(og);
	return 0;

}

int buildg_hic_cont(int which, char *fn, char *edge_fn, int min_wt, int use_sat, int norm, float min_mdw, int mlc, char *out_fn, int use_nw, int amode, int igm, int usep, int use_mst)
{
	int ret;
	if (which == 1) 
		ret = buildg_hic(fn, edge_fn, min_wt, use_sat, norm, min_mdw, mlc, out_fn, use_nw, amode, igm, usep, use_mst);
	else
		ret = buildg_hic2(fn, edge_fn, min_wt, use_sat, norm, min_mdw, mlc, out_fn, use_nw, amode, igm, usep, use_mst);
		
	return ret;
}

int main_bldg(int argc, char *argv[], int ishic)
{
	int c;
	uint32_t min_wt = 10; char *program = argv[0];
	char *sat_fn = 0, *ctg_fn = 0, *out_fn = 0;
	int use_sat = 0, mlc = 1;
	int norm = 0, amode = 0, use_nw = 1, use_mst = 0;
	float msn = .7, mdw = 0.95;
	int igm = 0, usep = 0;
	int md = 2; //new method
	--argc, ++argv;
	while (~(c=getopt(argc, argv, "w:ao:s:c:nem:f:k:hpg1B:"))) {
		switch (c) {
			case '1': 
				use_mst = 1;
				break;
			case 'w': 
				min_wt = atoi(optarg);
				break;
			case 'a': 
				amode = 1;
				break;
			case 'B': 
				md = atoi(optarg);
				break;
			case 'p': 
				usep = 1;
				break;
			case 'g': 
				igm = 1;
				break;
			case 's': 
				sat_fn = optarg;
				use_sat =  1;
				break;
			case 'c': 
				ctg_fn = optarg;
				break;
			case 'f': 
				mdw = atof(optarg);
				break;
			case 'k': 
				mlc = atoi(optarg);
				break;
			case 'n': 
				norm = 1;
				break;
			case 'e': 
				use_nw = 0;
				break;
			case 'o': 
				out_fn = optarg;
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s %s [<options>] <LINKS_MATRIX> \n", program, argv[0]);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -1    BOOL     use MST to construct scaffolding graph [FALSE, only for pin_hic]\n");
				fprintf(stderr, "         -a    BOOL     Hi-C scaffolding in accurate mode [FALSE]\n");
				fprintf(stderr, "         -w    INT      minimum weight for links [10]\n");
				fprintf(stderr, "         -k    INT      maximum linking candiates [1]\n");
				fprintf(stderr, "         -c    FILE     reference index file [nul]\n");
				fprintf(stderr, "         -n    BOOL     normalize weight [false]\n");
				fprintf(stderr, "         -p    BOOL     use product of length [False]\n");
				fprintf(stderr, "         -g    BOOL     ignore middle part of contigs [false]\n");
				fprintf(stderr, "         -e    BOOL     use normalized weight as edge weight [TRUE]\n");
				fprintf(stderr, "         -f    FLOAT    minimum weight difference [0.95]\n");
				fprintf(stderr, "         -s    FILE     sat file [nul]\n");
				fprintf(stderr, "         -o    FILE     output file [stdout]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 1 > argc) {
		fprintf(stderr,"[E::%s] require number of contig and links file!\n", __func__); goto help;
	}
	char *lnk_fn = argv[optind];
	fprintf(stderr, "Program starts\n");	
	if (!sat_fn) sat_fn = ctg_fn;
	int ret;
	if (!ishic)
		ret = buildg(sat_fn, lnk_fn, min_wt, use_sat, norm, mdw, mlc, out_fn);
	else 
		ret = buildg_hic_cont(md, sat_fn, lnk_fn, min_wt, use_sat, norm, mdw, mlc, out_fn, use_nw, amode, igm, usep, use_mst);
	
	fprintf(stderr, "Program ends\n");	
	return ret;	

}


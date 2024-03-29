/*
 * =====================================================================================
 *
 *       Filename:  make_brk.c
 *
 *    Description:	generate breaks for a scaffold
 *
 *        Version:  1.0
 *        Created:  14/11/2019 09:39:40
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D.G.), dfguan9@gmail.com
 *   Organization:  Harbin Institute of Technology
 *
 * =====================================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>	

#include "graph.h" 
#include "sdict.h"
#include "cdict.h"
#include "bed.h"
#include "kvec.h"
#include "utls.h"
#include "asset.h"
#include "bamlite.h"

typedef struct {
	/*int mq:15, rev:1, as:16;*/
	uint32_t s, ns;
   	uint32_t tid:31, rev:1;
	int qual;
}aln_inf_t;

typedef struct {
	float maxim, avg;
} cov_stat;

typedef struct {
	uint32_t cid;
	uint32_t c1s:31, c1rev:1;
	uint32_t c2s:31, c2rev:1;
}hit2_t;

typedef struct {
	uint32_t n, m;
	hit2_t *ary;
}hit2_ary_t;

void hit2_ary_push(hit2_ary_t *l, hit2_t *z)
{
	uint32_t max = -1;
	if (l->n >= l->m) {
		if (l->m > (max >> 1)) {
			fprintf(stderr, "Too many values here\n");
			exit(1);
		} else 
			l->m = l->m ? l->m << 1 : 16;
		l->ary = realloc(l->ary, sizeof(hit2_t) * l->m);//beware of overflow
	}
	l->ary[l->n++] = *z;
}

int cmp_hits2(const void *a, const void *b)
{
	hit2_t *m = (hit2_t *)a;
	hit2_t *n = (hit2_t *)b; //too many branches
	if (m->c1s > n->c1s) return 1;	
	else if (m->c1s < n->c1s) return -1;	
	else if (m->c1rev > n->c1rev) return 1;
	else if (m->c1rev < n->c1rev) return -1;
	else if (m->c2s > n->c2s) return 1;	
	else if (m->c2s < n->c2s) return -1;	
	else if (m->c2rev > n->c2rev) return 1;	
	else if (m->c2rev < n->c2rev) return -1;	
	else return 0;
}

int mb_col_contacts(hit2_ary_t *hit_ary, sdict_t *sd, ctg_pos_t *d)
{
	size_t i, j;
	/*sdict_t *use_sd = sd;*/
	hit2_t *hs = hit_ary->ary;
	size_t n = hit_ary->n;
	for (i = 0, j = 1; j <= n; ++j) {
		if (j == n || hs[i].cid != hs[j].cid || hs[i].c1s != hs[j].c1s || hs[i].c1rev != hs[j].c1rev || hs[i].c2s != hs[j].c2s || hs[i].c2rev != hs[j].c2rev) {
			/*if (hs[i].c2s > sd->seq[hs[i].cid].len) fprintf(stderr, "much larger\n");*/
			/*fprintf(stderr, "%s\t%u\t%u\n", sd->seq[hs[i].cid].name, hs[i].c1s, hs[i].c2s);*/
			pos_push(&d->ctg_pos[hs[i].cid], hs[i].c1s << 1);	
			pos_push(&d->ctg_pos[hs[i].cid], hs[i].c2s << 1 | 1);	
			i = j;	
		}
	}
	return 0;
}


int mb_col_hits2(aln_inf_t *f, int f_cnt, sdict_t *ctgs, sdict_t *scfs, hit2_ary_t *hit2_ary, int min_mq, char *qn)
{
	if (f[0].qual < min_mq || f[1].qual < min_mq) return 2;
	if (scfs->n_seq) {
			sd_seq_t *sq1 = &ctgs->seq[f[0].tid];
			sd_seq_t *sq2 = &ctgs->seq[f[1].tid];
			uint32_t ind1 = sq1->le; //maybe not well paired up
			uint32_t ind2 = sq2->le;
			uint32_t f0s = sq1->r_snp_n + (sq1->l_snp_n & 0x1 ? f[0].s : sq1->len - f[0].s);
			uint32_t f1s = sq2->r_snp_n + (sq2->l_snp_n & 0x1? f[1].s : sq2->len - f[1].s);
			fprintf(stdout, "%s\t%u\t%u\t%s/1\t%d\t%c\n", scfs->seq[ind1].name, f0s, f0s + 100, qn, f[0].qual, (f[0].rev ^ (sq1->l_snp_n & 0x1)) ? '-':'+');
			fprintf(stdout, "%s\t%u\t%u\t%s/2\t%d\t%c\n", scfs->seq[ind2].name, f1s, f1s + 100, qn, f[1].qual, (f[1].rev ^ (sq2->l_snp_n & 0x1)) ? '-':'+');
			return 0;	
	} else {
			uint32_t ind1 = f[0].tid;
			uint32_t ind2 = f[1].tid;
			uint32_t f0s = f[0].s;
			uint32_t f1s = f[1].s;
			fprintf(stdout, "%s\t%u\t%u\t%s/1\t%d\t%c\n", ctgs->seq[ind1].name, f0s, f0s + 100, qn, f[0].qual, (f[0].rev) ? '-':'+');
			fprintf(stdout, "%s\t%u\t%u\t%s/2\t%d\t%c\n", ctgs->seq[ind2].name, f1s, f1s + 100, qn, f[1].qual, (f[1].rev) ? '-':'+');
			return 0;	
	}
	return 1;
}

int mb_get_links_hic(char *links_fn, cdict2_t *cds, sdict_t *ctgs)
{	
	bed_file_t *bf = bed_open(links_fn);
	if (!bf) return 0;
	lnk_rec_t r;
	uint32_t line_n = 0;
	while (lnk_read(bf, &r) >= 0) {
		uint32_t ind1;
		ind1 = sd_get(ctgs, r.ctgn);		
		sd_put(ctgs, r.ctgn2);		
		/*uint32_t ind2 = sd_put2(ctgs, r.ctgn, 0, 0, 0, r.llen, r.rlen);		*/
		line_n += 1;
		cd2_add(&cds[ind1], r.is_l, r.ctgn2, r.is_l2, r.wt);	//this has been normalized	
	} 
	bed_close(bf);
	return 0;	
}

int mb_norm_links(cdict2_t *cds, sdict_t *ctgs)
{
	uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c ;
	for ( i = 0; i < n_cds; ++i) {
		/*char *name1 = ctgs->seq[i>>1].name;*/
		/*uint32_t snpn = i&1 ? ctgs->seq[i>>1].l_snp_n:ctgs->seq[i>>1].r_snp_n;*/
		uint32_t j;
		c = cds + i;
		float icnt;
        
		for (j = 0; j < c->n_cnt; ++j) {
            /*fprintf(stderr, "%s\n", c->cnts[j].name);*/
			char *name2 = c->cnts[j].name; 
            uint32_t ctg2_idx = sd_get(ctgs, name2);
			icnt = c->cnts[j].cnt[0] + c->cnts[j].cnt[1] + c->cnts[j].cnt[2] + c->cnts[j].cnt[3]; 
            c->cnts[j].ncnt = icnt / ctgs->seq[ctg2_idx].len;
            /*c->cnts[j].ncnt = (float) icnt / (ctgs->seq[ctg2_idx].l_snp_n + ctgs->seq[ctg2_idx].r_snp_n);*/
            /*uint32_t z;*/
            /*for ( z = 0; z < 4; ++z) c->cnts[j].fcnt[z] = (float) c->cnts[j].cnt[z]/(z >> 1 ? ctgs->seq[i].l_snp_n : ctgs->seq[i].r_snp_n) / ( z & 0x1 ? ctgs->seq[ctg2_idx].l_snp_n : ctgs->seq[ctg2_idx].r_snp_n);  */
			/*uint32_t snp2 = c->cnts[j].is_l ? ctgs->seq[sd_get(ctgs, name2)].l_snp_n:ctgs->seq[sd_get(ctgs,name2)].r_snp_n;*/
			/*fprintf(stderr, "%s\t%c\t%s\t%c\t%u\t%u\t%u\t%lf\n", name1, i&1?'+':'-', name2, c->cnts[j].is_l?'+':'-', icnt, snpn, snp2, 100000.0*(double)icnt/(snp2*snpn));*/
		}
	}
	return 0;
}


int deflimn(cdict2_t *cds, sdict_t *ctgs)
{
    uint32_t n_cds = ctgs->n_seq;
	uint32_t i;
	cdict2_t *c;
    kvec_t(uint32_t) v;
    kv_init(v);
	uint32_t dstr[8] = {0};
	int limn = 7;
	for ( i = 0; i < n_cds; ++i) {
        char *name = ctgs->seq[i].name;
        uint32_t scf_idx = ctgs->seq[i].le;
        uint32_t ctg_idx = i; 
		uint32_t j;
		c = cds + i;
        uint32_t fgt, flt, s_ok, p_ok;
        fgt = flt = s_ok = p_ok = 0;
        int ctgn = 0;
		for (j = 0; j < c->n_cnt; ++j) {
		   uint32_t ctg_idx2 =  sd_get(ctgs, c->cnts[j].name);	
            uint32_t scf_idx2 = ctgs->seq[ctg_idx2].le;
			/*if (ctgn <= 50) fprintf(stderr, "%u\t%s\t%u\t%s\t%f\t%f\t%f\t%f\n", ctg_idx, name, ctg_idx2, ctgs->seq[ctg_idx2].name, c->cnts[j].cnt[0], c->cnts[j].cnt[1], c->cnts[j].cnt[2], c->cnts[j].cnt[3]);*/
            if (scf_idx2 == scf_idx) {
				++ctgn;
				if (ctg_idx2 > ctg_idx) {
					++fgt;
					if (ctg_idx2 == ctg_idx + 1 && fgt < limn) s_ok = 1, ++dstr[fgt];
				} else {
					++flt;	
					if (ctg_idx2 + 1 == ctg_idx && flt < limn) p_ok = 1, ++dstr[flt];
				}
				/*if (ctg_idx2 > ctg_idx + 1) */
					/*fgt = 1;*/
				/*else if (ctg_idx2 + 1 == ctg_idx) */
					/*flt = 1;	*/
			}
			if ((p_ok && s_ok) || (fgt + 1 > limn  && flt + 1 > limn)) 
				break;
        }
		/*fprintf(stderr, "%s\t%d\n", name, ctgs->seq[ctg_idx].rs >> 1 & 1);*/
		/*if (!s_ok) fprintf(stderr, "miss successor add a break\n");	*/
		/*if (!p_ok) fprintf(stderr, "miss predecessor add a break\n");*/
		//1 tail 2 head 3 head + tail 0 middle
        if (!s_ok && (ctgs->seq[ctg_idx].rs & 1) == 0)  
			++dstr[limn];	
        if (!p_ok && (ctgs->seq[ctg_idx].rs >> 1 & 1) == 0)  
			++dstr[limn];	
        /*if (ctgs->seq[ctg_idx].rs & 1) */
            /*kv_push(uint32_t, v, ctg_idx << 1 | 1);  */
	}
	fprintf(stderr, "%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\n", dstr[1], dstr[2], dstr[2], dstr[3], dstr[4], dstr[5], dstr[6], dstr[7]);
}


int cal_cov_stat(cov_ary_t *ca, graph_t *g, sdict_t *ctgs, sdict_t *scfs)
{
	asm_t *as = &g->as.asms[g->as.casm];
	vertex_t *vt = g->vtx.vertices;
	path_t *pt = g->pt.paths;

	uint32_t n = as->n;
	uint32_t i;	
	uint32_t idx = 0;
	for ( i = 0; i < n; ++i) {
		uint32_t m = pt[as->pn[i] >> 1].n;
		cov_stat_t *ct = malloc(sizeof(cov_stat_t) * (2 * m - 1));
		uint32_t j;
		uint32_t scf_id = sd_get(scfs, pt[as->pn[i]>>1].name);
		for ( j = 0; j < 2 * m - 1; ++j ) {
			if ( j % 2) {
				uint32_t st = ctgs->seq[idx].r_snp_n + ctgs->seq[idx].len + 1; ;  
				uint32_t ed = st + GAP_SZ - 1;	
				cal_cov_stat4reg(&ca[scf_id], st, ed, &ct[j]);		
				fprintf(stderr, "%s\t%u\t%u\t%s\t%u\t%f\t%f\t%f\n", pt[as->pn[i]>>1].name, st, ed, "gap", GAP_SZ, ct[j].maxim, ct[j].avg, ct[j].rsd);
				++idx;
			} else {
				uint32_t st = ctgs->seq[idx].r_snp_n + 1;
				uint32_t ed = st + ctgs->seq[idx].len - 1;				
				cal_cov_stat4reg(&ca[scf_id], st, ed, &ct[j]);		
				fprintf(stderr, "%s\t%u\t%u\t%s\t%u\t%f\t%f\t%f\n", pt[as->pn[i]>>1].name, st, ed, ctgs->seq[idx].name, ctgs->seq[idx].len, ct[j].maxim, ct[j].avg, ct[j].rsd);
			}	
		}
		++idx;
		free(ct);
	}
	return 0;
}

int mb_init_scaffs(graph_t *g, sdict_t *ctgs, sdict_t *scfs)
{
	/*dump_sat(g);*/
	asm_t *as = &g->as.asms[g->as.casm];
	vertex_t *vt = g->vtx.vertices;
	path_t *pt = g->pt.paths;
	uint32_t n = as->n;
	uint32_t i;	
	/*for (i = 0; i < g->vtx.n; ++i) sd_put(ctgs, vt[i].name);	*/
	FILE *chrsf = fopen("chrom.sizes", "w");
		
	for ( i = 0; i < n; ++i) {
		uint32_t m = pt[as->pn[i] >> 1].n;
		/*fprintf(stderr, "%d\n", as->pn[i]); */
		uint32_t *p = pt[as->pn[i] >> 1].ns;
		uint32_t j, len = 0, len_ctg;
		//push scaffold name to scfs 
		/*int32_t scf_id = as->pn[i] >> 1;*/
		int32_t scf_id =sd_put2(scfs, pt[as->pn[i]>>1].name, 0, 0, 0, 0, 0);
		uint32_t sid;
		for ( j = 0; j < m; ++j ) { // pid, length,   
			len_ctg = vt[p[j]>>2].len;	
			sid = sd_put(ctgs, vt[p[j]>>2].name);
			//contig points to scaffold and set its start and direction
			/*fprintf(stderr, "sid: %u\t%s\t%s\n", sid, pt[as->pn[i]>>1].name, vt[p[j]>>2].name, p[j]);*/
			/*fprintf(stdout, "%s\t%s\t%u\n", pt[as->pn[i]>>1].name, vt[p[j]>>2].name, len);*/
			ctgs->seq[sid].len = len_ctg;
			/*ctgs->seq[sid].le = as->pn[i]>>1;*/
			ctgs->seq[sid].le = scf_id;
			ctgs->seq[sid].l_snp_n = (j << 1) | (p[j] & 0x1);
			ctgs->seq[sid].r_snp_n = len;
			ctgs->seq[sid].rs = 0;	
			if (!j) 
				ctgs->seq[sid].rs = 2;
			if (j == m - 1) 
				ctgs->seq[sid].rs |= 1; 
			len += len_ctg;
			/*else */
				/*len += len_ctg + 200;*/
		}
		//reset scaffold length le rs l_snp_n, r_snp_n
		fprintf(chrsf, "%s %u\n", pt[as->pn[i]>>1].name, len);
		sd_put4(scfs, pt[as->pn[i]>>1].name, len, len >> 1, (len >>1) + 1, len >> 1, len >> 1, pt[as->pn[i]>>1].is_circ);
		/*free(p);*/
	}
	fclose(chrsf);
	return 0;
}
int cmp_brks(const void *a, const void *b)
{
    uint32_t p = *(uint32_t *)a;
    uint32_t q = *(uint32_t *)b;
    if (p < q) return -1;
    else if (p == q) return 0;
    else return 1;
}

int cut_paths(graph_t *g, uint32_t *brks, uint32_t n_brks, sdict_t *ctgs, sdict_t *scfs) 
{
    qsort(brks, n_brks, sizeof(uint32_t), cmp_brks);  
    
    uint32_t i, j;
    /*for (i = 0; i < n_brks; ++i) {*/
        /*fprintf(stdout, "%s\t%d\t%d\n", ctgs->seq[brks[i]>>1].name, brks[i]>>1, brks[i]&1);*/
    /*}*/
    /*name_t nt, nt1;*/
	char *name, *name2;
    /*for ( i = 1, j = 0; i <= n_brks; ++i) {*/
        /*if (i == n_brks || brks[i] != brks[j]) {*/
            /*uint32_t ctg_idx = brks[j];*/
			/*fprintf(stdout, "%s\t%s\n", ctgs->seq[ctg_idx >> 1].name, ctgs->seq[(ctg_idx >> 1)+1].name);*/
			/*j = i;        */
        /*} */
    /*}*/
	uint32_t scf_id = -1;	
	kvec_t(uint32_t) pos;
	kv_init(pos);
	int brkn = 0;
	int adpn = 0;
	for ( i = 1, j = 0; i <= n_brks; ++i) {
		if (i == n_brks || brks[i] != brks[j]) {
			uint32_t ctg_idx = brks[j] >> 1;
			if (ctgs->seq[ctg_idx].le != scf_id) {
				if (~scf_id) brkn += pos.n, adpn += break_path(g, scf_id, pos.a, pos.n);
				kv_reset(pos);
				scf_id = ctgs->seq[ctg_idx].le;	
			}	
			kv_push(uint32_t, pos, ctgs->seq[ctg_idx].l_snp_n >> 1);
			j = i;	
		} 			
		/*uint32_t ctg_idx = brks[i];*/
			/*fprintf(stdout, "%u\t%s\t%s\n", ctgs->seq[ctg_idx>>1].le, ctgs->seq[ctg_idx >> 1].name, ctgs->seq[(ctg_idx >> 1)+1].name);*/
	}
	fprintf(stderr, "[M::%s] find %d breaks, add %d paths\n", __func__, brkn, adpn);
	kv_destroy(pos);
    return 0;
}

int mb_proc_bam(char *bam_fn, int min_mq, sdict_t *ctgs, sdict_t *scfs, hit2_ary_t *ha)
{
	
	bamFile fp;
	bam_header_t *h;
	bam1_t *b;
	fp = bam_open(bam_fn, "r"); //should check if bam is sorted
	if (fp == 0) {
		fprintf(stderr, "[E::%s] fail to open %s\n", __func__, bam_fn);
		return -1;
	}
	
	h = bam_header_read(fp);
	b = bam_init1();
	char **names = h->target_name;
	/*int i;*/
	/*for ( i = 0; i < h->n_targets; ++i) {*/
		/*char *name = h->target_name[i];*/
		/*uint32_t len = h->target_len[i];*/
		/*if (ws > len) ws = len;*/
		/*uint32_t le = (len - ws) >> 1;*/
		/*uint32_t rs = (len + ws) >> 1;*/
		/*uint32_t lenl, lenr;*/
		/*lenl = lenr = (len - ws) >> 1;*/
		/*sd_put2(ctgs, name, len, le, rs, lenl, lenr);*/
	/*}*/
	/*uint32_t cur_ws;*/
	/*for ( i = 0; i < h->n_targets; ++i) {*/
		/*char *name = h->target_name[i];*/
		/*uint32_t len = h->target_len[i];*/
		/*cur_ws = ws;*/
		/*if (len < (cur_ws << 1)) cur_ws = len >> 1;*/
		/*uint32_t le = cur_ws;*/
		/*uint32_t rs = len - cur_ws + 1;*/
		/*uint32_t lenl, lenr;*/
		/*lenl = lenr = cur_ws;*/
		/*sd_put2(ctgs, name, len, le, rs, lenl, lenr);*/
	/*}*/
	/*if (!ns->ct) { //not initiate yet*/
		/*init_gaps(gap_fn, ns, ctgs, max_ins_len);*/
	/*}*/

	char *cur_qn = 0;
	long bam_cnt = 0, interb_cnt = 0, intrab_cnt = 0;
	int is_set = 0;
	/*aln_inf_t aln[2];*/
	/*int aln_cnt;*/
	
	kvec_t(aln_inf_t) all;
	kv_init(all);
	kvec_t(aln_inf_t) five;
	kv_init(five);

	uint8_t rev;
	uint64_t rdp_counter  = 0;
	uint32_t rd1_cnt = 0, rd2_cnt = 0;
	uint32_t rd1_5cnt = 0, rd2_5cnt = 0;
	uint64_t used_rdp_counter = 0;
	/*fprintf(stderr, "Proc Bam %d\n", __LINE__);*/
	while (1) {
		//segment were mapped 
		if (bam_read1(fp, b) >= 0 ) {
			if (!cur_qn || strcmp(cur_qn, bam1_qname(b)) != 0) {
				if (rd1_cnt < 3 && rd2_cnt < 3 && rd1_5cnt == 1 && rd2_5cnt == 1 && !mb_col_hits2(five.a, five.n, ctgs, scfs, ha, min_mq, cur_qn)) ++used_rdp_counter;
				/*aln_cnt = 0;	*/
				/*rev = 0;*/
				/*is_set = 0;*/
				/*kv_reset(all);*/
				kv_reset(five);
				rd1_cnt = rd2_cnt = rd1_5cnt = rd2_5cnt = 0;
				if (cur_qn) ++rdp_counter, free(cur_qn); 
				cur_qn = strdup(bam1_qname(b));
			}
			b->core.flag & 0x40 ? ++rd1_cnt : ++rd2_cnt;	
			if (b->core.flag & 0x4) continue; //not aligned
			aln_inf_t tmp;
			rev = tmp.rev = !!(b->core.flag & 0x10);
			/*tmp.nrev = !!(b->core.flag & 0x20);*/
			//only collects five prime
			tmp.tid = sd_get(ctgs, names[b->core.tid]);
			/*tmp.ntid = b->core.mtid;*/
			tmp.s = b->core.pos + 1;
			tmp.qual = b->core.qual;
			/*tmp.ns = b->core.mpos + 1;*/
			/*kv_push(aln_inf_t, all, tmp);*/
				
			uint32_t *cigar = bam1_cigar(b);
			if ((!rev && bam_cigar_op(cigar[0]) == BAM_CMATCH) || (rev && bam_cigar_op(cigar[b->core.n_cigar-1]) == BAM_CMATCH)) {
				b->core.flag & 0x40 ? ++rd1_5cnt: ++rd2_5cnt; kv_push(aln_inf_t, five, tmp);
			}
			
			/*aln_cnt = (aln_cnt + 1 ) & 1;*/
			/*if ((++bam_cnt % 1000000) == 0) fprintf(stderr, "[M::%s] processing %ld bams\n", __func__, bam_cnt); */
		} else {
			if (rd1_cnt < 3 && rd2_cnt < 3 && rd1_5cnt == 1 && rd2_5cnt == 1 && !mb_col_hits2(five.a, five.n, ctgs, scfs, ha, min_mq, cur_qn)) ++used_rdp_counter;
			if (cur_qn) ++rdp_counter, free(cur_qn); 
			break;	
		}
	}
	fprintf(stderr, "[M::%s] finish processing %lld read pairs %lld (%.2f) passed\n", __func__, rdp_counter, used_rdp_counter, (double)used_rdp_counter/rdp_counter); 
	bam_destroy1(b);
	bam_header_destroy(h);
	bam_close(fp);
	kv_destroy(all);
	kv_destroy(five);
	return 0;
}
int col_ctgs(char *bam_fn, sdict_t *ctgs)
{
	bamFile fp;
	bam_header_t *h;
	bam1_t *b;
	fp = bam_open(bam_fn, "r"); //should check if bam is sorted
	if (fp == 0) {
		fprintf(stderr, "[E::%s] fail to open %s\n", __func__, bam_fn);
		return -1;
	}
	
	h = bam_header_read(fp);
	b = bam_init1();
	int i;
	for ( i = 0; i < h->n_targets; ++i) {
		char *name = h->target_name[i];
		uint32_t len = h->target_len[i];
		uint32_t le = len >> 1;
		uint32_t rs = (len >> 1) + 1;
		uint32_t lenl, lenr;
		lenl = lenr = len >> 1;
		sd_put4(ctgs, name, len, le, rs, lenl, lenr, 0);
	}
	bam_destroy1(b);
	bam_header_destroy(h);
	bam_close(fp);
	return 0;
}
int bam2hig_bed(char *sat_fn, char *bam_fn[], int n_bams, int min_mq, char *out_fn)
{
	
	sdict_t *ctgs = sd_init();	
	sdict_t *scfs = sd_init();
	int win_s = 0;
#ifdef VERBOSE
	fprintf(stderr, "[M::%s] initiate contigs\n", __func__);
#endif
	if (sat_fn && *sat_fn) {
		graph_t *og = load_sat(sat_fn);
		simp_graph(og);
		mb_init_scaffs(og, ctgs, scfs);	
		graph_destroy(og);
		if (!scfs->n_seq) {
			fprintf(stderr, "[E::%s] fail to collect scaffolds\n", __func__);	
			return 1;
		} 
	} else {
		col_ctgs(bam_fn[0], ctgs);	
	}
	if (!ctgs) {
		fprintf(stderr, "[E::%s] fail to collect contigs\n", __func__);	
		return 1;
	} 

#ifdef VERBOSE
	fprintf(stderr, "[M::%s] processing bam file\n", __func__);
#endif
	
	hit2_ary_t *hit2_ary = calloc(1, sizeof(hit2_ary_t));
	int i;	
	for ( i = 0; i < n_bams; ++i) {
		if (mb_proc_bam(bam_fn[i], min_mq, ctgs, scfs, hit2_ary)) {
			return 1;	
		}	
	}
	free(hit2_ary->ary); free(hit2_ary);
	sd_destroy(ctgs);
	sd_destroy(scfs);
	return 0;
}


int main(int argc, char *argv[])
{
	int c;
	char *sat_fn = 0, *links_fn = 0, *out_fn = 0;
	int limn = 2;
	char *program = argv[0];
	/*--argc, ++argv;*/
	while (~(c=getopt(argc, argv, "s:n:o:h"))) {
		switch (c) {
			case 's':
				sat_fn = optarg;
				break;	
			case 'n':
				limn = atoi(optarg);
				break;	
			case 'o':
				out_fn = (optarg);
				break;	
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
help:	
				fprintf(stderr, "\nUsage: %s [<options>] <BAMs> ...\n", program);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -s    FILE     SAT file [null]\n");
				fprintf(stderr, "         -n    INT      allow successor or predecessor be top N candidates [2]\n");
				fprintf(stderr, "         -o    FILE     output file [stdout]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		}		
	}
	if (optind + 1 > argc) {
		fprintf(stderr,"[E::%s] require bam file(s)!\n", __func__); goto help;
	}
	/*links_fn = argv[optind++];*/
	char **bam_fn = argv + optind;
	int n_bams = argc - optind;	
	int min_mq = 10;
	/*fprintf(stderr, "%d\t%s\t%s\n", n_bams, bam_fn[0], bam_fn[1]);*/
	fprintf(stderr, "Program starts\n");	
	bam2hig_bed(sat_fn, bam_fn, n_bams, min_mq, out_fn);
/*int mk_brks(char *sat_fn, char *links_fn, int limn, char *out_fn)*/
	/*mk_brks(sat_fn, links_fn, limn, out_fn);*/
	fprintf(stderr, "Program ends\n");	
	return 0;	
}


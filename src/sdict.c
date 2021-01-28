#include <stdio.h>
#include <string.h>
#include "sdict.h"

#include "khash.h"

KHASH_MAP_INIT_STR(str, uint32_t)
typedef khash_t(str) shash_t;

sdict_t *sd_init(void)
{
	sdict_t *d;
	d = (sdict_t*)calloc(1, sizeof(sdict_t));
	d->h = kh_init(str);
	return d;
}

void sd_destroy(sdict_t *d)
{
	uint32_t i;
	if (d == 0) return;
	if (d->h) kh_destroy(str, (shash_t*)d->h);
	for (i = 0; i < d->n_seq; ++i)
		free(d->seq[i].name);
	free(d->seq);
	free(d);
}

int32_t sd_put(sdict_t *d, const char *name)
{
	shash_t *h = (shash_t*)d->h;
	khint_t k;
	int absent;
	k = kh_put(str, h, name, &absent);
	if (absent) {
		sd_seq_t *s;
		if (d->n_seq == d->m_seq) {
			d->m_seq = d->m_seq? d->m_seq<<1 : 16;
			d->seq = (sd_seq_t*)realloc(d->seq, d->m_seq * sizeof(sd_seq_t));
		}
		s = &d->seq[d->n_seq];
		kh_key(h, k) = s->name = strdup(name);
		kh_val(h, k) = d->n_seq++;
	} // TODO: test if len is the same;
	return kh_val(h, k);
}

int32_t sd_put3(sdict_t *d, const char *name, uint32_t len, uint32_t le, uint32_t rs, uint32_t snp_n, uint32_t l_snp_n, uint32_t r_snp_n)
{
	shash_t *h = (shash_t*)d->h;
	khint_t k;
	int absent;
	k = kh_put(str, h, name, &absent);
	if (absent) {
		sd_seq_t *s;
		if (d->n_seq == d->m_seq) {
			d->m_seq = d->m_seq? d->m_seq<<1 : 16;
			d->seq = (sd_seq_t*)realloc(d->seq, d->m_seq * sizeof(sd_seq_t));
		}
		s = &d->seq[d->n_seq];
		s->le = le, s->rs = rs;
		s->snp_n = snp_n, s->l_snp_n = l_snp_n, s->r_snp_n = r_snp_n;
		s->len = len;
		kh_key(h, k) = s->name = strdup(name);
		kh_val(h, k) = d->n_seq++;
	} // TODO: test if len is the same;
	else {
		uint32_t ind = kh_val(h, k);
		sd_seq_t *s = &d->seq[ind];
		if (len) s->len = len;
		if (snp_n) s->snp_n = snp_n;
	   	if (l_snp_n) s->l_snp_n = l_snp_n;
		if (r_snp_n) s->r_snp_n = r_snp_n;	
		if (rs) s->rs = rs;
		if (le) s->le = le;	
	}
	return kh_val(h, k);
}
int32_t sd_put4(sdict_t *d, const char *name, uint32_t len, uint32_t le, uint32_t rs, uint32_t l_snp_n, uint32_t r_snp_n, uint32_t is_circ)
{
	shash_t *h = (shash_t*)d->h;
	khint_t k;
	int absent;
	k = kh_put(str, h, name, &absent);
	if (absent) {
		sd_seq_t *s;
		if (d->n_seq == d->m_seq) {
			d->m_seq = d->m_seq? d->m_seq<<1 : 16;
			d->seq = (sd_seq_t*)realloc(d->seq, d->m_seq * sizeof(sd_seq_t));
		}
		s = &d->seq[d->n_seq];
		s->le = le, s->rs = rs;
		s->l_snp_n = l_snp_n, s->r_snp_n = r_snp_n;
		s->len = len;
		kh_key(h, k) = s->name = strdup(name);
		kh_val(h, k) = d->n_seq++;
		s->is_circ = is_circ;
	} // TODO: test if len is the same;
	else {
		uint32_t ind = kh_val(h, k);
		sd_seq_t *s = &d->seq[ind];
		if (len) s->len = len;
	   	if (l_snp_n) s->l_snp_n = l_snp_n;
		if (r_snp_n) s->r_snp_n = r_snp_n;	
		if (rs) s->rs = rs;
		if (le) s->le = le;	

		s->is_circ = is_circ;
	}
	return kh_val(h, k);
}
int32_t sd_put2(sdict_t *d, const char *name, uint32_t len, uint32_t le, uint32_t rs, uint32_t l_snp_n, uint32_t r_snp_n)
{
	shash_t *h = (shash_t*)d->h;
	khint_t k;
	int absent;
	k = kh_put(str, h, name, &absent);
	if (absent) {
		sd_seq_t *s;
		if (d->n_seq == d->m_seq) {
			d->m_seq = d->m_seq? d->m_seq<<1 : 16;
			d->seq = (sd_seq_t*)realloc(d->seq, d->m_seq * sizeof(sd_seq_t));
		}
		s = &d->seq[d->n_seq];
		s->le = le, s->rs = rs;
		s->l_snp_n = l_snp_n, s->r_snp_n = r_snp_n;
		s->len = len;
		kh_key(h, k) = s->name = strdup(name);
		kh_val(h, k) = d->n_seq++;
	} // TODO: test if len is the same;
	else {
		uint32_t ind = kh_val(h, k);
		sd_seq_t *s = &d->seq[ind];
		if (len) s->len = len;
	   	if (l_snp_n) s->l_snp_n = l_snp_n;
		if (r_snp_n) s->r_snp_n = r_snp_n;	
		if (rs) s->rs = rs;
		if (le) s->le = le;	

	}
	return kh_val(h, k);
}

int32_t sd_get(const sdict_t *d, const char *name)
{
	shash_t *h = (shash_t*)d->h;
	khint_t k;
	k = kh_get(str, h, name);
	return k == kh_end(h)? -1 : kh_val(h, k);
}

void sd_hash(sdict_t *d)
{
	uint32_t i;
	shash_t *h;
	if (d->h) return;
	d->h = h = kh_init(str);
	for (i = 0; i < d->n_seq; ++i) {
		int absent;
		khint_t k;
		k = kh_put(str, h, d->seq[i].name, &absent);
		kh_val(h, k) = i;
	}
}

int32_t *sd_squeeze(sdict_t *d)
{
	int32_t *map, i, j;
	if (d->h) {
		kh_destroy(str, (shash_t*)d->h);
		d->h = 0;
	}
	map = (int32_t*)calloc(d->n_seq, 4);
	for (i = j = 0; i < d->n_seq; ++i) {
		d->seq[j] = d->seq[i], map[i] = j++;
	}
	d->n_seq = j;
	sd_hash(d);
	return map;
}

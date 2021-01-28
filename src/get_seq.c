/*
 * =====================================================================================
 *
 *       Filename:  get_seq.c
 *
 *    Description:  get path sequence from gfa
 *
 *        Version:  1.0
 *        Created:  20/11/2018 13:12:55
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Dengfeng Guan (D. Guan), dfguan@hit.edu.cn
 *   Organization:  Center for Bioinformatics, Harbin Institute of Technology
 *
 * =====================================================================================
 */
#include <stdlib.h>
#include <getopt.h>
#include <string.h>

#include "graph.h"

int get_seq(char *sat_fn, char *seq_fn, uint32_t min_l, char *out_fn)
{
	
#ifdef VERBOSE
		fprintf(stderr, "[M::%s] load scaffolding graph to memory \n", __func__);
#endif
	graph_t *g = load_sat(sat_fn);
#ifdef VERBOSE
		fprintf(stderr, "[M::%s] get contigs to memory \n", __func__);
#endif
	if (seq_fn) read_seq(g, seq_fn);
#ifdef VERBOSE
		fprintf(stderr, "[M::%s] get scaffolds \n", __func__);
#endif
	get_path(g, min_l, out_fn);
	graph_destroy(g);	
	return 0;
}

int main_get_seq(int argc, char *argv[])
{

	char *seq_fn = 0, *out_fn = 0;
	uint32_t min_l = 0;	
	
	char *program;
   	(program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
	
	--argc, ++argv;
	int c;
	while (~(c = getopt(argc, argv, "c:o:l:h"))) {
		switch (c) {
			case 'c':
				seq_fn = optarg;
				break;
			case 'o':
				out_fn = optarg;
				break;
			case 'l':
				min_l = atoi(optarg);
				break;
			default:
				if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
		help:	
				fprintf(stderr, "\nUsage: %s %s [<options>] <SAT> ...\n", program, argv[0]);
				fprintf(stderr, "Options:\n");
				fprintf(stderr, "         -c    STR      fasta file\n");
				fprintf(stderr, "         -o    STR      output file [stdout]\n");
				fprintf(stderr, "         -l    STR      minimum output scaffolds length [0]\n");
				fprintf(stderr, "         -h             help\n");
				return 1;	
		
		}
	}

	if (optind + 1 > argc) {
		fprintf(stderr, "[E::%s] Require a SAT file", __func__);
		goto help;
	}
	char *sat_fn = argv[optind];
	fprintf(stderr, "[M::%s] program starts\n", __func__);
	get_seq(sat_fn, seq_fn, min_l, out_fn);
	fprintf(stderr, "[M::%s] program ends\n", __func__);
	return 0;


}


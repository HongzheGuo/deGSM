#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ubwt.h"
#include "utils.h"

int ubwt2unipath_usage(void)
{
    err_printf("\n");
    err_printf("Usage:    ubwt unipath [option] <BWT-STR>\n\n");
    err_printf("Options:\n\n");
    err_printf("    -t    [INT] Number of threads. [1]\n");
    err_printf("    -f    [STR] Format of input bwt-str. [B]\n"); // TODO only allow binary format
    err_printf("                  \"B\": binary file, 4-bit per bp, 0/1/2/3/4:A/C/G/T/#(first 64-bit: length).\n");
    err_printf("                  \"P\": plain text.\n");
    err_printf("    -e    [STR] Edge sequence file in binary format. Required when output GFA format. [NULL]\n");
    err_printf("    -k    [INT] Length of k-mer. Required when output GFA format.\n");
    err_printf("    -a    [STR] Format of output file. [F]\n");
    err_printf("                  \"F\": FASTA format.\n");
    err_printf("                  \"G\": GFA format.\n");
    err_printf("    -o    [STR] Output file. [stdout]\n");
    //err_printf("    -b          Output binary file. [False].\n");
    //err_printf("    -c    [INT] Chunk size of each run when t>1. [%d]\n", CHUNK_SIZE);
    err_printf("\n");
    return 1; 
}

int ubwt2unipath(int argc, char *argv[])
{
    int c, t=1, input_b = 1, b_out = 0, g_out = 0, k=0, chunk_size = CHUNK_SIZE;
    char *edge_fn=NULL, *out_fn=NULL; FILE *out=stdout;
    while ((c = getopt(argc, argv, "t:f:e:k:o:a:bc:")) >= 0) {
        switch (c) {
            case 't': t = atoi(optarg); break;
            case 'f': if (optarg[0] == 'P' || optarg[0] == 'p') input_b = 0; break;
            case 'e': edge_fn = strdup(optarg); break;
            case 'k': k = atoi(optarg); break;
            case 'o': out_fn = strdup(optarg); break;
            case 'a': if (optarg[0] == 'G' || optarg[0] == 'g') g_out = 1; break;
            case 'b': b_out = 1; break;
            case 'c': chunk_size = atoi(optarg); break;
            default: return ubwt2unipath_usage();
        }
    }
    if (argc - 1 != optind) return ubwt2unipath_usage();

    if (g_out && b_out) err_fatal_simple("Error: argument -g: not allowed with argument -b.\n");
    if (g_out && (k == 0 || edge_fn == NULL)) err_fatal_simple("Error: arguments -e and -k are required when -g is used.\n");
    if (out_fn) {
        if (b_out) out = fopen(out_fn, "wb");
        else out = fopen(out_fn, "w");
    }

    // read ubwt str
    char *bwt_fn = strdup(argv[optind]); ubwt_count_t ubwt_l;
    uint8_t *ubwt_bstr = ubwt_read_bwt_str(bwt_fn, input_b, &ubwt_l);

    err_func_format_printf(__func__, "Constructing ubwt ...\n");
    // init ubwt
    ubwt_t *ubwt = (ubwt_t*)_err_malloc(sizeof(ubwt_t));
    ubwt_init(ubwt, ubwt_l);
    // claculate occ, c
    ubwt_count_t uni_c = ubwt_cal(ubwt, ubwt_bstr, ubwt_l);
    ubwt_update(ubwt);

    // read edge sequence
    uint8_t *edge=NULL; ubwt_count_t edge_l;
    if (g_out) { // generate GFA: unipath, link
        ubwt_gen_map(ubwt, uni_c, t, chunk_size);
        err_func_format_printf(__func__, "Constructing ubwt done.\n");
        edge = ubwt_read_edge(edge_fn, &edge_l);
        if (edge_l != ubwt_l) err_fatal_simple("Error: lengths of bwt string and edge sequence do not match.");
        ubwt_gen_gfa(ubwt, edge, uni_c, k, out, t, chunk_size);
    } else { // generate unipath
        err_func_format_printf(__func__, "Constructing ubwt done.\n");
        // TODO binary unipath output foramt
        ubwt_gen_unipath(ubwt, uni_c, out, b_out, t, chunk_size);
    }

    free(bwt_fn); if (edge_fn) free(edge_fn); if (edge) free(edge); if (out) fclose(out);
    ubwt_free(ubwt); free(ubwt_bstr);
    return 0;
}

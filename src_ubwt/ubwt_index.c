#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ubwt.h"
#include "utils.h"

int ubwt_index_usage(void)
{
    err_printf("\n");
    err_printf("Usage:   ubwt index [option] <BWT-STR>\n\n");
    err_printf("Options:\n\n");
    err_printf("    -t    [INT] Number of threads. [1].\n");
    err_printf("    -f    [STR] Format of input bwt-str. [B]\n"); 
    err_printf("                  \"B\": binary file, 4-bit per bp, 0/1/2/3/4:A/C/G/T/#(first 64-bit: length).\n");
    err_printf("                  \"P\": plain text.\n");
    err_printf("    -o    [STR] Prefix of bwt-index file. [<BWT-STR>]\n");
    fprintf(stderr, "\n");
    return 0;
}

void ubwt_dump(ubwt_t *ubwt, const char *prefix)
{
    char *fn = (char*)_err_calloc(strlen(prefix)+15, sizeof(char));
    // dump bwt
    strcpy(fn, prefix); strcat(fn, ".ubwti");
    err_func_format_printf(__func__, "Writing ubwt index to %s ...\n", fn);
    FILE *fp = xopen(fn, "wb");
    ubwt_count_t zero = 0;
    err_fwrite(&zero, sizeof(ubwt_count_t), 1, fp);
    err_fwrite(&ubwt->ubwt_l, sizeof(ubwt_count_t), 1, fp);
    err_fwrite(&ubwt->ubwt_size, sizeof(ubwt_count_t), 1, fp);
    err_fwrite(ubwt->C, sizeof(ubwt_count_t), _OCC_C, fp);
    err_fwrite(ubwt->ubwt, sizeof(ubwt_int_t), ubwt->ubwt_size, fp);
    err_fwrite(&ubwt->uni_c, sizeof(ubwt_count_t), 1, fp);
    err_fwrite(ubwt->ubwt_map, sizeof(ubwt_count_t), ubwt->uni_c, fp);
    err_fflush(fp);
    err_fclose(fp);
    
    free(fn);
    err_func_format_printf(__func__, "Writing ubwt index done.\n");
}

static void fread_fix(FILE *fp, ubwt_count_t size, void *a)
{
	const int bufsize = 0x1000000; // 16M block
	ubwt_count_t offset = 0;
	while (size) {
		int x = bufsize < size? bufsize : size;
		if ((x = err_fread_noeof(a + offset, 1, x, fp)) == 0) break;
		size -= x; offset += x;
	}
	//return offset;
}

ubwt_t *ubwt_restore_index(const char *prefix)
{
    char *fn = (char*)_err_calloc(strlen(prefix)+15, sizeof(char));
    ubwt_t *ubwt;
    // restore bwt
    strcpy(fn, prefix); strcat(fn, ".ubwti");
    err_func_format_printf(__func__, "Restoring ubwt index from %s ...\n", fn);
    FILE *fp = xopen(fn, "rb");
    ubwt = (ubwt_t*)_err_calloc(1, sizeof(ubwt_t));

    ubwt_count_t zero;
    err_fread_noeof(&zero, sizeof(ubwt_count_t), 1, fp);
    if (zero != 0) err_fatal_simple("Error: wrong ubwt index format.");

    err_fread_noeof(&ubwt->ubwt_l, sizeof(ubwt_count_t), 1, fp);
    err_fread_noeof(&ubwt->ubwt_size, sizeof(ubwt_count_t), 1, fp);
    err_fread_noeof(ubwt->C, sizeof(ubwt_count_t), _OCC_C, fp);
    ubwt->ubwt = (ubwt_int_t*)_err_calloc(ubwt->ubwt_size, sizeof(ubwt_int_t));
    fread_fix(fp, ubwt->ubwt_size * sizeof(ubwt_int_t), ubwt->ubwt);
    err_fread_noeof(&ubwt->uni_c, sizeof(ubwt_count_t), 1, fp);
    ubwt->ubwt_map = (ubwt_count_t*)_err_calloc(ubwt->uni_c, sizeof(ubwt_count_t));
    fread_fix(fp, ubwt->uni_c * sizeof(ubwt_count_t), ubwt->ubwt_map);
    err_fclose(fp); 

    // count table16
    ubwt_gen_bit_table16(ubwt);
    //debwt_gen_cnt_table8(db);

    free(fn);
    err_func_format_printf(__func__, "Restoring ubwt index done.\n");
    return ubwt;
}

int ubwt_index(int argc, char *argv[])
{
    int c, input_b=1, t=1; char *prefix=NULL;
    while ((c = getopt(argc, argv, "t:f:o:")) >= 0) {
        switch (c) {
            case 't': t = atoi(optarg); break;
            case 'f': if (optarg[0] == 'P') input_b = 0; break;
            case 'o': prefix = optarg; break;
            default: return ubwt_index_usage();
        }
    }
    if (optind + 1 > argc) return ubwt_index_usage();
    char *fn = strdup(argv[optind]); ubwt_count_t ubwt_l;
    if (prefix == NULL) prefix = fn;

    // read ubwt str
    uint8_t *ubwt_bstr = ubwt_read_bwt_str(fn, input_b, &ubwt_l);

    err_func_format_printf(__func__, "Constructing ubwt ...\n");
    // init ubwt
    ubwt_t *ubwt = (ubwt_t*)_err_malloc(sizeof(ubwt_t));
    ubwt_init(ubwt, ubwt_l);
    // claculate occ, c bwt_map
    ubwt->uni_c = ubwt_cal(ubwt, ubwt_bstr, ubwt_l);
    ubwt_update(ubwt);
    ubwt_gen_map(ubwt, ubwt->uni_c, t, CHUNK_SIZE);
    err_func_format_printf(__func__, "Constructing ubwt done.\n");

    // dump bwt index to disk
    ubwt_dump(ubwt, prefix); 

    ubwt_free(ubwt); free(ubwt_bstr);
    free(fn); 
    return 0;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ubwt2unipath.h"
#include "ubwt_index.h"
#include "ubwt_query.h"
#include "utils.h"

#define VERSION "1.0.0"
#define CONTACT "yangao07@hit.edu.cn"

int usage(void)
{
    err_printf("\n");
	err_printf("Program: ubwt (uni-path BWT utilities)\n");
	err_printf("Usage:   ubwt <command> [options]\n\n");
	err_printf("Commands: \n");
	err_printf("         unipath     generate uni-path sequence from BWT string\n");
    err_printf("         index       build BWT index for uni-path sequence\n");
	err_printf("         query       find exact match of query sequence with BWT index\n");
	err_printf("\n");
    err_printf("Version: %s\n", VERSION);
    err_printf("Contact: %s\n", CONTACT);
	err_printf("\n");
	return 1;
}

int main(int argc, char *argv[])
{
    if (argc < 2) return usage();
    if (strcmp(argv[1], "unipath") == 0)      return ubwt2unipath(argc-1, argv+1);
    else if (strcmp(argv[1], "index") == 0)   return ubwt_index(argc-1, argv+1);
	else if (strcmp(argv[1], "query") == 0)   return ubwt_query(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}

    return 0;
}

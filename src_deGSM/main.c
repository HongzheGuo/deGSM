#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <time.h>
#include <inttypes.h>
#include <string.h>
#include <getopt.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "main.h"
#include "bucket_sort.h"

#define VERSION "1.0"
#define	MAX_ROUTE_TMP	500
	
uint64_t trans[256];
int KMER_LENGTH_PlusTwo=32;//57
int KMER_LENGTH=30;//55
unsigned int memoryKmer=0;

char *get_bin_dir(char *bin);
void usage(void);

struct timeval t_start,t_end;

enum
{
	NO_COMPLEMENTARY_REVERSE
};

static const char *short_option = "t:k:m:l:u:o:q:gs";
static struct option long_option[] =
{
	{(char*)"noCom", no_argument,  0, NO_COMPLEMENTARY_REVERSE},
	{(char*)0, 0, 0, 0}
};

static int load_input_usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program:	deGSM: memory scalable construction of large scale de Bruijn Graph \n");
	fprintf(stderr, "Version:	%s\n", VERSION);
	fprintf(stderr, "Contact:	Hongzhe Guo <hzguo@hit.edu>\n\n");
	fprintf(stderr, "Usage:  	deGSM [options] <jellyfish_path> <output_file> <source_file1> [source_file2] ...\n\nOptions:\n");
	fprintf(stderr, "	-k INT	k-mer length [20-253], default %u\n", KMER_LENGTH);
	fprintf(stderr, "	-t INT	maximum thread number [1-32], default 8\n");
	fprintf(stderr, "	-m INT	max memory usage in De Bruijn Graph Building(G) [4-128], default 8G\n");
	fprintf(stderr, "	-l INT	Don't cal k-mer with count < lower-count\n");
	fprintf(stderr, "	-u INT	Don't cal k-mer with count > upper-count\n\n");
	fprintf(stderr, "	-q STR	the k-mer whose any base with quality under the character will be filtered out\n\n");
	fprintf(stderr, "	-g 	the format of input file is .gz format\n\n");
	fprintf(stderr, "	-s 	the format of input file is .sra format\n\n");
	fprintf(stderr, "	--noCom	Do not count the k-mers on complementary-reverse strand\n\n");
	//fprintf(stderr,"	Please make sure your sequence don't contain any uncertain characters like 'N'\n");
	fprintf(stderr, "	Input: sequence in fasta or fastq format\n");
	fprintf(stderr, "	Please refer to the following link for more detailed information about the options: https://github.com/hitbc/deGSM\n");
	fprintf(stderr, "\n");

	return 1;
}

int main(int argc, char *argv[])
{
	float cost_time = 0;
	long start = 0;
	long end = 0;
	uint64_t THREAD_NUM=8;
	int i;
	int tempM = 0;

	trans['A']=trans['a']=0;
	trans['C']=trans['c']=1;
	trans['G']=trans['g']=2;
	trans['T']=trans['t']=3;
	trans['#']=4;// patch for special characters
	trans['$']=5;

	char *jRoot=NULL;
	char *source = NULL;
	char real_path[1024];
        realpath(argv[0], real_path);
        char *bin = get_bin_dir(real_path);

	int c = 0;
	uint8_t com_flag = 0;
	uint8_t compress_flag = 0;
	
	memoryKmer = 8192;//1024 2048 4096

	flag_filter = 0;
	filter_min = 1;
	filter_max = 0Xffffff;

	char* tmp_route = calloc(500, 1);
	strcpy(tmp_route, bin);
	quality_flag = 0;
	quality_char = 'c';
	
	while((c = getopt_long(argc, argv, short_option, long_option, NULL)) != -1)
	{
		switch (c)
		{
		case 't':
			THREAD_NUM = atoi(optarg);
			if(THREAD_NUM==0) fprintf(stderr,"thread number must be a number!\n"), exit(1);
			break;
		case 'k':
			KMER_LENGTH = atoi(optarg);
			KMER_LENGTH_PlusTwo = KMER_LENGTH + 2;
			KMER_LENGTH_PlusOne = KMER_LENGTH + 1;
			//backup
			if(KMER_LENGTH<20||KMER_LENGTH_PlusTwo>253)
				fprintf(stderr, "-k: k-mer length (from 20 to 253, default 55)\n" ),exit(1);
			break;
		case 'm':
			tempM = atoi(optarg);
			if(tempM>0)
			{
				int j=0;
				while((optarg[j] < '9') && (optarg[j] > '0')) j++;

				if((optarg[j]=='G') || (optarg[j]=='g'))
					memoryKmer = tempM << 10;
				else if((optarg[j] == 'M') || (optarg[j] == 'm'))
					memoryKmer = tempM;
				else
					fprintf(stderr,"memory format: XXXS (XXX: number, S: G or M) (at least %dM)\n",12*(BUFFERSIZE>>20)),exit(1);

				if(memoryKmer < (12 * (BUFFERSIZE>>20)))  fprintf(stderr,"memory limitation at least %dM!\n", 12*(BUFFERSIZE>>20)),exit(1);
			}
			else fprintf(stderr,"memory limitation must > 0 (at least %dM)!\n", 12*(BUFFERSIZE>>20)),exit(1);

			break;
		case 'l':
			filter_min = atoi(optarg);
			if(filter_min < 2)
				fprintf(stderr, "-l: kmer lower-count should be more than one\n" ),exit(1);
			flag_filter = 1;
			break;
		case 'u':
			filter_max = atoi(optarg);
			if(filter_max < 2)
				fprintf(stderr, "-l: kmer lower-count should be more than one\n" ),exit(1);
			flag_filter = 1;
			break;
		case 'o':
			strcpy(tmp_route, optarg);
			break;
		case 'q':
			//strcpy(quality_char, optarg);
			quality_char = (*optarg);
			quality_flag = 1;
			break;
		case 'g':
			compress_flag = 1;
			break;	
		case 's':
			compress_flag = 2;
			break;	
		case NO_COMPLEMENTARY_REVERSE:
			com_flag = 1;
			break;
		default:
			return load_input_usage();
		}
	}

	if (argc - optind < 3) return load_input_usage();

	if(filter_min > filter_max)
		fprintf(stderr, "-l: kmer lower-count should not be more than upper-count\n" ),exit(1);

	filter_min--;
	filter_max++;

	jRoot = argv[optind];
	out_route = argv[optind + 1];

	source = argv[optind + 2];

	if(jRoot==NULL)
	{
		fprintf(stderr, "Must specify your jellyfish root to use deGSM\n" );
		exit(1);
	}

	char* timebuff = calloc(50, 1);

	////////////////////////////////////////////////////////////////////////////////
	char cmd[1024];
	fprintf(stderr,"---------------------------------------------------------------------------\n");
	fprintf(stderr,"run deBuild %s:\n",VERSION);
	fprintf(stderr,"sequence file or path: %s \n", source);
	fprintf(stderr,"output bwt file: %s\n",out_route);
	fprintf(stderr,"take %s as the tmp_path to store some temporary files, which will be removed in the end\n",tmp_route );
	fprintf(stderr,"use %lu-thread\n",THREAD_NUM);
	fprintf(stderr,"k-mer length: %d %d\n",KMER_LENGTH, KMER_LENGTH_PlusTwo);//KMER_LENGTH_PlusTwo
	fprintf(stderr,"memory size: %u\n", memoryKmer);
	
	fprintf(stderr,"try to run jellyfish in %s\n", jRoot);
	if(com_flag)
		fprintf(stderr,"Do not count the k-mers on complementary-reverse strand in Jellyfish\n");

	/////////////////////////////////////////jellyfish///////////////////////////////////////////////////
	
	sprintf(cmd, "bash %s/src_deGSM/kmercounting.sh %s %lu %s %s %d %u %u %c %u", bin, source, THREAD_NUM, tmp_route, jRoot, KMER_LENGTH_PlusTwo, com_flag, quality_flag, quality_char, compress_flag);

	//fprintf(stderr,"\ncmd: %s\n",cmd);
	
	strcpy(timebuff, "date +%H:%M:%S");
	system(timebuff);
	
	clock_t start_c = 0;
	clock_t start_c_tol = 0;
	clock_t end_c_tol = 0;
	clock_t end_c = 0;

	time_t mix_begin = 0;
	time_t mix_start = 0;
	time_t mix_end = 0;

	void *arg_collect[4];

#ifdef	JELLY_LFKMER
	if(system(cmd) != 0 )
	{
		fprintf(stderr, "\njellyfish exit abnormally.\nMaybe you need add jellyfish into LIBRARY PATH with export LD_LIBRARY_PATH=jellyfish_path/.libs\n");
		exit(1);
	}
	fprintf(stderr,"jellyfish done!\n");
	fprintf(stderr,"use deGSM construct unipaths string ...\n");
	
#endif
	
	FILE *fp_size = NULL;
	uint64_t filesize = 0;
	char path_tmp[MAX_ROUTE_TMP];
	strcpy(path_tmp, tmp_route);
	strcat(path_tmp, "/kmerInfo");
	
	fp_size = fopen(path_tmp, "rb");
	fseek(fp_size, 0L, SEEK_END);
	filesize = ftell(fp_size);
	if(fp_size)	fclose(fp_size);
	
	if(filesize < MIN_BLOCK_FILE)
	{
		THREAD_NUM = 1;
		fprintf(stderr, "\nrun in very small dataset\n");
	}
	
	
	
	arg_collect[0]=(void *)bin;
	arg_collect[1]=(void *)THREAD_NUM;
	arg_collect[2]=(void *)source;
	arg_collect[3]=(void *)tmp_route;

#ifdef	BEF_BRANCH
	gettimeofday(&t_start, NULL);
	start = ((long)t_start.tv_sec)*1000+(long)t_start.tv_usec/1000;

	bucket_sort((void *)arg_collect);

	strcpy(timebuff, "date +%H:%M:%S");
	system(timebuff);
	
	gettimeofday(&t_end, NULL);
	end = ((long)t_end.tv_sec)*1000+(long)t_end.tv_usec/1000;

	cost_time = end - start;
	fprintf(stderr, "bucketing and sorting time is %.2lf s\n", cost_time/1000);

#endif

	if(generateUnipath((void *)arg_collect)==1)
	{
		fprintf(stderr, "success generate unipaths!\n");
	}
	else
	{
		fprintf(stderr, "failed to generate unipaths!\n");
	}

	gettimeofday(&t_end, NULL);
	end = ((long)t_end.tv_sec)*1000+(long)t_end.tv_usec/1000;

	cost_time = end - start;
	fprintf(stderr, "deGSM De bruijn graph construction time is %.2lf s\n", cost_time/1000);

	strcpy(timebuff, "date +%H:%M:%S");
	system(timebuff);

	if(tmp_route)	free(tmp_route);
	
	return 0;
}

char *get_bin_dir(char *bin)
{
	char *end = strrchr(bin, '/');
	bin[end-bin] = '\0';
	return bin;
}

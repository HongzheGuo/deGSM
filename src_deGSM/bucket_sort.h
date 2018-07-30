#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include <sys/time.h>
#include <time.h>
#include <inttypes.h>
#include <string.h>
#include <math.h>
#include <pthread.h>
#include <stddef.h>

#define	MIN_BLOCK	100
#define	MIN_BLOCK_FILE	5000

#define MOD32 31
#define MOD64 63
#define BUFBITLEN 22
#define BUFFERSIZE (1<<BUFBITLEN)

#define PATH_LEN 1024
#define MAXU64 ((uint64_t)(-1))
#define ELIMINATE (((uint64_t)1<<((32-KMER_LENGTH)<<1))-1)
#define STKSIZ 40
#define CUTOFF 8   
//8
#define BUCKET_LENGTH  8//less than KMER_LENGTH_PlusTwo
#define BUCKET_CAPACITY (1<<(BUCKET_LENGTH<<1))

//(1 << 22)
#define BUFFERSIZENEW	(1 << 24)

#define BUFFERSIZENEWNEW	(1 << 28)

#define MAX_FA_LINE	500
#define MAX_FQ_LINE	1024
//2048
#define	MAX_FILE_NMAE	200

/*************************************************************/
#define	FIRST_LAST_KMER_EXTRACT
#define	JELLY_LFKMER
#define	BEF_BRANCH
#define	KEEP_TMP
#define	DEL_LFINFO
#define	FILTERFILE_STOP
#define	FILTERFIRSTMERGE_STOP	
/*************************************************************/

#define	GFA_COM

#define	TOMERGE_WRITE
#define	MEM_COPY_TEST
#define	BIT_NEW
#define	MULTI_COUNT
#define	BOUNDARY
#define	LAST_HASH_DIS

#define	FIRST_LAST_KMER


//need to be deleted
#define	ADD_TMP

#define	TWOBITTOTOFOURBIT
#define	COM_DEL

//debug
FILE* fp_linear;

#define copykmer(tmp_com, buff, start)\
	memcpy(tmp_com, buff + start, kmer_num_tmp);

//debug
#ifdef	DEBUG_FILTER0
FILE* debug_fp[32];
#endif

//#define	TIME_TEST
#define	MERGE_MODIFY

uint8_t	OneUnit;
uint64_t f_kmer;
uint64_t l_kmer;
uint64_t bfq_n;

typedef struct THREAD_ALI{
	uint8_t thread;
	uint8_t num;
} thread_data;

typedef struct THREAD_BRA{
	uint64_t* readp;
	uint8_t* readp_t;
	uint64_t* writep;
	int wta_n;
	uint8_t num;
} thread_branch;

#ifdef	BOUNDARY
uint64_t last_block_num;
#endif
uint64_t last_block_num_lfkmers;

FILE* fp_rk;
FILE* fp_fa;

uint8_t tempMove_new;
uint8_t tempMove_new2;
uint8_t tempMove_r_new2;
uint8_t tempMove_r_new3;
uint8_t temphash;
uint8_t first_bp_move;
uint8_t tempMove_new_r;
//uint8_t kmer_rem;
uint64_t kmer_rem;
uint8_t kmer_num_br;
uint8_t br_off;
uint8_t br_offm;
uint16_t kmer_num;
uint16_t kmer_num_new;
uint16_t kmer_num_filter;
uint64_t keyMove;
uint64_t KMask;
uint64_t KMasks;

uint8_t leftmove32;

char *out_route;

int KMER_LENGTH_PlusTwo;
int KMER_LENGTH_PlusOne;
int KMER_LENGTH;
unsigned int hashCount;
unsigned int memoryKmer;
uint64_t trans[256];

int16_t fa_line_len;
uint64_t line_c;

uint64_t BinarySearch(uint64_t mk, uint64_t *target, int64_t up);//Search the # and $ point
uint64_t convert(uint64_t mark);

char *getPath(char *bin, char *name);	
void decode_onum(uint64_t, int );

uint8_t lf_tage;
FILE* fp_lfk_l;
FILE* fp_lfk_r;
FILE* fp_lfk_b;
uint8_t** write_lf_buff;

uint8_t flag_filter;
uint32_t filter_min;
uint32_t filter_max;

uint8_t one_out;
uint8_t one_suf;
uint8_t one_sufp;

uint64_t firstls[12];
uint64_t type_one;
uint64_t type_two;
uint64_t sh_mask;

uint16_t kmer_num_filtermerge;
uint16_t kmer_num_fm;

char quality_char;
uint8_t quality_flag;

void getbranchkmer(uint64_t, char*);

#ifdef	BOUNDARY

#ifdef	GFA_COM
void multiTagcheck_addkmer(uint64_t , uint64_t , FILE* , FILE* , FILE*);
#else
void multiTagcheck_addkmer(uint64_t , uint64_t , FILE* , FILE* );
#endif

void multiMergeTagcheck_addkmer(uint64_t , uint64_t , FILE* , FILE* , FILE* );
#endif

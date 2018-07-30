
#define	SPLIT_MI
//(1 << 26)
#define	MAX_MWA_L	((uint64_t)1 << 25)
#define	MAX_ROUTE	500

uint64_t countRead;

inline int cmpMK2(uint64_t a, uint64_t b);
inline int cmpMK2_s(uint64_t* va, uint64_t* vb);
inline int cmpMK2_sm(uint64_t* va, uint64_t* vb);
inline int L_cmpMK2(uint64_t* , uint64_t* );
void decode(uint64_t obj);
int mergeK2(char *, char *, uint64_t , char* );

void *multiTagcheck(void *arg);
void *multiMergeTagcheck(void *arg);

void getbranchkmer(uint64_t, char*);
void last_bucket_sort(uint64_t, uint64_t );
void *last_multiCount(void* );
void *last_multiDistri(void* );
void *last_multiThreadSort(void* );
int L_cmpKmer(const void * , const void * );
int L_cmpKmers(const void * , const void * );

#define	MAX_KMER_SEQ_BIT	4

typedef struct kmer_seq{
	uint64_t seq[MAX_KMER_SEQ_BIT];
}kmer_s;

//k mer
#define getkmerno_psf(tmp_com, buff, start)\
	if(kmer_num > 1)\
	{\
		tmp_com[0] = (buff[start] << 2);\
		memcpy(tmp_com + 1, buff + start + 1, (kmer_num - 2) << 3);\
		tmp_com[kmer_num - 1] = buff[start + kmer_num - 1] & KMask;\
	}else	tmp_com[0] = ((buff[start] << 2) & KMask);

//k mer
#define getkmerno_ps(tmp_com, buff, start)\
	if(kmer_num > 1)\
	{\
		tmp_com[0] = (buff[start] << (tempMove_new + 2));\
		memcpy(tmp_com + 1, buff + start + 1, (kmer_num - 2) << 3);\
		tmp_com[kmer_num - 1] = buff[start + kmer_num - 1] & KMask;\
	}else	tmp_com[0] = ((buff[start] << 2) & KMask);
	
//k-1 mer
#define getkmerno_pssf(tmp_com, buff, start)\
	if(kmer_num > 1)\
	{\
		tmp_com[0] = (buff[start] << 2);\
		memcpy(tmp_com + 1, buff + start + 1, (kmer_num - 2) << 3);\
		tmp_com[kmer_num - 1] = buff[start + kmer_num - 1] & KMasks;\
	}else	tmp_com[0] = ((buff[start] << 2) & KMasks);

//k-1 mer
#define getkmerno_pss(tmp_com, buff, start)\
	if(kmer_num > 1)\
	{\
		tmp_com[0] = (buff[start] << (tempMove_new + 2));\
		memcpy(tmp_com + 1, buff + start + 1, (kmer_num - 2) << 3);\
		tmp_com[kmer_num - 1] = buff[start + kmer_num - 1] & KMasks;\
	}else	tmp_com[0] = ((buff[start] << 2) & KMasks);
		

//k mer + suffix 
#define getkmerno_p(tmp_com, buff, start)\
	tmp_com[0] = (buff[start] << (tempMove_new + 2));\
	memcpy(tmp_com + 1, buff + start + 1, (kmer_num - 1) << 3);\

//prefix + k mer
#define getkmerno_s(tmp_com, buff, start)\
	memcpy(tmp_com, buff + start, (kmer_num - 1) << 3);\
	tmp_com[kmer_num - 1] = buff[start + kmer_num - 1] & KMask;
	
	
//k mer tomerge
#define getkmerno_psm(tmp_com, buff, start)\
	tmp_com[0] = (buff[start] << 2);\
	memcpy(tmp_com + 1, buff + start + 1, (kmer_num - 1) << 3);	
	
//k-1 mer tomerge
#define getkmerno_pssm(tmp_com, buff, start)\
	tmp_com[0] = (buff[start] << 2);\
	memcpy(tmp_com + 1, buff + start + 1, (kmer_num - 2) << 3);\
	tmp_com[kmer_num - 1] = buff[start + kmer_num - 1] & KMasks;
	
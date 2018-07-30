//#include "collect#$.h"
//#include "INandOut.h"

#include "bucket_sort.h"
#include "geUnipath.h"

uint64_t *splitMK2=NULL;
uint64_t* bufMK2 = NULL;
uint64_t* bufMK2_m = NULL;
uint8_t* bufMK2_t = NULL;

pthread_rwlock_t alphabet[4];

uint64_t alphabet_num[5];
uint64_t alphabet_num_ori[5];
uint64_t alphabet_num_new[32][4];

FILE* fp_fy = NULL;
FILE* fp_ry = NULL;
FILE* fp_x = NULL;

uint64_t** write_filterin = NULL;
uint64_t** write_filterout = NULL;

uint64_t** write_kmermerge = NULL;
uint8_t** write_kmermerge_t = NULL;

#ifdef	GFA_COM
uint8_t** edges_t = NULL;
#endif

uint64_t** write_fy = NULL;
uint64_t** write_ry = NULL;
uint64_t** write_x = NULL;

uint32_t w_m_p[32];
uint32_t w_f_p[32];
uint32_t w_r_p[32];
uint32_t w_x_p[32];

uint32_t w_f_is[32];
uint32_t w_f_os[32];

uint32_t line_tol = 0;
uint64_t wta_n = 0;
uint64_t fy_n = 0;
uint64_t ry_n = 0;
uint64_t x_n = 0;

//maybe used by another file

uint64_t read_num = 0;
uint64_t* first_buffer = NULL;
uint64_t* second_buffer = NULL;
uint64_t** last_hashKmer_dis = NULL;
uint64_t** last_countKmer_dis = NULL;
uint64_t* last_hashKmer = NULL;
uint64_t l_tempMove = 44;
uint64_t l_maskBKT=(1<<(BUCKET_LENGTH<<1))-1;
uint64_t* last_countKmer = NULL;
uint64_t* l_segCount=NULL;

#ifdef	LAST_HASH_DIS
uint64_t* last_hashKmer_write = NULL;
#endif

uint8_t input_types[9];

struct timeval t_start1,t_end1;

uint64_t* countKmerf = NULL;
uint64_t** countKmer_disf = NULL;
uint64_t** hashKmer_disf = NULL;
uint64_t* hashKmerf = NULL;
uint64_t* segCountf=NULL;
uint64_t maskBKT;

void *multiCount_filter(void * );
void *multiDistri_filter(void * );
void *multiThreadSort_filter(void * );
void *multiTagcheck_filter(void * );
void *twobittofourbit(void * );
int cmpKmer_filter(const void * , const void * );
inline int cmpMK2_s_filter(uint64_t* , uint64_t* );
uint8_t filter_out_sort_tag(uint64_t , int, char* );
void multiTagcheck_addkmer_filter(uint64_t, uint64_t, uint64_t, uint64_t, uint64_t, FILE* );//, FILE*

#ifdef	GFA_COM
void multiTagcheck_filter_addkmer(uint64_t, uint64_t, FILE*, FILE*, FILE*);
#else
void multiTagcheck_filter_addkmer(uint64_t, uint64_t, FILE*, FILE*);
#endif

uint64_t BinarySearchf(uint64_t , uint64_t* , int64_t );

void Swap(uint64_t A[], int i, int j, uint8_t u)
{
	uint64_t tmp_com[kmer_num];

	memcpy(tmp_com, A + i*kmer_num, u);
	memcpy(A + i*kmer_num, A + j*kmer_num, u);
	memcpy(A + i*kmer_num, tmp_com, u);

	//uint64_t temp = A[i];
	//A[i] = A[j];
	//A[j] = temp;
}

void Swap_t(uint8_t A[], int i, int j)
{
	uint8_t temp = A[i];
	A[i] = A[j];
	A[j] = temp;
}

void BubbleSort(uint64_t* buf, uint8_t* buf_t, int n, uint8_t unit)
{
	int i = 0;
	int j = 0;

	for (j = 0; j < n - 1; j++)         // 每次最大元素就像气泡一样"浮"到数组的最后
	{
		for (i = 0; i < n - 1 - j; i++) // 依次比较相邻的两个元素,使较大的那个向后移
		{
			//if (buf[i] > buf[i + 1])        // 如果条件改成A[i] >= A[i + 1],则变为不稳定的排序算法
			if(cmpMK2_sm(buf + (i + 1) * kmer_num, buf + i * kmer_num))
			{
				Swap(buf, i, i + 1, unit);
				Swap_t(buf_t, i, i + 1);
			}
		}
	}

}

int generateUnipath(void *arg)
{
	void **argT=(void **)arg;
	char *bin=(char *)argT[0];
	uint64_t THREAD_NUM=(uint64_t)argT[1];
	char* tmp_route = (char *)argT[3];
	uint64_t i;

	///////////////////////////////////merge k+2-mer////////////////////////////////////////////////////
	splitMK2=(uint64_t *)calloc(THREAD_NUM+1,sizeof(uint64_t));

	mergeK2("/kmerInfo.", bin, THREAD_NUM, tmp_route);

	free(splitMK2);
	////////////////////////////////////////////////////////////////////////////////////////////////////
	return 1;
}

int mergeK2(char *prefix, char *bin, uint64_t THREAD_NUM, char* tmp_route)   //mergeK2("/kmerInfo.", bin, THREAD_NUM);
{
	float cost_time = 0;
	long start1 = 0;
	long end1 = 0;

	gettimeofday(&t_start1, NULL);
	start1 = ((long)t_start1.tv_sec)*1000+(long)t_start1.tv_usec/1000;

	clock_t start_c = 0, end_c = 0;
	time_t start_t = 0, end_t = 0;

#ifdef	BEF_BRANCH

	char cNum[4];
	char nameK2[40];
	char path_tmp[MAX_ROUTE];
	char path_whole_tmp[MAX_ROUTE];
	char path_whole_tmp_t[MAX_ROUTE];
	char path_edge_t[MAX_ROUTE];

	uint8_t wta_rem = 0;
	uint8_t kmer_num_tmp = 0;
	uint8_t kmer_num_tmp_new = 0;

	int temp = 0;
	int block_num = 0;
	int heapTail = 0;//create the heap
	int i = 0;
	int cdi1 = 0;
	int cdi2 = 0;
	int hone_num = 0;
	int bound_off = 0;
	int block_numf = 0;
	int check = 0;
	int heap_tmp = 0;
	int *heap=NULL;
	uint32_t num = 0;

	uint64_t anchor = 0;
	uint64_t tempK = 0;
	uint64_t bufferSize = 0;
	uint64_t pBuf_t[4];
	uint64_t pBufM = 0;
	uint64_t segment = 0;
	uint64_t pBufM_pre = 0;
	uint64_t bufferSizeM = 0;
#ifdef	SPLIT_MI
	uint64_t eliminate_mi = 0;
#endif
	uint64_t eliminate=~(ELIMINATE);//#define ELIMINATE (((uint64_t)1<<((32-KMER_LENGTH)<<1))-1) 29
	uint64_t tmp_comf[kmer_num];
	uint64_t tmp_coms[kmer_num];


	FILE** fpK2 = NULL;
	FILE* fp_tomerge = NULL;
	FILE* fp_tomerge_t = NULL;

#ifdef	GFA_COM
	FILE* fp_edges_t = NULL;
#endif

	FILE* fp_filterin = NULL;
	FILE* fp_kmermergewhole[4];
	FILE* fp_kmermergewhole_t[4];

	pthread_t* myThread = NULL;

	kmer_num_tmp = (kmer_num << 3);

	kmer_num_tmp_new = (kmer_num_filter << 3);

	wta_rem = ((KMER_LENGTH+1) & 0Xf);
	if(wta_rem)	wta_n = ((KMER_LENGTH+1) >> 4) + 1;
	else	wta_n = ((KMER_LENGTH+1) >> 4);

	memset(input_types, 0, 9);
	input_types[2] = 1;
	input_types[4] = 2;
	input_types[8] = 3;

	heap = (int *)calloc((hashCount > 4 ? hashCount : 4), sizeof(int));

	strcpy(path_tmp, tmp_route);
	strcat(path_tmp, "/kmerFY");
	fp_fy = fopen(path_tmp,"wb");
	if(fp_fy == NULL) fprintf(stderr,"Fail to open kmerFY\n"),exit(1);

	strcpy(path_tmp, tmp_route);
	strcat(path_tmp, "/kmerRY");
	fp_ry = fopen(path_tmp,"wb");
	if(fp_ry == NULL) fprintf(stderr,"Fail to open kmerRY\n"),exit(1);

	strcpy(path_tmp, tmp_route);
	strcat(path_tmp, "/kmerX");
	fp_x = fopen(path_tmp,"wb");
	if(fp_x == NULL) fprintf(stderr,"Fail to open kmerX\n"),exit(1);

	memset(alphabet_num_new, 0, 1024);

#ifdef	SPLIT_MI
	eliminate_mi = (eliminate << 2);
#endif

	fprintf(stderr, "wta_n: %"PRId64"\n", wta_n);

#ifdef	FILTERFILE_STOP

	if(flag_filter)
	{
		if(hashCount == 1)	bufferSize = (uint64_t )((double)((uint64_t)memoryKmer<<20) / (double)((THREAD_NUM * 9) * kmer_num_filtermerge));
		else	bufferSize = (uint64_t )((double)((uint64_t)memoryKmer<<20) / (double)((THREAD_NUM * 9 + hashCount) * kmer_num_filtermerge));
	}
	else	bufferSize = (uint64_t )((double)((uint64_t)memoryKmer<<20) / (double)(THREAD_NUM * (kmer_num_tmp + kmer_num_tmp_new + 1) + kmer_num_tmp_new * hashCount));

	fprintf(stderr, "bufferSize: %"PRId64" %u %u\n", bufferSize, filter_min, filter_max);//2651214

	bufferSizeM = bufferSize * THREAD_NUM;

	uint64_t readNum[hashCount];
	uint64_t pBuf[hashCount];
	uint64_t* bufK2[hashCount];
	char *pathK2[hashCount];

	fpK2 = calloc(hashCount, sizeof(FILE* ));


	if(flag_filter)
	{
		strcpy(path_tmp, tmp_route);
		strcat(path_tmp, "/kmerfilterin");

		fp_filterin = fopen(path_tmp, "wb");
		if(fp_filterin==NULL)	fprintf(stderr, "Fail to open file kmerfilterin"), exit(1);

		write_filterin = (uint64_t** )calloc(THREAD_NUM, sizeof(uint64_t* ));
		for(i=0; i<THREAD_NUM; i++)
			if((write_filterin[i] = (uint64_t* )calloc(bufferSize, kmer_num_filtermerge)) == NULL)
				fprintf(stderr, "Fail to allocate memory for write_filterin\n"), exit(1);

		write_filterout = (uint64_t** )calloc(THREAD_NUM, sizeof(uint64_t* ));
		for(i=0; i<THREAD_NUM; i++)
			if((write_filterout[i] = (uint64_t* )calloc(bufferSize << 3, kmer_num_filtermerge)) == NULL)
				fprintf(stderr, "Fail to allocate memory for write_filterout\n"), exit(1);
	}
	else
	{
		strcpy(path_whole_tmp, tmp_route);
		strcat(path_whole_tmp, "/kmerMergewhole");

		strcpy(path_whole_tmp_t, tmp_route);
		strcat(path_whole_tmp_t, "/kmerMergewhole_t");

		fp_tomerge = fopen(path_whole_tmp, "wb");
		if(fp_tomerge==NULL)	fprintf(stderr, "Fail to open file kmerMergewhole"), exit(1);

		fp_tomerge_t = fopen(path_whole_tmp_t, "wb");
		if(fp_tomerge_t==NULL)	fprintf(stderr, "Fail to open file kmerMergewhole_t"), exit(1);

#ifdef	GFA_COM
		strcpy(path_edge_t, tmp_route);
		strcat(path_edge_t, "/edges_t");
		fp_edges_t = fopen(path_edge_t, "wb");
		if(fp_edges_t==NULL)	fprintf(stderr, "Fail to open file edges_t"), exit(1);
#endif

		write_kmermerge = (uint64_t** )calloc(THREAD_NUM, sizeof(uint64_t* ));
		for(i=0; i<THREAD_NUM; i++)
			if((write_kmermerge[i] = (uint64_t* )calloc(bufferSize, kmer_num_tmp)) == NULL)
				fprintf(stderr, "Fail to allocate memory for write_kmermerge\n"), exit(1);

		write_kmermerge_t = (uint8_t** )calloc(THREAD_NUM, sizeof(uint8_t* ));
		for(i=0; i<THREAD_NUM; i++)
			if((write_kmermerge_t[i] = (uint8_t* )calloc(bufferSize, 1)) == NULL)
				fprintf(stderr, "Fail to allocate memory for write_kmermerge_t\n"), exit(1);

#ifdef	GFA_COM
		edges_t = (uint8_t** )calloc(THREAD_NUM, sizeof(uint8_t* ));
		for(i=0; i<THREAD_NUM; i++)
			if((edges_t[i] = (uint8_t* )calloc(bufferSize, 1)) == NULL)
				fprintf(stderr, "Fail to allocate memory for edges_t\n"), exit(1);
#endif

	}

	fprintf(stderr, "First %u-way merging and to identify nodes type using %"PRId64" threads\n", hashCount, THREAD_NUM);

	if(hashCount == 1)
	{
		strcpy(path_tmp, tmp_route);//bin
		strcat(path_tmp, prefix);
		strcat(path_tmp, "0");

		fpK2[0]=fopen(path_tmp,"rb");
		if(fpK2[0]==NULL) fprintf(stderr,"Fail to open kmerInfo.0\n"),exit(1);

		if(flag_filter)
		{
			fseek(fpK2[0], 0, SEEK_END);
			uint64_t fi_size=ftell (fpK2[0]);
			fseek(fpK2[0], 0, SEEK_SET);

			hone_num = (int)((double)(fi_size) / (double)(bufferSizeM * kmer_num_tmp_new));
			if((fi_size % (bufferSizeM * kmer_num_tmp_new)) == 0)	--hone_num;

			fprintf(stderr, "hone_num: %d\n", hone_num);
		}
	}
	else
	{

		for(num=0; num<hashCount; num++)
		{
			sprintf(cNum,"%u",num);
			strcpy(nameK2, tmp_route);
			strcat(nameK2, prefix);
			strcat(nameK2, cNum);

			bufK2[num]=(uint64_t* )calloc(bufferSize+1, kmer_num_tmp_new);

			fpK2[num]=fopen(nameK2,"rb");//pathK2[num]
			if(fpK2[num]==NULL) fprintf(stderr,"Fail to open kmerInfo.%u\n",num),exit(1);

			readNum[num]=fread(bufK2[num], kmer_num_tmp_new, bufferSize, fpK2[num]);

			pBuf[num]=0;
		}
		for(num=0; num<hashCount; num++)
		{
			heapTail=num;
			for(i=heapTail; i>0; i=(i-1)>>1)
			{
				if(cmpMK2_s(bufK2[num], bufK2[heap[(i-1)>>1]]))
				{
					heap[i]=heap[(i-1)>>1];
				}
				else break;
			}
			heap[i]=num;
		}
	}

	bufMK2 = (uint64_t* )calloc(bufferSizeM + 128, kmer_num_tmp_new);

	if(flag_filter)
	{
		segCountf = (uint64_t* )calloc(THREAD_NUM+1,sizeof(uint64_t));
		segCountf[THREAD_NUM]=BUCKET_CAPACITY;
	}

	while(heapTail>=0)
	{
		pBufM=0;

		if(hashCount == 1)
		{
			pBufM = fread(bufMK2 + bound_off * kmer_num_filter, kmer_num_tmp_new, bufferSizeM, fpK2[0]);
			if(pBufM == 0)	break;

			--hone_num;
			fprintf(stderr, "hashCount==1 %"PRId64"\n", pBufM);
		}
		else
		{
			while(pBufM<bufferSizeM&&heapTail>=0)   //get bufMK2
			{
				num=heap[0];

				memcpy(bufMK2 + (pBufM + bound_off) * kmer_num_filter, bufK2[num] + pBuf[num] * kmer_num_filter, kmer_num_tmp_new);

				pBufM++;
				++pBuf[num];

				if(pBuf[num]==readNum[num])
				{
					readNum[num] = fread(bufK2[num], kmer_num_tmp_new, bufferSize, fpK2[num]);

					pBuf[num]=0;

					if(readNum[num]==0)
					{
						fclose(fpK2[num]);
						free(bufK2[num]);

						num=heap[heapTail];
						heapTail--;
					}
				}
				for(i=0; (i<<1)+1<=heapTail; i=temp)
				{
					cdi1 = (i<<1)+2;
					cdi2 = cdi1-1;

					temp = (((cdi1<=heapTail) && cmpMK2_s(bufK2[heap[cdi1]] + pBuf[heap[cdi1]] * kmer_num_filter, bufK2[heap[cdi2]] + pBuf[heap[cdi2]] * kmer_num_filter)) ? cdi1:cdi2);

					heap_tmp = heap[temp];
					if(heap_tmp == num)	heap[i]	= num;
					else
					{
						if(cmpMK2_s(bufK2[heap_tmp] + pBuf[heap_tmp] * kmer_num_filter, bufK2[num] + pBuf[num] * kmer_num_filter))
						{
							heap[i]=heap_tmp;
						}
						else break;
					}
				}
				heap[i]=num;

				if(pBufM>=bufferSizeM-15)
				{
					if(OneUnit)
					{
						if(((bufMK2[(pBufM-1 + bound_off) * kmer_num_filter]<<4)&eliminate)!=((bufK2[heap[0]][pBuf[heap[0]] * kmer_num_filter]<<4)&eliminate))
							break;
					}
					else
					{
						getkmerno_ps(tmp_comf, bufMK2, (pBufM-1 + bound_off) * kmer_num_filter)
						getkmerno_ps(tmp_coms, bufK2[heap[0]], pBuf[heap[0]] * kmer_num_filter)
						if(memcmp(tmp_comf, tmp_coms, kmer_num_tmp))	break;
					}
				}
			}
		}

		//use bufMK2 to generate tags

		if(pBufM == 0)	break;

#ifdef	BOUNDARY
		//update pBufM
		pBufM += bound_off;
		//
#endif
		segment=pBufM/THREAD_NUM;//4194289

			
		splitMK2[0]=0; //lower bound (can reach)
		for(i=1; i<THREAD_NUM; i++)
		{
			splitMK2[i]=i*segment;

			//if(kmer_num == 1)
			if(OneUnit)
			{
				tempK = (bufMK2[splitMK2[i] * kmer_num_filter]<<4)&eliminate_mi;

				while(((bufMK2[(splitMK2[i]+1) * kmer_num_filter]<<4)&eliminate_mi)==tempK)
				{
					splitMK2[i]++;
				}
			}
			else
			{
				getkmerno_pss(tmp_comf, bufMK2, splitMK2[i] * kmer_num_filter)

				while(1)
				{
					getkmerno_pss(tmp_coms, bufMK2, (splitMK2[i] + 1) * kmer_num_filter)
					if(memcmp(tmp_comf, tmp_coms, kmer_num_tmp))	break;
					else	splitMK2[i]++;
				}
			}



			splitMK2[i]++;

		}
		splitMK2[THREAD_NUM]=pBufM;

		myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));

		thread_data* tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));
		for(i=0; i<THREAD_NUM; i++)
		{
			w_m_p[i] = 0;
			w_f_is[i] = 0;
			w_f_os[i] = 0;

			tt[i].thread = THREAD_NUM;
			tt[i].num = i;
			check=pthread_create(&myThread[i], NULL, multiTagcheck, tt + i);

			if(check)
			{
				fprintf(stderr, "Thread Num: %d, Error: pthread create return code: %d\n", i, check);
				exit(EXIT_FAILURE);
			}
		}

		for(i=0; i<THREAD_NUM; i++)
		{
			pthread_join(myThread[i], NULL);
		}
		free(myThread);
		free(tt);

		if(flag_filter)
		{
			if(OneUnit)
				for(i=0; i<THREAD_NUM; i++)
					fwrite(write_filterin[i], kmer_num << 4, w_f_is[i], fp_filterin);
			else
				for(i=0; i<THREAD_NUM; i++)
					fwrite(write_filterin[i], kmer_num_new << 3, w_f_is[i], fp_filterin);
		}
		else
		{
			for(i=0; i<THREAD_NUM; i++)
				fwrite(write_kmermerge[i], kmer_num_tmp, w_m_p[i], fp_tomerge);

			for(i=0; i<THREAD_NUM; i++)
				fwrite(write_kmermerge_t[i], 1, w_m_p[i], fp_tomerge_t);

#ifdef	GFA_COM
			for(i=0; i<THREAD_NUM; i++)
				fwrite(edges_t[i], 1, w_m_p[i], fp_edges_t);
#endif

		}

		//kmer-count filter to do kmers-filtered out sorting
		if(flag_filter)
		{
			if(hashCount == 1)
			{
				if(hone_num >= 0)
					if(!(check = filter_out_sort_tag(THREAD_NUM, block_numf, tmp_route)))	block_numf++;
			}
			else
			{
				if((pBufM < bufferSizeM + 128) && (heapTail >= 0))
					if(!(check = filter_out_sort_tag(THREAD_NUM, block_numf, tmp_route)))	block_numf++;
			}
		}

#ifdef	BOUNDARY
		//get num of last k-1 mers and copy these kmers to new bufMK2
		if(heapTail >= 0)
		{
			bound_off = pBufM - last_block_num;
			fprintf(stderr, "BOUNDARY: %d %"PRId64" %"PRId64"\n", bound_off, pBufM, last_block_num);
			memcpy(bufMK2, bufMK2 + (pBufM - bound_off) * kmer_num_filter, bound_off * kmer_num_tmp_new);
		}
#endif
		pBufM_pre = pBufM;
		//

		fprintf(stderr, "Merging block %d and multiTagcheck finish\n", block_num + 1);

		block_num++;
	}

	//to modify
	if(flag_filter)
		multiTagcheck_addkmer_filter(pBufM_pre, last_block_num, THREAD_NUM, w_f_is[THREAD_NUM - 1], w_f_os[THREAD_NUM - 1], fp_filterin);
#ifdef	GFA_COM
	else	multiTagcheck_addkmer(pBufM_pre, last_block_num, fp_tomerge, fp_tomerge_t, fp_edges_t);
#else
	else	multiTagcheck_addkmer(pBufM_pre, last_block_num, fp_tomerge, fp_tomerge_t);
#endif

	//kmer-count filter to do kmers-filtered out sorting
	if(flag_filter)
	{
		if(!(check = filter_out_sort_tag(THREAD_NUM, block_numf, tmp_route)))	block_numf++;
	}

	free(bufMK2);

	if(flag_filter)
	{
		//modify judge whether the filter-in is empty
		strcpy(nameK2, tmp_route);
		strcat(nameK2, "/kmerfilterin");
		FILE* fp_filterin_test = fopen(nameK2, "rb");
		fseek(fp_filterin_test, 0, SEEK_END);
		if(ftell(fp_filterin_test) == 0)
			fprintf(stderr, "\n\nAll kmers are filtered out, exit\n\n"), exit(1);
		if(fp_filterin_test) fclose(fp_filterin_test);
		//
		fclose(fp_filterin);

		for(i=0; i<THREAD_NUM; i++)
			if(write_filterin[i])	free(write_filterin[i]);
		if(write_filterin)	free(write_filterin);

		for(i=0; i<THREAD_NUM; i++)
			if(write_filterout[i])	free(write_filterout[i]);
		if(write_filterout)	free(write_filterout);

	}

	/**/
	for(i=0; i<hashCount; i++)
	{
		strcpy(path_tmp, "rm -f ");
		//strcat(path_tmp, bin);
		strcat(path_tmp, tmp_route);
		strcat(path_tmp, "/kmerInfo.");
		sprintf(cNum,"%d", i);
		strcat(path_tmp, cNum);
		system(path_tmp);
	}


	gettimeofday(&t_end1, NULL);
	end1 = ((long)t_end1.tv_sec)*1000+(long)t_end1.tv_usec/1000;

	cost_time = end1 - start1;
	fprintf(stderr, "First Merging time is %.2lf s\n", cost_time/1000);

#endif

#ifdef	FILTERFIRSTMERGE_STOP
	if(flag_filter)// && filter_out_f
	{
		gettimeofday(&t_start1, NULL);
		start1 = ((long)t_start1.tv_sec)*1000+(long)t_start1.tv_usec/1000;

		block_num = 0;

		fprintf(stderr, "Begin kmer-count filter files merging and do multi-Tagcheck: %d\n", block_numf);

		//cal bufferSize
		bufferSize = (uint64_t )((double)((uint64_t)memoryKmer<<20) / (double)(THREAD_NUM * (kmer_num_tmp + kmer_num_filtermerge + 1) + kmer_num_filtermerge * block_numf));

		bufferSizeM = bufferSize * THREAD_NUM;

		int block_filter = 1;
		strcpy(path_whole_tmp, tmp_route);
		strcat(path_whole_tmp, "/kmerMergewhole");

		strcpy(path_whole_tmp_t, tmp_route);
		strcat(path_whole_tmp_t, "/kmerMergewhole_t");

		FILE* fp_tomerge = fopen(path_whole_tmp, "wb");
		if(fp_tomerge==NULL)	fprintf(stderr, "Fail to open file kmerMergewhole"), exit(1);

		FILE* fp_tomerge_t = fopen(path_whole_tmp_t, "wb");
		if(fp_tomerge_t==NULL)	fprintf(stderr, "Fail to open file kmerMergewhole_t"), exit(1);

#ifdef	GFA_COM
		strcpy(path_edge_t, tmp_route);
		strcat(path_edge_t, "/edges_t");
		FILE* fp_edges_t = fopen(path_edge_t, "wb");
		if(fp_edges_t==NULL)	fprintf(stderr, "Fail to open file edges_t"), exit(1);
#endif

		write_kmermerge = (uint64_t** )calloc(THREAD_NUM, sizeof(uint64_t* ));
		for(i=0; i<THREAD_NUM; i++)
			if((write_kmermerge[i] = (uint64_t* )calloc(bufferSize, kmer_num_tmp)) == NULL)
				fprintf(stderr, "Fail to allocate memory for write_kmermerge\n"), exit(1);

		write_kmermerge_t = (uint8_t** )calloc(THREAD_NUM, sizeof(uint8_t* ));
		for(i=0; i<THREAD_NUM; i++)
			if((write_kmermerge_t[i] = (uint8_t* )calloc(bufferSize, 1)) == NULL)
				fprintf(stderr, "Fail to allocate memory for write_kmermerge_t\n"), exit(1);
#ifdef	GFA_COM
		edges_t = (uint8_t** )calloc(THREAD_NUM, sizeof(uint8_t* ));
		for(i=0; i<THREAD_NUM; i++)
			if((edges_t[i] = (uint8_t* )calloc(bufferSize, 1)) == NULL)
				fprintf(stderr, "Fail to allocate memory for edges_t\n"), exit(1);
#endif

		FILE** fpK2f = calloc(block_numf + 1, sizeof(FILE* ));
		uint64_t** bufK2f = calloc(block_numf + 1, sizeof(uint64_t* ));
		uint64_t* readNumf = calloc(block_numf + 1, sizeof(uint64_t ));
		uint64_t* pBuff = calloc(block_numf + 1, sizeof(uint64_t ));
		int* heapf = calloc(block_numf + 1, sizeof(int ));

		strcpy(path_tmp, tmp_route);
		strcat(path_tmp, "/kmerfilterin");
		fpK2f[0]=fopen(path_tmp,"rb");
		if(fpK2f[0]==NULL)
		{
			fprintf(stderr,"Fail to open kmerfilterin\n");
			exit(1);
		}

		if(block_numf)
		{
			strcpy(path_tmp, "ulimit -n 655360");
			system(path_tmp);

			bufK2f[0]=(uint64_t* )calloc(bufferSize+1, kmer_num_filtermerge);
			readNumf[0]=fread(bufK2f[0], kmer_num_filtermerge, bufferSize, fpK2f[0]);

			pBuff[0]=0;

			for(num=1; num<block_numf + 1; num++)
			{
				sprintf(cNum,"%u",num - 1);
				strcpy(nameK2, tmp_route);
				strcat(nameK2, "/kmerfilter.");
				strcat(nameK2, cNum);

				bufK2f[num]=(uint64_t* )calloc(bufferSize+1, kmer_num_filtermerge);

				fpK2f[num]=fopen(nameK2,"rb");
				if(fpK2f[num]==NULL) fprintf(stderr,"Fail to open kmerfilter.%u\n",num),exit(1);

				readNumf[num]=fread(bufK2f[num], kmer_num_filtermerge, bufferSize, fpK2f[num]);

				pBuff[num]=0;
			}

			for(num=0; num<block_numf + 1; num++)
			{
				heapTail=num;
				for(i=heapTail; i>0; i=(i-1)>>1)
				{
					if(cmpMK2_s_filter(bufK2f[num], bufK2f[heapf[(i-1)>>1]]))
					{
						heapf[i]=heapf[(i-1)>>1];
					}
					else break;
				}
				heapf[i]=num;
			}
		}

		bufMK2 = (uint64_t* )calloc(bufferSizeM + 128, kmer_num_filtermerge);

		bound_off = 0;
		while(heapTail>=0)
		{
			pBufM=0;

			if(block_numf == 0)
			{
				pBufM = fread(bufMK2 + bound_off * kmer_num_fm, kmer_num_filtermerge, bufferSizeM, fpK2f[0]);
				if(pBufM == 0)	break;

				fprintf(stderr, "block_numf==0 %"PRId64"\n", pBufM);
			}
			else
			{
				while(pBufM<bufferSizeM&&heapTail>=0)   //get bufMK2
				{
					num=heapf[0];

					memcpy(bufMK2 + (pBufM + bound_off) * kmer_num_fm, bufK2f[num] + pBuff[num] * kmer_num_fm, kmer_num_filtermerge);

					pBufM++;
					++pBuff[num];

					if(pBuff[num]==readNumf[num])
					{
						readNumf[num] = fread(bufK2f[num], kmer_num_filtermerge, bufferSize, fpK2f[num]);

						pBuff[num]=0;

						if(readNumf[num]==0)
						{
							fclose(fpK2f[num]);
							free(bufK2f[num]);

							num=heapf[heapTail];
							heapTail--;
						}
					}
					for(i=0; (i<<1)+1<=heapTail; i=temp)
					{
						cdi1 = (i<<1)+2;
						cdi2 = cdi1-1;

						temp = (((cdi1<=heapTail) && cmpMK2_s_filter(bufK2f[heapf[cdi1]] + pBuff[heapf[cdi1]] * kmer_num_fm, bufK2f[heapf[cdi2]] + pBuff[heapf[cdi2]] * kmer_num_fm)) ? cdi1:cdi2);

						if(cmpMK2_s_filter(bufK2f[heapf[temp]] + pBuff[heapf[temp]] * kmer_num_fm, bufK2f[num] + pBuff[num] * kmer_num_fm))
						{
							heapf[i]=heapf[temp];
						}
						else break;
					}
					heapf[i]=num;

					if(pBufM>=bufferSizeM-15)
					{
						if(OneUnit)
						{
							if(((bufMK2[(pBufM-1 + bound_off) * kmer_num_fm]<<2)&eliminate)!=((bufK2f[heapf[0]][pBuff[heapf[0]] * kmer_num_fm]<<2)&eliminate))
								break;
						}
						else
						{
							getkmerno_psf(tmp_comf, bufMK2, (pBufM-1 + bound_off) * kmer_num_fm)
							getkmerno_psf(tmp_coms, bufK2f[heapf[0]], pBuff[heapf[0]] * kmer_num_fm)
							if(memcmp(tmp_comf, tmp_coms, kmer_num_tmp))	break;
						}
					}
				}
			}

			//use bufMK2 to generate tags

			if(pBufM == 0)	break;

#ifdef	BOUNDARY
			//update pBufM
			pBufM += bound_off;
			//
#endif

			segment=pBufM/THREAD_NUM;//4194289

			splitMK2[0]=0; //lower bound (can reach)
			for(i=1; i<THREAD_NUM; i++)
			{
				splitMK2[i]=i*segment;

				if(OneUnit)
				{
					tempK = (bufMK2[splitMK2[i] * kmer_num_fm]<<2)&eliminate_mi;

					while(((bufMK2[(splitMK2[i]+1) * kmer_num_fm]<<2)&eliminate_mi)==tempK)
					{
						splitMK2[i]++;
					}
				}
				else
				{
					getkmerno_pssf(tmp_comf, bufMK2, splitMK2[i] * kmer_num_fm)

					while(1)
					{
						getkmerno_pssf(tmp_coms, bufMK2, (splitMK2[i] + 1) * kmer_num_fm)
						if(memcmp(tmp_comf, tmp_coms, kmer_num_tmp))	break;
						else	splitMK2[i]++;
					}
				}


				splitMK2[i]++;

			}
			splitMK2[THREAD_NUM]=pBufM;

			myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
			thread_data* tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));
			for(i=0; i<THREAD_NUM; i++)
			{
				w_m_p[i] = 0;

				tt[i].thread = THREAD_NUM;
				tt[i].num = i;
				check=pthread_create(&myThread[i], NULL, multiTagcheck_filter, tt + i);

				if(check)
				{
					fprintf(stderr, "Thread Num: %d, Error: pthread create return code: %d\n", i, check);
					exit(EXIT_FAILURE);
				}
			}

			for(i=0; i<THREAD_NUM; i++)
			{
				pthread_join(myThread[i], NULL);
			}
			free(myThread);
			free(tt);

			for(i=0; i<THREAD_NUM; i++)
				fwrite(write_kmermerge[i], kmer_num_tmp, w_m_p[i], fp_tomerge);

			for(i=0; i<THREAD_NUM; i++)
				fwrite(write_kmermerge_t[i], 1, w_m_p[i], fp_tomerge_t);
#ifdef	GFA_COM
			for(i=0; i<THREAD_NUM; i++)
				fwrite(edges_t[i], 1, w_m_p[i], fp_edges_t);
#endif

#ifdef	BOUNDARY
			//get num of last k-1 mers and copy these kmers to new bufMK2
			if(heapTail >= 0)
			{
				bound_off = pBufM - last_block_num;
				fprintf(stderr, "BOUNDARY: %d\n", bound_off);
				memcpy(bufMK2, bufMK2 + (pBufM - bound_off) * kmer_num_fm, bound_off * kmer_num_filtermerge);
			}
#endif
			pBufM_pre = pBufM;
			//


			fprintf(stderr, "Merging block %d kmer-count filter and multiTagcheck finish\n", block_filter);

			block_filter++;
		}
#ifdef	GFA_COM
		multiTagcheck_filter_addkmer(pBufM_pre, last_block_num, fp_tomerge, fp_tomerge_t, fp_edges_t);
#else
		multiTagcheck_filter_addkmer(pBufM_pre, last_block_num, fp_tomerge, fp_tomerge_t);
#endif
		if(fp_tomerge) fclose(fp_tomerge);
		if(fp_tomerge_t) fclose(fp_tomerge_t);

#ifdef	GFA_COM
		if(fp_edges_t) fclose(fp_edges_t);
#endif
		free(bufMK2);

		for(i=0; i<THREAD_NUM; i++)
			if(write_kmermerge[i])	free(write_kmermerge[i]);
		if(write_kmermerge)	free(write_kmermerge);

		for(i=0; i<THREAD_NUM; i++)
			if(write_kmermerge_t[i])	free(write_kmermerge_t[i]);
		if(write_kmermerge_t)	free(write_kmermerge_t);

#ifdef	GFA_COM
		for(i=0; i<THREAD_NUM; i++)
			if(edges_t[i])	free(edges_t[i]);
		if(edges_t)	free(edges_t);
#endif

		for(i=0; i<block_numf; i++)
		{
			strcpy(path_tmp, "rm -f ");
			strcat(path_tmp, tmp_route);
			strcat(path_tmp, "/kmerfilter.");
			sprintf(cNum,"%u",i);
			strcat(path_tmp, cNum);
			system(path_tmp);
		}

		strcpy(path_tmp, "rm -f ");
		strcat(path_tmp, tmp_route);
		strcat(path_tmp, "/kmerfilterin");
		system(path_tmp);

		gettimeofday(&t_end1, NULL);
		end1 = ((long)t_end1.tv_sec)*1000+(long)t_end1.tv_usec/1000;

		cost_time = end1 - start1;
		fprintf(stderr, "Filter Merging time is %.2lf s\n", cost_time/1000);

	}
	else
	{
		if(fp_tomerge)	fclose(fp_tomerge);
		if(fp_tomerge_t)	fclose(fp_tomerge_t);
#ifdef	GFA_COM
		if(fp_edges_t)	fclose(fp_edges_t);
#endif
		for(i=0; i<THREAD_NUM; i++)
			if(write_kmermerge[i])	free(write_kmermerge[i]);
		if(write_kmermerge)	free(write_kmermerge);

#ifdef	GFA_COM
		for(i=0; i<THREAD_NUM; i++)
			if(edges_t[i])	free(edges_t[i]);
		if(edges_t)	free(edges_t);
#endif

	}

	if(segCountf) free(segCountf);

#endif


	fprintf(stderr, "Begin second 4-way merging kmer file and do multi-Merge-Tagcheck\n");

	gettimeofday(&t_start1, NULL);
	start1 = ((long)t_start1.tv_sec)*1000+(long)t_start1.tv_usec/1000;

	//every kmerMerge. is in order and need to cat all kmerMerge. to find the first char A/C/G/T 's position and do 4-way merge

	bufferSize = (uint64_t )((double)((uint64_t)memoryKmer<<20) / (double)((kmer_num << 5) * (THREAD_NUM + 1) + THREAD_NUM));

	fprintf(stderr, "MergeTag: %"PRId64" %"PRId64"\n", ((uint64_t)memoryKmer<<20), (((kmer_num << 5) + 1) * (THREAD_NUM + 1) + 3));

	write_fy = (uint64_t** )calloc(THREAD_NUM, sizeof(uint64_t* ));
	for(i=0; i<THREAD_NUM; i++)
		if((write_fy[i] = (uint64_t* )calloc(bufferSize, kmer_num_tmp)) == NULL)
			fprintf(stderr, "Fail to allocate memory for write_fy\n"), exit(1);

	write_ry = (uint64_t** )calloc(THREAD_NUM, sizeof(uint64_t* ));
	for(i=0; i<THREAD_NUM; i++)
		if((write_ry[i] = (uint64_t* )calloc(bufferSize, kmer_num_tmp)) == NULL)
			fprintf(stderr, "Fail to allocate memory for write_ry\n"), exit(1);

	write_x = (uint64_t** )calloc(THREAD_NUM, sizeof(uint64_t* ));
	for(i=0; i<THREAD_NUM; i++)
		if((write_x[i] = (uint64_t* )calloc(bufferSize, kmer_num_tmp)) == NULL)
			fprintf(stderr, "Fail to allocate memory for write_x\n"), exit(1);


#ifdef	FILTERFIRSTMERGE_STOP
	memset(alphabet_num, 0, 40);//32

	for(i = 0; i < 4; i++)
		for(num = 0; num < THREAD_NUM; num++)
			alphabet_num[i] += alphabet_num_new[num][i];

	for(i = 4; i > 0; i--)
		alphabet_num[i] = alphabet_num[i - 1];
	alphabet_num[0] = 0;

#endif

	memcpy(alphabet_num_ori, alphabet_num, 40);
	if(!(alphabet_num_ori[1] && alphabet_num_ori[2] && alphabet_num_ori[3] && alphabet_num_ori[4]))
		fprintf(stderr, "err merging\n"), exit(1);

	for(i = 1; i < 5; i++)
		alphabet_num[i] += alphabet_num[i - 1];


	if(alphabet_num[4] < MIN_BLOCK)
	{
		fprintf(stderr, "process 4-way merging for small dataset\n");

		uint64_t readNum_m_s = 0;

		FILE* fp_kmermergewhole_s = fopen(path_whole_tmp, "rb");
		FILE* fp_kmermergewhole_t_s = fopen(path_whole_tmp_t, "rb");

		bufMK2_m = (uint64_t* )calloc(bufferSize + 64, kmer_num_tmp);
		bufMK2_t = (uint8_t* )calloc(bufferSize + 64, 1);

		fread(bufMK2_m, kmer_num_tmp, bufferSize, fp_kmermergewhole_s);
		readNum_m_s = fread(bufMK2_t, 1, bufferSize, fp_kmermergewhole_t_s);

		BubbleSort(bufMK2_m, bufMK2_t, readNum_m_s, kmer_num_tmp);

		splitMK2[0] = 0;
		splitMK2[1] = readNum_m_s;

		myThread = (pthread_t* )calloc(1, sizeof(pthread_t ));
		thread_data* tt = (thread_data* )calloc(1, sizeof(thread_data ));

		for(i=0; i<1; i++)
		{
			tt[i].thread = 1;
			tt[i].num = i;
			w_f_p[i] = 0;
			w_r_p[i] = 0;
			w_x_p[i] = 0;

			check=pthread_create(&myThread[i], NULL, multiMergeTagcheck, tt + i);

			if(check)
			{
				fprintf(stderr, "Thread Num: %d, Error: pthread create return code: %d\n",i,check);
				exit(EXIT_FAILURE);
			}
			//pthread_join(myThread[i],NULL);
		}

		for(i=0; i<1; i++)
		{
			pthread_join(myThread[i], NULL);
		}
		free(myThread);
		free(tt);

#ifdef	BOUNDARY
		multiMergeTagcheck_addkmer(readNum_m_s, last_block_num, fp_fy, fp_ry, fp_x);
#endif

		fwrite(write_fy[0], kmer_num_tmp, w_f_p[0], fp_fy);

		fwrite(write_ry[0], kmer_num_tmp, w_r_p[0], fp_ry);

		fwrite(write_x[0], kmer_num_tmp, w_x_p[0], fp_x);

		fy_n += w_f_p[0];
		ry_n += w_r_p[0];
		x_n += w_x_p[0];

		if(bufMK2_m)	free(bufMK2_m);
		if(bufMK2_t)	free(bufMK2_t);

		if(fp_kmermergewhole_s)	fclose(fp_kmermergewhole_s);
		if(fp_kmermergewhole_t_s)	fclose(fp_kmermergewhole_t_s);

	}
	else
	{
		long int off_tmp = 0;
		for(i = 0; i < 4; i++)
		{
			fp_kmermergewhole[i] = fopen(path_whole_tmp, "rb");
			if(fp_kmermergewhole[i] == NULL)	fprintf(stderr,"Fail to open file kmerMergewhole %d\n",i),exit(1);
			//fseek(fp_kmermergewhole[i], alphabet_num[i] * 8, SEEK_SET);
			off_tmp = alphabet_num[i] * ((uint64_t )kmer_num_tmp);
			fseek(fp_kmermergewhole[i], off_tmp, SEEK_SET);

			fp_kmermergewhole_t[i] = fopen(path_whole_tmp_t, "rb");
			if(fp_kmermergewhole_t[i] == NULL)	fprintf(stderr,"Fail to open file kmerMergewhole_t %d\n",i),exit(1);
			off_tmp = alphabet_num[i];
			fseek(fp_kmermergewhole_t[i], off_tmp, SEEK_SET);
			//fseek(fp_kmermergewhole_t[i], alphabet_num[i] * (kmer_num << 3), SEEK_SET);
		}

		//start merging the kmers to merge

		memset(heap, 0, 16);//hashCount

		uint64_t readcnt[4];
		uint64_t readtol[4];
		uint64_t readspa[4];
		uint64_t pBuf_m[4];
		uint64_t* bufK2_m[4];
		uint8_t* bufK2_t[4];
		uint64_t readNum_m[4];

		for(i=1; i<5; i++)
		{
			if(alphabet_num_ori[i] < bufferSize)	bufferSize = alphabet_num_ori[i];
		}

		bufferSizeM = bufferSize * THREAD_NUM;

		for(i=0; i<4; i++)
		{
			bufK2_m[i] = (uint64_t* )calloc(bufferSize+1, kmer_num_tmp);
			readtol[i] = (alphabet_num_ori[i+1] / bufferSize) + 1;
			readspa[i] = alphabet_num_ori[i+1] - ((readtol[i] - 1) * bufferSize);
			if(!readspa[i])	readtol[i]--;
		}

		bufMK2_m=(uint64_t* )calloc(bufferSizeM + 64, kmer_num_tmp);
		bufMK2_t=(uint8_t* )calloc(bufferSizeM + 64, 1);

		for(num=0; num<4; num++)
		{
			fread(bufK2_m[num], kmer_num_tmp, bufferSize, fp_kmermergewhole[num]);

			pBuf_m[num]=0;

			bufK2_t[num]=(uint8_t* )calloc(bufferSize+1, 1);
			readNum_m[num]=fread(bufK2_t[num], 1, bufferSize, fp_kmermergewhole_t[num]);
			readcnt[num] = 1;
		}

		for(num=0; num<4; num++)
		{
			heapTail=num;
			for(i=heapTail; i>0; i=(i-1)>>1)
			{
				if(cmpMK2_sm(bufK2_m[num], bufK2_m[heap[(i-1)>>1]]))
				{
					heap[i]=heap[(i-1)>>1];
				}
				else break;
			}

			heap[i]=num;
		}

		bound_off = 0;
		block_num = 0;

		while(heapTail>=0)
		{
			pBufM=0;

			while(pBufM<bufferSizeM&&heapTail>=0)   //get bufMK2
			{
				num=heap[0];

				memcpy(bufMK2_m + (pBufM + bound_off) * kmer_num, bufK2_m[num] + pBuf_m[num] * kmer_num, kmer_num_tmp);

				bufMK2_t[pBufM + bound_off]=bufK2_t[num][pBuf_m[num]];

				pBufM++;
				++pBuf_m[num];

				//deal with the heap
				if(pBuf_m[num]==readNum_m[num])
				{
					fread(bufK2_m[num], kmer_num_tmp, bufferSize,fp_kmermergewhole[num]);
					readNum_m[num] = fread(bufK2_t[num], 1, bufferSize, fp_kmermergewhole_t[num]);

					pBuf_m[num]=0;

					++readcnt[num];

					if(readcnt[num] == readtol[num])
					{
						readNum_m[num] = readspa[num];
					}
					else if(readcnt[num] == readtol[num] + 1)
					{

						num=heap[heapTail];
						heapTail--;
					}
				}

				for(i=0; (i<<1)+1<=heapTail; i=temp)
				{
					cdi1=(i<<1)+2;
					cdi2=cdi1-1;

					temp=(cdi1<=heapTail&&cmpMK2_sm(bufK2_m[heap[cdi1]] + pBuf_m[heap[cdi1]] * kmer_num, bufK2_m[heap[cdi2]] + pBuf_m[heap[cdi2]] * kmer_num))?cdi1:cdi2;

					heap_tmp = heap[temp];
					if(heap_tmp == num)	heap[i]=num;
					else
					{
						if(cmpMK2_sm(bufK2_m[heap_tmp] + pBuf_m[heap_tmp] * kmer_num, bufK2_m[num] + pBuf_m[num] * kmer_num))
						{
							heap[i]=heap_tmp;
						}
						else break;
					}
				}
				heap[i]=num;

				if(pBufM>=bufferSizeM-15)
				{
					if(kmer_num == 1)
					{
						if(((bufMK2_m[pBufM-1 + bound_off]<<2)&eliminate)!=((bufK2_m[heap[0]][pBuf_m[heap[0]]]<<2)&eliminate))
						{
							break;
						}
					}
					else
					{
						getkmerno_psm(tmp_comf, bufMK2_m, (pBufM-1 + bound_off) * kmer_num)
						getkmerno_psm(tmp_coms, bufK2_m[heap[0]], pBuf_m[heap[0]] * kmer_num)
						if(memcmp(tmp_comf, tmp_coms, kmer_num_tmp))	break;
					}
				}
			}

			//use bufMK2 to generate tags
#ifdef	BOUNDARY
			pBufM += bound_off;
#endif
			segment=pBufM/THREAD_NUM;

			splitMK2[0]=0; //lower bound (can reach)
			for(i=1; i<THREAD_NUM; i++)
			{
				splitMK2[i]=i*segment;

#ifndef	SPLIT_MI
				tempK = (bufMK2_m[splitMK2[i]]<<2)&eliminate;

				while(((bufMK2_m[splitMK2[i]+1]<<2)&eliminate)==tempK)
				{
					splitMK2[i]++;
				}
#else
				//
				//if(OneUnit)

				if(kmer_num == 1)
				{
					tempK = (bufMK2_m[splitMK2[i]]<<2);//&eliminate_mi

					while((bufMK2_m[splitMK2[i]+1]<<2)==tempK)   //&eliminate_mi
					{
						splitMK2[i]++;
					}
				}
				else
				{
					//getkmerno_pssm(tmp_comf, bufMK2_m, splitMK2[i] * kmer_num)
					getkmerno_psm(tmp_comf, bufMK2_m, splitMK2[i] * kmer_num)

					while(1)
					{
						//getkmerno_pssm(tmp_coms, bufMK2_m, (splitMK2[i] + 1) * kmer_num)
						getkmerno_psm(tmp_coms, bufMK2_m, (splitMK2[i] + 1) * kmer_num)

						if(memcmp(tmp_comf, tmp_coms, kmer_num_tmp))	break;
						else	splitMK2[i]++;
					}

				}


#endif
				splitMK2[i]++;

			}
			splitMK2[THREAD_NUM]=pBufM;

			myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
			thread_data* tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));

			for(i=0; i<THREAD_NUM; i++)
			{
				tt[i].thread = THREAD_NUM;
				tt[i].num = i;
				w_f_p[i] = 0;
				w_r_p[i] = 0;
				w_x_p[i] = 0;

				check=pthread_create(&myThread[i], NULL, multiMergeTagcheck, tt + i);

				if(check)
				{
					fprintf(stderr, "Thread Num: %d, Error: pthread create return code: %d\n",i,check);
					exit(EXIT_FAILURE);
				}
				//pthread_join(myThread[i],NULL);
			}

			for(i=0; i<THREAD_NUM; i++)
			{
				pthread_join(myThread[i], NULL);
			}
			free(myThread);
			free(tt);

			for(i=0; i<THREAD_NUM; i++)
				fwrite(write_fy[i], kmer_num_tmp, w_f_p[i], fp_fy);

			for(i=0; i<THREAD_NUM; i++)
				fwrite(write_ry[i], kmer_num_tmp, w_r_p[i], fp_ry);

			for(i=0; i<THREAD_NUM; i++)
				fwrite(write_x[i], kmer_num_tmp, w_x_p[i], fp_x);

			for(i=0; i<THREAD_NUM; i++)
				fy_n += w_f_p[i];
			for(i=0; i<THREAD_NUM; i++)
				ry_n += w_r_p[i];
			for(i=0; i<THREAD_NUM; i++)
				x_n += w_x_p[i];

			fprintf(stderr, "Block %d merging finish\n", block_num + 1);

#ifdef	BOUNDARY
			//get num of last k-1 mers and copy these kmers to new bufMK2
			if(heapTail >= 0)
			{
				bound_off = pBufM - last_block_num;
				//fprintf(stderr, "BOUNDARY merge: %d %d %d %d\n", bound_off, heapTail, pBufM, last_block_num);//-64 2 8 72
				memcpy(bufMK2_m, bufMK2_m + (pBufM - bound_off) * kmer_num, bound_off * kmer_num_tmp);
				memcpy(bufMK2_t, bufMK2_t + (pBufM - bound_off), bound_off);
			}
#endif

			pBufM_pre = pBufM;
			block_num++;
		}
#ifdef	BOUNDARY
		multiMergeTagcheck_addkmer(pBufM_pre, last_block_num, fp_fy, fp_ry, fp_x);
#endif


		for(i=0; i<4; i++)
			if(bufK2_m[i])	free(bufK2_m[i]);

		for(i=0; i<4; i++)
			if(bufK2_t[i])	free(bufK2_t[i]);

	}

	for(i=0; i<THREAD_NUM; i++)
		if(write_fy[i])	free(write_fy[i]);
	if(write_fy)	free(write_fy);

	for(i=0; i<THREAD_NUM; i++)
		if(write_ry[i])	free(write_ry[i]);
	if(write_ry)	free(write_ry);

	for(i=0; i<THREAD_NUM; i++)
		if(write_x[i])	free(write_x[i]);
	if(write_x)	free(write_x);

	fprintf(stderr, "\n\nstart-node end-node start-end-bode num: %"PRId64" %"PRId64" %"PRId64"\n\n", fy_n, ry_n, x_n);

	gettimeofday(&t_end1, NULL);
	end1 = ((long)t_end1.tv_sec)*1000+(long)t_end1.tv_usec/1000;

	cost_time = end1 - start1;
	fprintf(stderr, "4-way merging and multi-Merge-Tagcheck time is %.2lf s\n", cost_time/1000);

//

	FILE *fp_size = NULL;
	uint64_t filesize = 0;

	fp_size = fopen(path_whole_tmp_t, "r");
	fseek(fp_size, 0L, SEEK_END);
	filesize = ftell(fp_size);
	if(fp_size)	fclose(fp_size);

	fprintf(stderr, "total number of all k-mers: %"PRId64"\n", filesize);

	strcpy(path_tmp, tmp_route);
	strcat(path_tmp, "/kmerFY");
	fp_size = fopen(path_tmp, "r");
	fseek(fp_size, 0L, SEEK_END);
	filesize = ftell(fp_size);

	if(fp_size)	fclose(fp_size);

	strcpy(path_tmp, tmp_route);
	strcat(path_tmp, "/kmerRY");
	fp_size = fopen(path_tmp, "r");
	fseek(fp_size, 0L, SEEK_END);
	filesize = ftell(fp_size);
	if(fp_size)	fclose(fp_size);

	strcpy(path_tmp, tmp_route);
	strcat(path_tmp, "/kmerX");
	fp_size = fopen(path_tmp, "r");
	fseek(fp_size, 0L, SEEK_END);
	filesize = ftell(fp_size);
	if(fp_size)	fclose(fp_size);

	if(fp_fy)	fclose(fp_fy);
	if(fp_ry)	fclose(fp_ry);
	if(fp_x)	fclose(fp_x);

#endif

	getbranchkmer(THREAD_NUM, tmp_route);

#ifdef	BEF_BRANCH

	strcpy(path_tmp, "rm -f ");
	strcat(path_tmp, tmp_route);
	strcat(path_tmp, "/kmerMergewhole");
	system(path_tmp);

	strcpy(path_tmp, "rm -f ");
	strcat(path_tmp, tmp_route);
	strcat(path_tmp, "/kmerMergewhole_t");
	system(path_tmp);

	strcpy(path_tmp, "rm -f ");
	strcat(path_tmp, tmp_route);
	strcat(path_tmp, "/edges_t");
	system(path_tmp);
	
#endif

	return 1;
}
void *multiCount_filter(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	uint64_t num = argT->num;
	countKmer_disf[num] = (uint64_t* )calloc(BUCKET_CAPACITY+1,sizeof(uint64_t ));

	uint16_t kmer_num_count = kmer_num_filtermerge;
	uint64_t step = 0;
	uint64_t i = 0;

	uint64_t* read_p = NULL;

	if(w_f_os[num])
	{
		read_p = write_filterout[num];

		if(OneUnit)	step = 2;
		else	step = kmer_num_new;

		for(i = 0; i < w_f_os[num]; i++)
			countKmer_disf[num][(read_p[i * step] >> 46)&maskBKT]++;

		for(i=1; i<BUCKET_CAPACITY; i++)
			countKmer_disf[num][i]=countKmer_disf[num][i-1]+countKmer_disf[num][i];

		countKmer_disf[num][BUCKET_CAPACITY] = countKmer_disf[num][BUCKET_CAPACITY-1];

		hashKmer_disf[num] = (uint64_t* )calloc(countKmer_disf[num][BUCKET_CAPACITY], kmer_num_count);// * kmer_num, sizeof(uint64_t)
	}

	return (void*)NULL;
}

void *multiDistri_filter(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	uint64_t num=argT->num;
	uint64_t* read_p = NULL;
	uint64_t step = 0;
	uint64_t i = 0;
	uint64_t tempSeq = 0;
	uint16_t kmer_num_distri = kmer_num_filtermerge;

	uint64_t step_tmp = 0;
	if(w_f_os[num])
	{
		read_p = write_filterout[num];

		if(OneUnit)	step = 2;
		else	step = kmer_num_new;

		for(i = 0; i < w_f_os[num]; i++)
		{
			step_tmp = i * step;
			tempSeq = (read_p[step_tmp] >> 46)&maskBKT;
			countKmer_disf[num][tempSeq]--;
			memcpy(hashKmer_disf[num] + countKmer_disf[num][tempSeq] * step, read_p + step_tmp, kmer_num_distri);
		}
	}

	return (void *)NULL;
}

void *multiThreadSort_filter(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	uint64_t num=argT->num;
	if(countKmer_disf[num])	free(countKmer_disf[num]);


	uint64_t i,low,up;
	uint64_t kmer_copy_tmp = 0;

	kmer_copy_tmp = kmer_num_fm;
	uint64_t kmer_num_qsort = kmer_num_filtermerge;

	low=segCountf[num];
	up=segCountf[num+1];

	if(hashKmer_disf[num])	free(hashKmer_disf[num]);

	for(i=low; i<up; i++)
	{
		qsort(hashKmerf + countKmerf[i] * kmer_copy_tmp, countKmerf[i+1] - countKmerf[i], kmer_num_qsort, cmpKmer_filter);
	}

	return (void *)NULL;
}

int cmpKmer_filter(const void *a, const void *b)
{
	if(OneUnit)
	{
		uint64_t* vav = (uint64_t* )a;
		uint64_t* vbv = (uint64_t* )b;
		uint64_t va = *vav;
		uint64_t vb = *vbv;

		va <<= 2;
		vb <<= 2;

		if(va<vb) return -1;
		else if(va==vb) return 0;
		else return 1;
	}
	else
	{
		uint8_t i = 0;
		uint64_t va_t = 0, vb_t = 0;

		uint64_t* va = (uint64_t* )a;
		uint64_t* vb = (uint64_t* )b;

		va_t = (*va);
		vb_t = (*vb);

		va_t <<= 2;
		vb_t <<= 2;

		if(va_t < vb_t) return -1;
		else if(va_t > vb_t) return 1;
		else
		{
			if(kmer_num > 1)
			{
				for(i = 1; i< kmer_num - 1; i++)
				{
					if(va[i] < vb[i])	return -1;
					else if(va[i] > vb[i])	return 1;
				}

				if(kmer_num == kmer_num_new)
				{
					va_t = (va[kmer_num - 1] >> 8) << 8;
					vb_t = (vb[kmer_num - 1] >> 8) << 8;
					if(va_t < vb_t)	return -1;
					else if(va_t > vb_t)	return 1;
					else	return 0;
				}
				else
				{
					if(va[i] < vb[i])	return -1;
					else if(va[i] > vb[i])	return 1;
					else	return 0;
				}
			}
			else	return 0;
		}
	}
}

void *multiTagcheck_filter(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	int num = argT->num;
	int totalNum = argT->thread;

	uint64_t i = 0;

	uint64_t start=splitMK2[num], end=splitMK2[num+1];
	uint64_t bufMK1_len=end-start;

	if(bufMK1_len==0) return (void *)NULL;

	uint64_t eliminate=~(ELIMINATE);
	uint64_t eliminate_m = (eliminate << 2);
	uint64_t anchor = 0;
	uint64_t anchor_a = 0;
	uint64_t kmer_pre = 0;
	uint64_t kmer_pre_a = 0;

	uint8_t input = 0;
	uint8_t output = 0;
	uint8_t in_j = 0;
	uint8_t out_j = 0;
	uint8_t in_linear = 0;
	uint8_t out_linear = 0;
	uint8_t one_tmp = 1;
	uint8_t type_f = 0;
	uint8_t type_s = 0;
	uint8_t alph_tmp = 0;
	uint8_t alph_pre = 0;
	uint64_t char_tmp = 0;
	uint64_t w_m_p_t = 0;
	uint64_t anchor_km = 0;
	uint64_t anchor_km_pre = 0;
	uint64_t anchor_km_pre_i = 0;
	uint64_t j = 0;
	uint8_t kmer_tmp = 0;
	uint8_t pre_k = 0;
	uint64_t pre_i = 0;
	uint64_t pre_j = 0;
	uint8_t pre_flag = 0;
	uint8_t kmer_num_tmp = (kmer_num << 3);

	anchor_km_pre_i = start;

	uint64_t start_get = 0;
	uint8_t kmer_ins[24];
	uint8_t kmer_outs[24];
	uint64_t kmer_indexs[24];
	uint8_t kmer_i = 0;
	uint8_t kmer_k = 0;
	uint8_t aph_i = 0;
	uint8_t type_tmp = 0;
	uint8_t de_in = 0;
	uint8_t de_out = 0;
	uint8_t flag_filter_in = 0;

	if(OneUnit)
	{
		uint8_t kmer_ins[24];
		uint8_t kmer_outs[24];
		uint64_t kmer_indexs[24];
		uint8_t kmer_i = 0;
		uint8_t kmer_k = 0;
		uint64_t start_get_kmer = 0;
		uint8_t type_c = 0;

		for(i = start; i < end; i++)
		{
			start_get = i << 1;

			anchor_km = ((bufMK2[start_get] << 2)&eliminate_m);

			if((anchor_km != anchor_km_pre) && (i > start))
			{
				if(i - anchor_km_pre_i > 1)
				{
					pre_j = anchor_km_pre_i << 1;
					for(j = anchor_km_pre_i; j < i; j++)
					{
						start_get_kmer = j << 1;
						anchor_a = bufMK2[start_get_kmer];

						anchor = ((anchor_a<<2)&eliminate);

						if((anchor != kmer_pre) && (j > anchor_km_pre_i))
						{
							if(flag_filter_in)
							{
								if(in_j || out_j)
								{
									in_j ^= de_in;
									out_j ^= de_out;
								}
								kmer_ins[kmer_i] = in_j;
								kmer_outs[kmer_i] = out_j;
								kmer_indexs[kmer_i] = pre_j;

								kmer_i++;
							}

							in_j = 0;
							out_j = 0;
							de_in = 0;
							de_out = 0;
							flag_filter_in = 0;
							pre_j = start_get_kmer;
						}

						type_tmp = (anchor_a >> 62) & 0X3;

						if(type_tmp == 1)
						{
							de_out |= (one_tmp << bufMK2[start_get_kmer + 1]);
						}
						else if(type_tmp == 2)
						{
							de_in |= (one_tmp << bufMK2[start_get_kmer + 1]);
						}
						else
						{
							//debug
							if(flag_filter_in == 1)	printf("Err: %u %"PRId64"", type_tmp, anchor_a);

							in_j = (bufMK2[start_get_kmer + 1] >> 4) & 0Xf;
							out_j = bufMK2[start_get_kmer + 1] & 0Xf;

							flag_filter_in = 1;
						}

						kmer_pre = anchor;
					}

					if(flag_filter_in)
					{
						if(in_j || out_j)
						{
							in_j ^= de_in;
							out_j ^= de_out;
						}
						kmer_ins[kmer_i] = in_j;
						kmer_outs[kmer_i] = out_j;
						kmer_indexs[kmer_i] = pre_j;

						kmer_i++;
					}


				}
				else
				{
					start_get_kmer = anchor_km_pre_i * kmer_num_filter;
					kmer_pre_a = bufMK2[start_get_kmer];

					type_tmp = (kmer_pre_a >> 62) & 0X3;

					if(!type_tmp)
					{
						in_j = (bufMK2[start_get_kmer + 1] >> 4) & 0Xf;
						out_j = bufMK2[start_get_kmer + 1] & 0Xf;

						if(in_j || out_j)
						{
							in_j ^= de_in;
							out_j ^= de_out;
						}
						kmer_ins[kmer_i] = in_j;
						kmer_outs[kmer_i] = out_j;
						kmer_indexs[kmer_i] = start_get_kmer;

						kmer_i++;

					}
				}

				for(kmer_k = 0; kmer_k < kmer_i; kmer_k++)
				{
					kmer_pre_a = bufMK2[kmer_indexs[kmer_k]];
					in_j = kmer_ins[kmer_k];
					out_j = kmer_outs[kmer_k];

					if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
					else	in_linear = 0;

					if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
					else	out_linear = 0;

					if(in_linear && out_linear)			//linear
						type_f = 0;
					else if(in_linear && !out_linear)	//single-in multi-out
						type_f = 1;
					else if(!in_linear && out_linear)	//multi-in single-out
						type_f = 2;
					else								//multi-in multi-out
						type_f = 3;

					type_s = type_f;
					if((type_f == 0) || (type_f == 1))
					{
						pre_flag = 0;

						for(pre_k = 0; pre_k < kmer_k; pre_k++)
							if(kmer_ins[pre_k] & in_j)
							{
								pre_flag = 1;
								break;
							}

						for(pre_k = kmer_k + 1; (pre_k < kmer_i) && (!pre_flag); pre_k++)
							if(kmer_ins[pre_k] & in_j)
							{
								pre_flag = 1;
								break;
							}

						if(pre_flag)
						{
							if(type_f == 0)	type_s = 2;
							else	type_s = 3;
						}
					}

					kmer_tmp = ((kmer_pre_a >> 60) & 0Xf);
					write_kmermerge[num][w_m_p_t] = ((kmer_pre_a << 2) & eliminate);

					alph_tmp = (kmer_tmp & 0X3);
					alphabet_num_new[num][alph_tmp]++;

					type_s |= (input_types[in_j] << 2);

					type_s |= (out_j << 4);

					write_kmermerge_t[num][w_m_p_t] = type_s;
#ifdef	GFA_COM
					edges_t[num][w_m_p_t] = ((in_j << 4) | out_j);
#endif
					w_m_p_t++;
				}

				in_j = 0;
				out_j = 0;
				de_in = 0;
				de_out = 0;
				flag_filter_in = 0;

				kmer_i = 0;

				anchor_km_pre_i = i;
			}

			anchor_km_pre = anchor_km;
		}

		if(num != totalNum - 1)
		{
			in_j = 0;
			out_j = 0;
			de_in = 0;
			de_out = 0;
			flag_filter_in = 0;

			if(i - anchor_km_pre_i > 1)
			{
				pre_j = anchor_km_pre_i << 1;
				for(j = anchor_km_pre_i; j < i; j++)
				{
					start_get_kmer = j << 1;
					anchor_a = bufMK2[start_get_kmer];

					anchor = ((anchor_a<<2)&eliminate);

					if((anchor != kmer_pre) && (j > anchor_km_pre_i))
					{
						if(flag_filter_in)
						{
							if(in_j || out_j)
							{
								in_j ^= de_in;
								out_j ^= de_out;
							}
							kmer_ins[kmer_i] = in_j;
							kmer_outs[kmer_i] = out_j;
							kmer_indexs[kmer_i] = pre_j;

							kmer_i++;
						}

						in_j = 0;
						out_j = 0;
						de_in = 0;
						de_out = 0;
						flag_filter_in = 0;

						pre_j = start_get_kmer;
					}

					type_tmp = (anchor_a >> 62) & 0X3;
					if(type_tmp == 1)
					{
						de_out |= (one_tmp << bufMK2[start_get_kmer + 1]);
					}
					else if(type_tmp == 2)
					{
						de_in |= (one_tmp << bufMK2[start_get_kmer + 1]);
					}
					else
					{
						in_j = (bufMK2[start_get_kmer + 1] >> 4) & 0Xf;
						out_j = bufMK2[start_get_kmer + 1] & 0Xf;
						flag_filter_in = 1;
					}

					kmer_pre = anchor;
				}
				if(flag_filter_in)
				{
					if(in_j || out_j)
					{
						in_j ^= de_in;
						out_j ^= de_out;
					}
					kmer_ins[kmer_i] = in_j;
					kmer_outs[kmer_i] = out_j;
					kmer_indexs[kmer_i] = pre_j;

					kmer_i++;
				}


			}
			else
			{
				start_get_kmer = anchor_km_pre_i * kmer_num_filter;
				kmer_pre_a = bufMK2[start_get_kmer];

				type_tmp = (kmer_pre_a >> 62) & 0X3;

				if(!type_tmp)
				{
					in_j = (bufMK2[start_get_kmer + 1] >> 4) & 0Xf;
					out_j = bufMK2[start_get_kmer + 1] & 0Xf;

					if(in_j || out_j)
					{
						in_j ^= de_in;
						out_j ^= de_out;
					}
					kmer_ins[kmer_i] = in_j;
					kmer_outs[kmer_i] = out_j;
					kmer_indexs[kmer_i] = start_get_kmer;

					kmer_i++;

				}
			}

			for(kmer_k = 0; kmer_k < kmer_i; kmer_k++)
			{
				kmer_pre_a = bufMK2[kmer_indexs[kmer_k]];
				in_j = kmer_ins[kmer_k];
				out_j = kmer_outs[kmer_k];

				if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
				else	in_linear = 0;

				if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
				else	out_linear = 0;

				if(in_linear && out_linear)			//linear
					type_f = 0;
				else if(in_linear && !out_linear)	//single-in multi-out
					type_f = 1;
				else if(!in_linear && out_linear)	//multi-in single-out
					type_f = 2;
				else								//multi-in multi-out
					type_f = 3;

				type_s = type_f;
				if((type_f == 0) || (type_f == 1))
				{
					pre_flag = 0;

					for(pre_k = 0; pre_k < kmer_k; pre_k++)
						if(kmer_ins[pre_k] & in_j)
						{
							pre_flag = 1;
							break;
						}

					for(pre_k = kmer_k + 1; (pre_k < kmer_i) && (!pre_flag); pre_k++)
						if(kmer_ins[pre_k] & in_j)
						{
							pre_flag = 1;
							break;
						}

					if(pre_flag)
					{
						if(type_f == 0)	type_s = 2;
						else	type_s = 3;
					}
				}

				kmer_tmp = ((kmer_pre_a >> 60) & 0Xf);
				write_kmermerge[num][w_m_p_t] = ((kmer_pre_a << 2) & eliminate);

				alph_tmp = (kmer_tmp & 0X3);
				alphabet_num_new[num][alph_tmp]++;

				type_s |= (input_types[in_j] << 2);

				type_s |= (out_j << 4);

				write_kmermerge_t[num][w_m_p_t] = type_s;
#ifdef	GFA_COM
				edges_t[num][w_m_p_t] = ((in_j << 4) | out_j);
#endif
				w_m_p_t++;
			}

		}
		else	last_block_num = anchor_km_pre_i;
	}
	else
	{
		uint64_t anchor_kms[kmer_num];
		uint64_t anchor_kms_pre[kmer_num];
		uint64_t kmer_pre_as[kmer_num_new];

		uint8_t kmer_num_tmp_new = (kmer_num_new << 3);
		uint8_t type_c = 0;
		int16_t tra_i = 0;
		uint8_t shc = 0;
		uint8_t shc_pre = 0;
		uint8_t shtype = 0;
		uint8_t add_edge = 0;
		uint64_t w_f_o = 0;

		for(i = start; i < end; i++)
		{
			start_get = (i * kmer_num_fm);

			getkmerno_pssf(anchor_kms, bufMK2, start_get)

			if(memcmp(anchor_kms, anchor_kms_pre, kmer_num_tmp) && (i > start))
			{
				if(i - anchor_km_pre_i > 1)
				{
					pre_j = anchor_km_pre_i * kmer_num_fm;
					for(j = anchor_km_pre_i; j < i; j++)
					{
						char_tmp = j * kmer_num_fm;

						anchor = ((bufMK2[char_tmp + kmer_num - 1] << leftmove32) & KMask);

						if((anchor != kmer_pre) && (j > anchor_km_pre_i))
						{
							if(flag_filter_in)
							{
								if(in_j || out_j)
								{
									in_j ^= de_in;
									out_j ^= de_out;
								}
								kmer_ins[kmer_i] = in_j;
								kmer_outs[kmer_i] = out_j;
								kmer_indexs[kmer_i] = pre_j;

								kmer_i++;
							}

							in_j = 0;
							out_j = 0;
							de_in = 0;
							de_out = 0;
							flag_filter_in = 0;

							pre_j = char_tmp;//j
						}

						type_tmp = (bufMK2[char_tmp] >> 62) & 0X3;

						if(type_tmp == 1)
						{
							de_out |= (one_tmp << (bufMK2[char_tmp + kmer_num_fm - 1] & 0X3));
						}
						else if(type_tmp == 2)
						{
							de_in |= (one_tmp << (bufMK2[char_tmp + kmer_num_fm - 1] & 0X3));
						}
						else
						{
							in_j = (bufMK2[char_tmp + kmer_num_fm - 1] >> 4) & 0Xf;
							out_j = bufMK2[char_tmp + kmer_num_fm - 1] & 0Xf;
							flag_filter_in = 1;
						}

						kmer_pre = anchor;
					}

					if(flag_filter_in)
					{
						if(in_j || out_j)
						{
							in_j ^= de_in;
							out_j ^= de_out;
						}
						kmer_ins[kmer_i] = in_j;
						kmer_outs[kmer_i] = out_j;
						kmer_indexs[kmer_i] = pre_j;

						kmer_i++;
					}


				}
				else
				{
					char_tmp = anchor_km_pre_i * kmer_num_fm;

					memcpy(kmer_pre_as, bufMK2 + char_tmp, kmer_num_tmp_new);

					type_tmp = (kmer_pre_as[0] >> 62) & 0X3;

					if(!type_tmp)
					{
						in_j = (kmer_pre_as[kmer_num_fm - 1] >> 4) & 0Xf;
						out_j = kmer_pre_as[kmer_num_fm - 1] & 0Xf;

						if(in_j || out_j)
						{
							in_j ^= de_in;
							out_j ^= de_out;
						}
						kmer_ins[kmer_i] = in_j;
						kmer_outs[kmer_i] = out_j;
						kmer_indexs[kmer_i] = char_tmp;

						kmer_i++;

					}
				}

				for(kmer_k = 0; kmer_k < kmer_i; kmer_k++)
				{
					memcpy(kmer_pre_as, bufMK2 + kmer_indexs[kmer_k], kmer_num_tmp);

					in_j = kmer_ins[kmer_k];
					out_j = kmer_outs[kmer_k];

					if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
					else	in_linear = 0;

					if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
					else	out_linear = 0;

					if(in_linear && out_linear)			//linear
						type_f = 0;
					else if(in_linear && !out_linear)	//single-in multi-out
						type_f = 1;
					else if(!in_linear && out_linear)	//multi-in single-out
						type_f = 2;
					else								//multi-in multi-out
						type_f = 3;

					type_s = type_f;
					if((type_f == 0) || (type_f == 1))
					{
						pre_flag = 0;

						for(pre_k = 0; pre_k < kmer_k; pre_k++)
							if(kmer_ins[pre_k] & in_j)
							{
								pre_flag = 1;
								break;
							}

						for(pre_k = kmer_k + 1; (pre_k < kmer_i) && (!pre_flag); pre_k++)
							if(kmer_ins[pre_k] & in_j)
							{
								pre_flag = 1;
								break;
							}

						if(pre_flag)
						{
							if(type_f == 0)	type_s = 2;
							else	type_s = 3;
						}
					}

					kmer_tmp = ((kmer_pre_as[0] >> 60) & 0Xf);

					kmer_pre_as[0] <<= 2;
					kmer_pre_as[kmer_num - 1] &= KMask;
					memcpy(write_kmermerge[num] + (w_m_p_t * kmer_num), kmer_pre_as, kmer_num_tmp);

					alph_tmp = (kmer_tmp & 0X3);
					alphabet_num_new[num][alph_tmp]++;

					type_s |= (input_types[in_j] << 2);

					type_s |= (out_j << 4);

					write_kmermerge_t[num][w_m_p_t] = type_s;
#ifdef	GFA_COM
					edges_t[num][w_m_p_t] = ((in_j << 4) | out_j);
#endif
					w_m_p_t++;
				}


				in_j = 0;
				out_j = 0;
				de_in = 0;
				de_out = 0;
				kmer_i = 0;
				flag_filter_in = 0;

				anchor_km_pre_i = i;
			}

			memcpy(anchor_kms_pre, anchor_kms, kmer_num_tmp);
		}

		if(num != totalNum - 1)
		{
			in_j = 0;
			out_j = 0;
			de_in = 0;
			de_out = 0;
			flag_filter_in = 0;

			if(i - anchor_km_pre_i > 1)
			{
				pre_j = anchor_km_pre_i * kmer_num_fm;
				for(j = anchor_km_pre_i; j < i; j++)
				{
					char_tmp = j * kmer_num_fm;

					anchor = ((bufMK2[char_tmp + kmer_num - 1] << leftmove32) & KMask);

					if((anchor != kmer_pre) && (j > anchor_km_pre_i))
					{
						if(flag_filter_in)
						{
							if(in_j || out_j)
							{
								in_j ^= de_in;
								out_j ^= de_out;
							}
							kmer_ins[kmer_i] = in_j;
							kmer_outs[kmer_i] = out_j;
							kmer_indexs[kmer_i] = pre_j;

							kmer_i++;
						}

						in_j = 0;
						out_j = 0;
						de_in = 0;
						de_out = 0;
						flag_filter_in = 0;

						pre_j = char_tmp;//j
					}

					type_tmp = (bufMK2[char_tmp] >> 62) & 0X3;
					if(type_tmp == 1)
					{
						de_out |= (one_tmp << (bufMK2[char_tmp + kmer_num_fm - 1] & 0X3));
					}
					else if(type_tmp == 2)
					{
						de_in |= (one_tmp << (bufMK2[char_tmp + kmer_num_fm - 1] & 0X3));
					}
					else
					{
						in_j = (bufMK2[char_tmp + kmer_num_fm - 1] >> 4) & 0Xf;
						out_j = bufMK2[char_tmp + kmer_num_fm - 1] & 0Xf;
						flag_filter_in = 1;
					}

					kmer_pre = anchor;
				}

				if(flag_filter_in)
				{
					if(in_j || out_j)
					{
						in_j ^= de_in;
						out_j ^= de_out;
					}
					kmer_ins[kmer_i] = in_j;
					kmer_outs[kmer_i] = out_j;
					kmer_indexs[kmer_i] = pre_j;

					kmer_i++;
				}


			}
			else
			{
				char_tmp = anchor_km_pre_i * kmer_num_fm;

				memcpy(kmer_pre_as, bufMK2 + char_tmp, kmer_num_tmp_new);

				type_tmp = (kmer_pre_as[0] >> 62) & 0X3;

				if(!type_tmp)
				{
					in_j = (kmer_pre_as[kmer_num_fm - 1] >> 4) & 0Xf;
					out_j = kmer_pre_as[kmer_num_fm - 1] & 0Xf;

					if(in_j || out_j)
					{
						in_j ^= de_in;
						out_j ^= de_out;
					}
					kmer_ins[kmer_i] = in_j;
					kmer_outs[kmer_i] = out_j;
					kmer_indexs[kmer_i] = char_tmp;

					kmer_i++;

				}
			}

			for(kmer_k = 0; kmer_k < kmer_i; kmer_k++)
			{
				memcpy(kmer_pre_as, bufMK2 + kmer_indexs[kmer_k], kmer_num_tmp);

				in_j = kmer_ins[kmer_k];
				out_j = kmer_outs[kmer_k];

				if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
				else	in_linear = 0;

				if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
				else	out_linear = 0;

				if(in_linear && out_linear)			//linear
					type_f = 0;
				else if(in_linear && !out_linear)	//single-in multi-out
					type_f = 1;
				else if(!in_linear && out_linear)	//multi-in single-out
					type_f = 2;
				else								//multi-in multi-out
					type_f = 3;

				type_s = type_f;
				if((type_f == 0) || (type_f == 1))
				{
					pre_flag = 0;

					for(pre_k = 0; pre_k < kmer_k; pre_k++)
						if(kmer_ins[pre_k] & in_j)
						{
							pre_flag = 1;
							break;
						}

					for(pre_k = kmer_k + 1; (pre_k < kmer_i) && (!pre_flag); pre_k++)
						if(kmer_ins[pre_k] & in_j)
						{
							pre_flag = 1;
							break;
						}

					if(pre_flag)
					{
						if(type_f == 0)	type_s = 2;
						else	type_s = 3;
					}
				}

				kmer_tmp = ((kmer_pre_as[0] >> 60) & 0Xf);

				kmer_pre_as[0] <<= 2;
				kmer_pre_as[kmer_num - 1] &= KMask;
				memcpy(write_kmermerge[num] + (w_m_p_t * kmer_num), kmer_pre_as, kmer_num_tmp);

				alph_tmp = (kmer_tmp & 0X3);
				alphabet_num_new[num][alph_tmp]++;


				//type_s |= ((kmer_tmp >> 2) << 2);
				type_s |= (input_types[in_j] << 2);

				//printf("%u\n", input_types[in_j]);

				type_s |= (out_j << 4);

				write_kmermerge_t[num][w_m_p_t] = type_s;
#ifdef	GFA_COM
				edges_t[num][w_m_p_t] = ((in_j << 4) | out_j);
#endif
				w_m_p_t++;
			}

		}
		else	last_block_num = anchor_km_pre_i;
	}

	w_m_p[num] = w_m_p_t;

	return (void *)NULL;
}

#ifdef	GFA_COM
void multiTagcheck_filter_addkmer(uint64_t i, uint64_t anchor_km_pre_i, FILE* fp_tomerge, FILE* fp_tomerge_t, FILE* fp_edges_t)
#else
void multiTagcheck_filter_addkmer(uint64_t i, uint64_t anchor_km_pre_i, FILE* fp_tomerge, FILE* fp_tomerge_t)
#endif
{
	uint8_t in_j = 0;
	uint8_t out_j = 0;
	uint8_t in_linear = 0;
	uint8_t out_linear = 0;
	uint8_t one_tmp = 1;
	uint8_t type_f = 0;
	uint8_t type_s = 0;

#ifdef	GFA_COM
	uint8_t type_e = 0;
#endif

	uint8_t input = 0;
	uint8_t output = 0;
	uint8_t pre_k = 0;
	uint8_t pre_tmp = 0;
	uint8_t pre_flag = 0;
	uint8_t kmer_tmp = 0;
	uint8_t fl_f = 0;
	uint8_t kmer_num_tmp = (kmer_num << 3);

	uint64_t j = 0;
	uint64_t pre_j = 0;
	uint64_t start_get = 0;
	uint64_t anchor = 0;
	uint64_t anchor_a = 0;
	uint64_t kmer_pre = 0;
	uint64_t kmer_pre_a = 0;
	uint64_t write_tmp = 0;
	uint8_t kmer_ins[24];
	uint8_t kmer_outs[24];
	uint64_t kmer_indexs[24];
	uint8_t kmer_i = 0;
	uint8_t kmer_k = 0;
	uint8_t de_in = 0;
	uint8_t de_out = 0;
	uint8_t type_tmp = 0;
	uint8_t alph_tmp = 0;
	uint8_t flag_filter_in = 0;

	if(OneUnit)
	{
		uint64_t eliminate=~(ELIMINATE);
		uint64_t start_get_kmer = 0;

		if(i - anchor_km_pre_i > 1)
		{
			pre_j = anchor_km_pre_i << 1;
			for(j = anchor_km_pre_i; j < i; j++)
			{
				start_get_kmer = j << 1;
				anchor_a = bufMK2[start_get_kmer];

				anchor = ((anchor_a<<2)&eliminate);

				if((anchor != kmer_pre) && (j > anchor_km_pre_i))
				{
					if(flag_filter_in)
					{
						if(in_j || out_j)
						{
							in_j ^= de_in;
							out_j ^= de_out;
						}
						kmer_ins[kmer_i] = in_j;
						kmer_outs[kmer_i] = out_j;
						kmer_indexs[kmer_i] = pre_j;

						kmer_i++;
					}

					in_j = 0;
					out_j = 0;
					de_in = 0;
					de_out = 0;
					flag_filter_in = 0;

					pre_j = start_get_kmer;
				}

				type_tmp = (anchor_a >> 62) & 0X3;
				if(type_tmp == 1)
				{
					de_out |= (one_tmp << bufMK2[start_get_kmer + 1]);
				}
				else if(type_tmp == 2)
				{
					de_in |= (one_tmp << bufMK2[start_get_kmer + 1]);
				}
				else
				{
					//debug
					if(flag_filter_in == 1)	printf("Err: %u %"PRId64"", type_tmp, anchor_a);

					in_j = (bufMK2[start_get_kmer + 1] >> 4) & 0Xf;
					out_j = bufMK2[start_get_kmer + 1] & 0Xf;
					flag_filter_in = 1;
				}

				kmer_pre = anchor;
			}

			if(flag_filter_in)
			{
				if(in_j || out_j)
				{
					in_j ^= de_in;
					out_j ^= de_out;
				}
				kmer_ins[kmer_i] = in_j;
				kmer_outs[kmer_i] = out_j;
				kmer_indexs[kmer_i] = pre_j;

				kmer_i++;
			}


		}
		else
		{
			start_get_kmer = anchor_km_pre_i * kmer_num_filter;
			kmer_pre_a = bufMK2[start_get_kmer];

			type_tmp = (kmer_pre_a >> 62) & 0X3;

			if(!type_tmp)
			{
				in_j = (bufMK2[start_get_kmer + 1] >> 4) & 0Xf;
				out_j = bufMK2[start_get_kmer + 1] & 0Xf;

				if(in_j || out_j)
				{
					in_j ^= de_in;
					out_j ^= de_out;
				}
				kmer_ins[kmer_i] = in_j;
				kmer_outs[kmer_i] = out_j;
				kmer_indexs[kmer_i] = start_get_kmer;

				kmer_i++;

			}
		}

		for(kmer_k = 0; kmer_k < kmer_i; kmer_k++)
		{
			kmer_pre_a = bufMK2[kmer_indexs[kmer_k]];
			in_j = kmer_ins[kmer_k];
			out_j = kmer_outs[kmer_k];

			if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
			else	in_linear = 0;

			if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
			else	out_linear = 0;

			if(in_linear && out_linear)			//linear
				type_f = 0;
			else if(in_linear && !out_linear)	//single-in multi-out
				type_f = 1;
			else if(!in_linear && out_linear)	//multi-in single-out
				type_f = 2;
			else								//multi-in multi-out
				type_f = 3;

			type_s = type_f;
			if((type_f == 0) || (type_f == 1))
			{
				pre_flag = 0;

				for(pre_k = 0; pre_k < kmer_k; pre_k++)
					if(kmer_ins[pre_k] & in_j)
					{
						pre_flag = 1;
						break;
					}

				for(pre_k = kmer_k + 1; (pre_k < kmer_i) && (!pre_flag); pre_k++)
					if(kmer_ins[pre_k] & in_j)
					{
						pre_flag = 1;
						break;
					}

				if(pre_flag)
				{
					if(type_f == 0)	type_s = 2;
					else	type_s = 3;
				}
			}

			kmer_tmp = ((kmer_pre_a >> 60) & 0Xf);

			write_tmp = (kmer_pre_a << 2) & eliminate;

			alph_tmp = (kmer_tmp & 0X3);
			alphabet_num_new[0][alph_tmp]++;

			type_s |= (input_types[in_j] << 2);

			type_s |= (out_j << 4);

			fwrite(&write_tmp, 8, 1, fp_tomerge);
			fwrite(&type_s, 1, 1, fp_tomerge_t);
#ifdef	GFA_COM
			type_e = ((in_j << 4)|out_j);
			fwrite(&type_e, 1, 1, fp_edges_t);
#endif
		}
	}
	else
	{
		uint64_t anchor_kms[kmer_num];
		uint64_t anchor_kms_pre[kmer_num];
		uint64_t kmer_pre_as[kmer_num_new];
		uint64_t char_tmp = 0;
		uint8_t kmer_num_tmp_new = (kmer_num_new << 3);

		if(i - anchor_km_pre_i > 1)
		{
			pre_j = anchor_km_pre_i * kmer_num_fm;
			for(j = anchor_km_pre_i; j < i; j++)
			{
				char_tmp = j * kmer_num_fm;

				anchor = ((bufMK2[char_tmp + kmer_num - 1] << leftmove32) & KMask);

				if((anchor != kmer_pre) && (j > anchor_km_pre_i))
				{
					if(flag_filter_in)
					{
						if(in_j || out_j)
						{
							in_j ^= de_in;
							out_j ^= de_out;
						}
						kmer_ins[kmer_i] = in_j;
						kmer_outs[kmer_i] = out_j;
						kmer_indexs[kmer_i] = pre_j;

						kmer_i++;
					}

					in_j = 0;
					out_j = 0;
					de_in = 0;
					de_out = 0;
					flag_filter_in = 0;

					pre_j = char_tmp;//j
				}

				type_tmp = (bufMK2[char_tmp] >> 62) & 0X3;
				if(type_tmp == 1)
				{
					de_out |= (one_tmp << (bufMK2[char_tmp + kmer_num_fm - 1] & 0X3));
				}
				else if(type_tmp == 2)
				{
					de_in |= (one_tmp << (bufMK2[char_tmp + kmer_num_fm - 1] & 0X3));
				}
				else
				{
					in_j = (bufMK2[char_tmp + kmer_num_fm - 1] >> 4) & 0Xf;
					out_j = bufMK2[char_tmp + kmer_num_fm - 1] & 0Xf;
					flag_filter_in = 1;
				}

				kmer_pre = anchor;
			}

			if(flag_filter_in)
			{
				if(in_j || out_j)
				{
					in_j ^= de_in;
					out_j ^= de_out;
				}
				kmer_ins[kmer_i] = in_j;
				kmer_outs[kmer_i] = out_j;
				kmer_indexs[kmer_i] = pre_j;

				kmer_i++;
			}


		}
		else
		{
			char_tmp = anchor_km_pre_i * kmer_num_fm;

			memcpy(kmer_pre_as, bufMK2 + char_tmp, kmer_num_tmp_new);

			type_tmp = (kmer_pre_as[0] >> 62) & 0X3;

			if(!type_tmp)
			{
				in_j = (kmer_pre_as[kmer_num_fm - 1] >> 4) & 0Xf;
				out_j = kmer_pre_as[kmer_num_fm - 1] & 0Xf;

				if(in_j || out_j)
				{
					in_j ^= de_in;
					out_j ^= de_out;
				}
				kmer_ins[kmer_i] = in_j;
				kmer_outs[kmer_i] = out_j;
				kmer_indexs[kmer_i] = char_tmp;

				kmer_i++;

			}
		}

		for(kmer_k = 0; kmer_k < kmer_i; kmer_k++)
		{
			memcpy(kmer_pre_as, bufMK2 + kmer_indexs[kmer_k], kmer_num_tmp);

			in_j = kmer_ins[kmer_k];
			out_j = kmer_outs[kmer_k];

			if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
			else	in_linear = 0;

			if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
			else	out_linear = 0;

			if(in_linear && out_linear)			//linear
				type_f = 0;
			else if(in_linear && !out_linear)	//single-in multi-out
				type_f = 1;
			else if(!in_linear && out_linear)	//multi-in single-out
				type_f = 2;
			else								//multi-in multi-out
				type_f = 3;

			type_s = type_f;
			if((type_f == 0) || (type_f == 1))
			{
				pre_flag = 0;

				for(pre_k = 0; pre_k < kmer_k; pre_k++)
					if(kmer_ins[pre_k] & in_j)
					{
						pre_flag = 1;
						break;
					}

				for(pre_k = kmer_k + 1; (pre_k < kmer_i) && (!pre_flag); pre_k++)
					if(kmer_ins[pre_k] & in_j)
					{
						pre_flag = 1;
						break;
					}

				if(pre_flag)
				{
					if(type_f == 0)	type_s = 2;
					else	type_s = 3;
				}
			}

			kmer_tmp = ((kmer_pre_as[0] >> 60) & 0Xf);

			kmer_pre_as[0] <<= 2;
			kmer_pre_as[kmer_num - 1] &= KMask;

			alph_tmp = (kmer_tmp & 0X3);
			alphabet_num_new[0][alph_tmp]++;

			type_s |= (input_types[in_j] << 2);

			type_s |= (out_j << 4);

			fwrite(kmer_pre_as, kmer_num_tmp, 1, fp_tomerge);
			fwrite(&type_s, 1, 1, fp_tomerge_t);
#ifdef	GFA_COM
			type_e = ((in_j << 4)|out_j);
			fwrite(&type_e, 1, 1, fp_edges_t);
#endif
		}
	}
}

void *multiTagcheck(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	int num = argT->num;
	int totalNum = argT->thread;

	uint64_t i = 0;

	uint64_t start=splitMK2[num], end=splitMK2[num+1];
	uint64_t bufMK1_len=end-start;

	if(bufMK1_len==0) return (void *)NULL;

	uint64_t eliminate=~(ELIMINATE);
	uint64_t eliminate_m = (eliminate << 2);
	uint64_t eliminate_a = ~(((uint64_t)1<<((31-KMER_LENGTH)<<1))-1);//select KMER_LENGTH+1 mers
	uint64_t anchor = 0;
	uint64_t anchor_a = 0;

	uint64_t kmer_pre = 0;
	uint64_t kmer_pre_a = 0;
	uint64_t kmer_pre_p = 0;
	uint8_t input = 0;
	uint8_t output = 0;
	uint8_t in_j = 0;
	uint8_t out_j = 0;
	uint8_t in_linear = 0;
	uint8_t out_linear = 0;
	uint8_t one_tmp = 1;
	uint8_t type_f = 0;
	uint8_t type_s = 0;

	uint64_t well_tmp = 4;
	uint8_t alph_tmp = 0;
	uint8_t alph_pre = 0;
	uint64_t char_tmp = 0;
	uint64_t w_m_p_t = 0;
	uint64_t anchor_km = 0;
	uint64_t anchor_km_pre = 0;
	uint64_t anchor_km_pre_i = 0;
	uint64_t j = 0;

	uint8_t kmer_tmp = 0;
	uint8_t pres[128];//64 512 1024
	uint8_t pres_i = 0;
	uint8_t pre_k = 0;
	uint64_t pre_i = 0;
	uint64_t pre_j = 0;
	uint8_t pre_tmp = 0;

	uint8_t pre_flag = 0;
	uint8_t kmer_num_tmp = (kmer_num << 3);

	anchor_km_pre_i = start;

	uint8_t fl_f = 0;

	uint64_t start_get = 0;
	uint64_t w_f_i = 0;
	uint64_t w_f_o = 0;

	if(OneUnit)
	{
		if(flag_filter)
		{
			uint8_t kmer_ins[24];
			uint8_t kmer_outs[24];
			uint8_t kmer_types[24];
			uint64_t kmer_indexs[24];
			uint32_t kmer_n = 0;
			uint8_t kmer_i = 0;
			uint8_t kmer_k = 0;
			//uint8_t multi_kmer = 0;
			uint64_t start_get_kmer = 0;
			uint8_t type_c = 0;
			uint8_t add_edge = 0;
			uint8_t aph_i = 0;

			for(i = start; i < end; i++)
			{
				start_get = i * kmer_num_filter;

				anchor_km = ((bufMK2[start_get] << 4)&eliminate_m);

				if((anchor_km != anchor_km_pre) && (i > start))
				{
					if(i - anchor_km_pre_i > 1)
					{
						pre_j = anchor_km_pre_i * kmer_num_filter;
						for(j = anchor_km_pre_i; j < i; j++)
						{
							start_get_kmer = j * kmer_num_filter;
							anchor_a = bufMK2[start_get_kmer];

							anchor = ((anchor_a<<4)&eliminate);

							if((anchor != kmer_pre) && (j > anchor_km_pre_i))
							{
								kmer_ins[kmer_i] = in_j;
								kmer_outs[kmer_i] = out_j;
								kmer_indexs[kmer_i] = pre_j;

								if((kmer_n > filter_min) && (kmer_n < filter_max))
									kmer_types[kmer_i] = 0;
								else	kmer_types[kmer_i] = 1;

								kmer_i++;

								in_j = 0;
								out_j = 0;
								pre_j = start_get_kmer;//j

								kmer_n = 0;
							}

							input = ((anchor_a >> 60) & 0Xf);
							output = ((anchor_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);

							fl_f = (input >> 2);
							input &= 0X3;

							if(fl_f == 1)
							{
								out_j |= (one_tmp << output);
							}
							else if(fl_f == 2)
							{
								in_j |= (one_tmp << input);
							}
							else if(fl_f == 0)
							{
								in_j |= (one_tmp << input);
								out_j |= (one_tmp << output);
							}

							kmer_pre = anchor;

							kmer_n += bufMK2[start_get_kmer + 1];

						}

						kmer_ins[kmer_i] = in_j;
						kmer_outs[kmer_i] = out_j;
						kmer_indexs[kmer_i] = pre_j;

						if((kmer_n > filter_min) && (kmer_n < filter_max))
							kmer_types[kmer_i] = 0;
						else	kmer_types[kmer_i] = 1;

						kmer_i++;
					}
					else
					{
						start_get_kmer = anchor_km_pre_i * kmer_num_filter;
						kmer_pre_a = bufMK2[start_get_kmer];
						type_s = 0;

						fl_f = (kmer_pre_a >> 60);

						input = (fl_f & 0X3);

						fl_f >>= 2;

						if(fl_f == 1)
						{
							type_s = 2;

							in_j = 0;
							output = ((kmer_pre_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);
							out_j = (one_tmp << output);
						}
						else if(fl_f == 2)
						{
							type_s = 1;

							in_j = (one_tmp << input);
							out_j = 0;
						}
						else if(fl_f == 0)
						{
							in_j = (one_tmp << input);
							output = ((kmer_pre_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);
							out_j = (one_tmp << output);
						}
						else
						{
							in_j = 0;
							out_j = 0;
						}

						kmer_n = bufMK2[start_get_kmer + 1];

						kmer_ins[kmer_i] = in_j;
						kmer_outs[kmer_i] = out_j;
						kmer_indexs[kmer_i] = start_get_kmer;

						if((kmer_n > filter_min) && (kmer_n < filter_max))
							kmer_types[kmer_i] = 0;
						else	kmer_types[kmer_i] = 1;

						kmer_i++;
					}

					for(kmer_k = 0; kmer_k < kmer_i; kmer_k++)
					{
						kmer_pre_a = bufMK2[kmer_indexs[kmer_k]];
						in_j = kmer_ins[kmer_k];
						out_j = kmer_outs[kmer_k];
						type_c = kmer_types[kmer_k];

						if(type_c)
						{
							if(in_j)
							{
								add_edge = (kmer_pre_a >> one_suf) & 0X3;//last char

								kmer_pre_a >>= one_sufp;
								kmer_pre_a <<= (one_suf + 6);
								kmer_pre_a >>= 4;

								for(aph_i = 0; aph_i < 4; aph_i++)
								{
									if(!((in_j >> aph_i) & 0X1))	continue;

									kmer_pre_a <<= 4;
									kmer_pre_a >>= 4;
									kmer_pre_a |= firstls[aph_i + 4];
									kmer_pre_a |= type_one;

									write_filterout[num][(w_f_o << 1)] = kmer_pre_a;

									write_filterout[num][(w_f_o << 1) + 1] = add_edge;//suf
									w_f_o++;
								}

								kmer_pre_a = bufMK2[kmer_indexs[kmer_k]];
							}

							if(out_j)
							{
								add_edge = (kmer_pre_a >> 58) & 0X3;//first char
								kmer_pre_a <<= 6;
								kmer_pre_a >>= 2;
								kmer_pre_a |= type_two;

								for(aph_i = 0; aph_i < 4; aph_i++)
								{
									if(!((out_j >> aph_i) & 0X1))	continue;

									kmer_pre_a >>= (one_sufp + 2);
									kmer_pre_a <<= 2;
									kmer_pre_a |= aph_i;
									kmer_pre_a <<= one_sufp;

									write_filterout[num][(w_f_o << 1)] = kmer_pre_a;

									write_filterout[num][(w_f_o << 1) + 1] = add_edge;//pre
									w_f_o++;
								}
							}
						}
						else
						{
							kmer_pre_a <<= 4;
							kmer_pre_a &= eliminate;
							kmer_pre_a >>= 2;
							write_filterin[num][(w_f_i << 1)] = kmer_pre_a;

							add_edge = ((in_j << 4) | out_j);
							write_filterin[num][(w_f_i << 1) + 1] = add_edge;

							w_f_i++;
						}
					}

					in_j = 0;
					out_j = 0;
					kmer_n = 0;
					kmer_i = 0;

					anchor_km_pre_i = i;
				}

				anchor_km_pre = anchor_km;
			}

			if(num != totalNum - 1)
			{
				in_j = 0;
				out_j = 0;
				kmer_n = 0;

				if(i - anchor_km_pre_i > 1)
				{
					pre_j = anchor_km_pre_i * kmer_num_filter;
					for(j = anchor_km_pre_i; j < i; j++)
					{
						start_get_kmer = j * kmer_num_filter;
						anchor_a = bufMK2[start_get_kmer];

						anchor = ((anchor_a<<4)&eliminate);

						if((anchor != kmer_pre) && (j > anchor_km_pre_i))
						{
							kmer_ins[kmer_i] = in_j;
							kmer_outs[kmer_i] = out_j;
							kmer_indexs[kmer_i] = pre_j;

							if((kmer_n > filter_min) && (kmer_n < filter_max))
								kmer_types[kmer_i] = 0;
							else	kmer_types[kmer_i] = 1;

							kmer_i++;

							in_j = 0;
							out_j = 0;
							pre_j = start_get_kmer;//j

							kmer_n = 0;
						}

						input = ((anchor_a >> 60) & 0Xf);
						output = ((anchor_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);

						fl_f = (input >> 2);
						input &= 0X3;

						if(fl_f == 1)
						{
							out_j |= (one_tmp << output);
						}
						else if(fl_f == 2)
						{
							in_j |= (one_tmp << input);
						}
						else if(fl_f == 0)
						{
							in_j |= (one_tmp << input);
							out_j |= (one_tmp << output);
						}

						kmer_pre = anchor;

						kmer_n += bufMK2[start_get_kmer + 1];
					}

					kmer_ins[kmer_i] = in_j;
					kmer_outs[kmer_i] = out_j;
					kmer_indexs[kmer_i] = pre_j;

					if((kmer_n > filter_min) && (kmer_n < filter_max))
						kmer_types[kmer_i] = 0;
					else	kmer_types[kmer_i] = 1;

					kmer_i++;
				}
				else
				{
					start_get_kmer = anchor_km_pre_i * kmer_num_filter;
					kmer_pre_a = bufMK2[start_get_kmer];
					type_s = 0;

					fl_f = (kmer_pre_a >> 60);

					input = (fl_f & 0X3);

					fl_f >>= 2;

					if(fl_f == 1)
					{
						type_s = 2;

						in_j = 0;
						output = ((kmer_pre_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);
						out_j = (one_tmp << output);
					}
					else if(fl_f == 2)
					{
						type_s = 1;

						in_j = (one_tmp << input);
						out_j = 0;
					}
					else if(fl_f == 0)
					{
						in_j = (one_tmp << input);
						output = ((kmer_pre_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);
						out_j = (one_tmp << output);
					}
					else
					{
						in_j = 0;
						out_j = 0;
					}

					kmer_n = bufMK2[start_get_kmer + 1];

					kmer_ins[kmer_i] = in_j;
					kmer_outs[kmer_i] = out_j;
					kmer_indexs[kmer_i] = start_get_kmer;

					if((kmer_n > filter_min) && (kmer_n < filter_max))
						kmer_types[kmer_i] = 0;
					else	kmer_types[kmer_i] = 1;

					kmer_i++;
				}

				for(kmer_k = 0; kmer_k < kmer_i; kmer_k++)
				{
					kmer_pre_a = bufMK2[kmer_indexs[kmer_k]];
					in_j = kmer_ins[kmer_k];
					out_j = kmer_outs[kmer_k];
					type_c = kmer_types[kmer_k];

					if(type_c)
					{
						if(in_j)
						{
							add_edge = (kmer_pre_a >> one_suf) & 0X3;//last char

							kmer_pre_a >>= one_sufp;
							kmer_pre_a <<= (one_suf + 6);
							kmer_pre_a >>= 4;

							for(aph_i = 0; aph_i < 4; aph_i++)
							{
								if(!((in_j >> aph_i) & 0X1))	continue;

								kmer_pre_a <<= 4;
								kmer_pre_a >>= 4;
								kmer_pre_a |= firstls[aph_i + 4];
								kmer_pre_a |= type_one;

								write_filterout[num][(w_f_o << 1)] = kmer_pre_a;

								write_filterout[num][(w_f_o << 1) + 1] = add_edge;//suf
								w_f_o++;
							}

							kmer_pre_a = bufMK2[kmer_indexs[kmer_k]];
						}

						if(out_j)
						{
							add_edge = (kmer_pre_a >> 58) & 0X3;//first char
							kmer_pre_a <<= 6;
							//kmer_pre_a &= eliminate_m;
							kmer_pre_a >>= 2;
							kmer_pre_a |= type_two;

							for(aph_i = 0; aph_i < 4; aph_i++)
							{
								if(!((out_j >> aph_i) & 0X1))	continue;

								kmer_pre_a >>= (one_sufp + 2);
								kmer_pre_a <<= 2;
								kmer_pre_a |= aph_i;
								kmer_pre_a <<= one_sufp;

								write_filterout[num][(w_f_o << 1)] = kmer_pre_a;

								write_filterout[num][(w_f_o << 1) + 1] = add_edge;//pre
								w_f_o++;
							}
						}
					}
					else
					{
						kmer_pre_a <<= 4;
						kmer_pre_a &= eliminate;
						kmer_pre_a >>= 2;
						write_filterin[num][(w_f_i << 1)] = kmer_pre_a;

						add_edge = ((in_j << 4) | out_j);
						write_filterin[num][(w_f_i << 1) + 1] = add_edge;

						w_f_i++;
					}
				}

			}
			else	last_block_num = anchor_km_pre_i;
		}
		else
		{
			for(i = start; i < end; i++)
			{
#ifdef	FIRST_LAST_KMER
				anchor_km = ((bufMK2[i] << 4)&eliminate_m);
#else
				anchor_km = ((bufMK2[i] << 2)&eliminate_m);
#endif
				if((anchor_km != anchor_km_pre) && (i > start))
				{
					if(i - anchor_km_pre_i > 1)
					{
						pre_j = anchor_km_pre_i;
						for(j = anchor_km_pre_i; j < i; j++)
						{
							anchor_a = bufMK2[j];
#ifdef	FIRST_LAST_KMER
							anchor = ((anchor_a<<4)&eliminate);
#else
							anchor = ((anchor_a<<2)&eliminate);
#endif

							if((anchor != kmer_pre) && (j > anchor_km_pre_i))
							{
								if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
								else	in_linear = 0;

								if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
								else	out_linear = 0;

								if(in_linear && out_linear)			//linear
									type_f = 0;
								else if(in_linear && !out_linear)	//single-in multi-out
									type_f = 1;
								else if(!in_linear && out_linear)	//multi-in single-out
									type_f = 2;
								else								//multi-in multi-out
									type_f = 3;

								type_s = type_f;
								if((type_f == 0) || (type_f == 1))
								{
									pre_flag = 0;
									pre_tmp = __builtin_ffs(in_j) - 1;
									for(pre_k = 0; pre_k < pre_j - anchor_km_pre_i; pre_k++)
										if(pres[pre_k] == pre_tmp)
										{
											pre_flag = 1;
											break;
										}

									for(pre_k = j - anchor_km_pre_i; (pre_k < i - anchor_km_pre_i) && (!pre_flag); pre_k++)
										if(pres[pre_k] == pre_tmp)
										{
											pre_flag = 1;
											break;
										}

									if(pre_flag)
									{
										if(type_f == 0)	type_s = 2;
										else	type_s = 3;
									}
								}

#ifdef	FIRST_LAST_KMER
								kmer_tmp = ((kmer_pre_a >> 58) & 0Xf);
								write_kmermerge[num][w_m_p_t] = ((kmer_pre_a << 4) & eliminate);
#else
								kmer_tmp = (kmer_pre_a >> 60);
								write_kmermerge[num][w_m_p_t] = ((kmer_pre_a << 2) & eliminate);
#endif
								alph_tmp = (kmer_tmp & 0X3);
								alphabet_num_new[num][alph_tmp]++;

								type_s |= (input_types[in_j] << 2);

								type_s |= (out_j << 4);

								write_kmermerge_t[num][w_m_p_t] = type_s;
#ifdef	GFA_COM
								edges_t[num][w_m_p_t] = ((in_j << 4) | out_j);
#endif
								w_m_p_t++;

								in_j = 0;
								out_j = 0;

								pre_j = j;
							}

#ifdef	FIRST_LAST_KMER
							input = ((anchor_a >> 60) & 0Xf);
							output = ((anchor_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);

							fl_f = (input >> 2);
							input &= 0X3;

							if(fl_f == 1)
							{
								out_j |= (one_tmp << output);
							}
							else if(fl_f == 2)
							{
								in_j |= (one_tmp << input);

							}
							else
							{
								in_j |= (one_tmp << input);
								out_j |= (one_tmp << output);
							}
#else
							input = ((anchor_a >> 62) & 0X3);
							output = ((anchor_a >> (64 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);

							in_j |= (one_tmp << input);
							out_j |= (one_tmp << output);
#endif

							kmer_pre = anchor;
							kmer_pre_a = anchor_a;
						}

						if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
						else	in_linear = 0;

						if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
						else	out_linear = 0;

						if(in_linear && out_linear)			//linear
							type_f = 0;
						else if(in_linear && !out_linear)	//single-in multi-out
							type_f = 1;
						else if(!in_linear && out_linear)	//multi-in single-out
							type_f = 2;
						else								//multi-in multi-out
							type_f = 3;

						type_s = type_f;
						if((type_f == 0) || (type_f == 1))
						{
							pre_flag = 0;
							pre_tmp = __builtin_ffs(in_j) - 1;
							for(pre_k = 0; pre_k < pre_j - anchor_km_pre_i; pre_k++)
								if(pres[pre_k] == pre_tmp)
								{
									pre_flag = 1;
									break;
								}

							for(pre_k = j - anchor_km_pre_i; (pre_k < i - anchor_km_pre_i) && (!pre_flag); pre_k++)
								if(pres[pre_k] == pre_tmp)
								{
									pre_flag = 1;
									break;
								}

							if(pre_flag)
							{
								if(type_f == 0)	type_s = 2;
								else	type_s = 3;
							}
						}
					}
					else
					{
						kmer_pre_a = bufMK2[anchor_km_pre_i];
						type_s = 0;

#ifdef	FIRST_LAST_KMER
						fl_f = (kmer_pre_a >> 60);

						input = (fl_f & 0X3);
						in_j = (one_tmp << input);

						fl_f >>= 2;

						if(fl_f == 1)
						{
							type_s = 2;
							output = ((kmer_pre_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);
							out_j = (one_tmp << output);

						}
						else if(fl_f == 2)
						{
							type_s = 1;
							out_j = 0;

						}
						else
						{
							output = ((kmer_pre_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);
							out_j = (one_tmp << output);
						}

#else
						output = ((kmer_pre_a >> (64 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);

						input = (kmer_pre_a >> 62);
						in_j = (one_tmp << input);

						out_j = (one_tmp << output);
#endif

					}

#ifdef	FIRST_LAST_KMER
					kmer_tmp = ((kmer_pre_a >> 58) & 0Xf);
					write_kmermerge[num][w_m_p_t] = ((kmer_pre_a << 4) & eliminate);
#else
					kmer_tmp = (kmer_pre_a >> 60);
					write_kmermerge[num][w_m_p_t] = ((kmer_pre_a << 2) & eliminate);
#endif

					alph_tmp = (kmer_tmp & 0X3);
					alphabet_num_new[num][alph_tmp]++;

					type_s |= (input_types[in_j] << 2);
					type_s |= (out_j << 4);

					write_kmermerge_t[num][w_m_p_t] = type_s;
#ifdef	GFA_COM
					edges_t[num][w_m_p_t] = ((in_j << 4) | out_j);
#endif
					w_m_p_t++;

					in_j = 0;
					out_j = 0;

					anchor_km_pre_i = i;

					pres_i = 0;
				}

				anchor_km_pre = anchor_km;

#ifdef	FIRST_LAST_KMER
				kmer_tmp = (bufMK2[i] >> 60);
				if((kmer_tmp >> 2) != 1)
				{
					pres[pres_i] = (kmer_tmp & 0X3);
				}
				else
				{
					pres[pres_i] = 4;
				}
				pres_i++;
#else
				pres[pres_i] = ((bufMK2[i] >> 62) & 0X3);
				pres_i++;
#endif

			}

			if(num != totalNum - 1)
			{
				if(i - anchor_km_pre_i > 1)
				{
					pre_j = anchor_km_pre_i;
					for(j = anchor_km_pre_i; j < i; j++)
					{
						anchor_a = bufMK2[j];
#ifdef	FIRST_LAST_KMER
						anchor = ((anchor_a<<4)&eliminate);
#else
						anchor = ((anchor_a<<2)&eliminate);
#endif

						if((anchor != kmer_pre) && (j > anchor_km_pre_i))
						{
							if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
							else	in_linear = 0;

							if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
							else	out_linear = 0;

							if(in_linear && out_linear)			//linear
								type_f = 0;
							else if(in_linear && !out_linear)	//single-in multi-out
								type_f = 1;
							else if(!in_linear && out_linear)	//multi-in single-out
								type_f = 2;
							else								//multi-in multi-out
								type_f = 3;

							type_s = type_f;
							if((type_f == 0) || (type_f == 1))
							{
								pre_flag = 0;
								pre_tmp = __builtin_ffs(in_j) - 1;
								for(pre_k = 0; pre_k < pre_j - anchor_km_pre_i; pre_k++)
									if(pres[pre_k] == pre_tmp)
									{
										pre_flag = 1;
										break;
									}

								for(pre_k = j - anchor_km_pre_i; (pre_k < i - anchor_km_pre_i) && (!pre_flag); pre_k++)
									if(pres[pre_k] == pre_tmp)
									{
										pre_flag = 1;
										break;
									}

								if(pre_flag)
								{
									if(type_f == 0)	type_s = 2;
									else	type_s = 3;
								}
							}
#ifdef	FIRST_LAST_KMER
							kmer_tmp = ((kmer_pre_a >> 58) & 0Xf);
							write_kmermerge[num][w_m_p_t] = ((kmer_pre_a << 4) & eliminate);
#else
							kmer_tmp = (kmer_pre_a >> 60);
							write_kmermerge[num][w_m_p_t] = ((kmer_pre_a << 2) & eliminate);
#endif
							alph_tmp = (kmer_tmp & 0X3);
							alphabet_num_new[num][alph_tmp]++;

							//type_s |= ((kmer_tmp >> 2) << 2);
							type_s |= (input_types[in_j] << 2);
							type_s |= (out_j << 4);

							write_kmermerge_t[num][w_m_p_t] = type_s;
#ifdef	GFA_COM
							edges_t[num][w_m_p_t] = ((in_j << 4) | out_j);
#endif
							w_m_p_t++;

							in_j = 0;
							out_j = 0;

							pre_j = j;
						}
#ifdef	FIRST_LAST_KMER
						input = ((anchor_a >> 60) & 0Xf);
						output = ((anchor_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);

						fl_f = (input >> 2);
						input &= 0X3;

						if(fl_f == 1)
						{
							out_j |= (one_tmp << output);

						}
						else if(fl_f == 2)
						{
							in_j |= (one_tmp << input);

						}
						else
						{
							in_j |= (one_tmp << input);
							out_j |= (one_tmp << output);
						}
#else
						input = ((anchor_a >> 62) & 0X3);
						output = ((anchor_a >> (64 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);

						in_j |= (one_tmp << input);
						out_j |= (one_tmp << output);
#endif

						kmer_pre = anchor;
						kmer_pre_a = anchor_a;
					}

					if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
					else	in_linear = 0;

					if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
					else	out_linear = 0;

					if(in_linear && out_linear)			//linear
						type_f = 0;
					else if(in_linear && !out_linear)	//single-in multi-out
						type_f = 1;
					else if(!in_linear && out_linear)	//multi-in single-out
						type_f = 2;
					else								//multi-in multi-out
						type_f = 3;

					type_s = type_f;
					if((type_f == 0) || (type_f == 1))
					{
						pre_flag = 0;
						pre_tmp = __builtin_ffs(in_j) - 1;
						for(pre_k = 0; pre_k < pre_j - anchor_km_pre_i; pre_k++)
							if(pres[pre_k] == pre_tmp)
							{
								pre_flag = 1;
								break;
							}

						for(pre_k = j - anchor_km_pre_i; (pre_k < i - anchor_km_pre_i) && (!pre_flag); pre_k++)
							if(pres[pre_k] == pre_tmp)
							{
								pre_flag = 1;
								break;
							}

						if(pre_flag)
						{
							if(type_f == 0)	type_s = 2;
							else	type_s = 3;
						}
					}
				}
				else
				{
					kmer_pre_a = bufMK2[anchor_km_pre_i];
					type_s = 0;

#ifdef	FIRST_LAST_KMER
					fl_f = (kmer_pre_a >> 60);

					input = fl_f & 0X3;
					in_j = (one_tmp << input);

					fl_f >>= 2;

					if(fl_f == 1)
					{
						type_s = 2;

						output = ((kmer_pre_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);
						out_j = (one_tmp << output);
					}
					else if(fl_f == 2)
					{
						type_s = 1;
						out_j = 0;
					}
					else
					{
						output = ((kmer_pre_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);
						out_j = (one_tmp << output);
					}

#else
					output = ((kmer_pre_a >> (64 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);

					input = (kmer_pre_a >> 62);
					in_j = (one_tmp << input);

					out_j = (one_tmp << output);
#endif

				}

#ifdef	FIRST_LAST_KMER
				kmer_tmp = ((kmer_pre_a >> 58) & 0Xf);
				write_kmermerge[num][w_m_p_t] = ((kmer_pre_a << 4) & eliminate);
#else
				kmer_tmp = (kmer_pre_a >> 60);
				write_kmermerge[num][w_m_p_t] = ((kmer_pre_a << 2) & eliminate);
#endif

				alph_tmp = (kmer_tmp & 0X3);
				alphabet_num_new[num][alph_tmp]++;

				type_s |= (input_types[in_j] << 2);
				type_s |= (out_j << 4);

				write_kmermerge_t[num][w_m_p_t] = type_s;
#ifdef	GFA_COM
				edges_t[num][w_m_p_t] = ((in_j << 4) | out_j);
#endif
				w_m_p_t++;

			}
			else	last_block_num = anchor_km_pre_i;

		}
	}
	else
	{
		uint64_t anchor_kms[kmer_num];
		uint64_t anchor_kms_pre[kmer_num];
		uint64_t kmer_pre_as[kmer_num_new];

		uint8_t kmer_num_tmp_new = (kmer_num_new << 3);


		if(flag_filter)
		{
			uint8_t kmer_ins[24];
			uint8_t kmer_outs[24];
			uint8_t kmer_types[24];
			uint64_t kmer_indexs[24];
			uint32_t kmer_n = 0;
			uint8_t kmer_i = 0;
			uint8_t kmer_k = 0;
			uint8_t aph_i = 0;
			uint8_t type_c = 0;
			int16_t tra_i = 0;
			uint8_t shc = 0;
			uint8_t shc_pre = 0;
			uint8_t shtype = 0;
			uint8_t add_edge = 0;

			for(i = start; i < end; i++)
			{
				start_get = (i * kmer_num_filter);

				getkmerno_pss(anchor_kms, bufMK2, start_get)

				if(memcmp(anchor_kms, anchor_kms_pre, kmer_num_tmp) && (i > start))
				{
					if(i - anchor_km_pre_i > 1)
					{
						pre_j = anchor_km_pre_i * kmer_num_filter;
						for(j = anchor_km_pre_i; j < i; j++)
						{
							char_tmp = j * kmer_num_filter;

							anchor = ((bufMK2[char_tmp + kmer_num - 1] << leftmove32) & KMask);

							if((anchor != kmer_pre) && (j > anchor_km_pre_i))
							{
								kmer_ins[kmer_i] = in_j;
								kmer_outs[kmer_i] = out_j;
								kmer_indexs[kmer_i] = pre_j;

								if((kmer_n > filter_min) && (kmer_n < filter_max))
									kmer_types[kmer_i] = 0;
								else	kmer_types[kmer_i] = 1;

								kmer_i++;

								in_j = 0;
								out_j = 0;
								kmer_n = 0;

								pre_j = char_tmp;//j
							}

							input = ((bufMK2[char_tmp] >> tempMove_r_new3) & 0X3);
							output = ((bufMK2[char_tmp + kmer_num - 1] >> keyMove) & 0X3);

							fl_f = (bufMK2[char_tmp + kmer_num_new - 1] & 0X3);

							if(fl_f == 1)
							{
								out_j |= (one_tmp << output);
							}
							else if(fl_f == 2)
							{
								in_j |= (one_tmp << input);
							}
							else if(fl_f == 0)
							{
								in_j |= (one_tmp << input);
								out_j |= (one_tmp << output);
							}

							kmer_pre = anchor;

							kmer_n += (bufMK2[char_tmp + kmer_num_filter - 1] >> 32);
						}

						kmer_ins[kmer_i] = in_j;
						kmer_outs[kmer_i] = out_j;
						kmer_indexs[kmer_i] = pre_j;

						if((kmer_n > filter_min) && (kmer_n < filter_max))
							kmer_types[kmer_i] = 0;
						else	kmer_types[kmer_i] = 1;

						kmer_i++;
					}
					else
					{
						char_tmp = anchor_km_pre_i * kmer_num_filter;

						memcpy(kmer_pre_as, bufMK2 + char_tmp, kmer_num_tmp_new);

						fl_f = (kmer_pre_as[kmer_num_new - 1] & 0X3);

						if(fl_f == 1)
						{
							in_j = 0;
							output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
							out_j = (one_tmp << output);
						}
						else if(fl_f == 2)
						{
							input = ((kmer_pre_as[0] >> tempMove_r_new3) & 0X3);
							in_j = (one_tmp << input);
							out_j = 0;
						}
						else if(fl_f == 0)
						{
							input = ((kmer_pre_as[0] >> tempMove_r_new3) & 0X3);
							in_j = (one_tmp << input);
							output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
							out_j = (one_tmp << output);
						}
						else
						{
							in_j = 0;
							out_j = 0;
						}

						kmer_n = (bufMK2[char_tmp + kmer_num_filter - 1] >> 32);

						kmer_ins[kmer_i] = in_j;
						kmer_outs[kmer_i] = out_j;
						kmer_indexs[kmer_i] = char_tmp;

						if((kmer_n > filter_min) && (kmer_n < filter_max))
							kmer_types[kmer_i] = 0;
						else	kmer_types[kmer_i] = 1;

						kmer_i++;
					}

					for(kmer_k = 0; kmer_k < kmer_i; kmer_k++)
					{
						kmer_pre_as[kmer_num_new - 1] = 0;
						memcpy(kmer_pre_as, bufMK2 + kmer_indexs[kmer_k], kmer_num_tmp);

						in_j = kmer_ins[kmer_k];
						out_j = kmer_outs[kmer_k];
						type_c = kmer_types[kmer_k];

						if(type_c)
						{
							if(in_j)//1 shift right
							{
								add_edge = (kmer_pre_as[kmer_num - 1] >> (keyMove + 2)) & 0X3;

								if(kmer_num > 1)
								{
									shc_pre = 0;
									for(tra_i = 0; tra_i < kmer_num; tra_i++)
									{
										shc = kmer_pre_as[tra_i] & 0X3;
										kmer_pre_as[tra_i] >>= 2;
										kmer_pre_as[tra_i] |= firstls[shc_pre];

										shc_pre = shc;
									}
									kmer_pre_as[kmer_num - 1] &= KMask;

									kmer_pre_as[0] <<= (tempMove_new + 4);
									kmer_pre_as[0] >>= 4;

									kmer_pre_as[kmer_num_new - 1] |= add_edge;

									for(aph_i = 0; aph_i < 4; aph_i++)
									{
										if(!((in_j >> aph_i) & 0X1))	continue;

										kmer_pre_as[0] <<= 4;
										kmer_pre_as[0] >>= 4;
										kmer_pre_as[0] |= firstls[aph_i + 4];
										kmer_pre_as[0] |= type_one;//1

										memcpy(write_filterout[num] + (w_f_o * kmer_num_new), kmer_pre_as, kmer_num_tmp_new);
										w_f_o++;
									}
								}
								else
								{
									kmer_pre_as[0] >>= 4;
									kmer_pre_as[0] <<= 2;

									kmer_pre_as[kmer_num_new - 1] |= add_edge;

									for(aph_i = 0; aph_i < 4; aph_i++)
									{
										if(!((in_j >> aph_i) & 0X1))	continue;

										kmer_pre_as[0] <<= 4;
										kmer_pre_as[0] >>= 4;
										kmer_pre_as[0] |= firstls[aph_i + 4];
										kmer_pre_as[0] |= type_one;//1

										memcpy(write_filterout[num] + (w_f_o * kmer_num_new), kmer_pre_as, kmer_num_tmp_new);
										w_f_o++;
									}
								}

								kmer_pre_as[kmer_num_new - 1] = 0;
								memcpy(kmer_pre_as, bufMK2 + kmer_indexs[kmer_k], kmer_num_tmp);
							}

							if(out_j)//2 shift left
							{
								add_edge = ((kmer_pre_as[0] >> first_bp_move) & 0X3);

								if(kmer_num > 1)
								{
									shc_pre = 0;//firstls[aph_i + 8]
									kmer_pre_as[kmer_num - 1] &= KMask;
									for(tra_i = kmer_num - 1; tra_i >= 0; tra_i--)
									{
										shc = (kmer_pre_as[tra_i] >> 62) & 0X3;
										kmer_pre_as[tra_i] <<= 2;
										kmer_pre_as[tra_i] |= shc_pre;

										shc_pre = shc;
									}

									kmer_pre_as[0] <<= (tempMove_new + 2);
									kmer_pre_as[0] >>= 2;
									kmer_pre_as[0] |= type_two;//2



									for(aph_i = 0; aph_i < 4; aph_i++)
									{
										if(!((out_j >> aph_i) & 0X1))	continue;

										kmer_pre_as[kmer_num - 1] >>= (keyMove + 4);
										kmer_pre_as[kmer_num - 1] <<= 2;
										kmer_pre_as[kmer_num - 1] |= aph_i;
										kmer_pre_as[kmer_num - 1] <<= (keyMove + 2);

										kmer_pre_as[kmer_num_new - 1] |= add_edge;

										memcpy(write_filterout[num] + (w_f_o * kmer_num_new), kmer_pre_as, kmer_num_tmp_new);
										w_f_o++;
									}
								}
								else
								{
									kmer_pre_as[0] <<= 4;
									kmer_pre_as[0] >>= 2;
									kmer_pre_as[0] |= type_two;//2

									kmer_pre_as[kmer_num_new - 1] |= add_edge;

									for(aph_i = 0; aph_i < 4; aph_i++)
									{
										if(!((out_j >> aph_i) & 0X1))	continue;

										kmer_pre_as[0] >>= 4;
										kmer_pre_as[0] <<= 2;
										kmer_pre_as[0] |= aph_i;
										kmer_pre_as[0] <<= 2;

										memcpy(write_filterout[num] + (w_f_o * kmer_num_new), kmer_pre_as, kmer_num_tmp_new);
										w_f_o++;
									}

								}
							}
						}
						else
						{
							if(kmer_num > 1)
							{
								kmer_pre_as[0] <<= (tempMove_new + 2);
								kmer_pre_as[0] >>= 2;
								kmer_pre_as[kmer_num - 1] &= KMask;
							}
							else
							{
								kmer_pre_as[0] <<= 2;
								kmer_pre_as[0] >>= 2;
								kmer_pre_as[0] >>= 2;
								kmer_pre_as[0] <<= 2;
							}

							add_edge = ((in_j << 4) | out_j);
							kmer_pre_as[kmer_num_new - 1] |= add_edge;
							memcpy(write_filterin[num] + (w_f_i * kmer_num_new), kmer_pre_as, kmer_num_tmp_new);

							w_f_i++;
						}
					}

					in_j = 0;
					out_j = 0;
					kmer_n = 0;
					kmer_i = 0;

					anchor_km_pre_i = i;
				}

				memcpy(anchor_kms_pre, anchor_kms, kmer_num_tmp);
			}

			if(num != totalNum - 1)
			{
				in_j = 0;
				out_j = 0;
				kmer_n = 0;

				if(i - anchor_km_pre_i > 1)
				{
					pre_j = anchor_km_pre_i * kmer_num_filter;
					for(j = anchor_km_pre_i; j < i; j++)
					{
						char_tmp = j * kmer_num_filter;

						anchor = ((bufMK2[char_tmp + kmer_num - 1] << leftmove32) & KMask);

						if((anchor != kmer_pre) && (j > anchor_km_pre_i))
						{
							kmer_ins[kmer_i] = in_j;
							kmer_outs[kmer_i] = out_j;
							kmer_indexs[kmer_i] = pre_j;

							if((kmer_n > filter_min) && (kmer_n < filter_max))
								kmer_types[kmer_i] = 0;
							else	kmer_types[kmer_i] = 1;

							kmer_i++;

							in_j = 0;
							out_j = 0;
							kmer_n = 0;

							pre_j = char_tmp;//j
						}

						input = ((bufMK2[char_tmp] >> tempMove_r_new3) & 0X3);
						output = ((bufMK2[char_tmp + kmer_num - 1] >> keyMove) & 0X3);

						fl_f = (bufMK2[char_tmp + kmer_num_new - 1] & 0X3);

						if(fl_f == 1)
						{
							out_j |= (one_tmp << output);
						}
						else if(fl_f == 2)
						{
							in_j |= (one_tmp << input);
						}
						else if(fl_f == 0)
						{
							in_j |= (one_tmp << input);
							out_j |= (one_tmp << output);
						}

						kmer_pre = anchor;

						kmer_n += (bufMK2[char_tmp + kmer_num_filter - 1] >> 32);
					}

					kmer_ins[kmer_i] = in_j;
					kmer_outs[kmer_i] = out_j;
					kmer_indexs[kmer_i] = pre_j;

					if((kmer_n > filter_min) && (kmer_n < filter_max))
						kmer_types[kmer_i] = 0;
					else	kmer_types[kmer_i] = 1;

					kmer_i++;
				}
				else
				{
					//multi_kmer = 0;

					char_tmp = anchor_km_pre_i * kmer_num_filter;

					memcpy(kmer_pre_as, bufMK2 + char_tmp, kmer_num_tmp_new);

					//type_s = 0;

					fl_f = (kmer_pre_as[kmer_num_new - 1] & 0X3);

					if(fl_f == 1)
					{
						//type_s = 2;

						in_j = 0;
						output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
						out_j = (one_tmp << output);
					}
					else if(fl_f == 2)
					{
						//type_s = 1;

						input = ((kmer_pre_as[0] >> tempMove_r_new3) & 0X3);
						in_j = (one_tmp << input);
						out_j = 0;
					}
					else if(fl_f == 0)
					{
						input = ((kmer_pre_as[0] >> tempMove_r_new3) & 0X3);
						in_j = (one_tmp << input);
						output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
						out_j = (one_tmp << output);
					}
					else
					{
						in_j = 0;
						out_j = 0;
					}

					kmer_n = (bufMK2[char_tmp + kmer_num_filter - 1] >> 32);

					kmer_ins[kmer_i] = in_j;
					kmer_outs[kmer_i] = out_j;
					kmer_indexs[kmer_i] = char_tmp;

					if((kmer_n > filter_min) && (kmer_n < filter_max))
						kmer_types[kmer_i] = 0;
					else	kmer_types[kmer_i] = 1;

					kmer_i++;
				}

				for(kmer_k = 0; kmer_k < kmer_i; kmer_k++)
				{
					kmer_pre_as[kmer_num_new - 1] = 0;
					memcpy(kmer_pre_as, bufMK2 + kmer_indexs[kmer_k], kmer_num_tmp);

					in_j = kmer_ins[kmer_k];
					out_j = kmer_outs[kmer_k];
					type_c = kmer_types[kmer_k];

					if(type_c)
					{
						if(in_j)//1 shift right
						{
							add_edge = (kmer_pre_as[kmer_num - 1] >> (keyMove + 2)) & 0X3;

							if(kmer_num > 1)
							{
								shc_pre = 0;
								for(tra_i = 0; tra_i < kmer_num; tra_i++)
								{
									shc = kmer_pre_as[tra_i] & 0X3;
									kmer_pre_as[tra_i] >>= 2;
									kmer_pre_as[tra_i] |= firstls[shc_pre];

									shc_pre = shc;
								}
								kmer_pre_as[kmer_num - 1] &= KMask;

								kmer_pre_as[0] <<= (tempMove_new + 4);
								kmer_pre_as[0] >>= 4;

								kmer_pre_as[kmer_num_new - 1] |= add_edge;

								for(aph_i = 0; aph_i < 4; aph_i++)
								{
									if(!((in_j >> aph_i) & 0X1))	continue;

									kmer_pre_as[0] <<= 4;
									kmer_pre_as[0] >>= 4;
									kmer_pre_as[0] |= firstls[aph_i + 4];
									kmer_pre_as[0] |= type_one;//1

									memcpy(write_filterout[num] + (w_f_o * kmer_num_new), kmer_pre_as, kmer_num_tmp_new);
									w_f_o++;
								}
							}
							else
							{
								kmer_pre_as[0] >>= 4;
								kmer_pre_as[0] <<= 2;

								kmer_pre_as[kmer_num_new - 1] |= add_edge;

								for(aph_i = 0; aph_i < 4; aph_i++)
								{
									if(!((in_j >> aph_i) & 0X1))	continue;

									kmer_pre_as[0] <<= 4;
									kmer_pre_as[0] >>= 4;
									kmer_pre_as[0] |= firstls[aph_i + 4];
									kmer_pre_as[0] |= type_one;//1

									memcpy(write_filterout[num] + (w_f_o * kmer_num_new), kmer_pre_as, kmer_num_tmp_new);
									w_f_o++;
								}
							}

							kmer_pre_as[kmer_num_new - 1] = 0;
							memcpy(kmer_pre_as, bufMK2 + kmer_indexs[kmer_k], kmer_num_tmp);
						}

						if(out_j)//2 shift left
						{
							add_edge = ((kmer_pre_as[0] >> first_bp_move) & 0X3);

							if(kmer_num > 1)
							{
								shc_pre = 0;//firstls[aph_i + 8]
								kmer_pre_as[kmer_num - 1] &= KMask;
								for(tra_i = kmer_num - 1; tra_i >= 0; tra_i--)
								{
									shc = (kmer_pre_as[tra_i] >> 62) & 0X3;
									kmer_pre_as[tra_i] <<= 2;
									kmer_pre_as[tra_i] |= shc_pre;

									shc_pre = shc;
								}

								kmer_pre_as[0] <<= (tempMove_new + 2);
								kmer_pre_as[0] >>= 2;
								kmer_pre_as[0] |= type_two;//2



								for(aph_i = 0; aph_i < 4; aph_i++)
								{
									if(!((out_j >> aph_i) & 0X1))	continue;

									kmer_pre_as[kmer_num - 1] >>= (keyMove + 4);
									kmer_pre_as[kmer_num - 1] <<= 2;
									kmer_pre_as[kmer_num - 1] |= aph_i;
									kmer_pre_as[kmer_num - 1] <<= (keyMove + 2);

									kmer_pre_as[kmer_num_new - 1] |= add_edge;

									memcpy(write_filterout[num] + (w_f_o * kmer_num_new), kmer_pre_as, kmer_num_tmp_new);
									w_f_o++;
								}
							}
							else
							{
								kmer_pre_as[0] <<= 4;
								kmer_pre_as[0] >>= 2;
								kmer_pre_as[0] |= type_two;//2

								kmer_pre_as[kmer_num_new - 1] |= add_edge;

								for(aph_i = 0; aph_i < 4; aph_i++)
								{
									if(!((out_j >> aph_i) & 0X1))	continue;

									kmer_pre_as[0] >>= 4;
									kmer_pre_as[0] <<= 2;
									kmer_pre_as[0] |= aph_i;
									kmer_pre_as[0] <<= 2;

									memcpy(write_filterout[num] + (w_f_o * kmer_num_new), kmer_pre_as, kmer_num_tmp_new);
									w_f_o++;
								}

							}
						}
					}
					else
					{
						if(kmer_num > 1)
						{
							kmer_pre_as[0] <<= (tempMove_new + 2);
							kmer_pre_as[0] >>= 2;
							kmer_pre_as[kmer_num - 1] &= KMask;
						}
						else
						{
							kmer_pre_as[0] <<= 2;
							kmer_pre_as[0] >>= 2;
							kmer_pre_as[0] >>= 2;
							kmer_pre_as[0] <<= 2;
						}

						add_edge = ((in_j << 4) | out_j);
						kmer_pre_as[kmer_num_new - 1] |= add_edge;
						memcpy(write_filterin[num] + (w_f_i * kmer_num_new), kmer_pre_as, kmer_num_tmp_new);

						w_f_i++;
					}
				}
			}
			else	last_block_num = anchor_km_pre_i;
		}
		else
		{
			for(i = start; i < end; i++)
			{
				start_get = (i * kmer_num_new);

				getkmerno_pss(anchor_kms, bufMK2, start_get)

				if(memcmp(anchor_kms, anchor_kms_pre, kmer_num_tmp) && (i > start))
				{
					if(i - anchor_km_pre_i > 1)
					{
						pre_j = anchor_km_pre_i;
						for(j = anchor_km_pre_i; j < i; j++)
						{
							anchor = ((bufMK2[(j * kmer_num_new) + kmer_num - 1] << leftmove32) & KMask);

							if((anchor != kmer_pre) && (j > anchor_km_pre_i))
							{
								if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
								else	in_linear = 0;

								if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
								else	out_linear = 0;

								if(in_linear && out_linear)			//linear
									type_f = 0;
								else if(in_linear && !out_linear)	//single-in multi-out
									type_f = 1;
								else if(!in_linear && out_linear)	//multi-in single-out
									type_f = 2;
								else								//multi-in multi-out
									type_f = 3;

								type_s = type_f;
								if((type_f == 0) || (type_f == 1))
								{
									pre_flag = 0;
									pre_tmp = __builtin_ffs(in_j) - 1;
									for(pre_k = 0; pre_k < pre_j - anchor_km_pre_i; pre_k++)
										if(pres[pre_k] == pre_tmp)
										{
											pre_flag = 1;
											break;
										}

									for(pre_k = j - anchor_km_pre_i; (pre_k < i - anchor_km_pre_i) && (!pre_flag); pre_k++)
										if(pres[pre_k] == pre_tmp)
										{
											pre_flag = 1;
											break;
										}

									if(pre_flag)
									{
										if(type_f == 0)	type_s = 2;
										else	type_s = 3;
									}
								}

								kmer_tmp = (kmer_pre_as[0] >> first_bp_move);

								alph_tmp = (kmer_tmp & 0X3);
								alphabet_num_new[num][alph_tmp]++;

								type_s |= (input_types[in_j] << 2);
								type_s |= (out_j << 4);

								kmer_pre_as[0] <<= (tempMove_new + 2);
								kmer_pre_as[kmer_num - 1] &= KMask;
								memcpy(write_kmermerge[num] + (w_m_p_t * kmer_num), kmer_pre_as, kmer_num_tmp);

								write_kmermerge_t[num][w_m_p_t] = type_s;
#ifdef	GFA_COM
								edges_t[num][w_m_p_t] = ((in_j << 4) | out_j);
#endif
								w_m_p_t++;

								in_j = 0;
								out_j = 0;

								pre_j = j;
							}

							char_tmp = j * kmer_num_new;

							input = ((bufMK2[char_tmp] >> tempMove_r_new3) & 0X3);
							output = ((bufMK2[char_tmp + kmer_num - 1] >> keyMove) & 0X3);
#ifdef	FIRST_LAST_KMER
							fl_f = (bufMK2[char_tmp + kmer_num_new - 1] & 0X3);

							if(fl_f == 1)
							{
								out_j |= (one_tmp << output);
							}
							else if(fl_f == 2)
							{
								in_j |= (one_tmp << input);
							}
							else
							{
								in_j |= (one_tmp << input);
								out_j |= (one_tmp << output);
							}
#else
							in_j |= (one_tmp << input);
							out_j |= (one_tmp << output);
#endif

							kmer_pre = anchor;

							memcpy(kmer_pre_as, bufMK2 + (j * kmer_num_new), kmer_num_tmp);
						}

						if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
						else	in_linear = 0;

						if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
						else	out_linear = 0;

						if(in_linear && out_linear)			//linear
							type_f = 0;
						else if(in_linear && !out_linear)	//single-in multi-out
							type_f = 1;
						else if(!in_linear && out_linear)	//multi-in single-out
							type_f = 2;
						else								//multi-in multi-out
							type_f = 3;

						type_s = type_f;
						if((type_f == 0) || (type_f == 1))
						{
							pre_flag = 0;
							pre_tmp = __builtin_ffs(in_j) - 1;
							for(pre_k = 0; pre_k < pre_j - anchor_km_pre_i; pre_k++)
								if(pres[pre_k] == pre_tmp)
								{
									pre_flag = 1;
									break;
								}

							for(pre_k = j - anchor_km_pre_i; (pre_k < i - anchor_km_pre_i) && (!pre_flag); pre_k++)
								if(pres[pre_k] == pre_tmp)
								{
									pre_flag = 1;
									break;
								}

							if(pre_flag)
							{
								if(type_f == 0)	type_s = 2;
								else	type_s = 3;
							}
						}
					}
					else
					{
//fprintf(stderr, "1.3\n");
						//kmer_pre_a = bufMK2[anchor_km_pre_i];
						memcpy(kmer_pre_as, bufMK2 + (anchor_km_pre_i * kmer_num_new), kmer_num_tmp_new);

						type_s = 0;

#ifdef	FIRST_LAST_KMER
						fl_f = (kmer_pre_as[kmer_num_new - 1] & 0X3);

						input = ((kmer_pre_as[0] >> tempMove_r_new3) & 0X3);
						in_j = (one_tmp << input);

						if(fl_f == 1)
						{
							type_s = 2;
							output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
							out_j = (one_tmp << output);
						}
						else if(fl_f == 2)
						{
							type_s = 1;
							out_j = 0;
						}
						else
						{
							output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
							out_j = (one_tmp << output);
						}
#else
						input = ((kmer_pre_as[0] >> tempMove_r_new3) & 0X3);
						in_j = (one_tmp << input);

						output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
						out_j = (one_tmp << output);
#endif
					}

					kmer_tmp = (kmer_pre_as[0] >> first_bp_move);
					alph_tmp = (kmer_tmp & 0X3);

					alphabet_num_new[num][alph_tmp]++;

					type_s |= (input_types[in_j] << 2);
					type_s |= (out_j << 4);

					kmer_pre_as[0] <<= (tempMove_new + 2);
					kmer_pre_as[kmer_num - 1] &= KMask;
					memcpy(write_kmermerge[num] + w_m_p_t * kmer_num, kmer_pre_as, kmer_num_tmp);

					write_kmermerge_t[num][w_m_p_t] = type_s;
#ifdef	GFA_COM
					edges_t[num][w_m_p_t] = ((in_j << 4) | out_j);
#endif
					w_m_p_t++;

					in_j = 0;
					out_j = 0;

					anchor_km_pre_i = i;
					pres_i = 0;
				}

				memcpy(anchor_kms_pre, anchor_kms, kmer_num_tmp);

				if((bufMK2[start_get + kmer_num_new - 1] & 0X3) != 1)
				{
					pres[pres_i] = ((bufMK2[start_get] >> tempMove_r_new3) & 0X3);
				}
				else
				{
					pres[pres_i] = 4;
				}
				pres_i++;

			}

			if(num != totalNum - 1)
			{

				if(i - anchor_km_pre_i > 1)
				{
					pre_j = anchor_km_pre_i;
					for(j = anchor_km_pre_i; j < i; j++)
					{
						anchor = ((bufMK2[(j * kmer_num_new) + kmer_num - 1] << leftmove32) & KMask);

						if((anchor != kmer_pre) && (j > anchor_km_pre_i))
						{
							if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
							else	in_linear = 0;

							if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
							else	out_linear = 0;

							if(in_linear && out_linear)			//linear
								type_f = 0;
							else if(in_linear && !out_linear)	//single-in multi-out
								type_f = 1;
							else if(!in_linear && out_linear)	//multi-in single-out
								type_f = 2;
							else								//multi-in multi-out
								type_f = 3;

							type_s = type_f;
							if((type_f == 0) || (type_f == 1))
							{
								pre_flag = 0;
								pre_tmp = __builtin_ffs(in_j) - 1;
								for(pre_k = 0; pre_k < pre_j - anchor_km_pre_i; pre_k++)
									if(pres[pre_k] == pre_tmp)
									{
										pre_flag = 1;
										break;
									}

								for(pre_k = j - anchor_km_pre_i; (pre_k < i - anchor_km_pre_i) && (!pre_flag); pre_k++)
									if(pres[pre_k] == pre_tmp)
									{
										pre_flag = 1;
										break;
									}

								//if(input_km_cnt[__builtin_ffs(in_j) - 1] > 1)
								if(pre_flag)
								{
									if(type_f == 0)	type_s = 2;
									else	type_s = 3;
								}
							}

							kmer_tmp = (kmer_pre_as[0] >> first_bp_move);

							alph_tmp = (kmer_tmp & 0X3);
							alphabet_num_new[num][alph_tmp]++;

							type_s |= (input_types[in_j] << 2);
							type_s |= (out_j << 4);

							kmer_pre_as[0] <<= (tempMove_new + 2);
							kmer_pre_as[kmer_num - 1] &= KMask;
							memcpy(write_kmermerge[num] + (w_m_p_t * kmer_num), kmer_pre_as, kmer_num_tmp);

							write_kmermerge_t[num][w_m_p_t] = type_s;
#ifdef	GFA_COM
							edges_t[num][w_m_p_t] = ((in_j << 4) | out_j);
#endif
							w_m_p_t++;

							in_j = 0;
							out_j = 0;

							pre_j = j;
						}

						char_tmp = j * kmer_num_new;
						input = ((bufMK2[char_tmp] >> tempMove_r_new3) & 0X3);
						output = ((bufMK2[char_tmp + kmer_num - 1] >> keyMove) & 0X3);

#ifdef	FIRST_LAST_KMER

						fl_f = (bufMK2[char_tmp + kmer_num_new - 1] & 0X3);

						if(fl_f == 1)
						{
							out_j |= (one_tmp << output);
						}
						else if(fl_f == 2)
						{
							in_j |= (one_tmp << input);
						}
						else
						{
							in_j |= (one_tmp << input);
							out_j |= (one_tmp << output);
						}

#else
						in_j |= (one_tmp << input);
						out_j |= (one_tmp << output);
#endif

						kmer_pre = anchor;

						memcpy(kmer_pre_as, bufMK2 + (j * kmer_num_new), kmer_num_tmp);
					}

					if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
					else	in_linear = 0;

					if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
					else	out_linear = 0;

					if(in_linear && out_linear)			//linear
						type_f = 0;
					else if(in_linear && !out_linear)	//single-in multi-out
						type_f = 1;
					else if(!in_linear && out_linear)	//multi-in single-out
						type_f = 2;
					else								//multi-in multi-out
						type_f = 3;

					type_s = type_f;
					if((type_f == 0) || (type_f == 1))
					{
						pre_flag = 0;
						pre_tmp = __builtin_ffs(in_j) - 1;
						for(pre_k = 0; pre_k < pre_j - anchor_km_pre_i; pre_k++)
							if(pres[pre_k] == pre_tmp)
							{
								pre_flag = 1;
								break;
							}

						for(pre_k = j - anchor_km_pre_i; (pre_k < i - anchor_km_pre_i) && (!pre_flag); pre_k++)
							if(pres[pre_k] == pre_tmp)
							{
								pre_flag = 1;
								break;
							}

						if(pre_flag)
						{
							if(type_f == 0)	type_s = 2;
							else	type_s = 3;
						}
					}
				}
				else
				{
					memcpy(kmer_pre_as, bufMK2 + (anchor_km_pre_i * kmer_num_new), kmer_num_tmp_new);

					type_s = 0;
#ifdef	FIRST_LAST_KMER
					input = ((kmer_pre_as[0] >> tempMove_r_new3) & 0X3);
					in_j = (one_tmp << input);

					fl_f = (kmer_pre_as[kmer_num_new - 1] & 0X3);
					if(fl_f == 1)
					{
						type_s = 2;
						output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
						out_j = (one_tmp << output);
					}
					else if(fl_f == 2)
					{
						type_s = 1;
						out_j = 0;
					}
					else
					{
						output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
						out_j = (one_tmp << output);
					}
#else
					input = ((kmer_pre_as[0] >> tempMove_r_new3) & 0X3);
					in_j = (one_tmp << input);

					output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
					out_j = (one_tmp << output);
#endif
				}

				kmer_tmp = (kmer_pre_as[0] >> first_bp_move);
				alph_tmp = (kmer_tmp & 0X3);

				alphabet_num_new[num][alph_tmp]++;

				type_s |= (input_types[in_j] << 2);
				type_s |= (out_j << 4);

				kmer_pre_as[0] <<= (tempMove_new + 2);
				kmer_pre_as[kmer_num - 1] &= KMask;
				memcpy(write_kmermerge[num] + w_m_p_t * kmer_num, kmer_pre_as, kmer_num_tmp);

				write_kmermerge_t[num][w_m_p_t] = type_s;
#ifdef	GFA_COM
				edges_t[num][w_m_p_t] = ((in_j << 4) | out_j);
#endif
				w_m_p_t++;

			}
			else	last_block_num = anchor_km_pre_i;


		}

	}

	if(flag_filter)
	{
		w_f_is[num] = w_f_i;
		w_f_os[num] = w_f_o;
	}
	w_m_p[num] = w_m_p_t;

	return (void *)NULL;
}

uint8_t filter_out_sort_tag(uint64_t THREAD_NUM, int block_numf, char* tmp_route)
{
	int i = 0;
	int check = 0;
	uint64_t i_in = 0;
	uint64_t totalKmerNum = 0;
	uint64_t locateNum = 0;
	uint64_t locate = 0;
	uint64_t hash_c = 0;
	uint64_t hash_dis = 0;
	char nameK2[40];
	char cNum[10];

	pthread_t* myThread = NULL;
	thread_data* tt = NULL;

	hashKmer_disf = (uint64_t** )calloc(THREAD_NUM, sizeof(uint64_t* ));

	countKmer_disf = (uint64_t** )calloc(THREAD_NUM, sizeof(uint64_t* ));

	maskBKT=((int)1<<(BUCKET_LENGTH<<1))-1;

	countKmerf=(uint64_t *)calloc(BUCKET_CAPACITY+1,sizeof(uint64_t));

	myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
	tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));
	for(i=0; i<THREAD_NUM; i++)
	{
		tt[i].num = i;
		tt[i].thread = THREAD_NUM;

		check=pthread_create(&myThread[i], NULL, multiCount_filter, tt + i);
		if(check)
		{
			fprintf(stderr,"threadNum: %d, Error - pthread_create() return code: %d\n",i,check);
			exit(EXIT_FAILURE);
		}
	}
	for(i=0; i<THREAD_NUM; i++)
	{
		pthread_join(myThread[i], NULL);
	}
	free(myThread);
	free(tt);


	fprintf(stderr, "Block %d kmer-count filter and counting finish\n", block_numf + 1);

	myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
	tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));
	for(i=0; i<THREAD_NUM; i++)
	{
		tt[i].thread = THREAD_NUM;
		tt[i].num = i;

		check=pthread_create(&myThread[i], NULL, multiDistri_filter, tt + i);

		if(check)
		{
			fprintf(stderr,"threadNum: %d, Error - pthread_create() return code: %d\n",i,check);
			exit(EXIT_FAILURE);
		}
	}
	for(i=0; i<THREAD_NUM; i++)
	{
		pthread_join( myThread[i], NULL);
	}
	free(myThread);
	free(tt);

	fprintf(stderr, "Block %d kmer-count filter distributing finish\n", block_numf + 1);

	countKmerf[0] = 0;

	for(i=1; i<BUCKET_CAPACITY+1; i++)
	{
		for(i_in=0; i_in<THREAD_NUM; i_in++)
		{
			countKmerf[i] += (countKmer_disf[i_in][i] - countKmer_disf[i_in][i - 1]);
		}
		countKmerf[i] += countKmerf[i - 1];
	}
	totalKmerNum = countKmerf[BUCKET_CAPACITY];
	hashKmerf = (uint64_t* )calloc(totalKmerNum, kmer_num_filtermerge);

	uint64_t kmer_num_cp = kmer_num_fm;

	for(i=0; i<BUCKET_CAPACITY; i++)
	{
		for(i_in=0, hash_c = 0; i_in<THREAD_NUM; i_in++)
		{
			hash_dis = countKmer_disf[i_in][i+1] - countKmer_disf[i_in][i];

			if(hash_dis)
				memcpy(hashKmerf + (countKmerf[i] + hash_c) * kmer_num_cp, hashKmer_disf[i_in] + countKmer_disf[i_in][i] * kmer_num_cp, hash_dis * kmer_num_filtermerge);//(hash_dis * kmer_num) << 3

			hash_c += hash_dis;
		}
	}

	if(!totalKmerNum)	return 1;

	segCountf[0]=0;
	for(i=1; i<THREAD_NUM; i++)
	{
		locateNum = i*(totalKmerNum/THREAD_NUM);
		locate = BinarySearchf(locateNum,countKmerf,BUCKET_CAPACITY-1);
		segCountf[i] = locate;
	}

	myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
	tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));
	for(i=0; i<THREAD_NUM; i++)
	{
		tt[i].thread = THREAD_NUM;
		tt[i].num = i;
		check=pthread_create(&myThread[i], NULL, multiThreadSort_filter, tt + i);

		if(check)
		{
			fprintf(stderr,"threadNum: %d, Error - pthread_create() return code: %d\n",i,check);
			exit(EXIT_FAILURE);
		}
	}

	for(i=0; i<THREAD_NUM; i++)
	{
		pthread_join( myThread[i], NULL);
	}
	free(myThread);
	free(tt);

	fprintf(stderr, "Block %d kmer-count filter sorting finish\n", block_numf + 1);

	strcpy(nameK2, tmp_route);
	strcat(nameK2, "/kmerfilter.");
	sprintf(cNum,"%u", block_numf);
	strcat(nameK2, cNum);

	FILE* fpKmer_w = fopen(nameK2,"wb");

	fwrite(hashKmerf, kmer_num_filtermerge, totalKmerNum, fpKmer_w);

	fclose(fpKmer_w);
	free(hashKmerf);
	free(countKmerf);

	return 0;
}

uint64_t BinarySearchf(uint64_t mk, uint64_t* target, int64_t up)
{
	int64_t low=0,mid=0;
	while(low<=up)
	{
		mid=(low+up)>>1;
		if(mk<target[mid])   up=mid-1;
		else if(mk>target[mid])   low=mid+1;
		if(mk==target[mid])    return mid;
	}
	return low;
}
#ifdef	GFA_COM
void multiTagcheck_addkmer(uint64_t i, uint64_t anchor_km_pre_i, FILE* fp_merge, FILE* fp_merge_t, FILE* fp_edges_t)
#else
void multiTagcheck_addkmer(uint64_t i, uint64_t anchor_km_pre_i, FILE* fp_merge, FILE* fp_merge_t)
#endif
{
	uint8_t in_j = 0;
	uint8_t out_j = 0;
	uint8_t in_linear = 0;
	uint8_t out_linear = 0;
	uint8_t one_tmp = 1;
	uint8_t type_f = 0;
	uint8_t type_s = 0;
#ifdef	GFA_COM
	uint8_t type_e = 0;
#endif
	uint8_t input = 0;
	uint8_t output = 0;
	uint8_t pre_k = 0;
	uint8_t pre_tmp = 0;
	uint8_t pre_flag = 0;
	uint8_t kmer_tmp = 0;
	uint8_t pres[128];//64
	uint8_t pres_i = 0;
	uint8_t alph_tmp = 0;
	uint8_t fl_f = 0;
	uint8_t kmer_num_tmp = (kmer_num << 3);

	uint64_t j = 0;
	uint64_t pre_j = 0;
	uint64_t start_get = 0;
	uint64_t anchor = 0;
	uint64_t anchor_a = 0;
	uint64_t kmer_pre = 0;
	uint64_t kmer_pre_a = 0;
	uint64_t write_tmp = 0;

	if(OneUnit)
	{
		for(j = anchor_km_pre_i; j < i; j++)
		{
			kmer_tmp = (bufMK2[j] >> 60);
			if((kmer_tmp >> 2) != 1)
			{
				pres[pres_i] = (kmer_tmp & 0X3);
			}
			else
			{
				pres[pres_i] = 4;
			}
			pres_i++;
		}
	}
	else
	{
		for(j = anchor_km_pre_i; j < i; j++)
		{
			start_get = (j * kmer_num_new);
			if((bufMK2[start_get + kmer_num_new - 1] & 0X3) != 1)
			{
				pres[pres_i] = ((bufMK2[start_get] >> tempMove_r_new3) & 0X3);
			}
			else
			{
				pres[pres_i] = 4;
			}
			pres_i++;
		}
	}

	printf("BON %u %u %"PRId64" %"PRId64"\n", pres_i, pres[pres_i], anchor_km_pre_i, i);

	if(OneUnit)
	{
		uint64_t eliminate=~(ELIMINATE);

		if(i - anchor_km_pre_i > 1)
		{
			pre_j = anchor_km_pre_i;
			for(j = anchor_km_pre_i; j < i; j++)
			{
				anchor_a = bufMK2[j];
#ifdef	FIRST_LAST_KMER
				anchor = ((anchor_a<<4)&eliminate);
#else
				anchor = ((anchor_a<<2)&eliminate);
#endif

				if((anchor != kmer_pre) && (j > anchor_km_pre_i))
				{
					if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
					else	in_linear = 0;

					if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
					else	out_linear = 0;

					if(in_linear && out_linear)			//linear
						type_f = 0;
					else if(in_linear && !out_linear)	//single-in multi-out
						type_f = 1;
					else if(!in_linear && out_linear)	//multi-in single-out
						type_f = 2;
					else								//multi-in multi-out
						type_f = 3;

					type_s = type_f;
					if((type_f == 0) || (type_f == 1))
					{
						pre_flag = 0;
						pre_tmp = __builtin_ffs(in_j) - 1;
						for(pre_k = 0; pre_k < pre_j - anchor_km_pre_i; pre_k++)
							if(pres[pre_k] == pre_tmp)
							{
								pre_flag = 1;
								break;
							}

						for(pre_k = j - anchor_km_pre_i; (pre_k < i - anchor_km_pre_i) && (!pre_flag); pre_k++)
							if(pres[pre_k] == pre_tmp)
							{
								pre_flag = 1;
								break;
							}

						if(pre_flag)
						{
							if(type_f == 0)	type_s = 2;
							else	type_s = 3;
						}
					}

					kmer_tmp = ((kmer_pre_a >> 58) & 0Xf);

					write_tmp = ((kmer_pre_a << 4) & eliminate);
					fwrite(&write_tmp, 8, 1, fp_merge);

					alph_tmp = (kmer_tmp & 0X3);
					alphabet_num_new[0][alph_tmp]++;

					type_s |= (input_types[in_j] << 2);
					type_s |= (out_j << 4);

					fwrite(&type_s, 1, 1, fp_merge_t);
#ifdef	GFA_COM
					type_e = ((in_j << 4)|out_j);
					fwrite(&type_e, 1, 1, fp_edges_t);
#endif
					in_j = 0;
					out_j = 0;

					pre_j = j;
				}
#ifdef	FIRST_LAST_KMER
				input = ((anchor_a >> 60) & 0Xf);
				output = ((anchor_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);

				fl_f = (input >> 2);
				input &= 0X3;

				if(fl_f == 1)
				{
					out_j |= (one_tmp << output);
				}
				else if(fl_f == 2)
				{
					in_j |= (one_tmp << input);
				}
				else
				{
					in_j |= (one_tmp << input);
					out_j |= (one_tmp << output);
				}
#else
				input = ((anchor_a >> 62) & 0X3);
				output = ((anchor_a >> (64 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);

				in_j |= (one_tmp << input);
				out_j |= (one_tmp << output);
#endif

				kmer_pre = anchor;
				kmer_pre_a = anchor_a;
			}

			if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
			else	in_linear = 0;

			if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
			else	out_linear = 0;

			if(in_linear && out_linear)			//linear
				type_f = 0;
			else if(in_linear && !out_linear)	//single-in multi-out
				type_f = 1;
			else if(!in_linear && out_linear)	//multi-in single-out
				type_f = 2;
			else								//multi-in multi-out
				type_f = 3;

			type_s = type_f;
			if((type_f == 0) || (type_f == 1))
			{
				pre_flag = 0;
				pre_tmp = __builtin_ffs(in_j) - 1;
				for(pre_k = 0; pre_k < pre_j - anchor_km_pre_i; pre_k++)
					if(pres[pre_k] == pre_tmp)
					{
						pre_flag = 1;
						break;
					}

				for(pre_k = j - anchor_km_pre_i; (pre_k < i - anchor_km_pre_i) && (!pre_flag); pre_k++)
					if(pres[pre_k] == pre_tmp)
					{
						pre_flag = 1;
						break;
					}

				if(pre_flag)
				{
					if(type_f == 0)	type_s = 2;
					else	type_s = 3;
				}
			}
		}
		else
		{
			kmer_pre_a = bufMK2[anchor_km_pre_i];
			type_s = 0;

#ifdef	FIRST_LAST_KMER
			fl_f = (kmer_pre_a >> 60);

			input = fl_f & 0X3;
			in_j = (one_tmp << input);

			fl_f >>= 2;

			if(fl_f == 1)
			{
				type_s = 2;

				output = ((kmer_pre_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);
				out_j = (one_tmp << output);
			}
			else if(fl_f == 2)
			{
				type_s = 1;
				out_j = 0;
			}
			else
			{
				output = ((kmer_pre_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);
				out_j = (one_tmp << output);
			}

#else
			output = ((kmer_pre_a >> (64 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);

			input = (kmer_pre_a >> 62);
			in_j = (one_tmp << input);

			out_j = (one_tmp << output);
#endif

		}

		kmer_tmp = ((kmer_pre_a >> 58) & 0Xf);
		write_tmp = ((kmer_pre_a << 4) & eliminate);
		fwrite(&write_tmp, 8, 1, fp_merge);

		alph_tmp = (kmer_tmp & 0X3);
		alphabet_num_new[0][alph_tmp]++;

		type_s |= (input_types[in_j] << 2);
		type_s |= (out_j << 4);

		fwrite(&type_s, 1, 1, fp_merge_t);
#ifdef	GFA_COM
		type_e = ((in_j << 4)|out_j);
		fwrite(&type_e, 1, 1, fp_edges_t);
#endif

	}
	else
	{
		uint64_t char_tmp = 0;
		uint64_t anchor_kms[kmer_num];
		uint64_t anchor_kms_pre[kmer_num];
		uint64_t kmer_pre_as[kmer_num_new];

		uint8_t kmer_num_tmp_new = (kmer_num_new << 3);

		if(i - anchor_km_pre_i > 1)
		{
			pre_j = anchor_km_pre_i;
			for(j = anchor_km_pre_i; j < i; j++)
			{
				anchor = ((bufMK2[(j * kmer_num_new) + kmer_num - 1] << leftmove32) & KMask);

				if((anchor != kmer_pre) && (j > anchor_km_pre_i))
				{
					if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
					else	in_linear = 0;

					if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
					else	out_linear = 0;

					if(in_linear && out_linear)			//linear
						type_f = 0;
					else if(in_linear && !out_linear)	//single-in multi-out
						type_f = 1;
					else if(!in_linear && out_linear)	//multi-in single-out
						type_f = 2;
					else								//multi-in multi-out
						type_f = 3;

					type_s = type_f;
					if((type_f == 0) || (type_f == 1))
					{
						pre_flag = 0;
						pre_tmp = __builtin_ffs(in_j) - 1;
						for(pre_k = 0; pre_k < pre_j - anchor_km_pre_i; pre_k++)
							if(pres[pre_k] == pre_tmp)
							{
								pre_flag = 1;
								break;
							}

						for(pre_k = j - anchor_km_pre_i; (pre_k < i - anchor_km_pre_i) && (!pre_flag); pre_k++)
							if(pres[pre_k] == pre_tmp)
							{
								pre_flag = 1;
								break;
							}

						if(pre_flag)
						{
							if(type_f == 0)	type_s = 2;
							else	type_s = 3;
						}
					}

					kmer_tmp = (kmer_pre_as[0] >> first_bp_move);

					alph_tmp = (kmer_tmp & 0X3);
					alphabet_num_new[0][alph_tmp]++;

					type_s |= (input_types[in_j] << 2);
					type_s |= (out_j << 4);

					kmer_pre_as[0] <<= (tempMove_new + 2);
					kmer_pre_as[kmer_num - 1] &= KMask;

					fwrite(kmer_pre_as, kmer_num_tmp, 1, fp_merge);

					fwrite(&type_s, 1, 1, fp_merge_t);
#ifdef	GFA_COM
					type_e = ((in_j << 4)|out_j);
					fwrite(&type_e, 1, 1, fp_edges_t);
#endif
					in_j = 0;
					out_j = 0;

					pre_j = j;
				}

				char_tmp = j * kmer_num_new;
				input = ((bufMK2[char_tmp] >> tempMove_r_new3) & 0X3);
				output = ((bufMK2[char_tmp + kmer_num - 1] >> keyMove) & 0X3);

#ifdef	FIRST_LAST_KMER

				fl_f = (bufMK2[char_tmp + kmer_num_new - 1] & 0X3);

				if(fl_f == 1)
				{
					out_j |= (one_tmp << output);
				}
				else if(fl_f == 2)
				{
					in_j |= (one_tmp << input);
				}
				else
				{
					in_j |= (one_tmp << input);
					out_j |= (one_tmp << output);
				}

#else

				in_j |= (one_tmp << input);
				out_j |= (one_tmp << output);
#endif

				kmer_pre = anchor;

				memcpy(kmer_pre_as, bufMK2 + (j * kmer_num_new), kmer_num_tmp);
			}

			if((in_j == 8) || (in_j == 4) || (in_j == 2) || (in_j == 1))	in_linear = 1;
			else	in_linear = 0;

			if((out_j == 8) || (out_j == 4) || (out_j == 2) || (out_j == 1))	out_linear = 1;
			else	out_linear = 0;

			if(in_linear && out_linear)			//linear
				type_f = 0;
			else if(in_linear && !out_linear)	//single-in multi-out
				type_f = 1;
			else if(!in_linear && out_linear)	//multi-in single-out
				type_f = 2;
			else								//multi-in multi-out
				type_f = 3;

			type_s = type_f;
			if((type_f == 0) || (type_f == 1))
			{
				pre_flag = 0;
				pre_tmp = __builtin_ffs(in_j) - 1;
				for(pre_k = 0; pre_k < pre_j - anchor_km_pre_i; pre_k++)
					if(pres[pre_k] == pre_tmp)
					{
						pre_flag = 1;
						break;
					}

				for(pre_k = j - anchor_km_pre_i; (pre_k < i - anchor_km_pre_i) && (!pre_flag); pre_k++)
					if(pres[pre_k] == pre_tmp)
					{
						pre_flag = 1;
						break;
					}

				if(pre_flag)
				{
					if(type_f == 0)	type_s = 2;
					else	type_s = 3;
				}
			}
		}
		else
		{
			memcpy(kmer_pre_as, bufMK2 + (anchor_km_pre_i * kmer_num_new), kmer_num_tmp_new);

			type_s = 0;
#ifdef	FIRST_LAST_KMER
			input = ((kmer_pre_as[0] >> tempMove_r_new3) & 0X3);
			in_j = (one_tmp << input);

			fl_f = (kmer_pre_as[kmer_num_new - 1] & 0X3);
			if(fl_f == 1)
			{
				type_s = 2;
				output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
				out_j = (one_tmp << output);
			}
			else if(fl_f == 2)
			{
				type_s = 1;
				out_j = 0;
			}
			else
			{
				output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
				out_j = (one_tmp << output);
			}
#else
			input = ((kmer_pre_as[0] >> tempMove_r_new3) & 0X3);
			in_j = (one_tmp << input);

			output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
			out_j = (one_tmp << output);
#endif

		}

		kmer_tmp = (kmer_pre_as[0] >> first_bp_move);
		alph_tmp = (kmer_tmp & 0X3);

		alphabet_num_new[0][alph_tmp]++;

		type_s |= (input_types[in_j] << 2);
		type_s |= (out_j << 4);

		kmer_pre_as[0] <<= (tempMove_new + 2);
		kmer_pre_as[kmer_num - 1] &= KMask;

		fwrite(kmer_pre_as, kmer_num_tmp, 1, fp_merge);

		fwrite(&type_s, 1, 1, fp_merge_t);

#ifdef	GFA_COM
		type_e = ((in_j << 4)|out_j);
		fwrite(&type_e, 1, 1, fp_edges_t);
#endif

	}
}

void multiTagcheck_addkmer_filter(uint64_t i, uint64_t anchor_km_pre_i, uint64_t THREAD_NUM, uint64_t w_f_i_v, uint64_t w_f_o_v, FILE* fp_filterin)
{

	uint8_t in_j = 0;
	uint8_t out_j = 0;
	uint8_t in_linear = 0;
	uint8_t out_linear = 0;
	uint8_t one_tmp = 1;
	uint8_t type_f = 0;
	uint8_t input = 0;
	uint8_t output = 0;
	uint8_t pre_k = 0;
	uint8_t pre_tmp = 0;
	uint8_t pre_flag = 0;
	uint8_t kmer_tmp = 0;
	uint8_t fl_f = 0;
	uint8_t kmer_num_tmp = (kmer_num << 3);

	uint64_t j = 0;
	uint64_t pre_j = 0;
	uint64_t start_get = 0;
	uint64_t anchor = 0;
	uint64_t anchor_a = 0;
	uint64_t kmer_pre = 0;
	uint64_t kmer_pre_a = 0;
	uint64_t write_tmp = 0;
	uint8_t num = 0;
	uint64_t w_f_i = w_f_i_v;
	uint64_t w_f_o = w_f_o_v;

	num = THREAD_NUM - 1;

	if(OneUnit)
	{
		uint64_t eliminate=~(ELIMINATE);
		uint8_t kmer_ins[24];
		uint8_t kmer_outs[24];
		uint8_t kmer_types[24];
		uint64_t kmer_indexs[24];
		uint32_t kmer_n = 0;
		uint8_t kmer_i = 0;
		uint8_t kmer_k = 0;
		uint64_t start_get_kmer = 0;
		uint8_t type_c = 0;
		uint8_t add_edge = 0;
		uint8_t aph_i = 0;

		if(i - anchor_km_pre_i > 1)
		{
			pre_j = anchor_km_pre_i * kmer_num_filter;
			for(j = anchor_km_pre_i; j < i; j++)
			{
				start_get_kmer = j * kmer_num_filter;
				anchor_a = bufMK2[start_get_kmer];

				anchor = ((anchor_a<<4)&eliminate);

				if((anchor != kmer_pre) && (j > anchor_km_pre_i))
				{
					kmer_ins[kmer_i] = in_j;
					kmer_outs[kmer_i] = out_j;
					kmer_indexs[kmer_i] = pre_j;

					if((kmer_n > filter_min) && (kmer_n < filter_max))
						kmer_types[kmer_i] = 0;
					else	kmer_types[kmer_i] = 1;

					kmer_i++;

					in_j = 0;
					out_j = 0;
					pre_j = start_get_kmer;//j

					kmer_n = 0;
				}

				input = ((anchor_a >> 60) & 0Xf);
				output = ((anchor_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);

				fl_f = (input >> 2);
				input &= 0X3;

				if(fl_f == 1)
				{
					out_j |= (one_tmp << output);
				}
				else if(fl_f == 2)
				{
					in_j |= (one_tmp << input);
				}
				else if(fl_f == 0)
				{
					in_j |= (one_tmp << input);
					out_j |= (one_tmp << output);
				}

				kmer_pre = anchor;

				kmer_n += bufMK2[start_get_kmer + 1];
			}

			kmer_ins[kmer_i] = in_j;
			kmer_outs[kmer_i] = out_j;
			kmer_indexs[kmer_i] = pre_j;

			if((kmer_n > filter_min) && (kmer_n < filter_max))
				kmer_types[kmer_i] = 0;
			else	kmer_types[kmer_i] = 1;

			kmer_i++;
		}
		else
		{
			start_get_kmer = anchor_km_pre_i * kmer_num_filter;
			kmer_pre_a = bufMK2[start_get_kmer];

			fl_f = (kmer_pre_a >> 60);

			input = (fl_f & 0X3);

			fl_f >>= 2;

			if(fl_f == 1)
			{
				in_j = 0;
				output = ((kmer_pre_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);
				out_j = (one_tmp << output);
			}
			else if(fl_f == 2)
			{
				in_j = (one_tmp << input);
				out_j = 0;
			}
			else if(fl_f == 0)
			{
				in_j = (one_tmp << input);
				output = ((kmer_pre_a >> (62 - (KMER_LENGTH_PlusTwo << 1))) & 0X3);
				out_j = (one_tmp << output);
			}
			else
			{
				in_j = 0;
				out_j = 0;
			}

			kmer_n = bufMK2[start_get_kmer + 1];

			kmer_ins[kmer_i] = in_j;
			kmer_outs[kmer_i] = out_j;
			kmer_indexs[kmer_i] = start_get_kmer;

			if((kmer_n > filter_min) && (kmer_n < filter_max))
				kmer_types[kmer_i] = 0;
			else	kmer_types[kmer_i] = 1;

			kmer_i++;
		}

		for(kmer_k = 0; kmer_k < kmer_i; kmer_k++)
		{
			kmer_pre_a = bufMK2[kmer_indexs[kmer_k]];
			in_j = kmer_ins[kmer_k];
			out_j = kmer_outs[kmer_k];
			type_c = kmer_types[kmer_k];

			if(type_c)
			{
				if(in_j)
				{
					add_edge = (kmer_pre_a >> one_suf) & 0X3;//last char

					kmer_pre_a >>= one_sufp;
					kmer_pre_a <<= (one_suf + 6);
					kmer_pre_a >>= 4;

					for(aph_i = 0; aph_i < 4; aph_i++)
					{
						if(!((in_j >> aph_i) & 0X1))	continue;

						kmer_pre_a <<= 4;
						kmer_pre_a >>= 4;
						kmer_pre_a |= firstls[aph_i + 4];
						kmer_pre_a |= type_one;

						write_filterout[num][(w_f_o << 1)] = kmer_pre_a;

						write_filterout[num][(w_f_o << 1) + 1] = add_edge;//suf
						w_f_o++;
					}

					kmer_pre_a = bufMK2[kmer_indexs[kmer_k]];
				}

				if(out_j)
				{
					add_edge = (kmer_pre_a >> 58) & 0X3;//first char
					kmer_pre_a <<= 6;

					kmer_pre_a >>= 2;
					kmer_pre_a |= type_two;

					for(aph_i = 0; aph_i < 4; aph_i++)
					{
						if(!((out_j >> aph_i) & 0X1))	continue;

						kmer_pre_a >>= (one_sufp + 2);
						kmer_pre_a <<= 2;
						kmer_pre_a |= aph_i;
						kmer_pre_a <<= one_sufp;

						write_filterout[num][(w_f_o << 1)] = kmer_pre_a;

						write_filterout[num][(w_f_o << 1) + 1] = add_edge;//pre
						w_f_o++;
					}
				}
			}
			else
			{
				kmer_pre_a <<= 4;
				kmer_pre_a &= eliminate;
				kmer_pre_a >>= 2;

				add_edge = ((in_j << 4) | out_j);

				fwrite(&kmer_pre_a, 8, 1, fp_filterin);
				fwrite(&add_edge, 8, 1, fp_filterin);

			}
		}

	}
	else
	{
		uint64_t char_tmp = 0;
		uint64_t anchor_kms[kmer_num];
		uint64_t anchor_kms_pre[kmer_num];
		uint64_t kmer_pre_as[kmer_num_new];

		uint8_t kmer_num_tmp_new = (kmer_num_new << 3);

		uint8_t kmer_ins[24];
		uint8_t kmer_outs[24];
		uint8_t kmer_types[24];
		uint64_t kmer_indexs[24];
		uint32_t kmer_n = 0;
		uint8_t kmer_i = 0;
		uint8_t kmer_k = 0;
		uint8_t aph_i = 0;
		uint8_t type_c = 0;
		int16_t tra_i = 0;
		uint8_t shc = 0;
		uint8_t shc_pre = 0;
		uint8_t shtype = 0;
		uint8_t add_edge = 0;

		if(i - anchor_km_pre_i > 1)
		{
			pre_j = anchor_km_pre_i * kmer_num_filter;
			for(j = anchor_km_pre_i; j < i; j++)
			{
				char_tmp = j * kmer_num_filter;

				anchor = ((bufMK2[char_tmp + kmer_num - 1] << leftmove32) & KMask);

				if((anchor != kmer_pre) && (j > anchor_km_pre_i))
				{
					kmer_ins[kmer_i] = in_j;
					kmer_outs[kmer_i] = out_j;
					kmer_indexs[kmer_i] = pre_j;

					if((kmer_n > filter_min) && (kmer_n < filter_max))
						kmer_types[kmer_i] = 0;
					else	kmer_types[kmer_i] = 1;

					kmer_i++;

					in_j = 0;
					out_j = 0;
					kmer_n = 0;

					pre_j = char_tmp;//j
				}

				input = ((bufMK2[char_tmp] >> tempMove_r_new3) & 0X3);
				output = ((bufMK2[char_tmp + kmer_num - 1] >> keyMove) & 0X3);

				fl_f = (bufMK2[char_tmp + kmer_num_new - 1] & 0X3);

				if(fl_f == 1)
				{
					out_j |= (one_tmp << output);
				}
				else if(fl_f == 2)
				{
					in_j |= (one_tmp << input);
				}
				else if(fl_f == 0)
				{
					in_j |= (one_tmp << input);
					out_j |= (one_tmp << output);
				}

				kmer_pre = anchor;

				kmer_n += (bufMK2[char_tmp + kmer_num_filter - 1] >> 32);
			}

			kmer_ins[kmer_i] = in_j;
			kmer_outs[kmer_i] = out_j;
			kmer_indexs[kmer_i] = pre_j;

			if((kmer_n > filter_min) && (kmer_n < filter_max))
				kmer_types[kmer_i] = 0;
			else	kmer_types[kmer_i] = 1;

			kmer_i++;
		}
		else
		{
			char_tmp = anchor_km_pre_i * kmer_num_filter;

			memcpy(kmer_pre_as, bufMK2 + char_tmp, kmer_num_tmp_new);

			fl_f = (kmer_pre_as[kmer_num_new - 1] & 0X3);

			if(fl_f == 1)
			{
				in_j = 0;
				output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
				out_j = (one_tmp << output);
			}
			else if(fl_f == 2)
			{
				input = ((kmer_pre_as[0] >> tempMove_r_new3) & 0X3);
				in_j = (one_tmp << input);
				out_j = 0;
			}
			else if(fl_f == 0)
			{
				input = ((kmer_pre_as[0] >> tempMove_r_new3) & 0X3);
				in_j = (one_tmp << input);
				output = ((kmer_pre_as[kmer_num - 1] >> keyMove) & 0X3);
				out_j = (one_tmp << output);
			}
			else
			{
				in_j = 0;
				out_j = 0;
			}

			kmer_n = (bufMK2[char_tmp + kmer_num_filter - 1] >> 32);

			kmer_ins[kmer_i] = in_j;
			kmer_outs[kmer_i] = out_j;
			kmer_indexs[kmer_i] = char_tmp;

			if((kmer_n > filter_min) && (kmer_n < filter_max))
				kmer_types[kmer_i] = 0;
			else	kmer_types[kmer_i] = 1;

			kmer_i++;
		}

		for(kmer_k = 0; kmer_k < kmer_i; kmer_k++)
		{
			kmer_pre_as[kmer_num_new - 1] = 0;
			memcpy(kmer_pre_as, bufMK2 + kmer_indexs[kmer_k], kmer_num_tmp);

			in_j = kmer_ins[kmer_k];
			out_j = kmer_outs[kmer_k];
			type_c = kmer_types[kmer_k];

			if(type_c)
			{
				if(in_j)//1 shift right
				{
					add_edge = (kmer_pre_as[kmer_num - 1] >> (keyMove + 2)) & 0X3;

					if(kmer_num > 1)
					{
						shc_pre = 0;
						for(tra_i = 0; tra_i < kmer_num; tra_i++)
						{
							shc = kmer_pre_as[tra_i] & 0X3;
							kmer_pre_as[tra_i] >>= 2;
							kmer_pre_as[tra_i] |= firstls[shc_pre];

							shc_pre = shc;
						}
						kmer_pre_as[kmer_num - 1] &= KMask;

						kmer_pre_as[0] <<= (tempMove_new + 4);
						kmer_pre_as[0] >>= 4;

						kmer_pre_as[kmer_num_new - 1] |= add_edge;

						for(aph_i = 0; aph_i < 4; aph_i++)
						{
							if(!((in_j >> aph_i) & 0X1))	continue;

							kmer_pre_as[0] <<= 4;
							kmer_pre_as[0] >>= 4;
							kmer_pre_as[0] |= firstls[aph_i + 4];
							kmer_pre_as[0] |= type_one;//1

							memcpy(write_filterout[num] + (w_f_o * kmer_num_new), kmer_pre_as, kmer_num_tmp_new);
							w_f_o++;
						}
					}
					else
					{
						kmer_pre_as[0] >>= 4;
						kmer_pre_as[0] <<= 2;

						kmer_pre_as[kmer_num_new - 1] |= add_edge;

						for(aph_i = 0; aph_i < 4; aph_i++)
						{
							if(!((in_j >> aph_i) & 0X1))	continue;

							kmer_pre_as[0] <<= 4;
							kmer_pre_as[0] >>= 4;
							kmer_pre_as[0] |= firstls[aph_i + 4];
							kmer_pre_as[0] |= type_one;//1

							memcpy(write_filterout[num] + (w_f_o * kmer_num_new), kmer_pre_as, kmer_num_tmp_new);
							w_f_o++;
						}
					}

					kmer_pre_as[kmer_num_new - 1] = 0;
					memcpy(kmer_pre_as, bufMK2 + kmer_indexs[kmer_k], kmer_num_tmp);
				}

				if(out_j)//2 shift left
				{
					add_edge = ((kmer_pre_as[0] >> first_bp_move) & 0X3);

					if(kmer_num > 1)
					{
						shc_pre = 0;//firstls[aph_i + 8]
						kmer_pre_as[kmer_num - 1] &= KMask;
						for(tra_i = kmer_num - 1; tra_i >= 0; tra_i--)
						{
							shc = (kmer_pre_as[tra_i] >> 62) & 0X3;
							kmer_pre_as[tra_i] <<= 2;
							kmer_pre_as[tra_i] |= shc_pre;

							shc_pre = shc;
						}

						kmer_pre_as[0] <<= (tempMove_new + 2);
						kmer_pre_as[0] >>= 2;
						kmer_pre_as[0] |= type_two;//2

						for(aph_i = 0; aph_i < 4; aph_i++)
						{
							if(!((out_j >> aph_i) & 0X1))	continue;

							kmer_pre_as[kmer_num - 1] >>= (keyMove + 4);
							kmer_pre_as[kmer_num - 1] <<= 2;
							kmer_pre_as[kmer_num - 1] |= aph_i;
							kmer_pre_as[kmer_num - 1] <<= (keyMove + 2);

							kmer_pre_as[kmer_num_new - 1] |= add_edge;

							memcpy(write_filterout[num] + (w_f_o * kmer_num_new), kmer_pre_as, kmer_num_tmp_new);
							w_f_o++;
						}
					}
					else
					{
						kmer_pre_as[0] <<= 4;
						kmer_pre_as[0] >>= 2;
						kmer_pre_as[0] |= type_two;//2

						kmer_pre_as[kmer_num_new - 1] |= add_edge;

						for(aph_i = 0; aph_i < 4; aph_i++)
						{
							if(!((out_j >> aph_i) & 0X1))	continue;

							kmer_pre_as[0] >>= 4;
							kmer_pre_as[0] <<= 2;
							kmer_pre_as[0] |= aph_i;
							kmer_pre_as[0] <<= 2;

							memcpy(write_filterout[num] + (w_f_o * kmer_num_new), kmer_pre_as, kmer_num_tmp_new);
							w_f_o++;
						}

					}
				}
			}
			else
			{
				if(kmer_num > 1)
				{
					kmer_pre_as[0] <<= (tempMove_new + 2);
					kmer_pre_as[0] >>= 2;
					kmer_pre_as[kmer_num - 1] &= KMask;
				}
				else
				{
					kmer_pre_as[0] <<= 2;
					kmer_pre_as[0] >>= 2;
					kmer_pre_as[0] >>= 2;
					kmer_pre_as[0] <<= 2;
				}

				add_edge = ((in_j << 4) | out_j);
				kmer_pre_as[kmer_num_new - 1] |= add_edge;

				fwrite(kmer_pre_as, kmer_num_tmp_new, 1, fp_filterin);
			}
		}
	}

	w_f_os[num] = w_f_o;
}


void *multiMergeTagcheck(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	int num = argT->num;
	int totalNum = argT->thread;

	uint64_t i = 0;
	uint64_t i_i = 0;
	uint64_t pre_i = 0;
	uint64_t anchor = 0;
	uint64_t anchor_a_pre = 0;
	uint64_t kmer_pre = 0;
	uint64_t start=splitMK2[num], end=splitMK2[num+1];
	uint64_t bufMK1_len=end-start;
	uint64_t buf_t = 0;
	uint8_t buf_t_tmp = 0;
	uint8_t kmer_num_tmp = (kmer_num << 3);

	uint64_t f_num = 0;
	uint64_t r_num = 0;
	uint64_t x_num = 0;

	if(bufMK1_len==0) return (void *)NULL;

	pre_i = start;

	uint8_t pres[4];
	uint8_t pres_i = 0;
	uint8_t pre_k = 0;
	uint8_t pre_flag = 0;
	uint8_t buf_t_s = 0;

	if(kmer_num == 1)
	{
		for(i = start; i < end; i++)
		{
			anchor = (bufMK2_m[i]<<2);

			if((anchor != kmer_pre) && (i > start))
			{
				if(i - pre_i > 1)
				{
					for(i_i = pre_i; i_i < i; i_i++)   //0,2
					{
						buf_t = bufMK2_t[i_i];
						buf_t_s = (buf_t >> 4);
						buf_t &= 0Xf;

						pre_flag = 0;
						for(pre_k = 0; pre_k < pres_i; pre_k++)
							if((pres[pre_k] & buf_t_s) && (pre_k != (i_i - pre_i)))
								pre_flag = 1;

						anchor_a_pre = bufMK2_m[i_i];
						anchor_a_pre >>= 2;

						if(buf_t >> 2)
							anchor_a_pre |= ((buf_t >> 2) << 62);

						buf_t_tmp = (buf_t&0X3);

						if(pre_flag)
						{
							if(buf_t_tmp < 2)  	//(buf_t_tmp == 0) || (buf_t_tmp == 1) 0->1, end, forward Y
							{
								write_fy[num][f_num] =  anchor_a_pre;
								f_num++;
							}
							else  				//2->3, start + end, x
							{
								//memcpy(write_x[num] + w_x_p[num], &anchor_a_pre, 8);
								write_x[num][x_num] = anchor_a_pre;
								x_num++;
							}
						}
						else
						{
							if(buf_t_tmp == 2)  	//start, reverse Y
							{
								write_ry[num][r_num] = anchor_a_pre;
								r_num++;
							}
							else if(buf_t_tmp == 1)  	//0->1, end, forward Y
							{
								write_fy[num][f_num] =  anchor_a_pre;
								f_num++;
							}
							else if(buf_t_tmp == 3)  	//2->3, start + end, x
							{
								write_x[num][x_num] = anchor_a_pre;
								x_num++;
							}
						}
					}
				}
				else
				{
					buf_t = bufMK2_t[pre_i];
					buf_t &= 0Xf;

					anchor_a_pre = bufMK2_m[pre_i];
					anchor_a_pre >>= 2;
					if(buf_t >> 2)
						anchor_a_pre |= ((buf_t >> 2) << 62);

					buf_t_tmp = (buf_t&0X3);
					if(buf_t_tmp == 2)  	//start, reverse Y
					{
						write_ry[num][r_num] = anchor_a_pre;
						r_num++;
					}
					else if(buf_t_tmp == 1)  	//0->1, end, forward Y
					{
						write_fy[num][f_num] =  anchor_a_pre;
						f_num++;
					}
					else if(buf_t_tmp == 3)  	//2->3, start + end, x
					{
						write_x[num][x_num] = anchor_a_pre;
						x_num++;
					}
				}
				pre_i = i;
				pres_i = 0;
			}

			kmer_pre = anchor;

			pres[pres_i] = (bufMK2_t[i] >> 4);
			pres_i++;
		}

#ifdef	BOUNDARY
		if(num != totalNum - 1)
		{
#endif
			if(i - pre_i > 1)
			{
				for(i_i = pre_i; i_i < i; i_i++)   //0,2
				{
					buf_t = bufMK2_t[i_i];
					buf_t_s = (buf_t >> 4);
					buf_t &= 0Xf;

					pre_flag = 0;
					for(pre_k = 0; pre_k < pres_i; pre_k++)
						if((pres[pre_k] & buf_t_s) && (pre_k != (i_i - pre_i)))
							pre_flag = 1;

					anchor_a_pre = bufMK2_m[i_i];
					anchor_a_pre >>= 2;

					if(buf_t >> 2)
						anchor_a_pre |= ((buf_t >> 2) << 62);

					buf_t_tmp = (buf_t&0X3);

					if(pre_flag)
					{
						if(buf_t_tmp < 2)  	//(buf_t_tmp == 0) || (buf_t_tmp == 1) 0->1, end, forward Y
						{
							write_fy[num][f_num] =  anchor_a_pre;
							f_num++;
						}
						else  				//2->3, start + end, x
						{
							write_x[num][x_num] = anchor_a_pre;
							x_num++;
						}
					}
					else
					{
						if(buf_t_tmp == 2)  	//start, reverse Y
						{
							write_ry[num][r_num] = anchor_a_pre;
							r_num++;
						}
						else if(buf_t_tmp == 1)  	//0->1, end, forward Y
						{
							write_fy[num][f_num] =  anchor_a_pre;
							f_num++;
						}
						else if(buf_t_tmp == 3)  	//2->3, start + end, x
						{
							write_x[num][x_num] = anchor_a_pre;
							x_num++;
						}
					}
				}
			}
			else
			{
				buf_t = bufMK2_t[pre_i];
				buf_t &= 0Xf;

				anchor_a_pre = bufMK2_m[pre_i];
				anchor_a_pre >>= 2;
				if(buf_t >> 2)
					anchor_a_pre |= ((buf_t >> 2) << 62);

				buf_t_tmp = (buf_t&0X3);
				if(buf_t_tmp == 2)  	//start, reverse Y
				{
					write_ry[num][r_num] = anchor_a_pre;
					r_num++;
				}
				else if(buf_t_tmp == 1)  	//0->1, end, forward Y
				{
					write_fy[num][f_num] =  anchor_a_pre;
					f_num++;
				}
				else if(buf_t_tmp == 3)  	//2->3, start + end, x
				{
					write_x[num][x_num] = anchor_a_pre;
					x_num++;
				}
			}

#ifdef	BOUNDARY
		}
		else	last_block_num = pre_i;
#endif
	}
	else
	{
		uint64_t anchor_com[kmer_num];
		uint64_t anchor_com_pre[kmer_num];
		uint64_t anchor_com_a_pre[kmer_num];

		for(i = start; i < end; i++)
		{
			//anchor = (bufMK2_m[i]<<2);
			getkmerno_psm(anchor_com, bufMK2_m, i * kmer_num)

			if(memcmp(anchor_com, anchor_com_pre, kmer_num_tmp) && (i > start))
			{
				if(i - pre_i > 1)
				{
					for(i_i = pre_i; i_i < i; i_i++)   //0,2
					{
						buf_t = bufMK2_t[i_i];
						buf_t_s = (buf_t >> 4);
						buf_t &= 0Xf;

						pre_flag = 0;
						for(pre_k = 0; pre_k < pres_i; pre_k++)
							if((pres[pre_k] & buf_t_s) && (pre_k != (i_i - pre_i)))
								pre_flag = 1;

						getkmerno_s(anchor_com_a_pre, bufMK2_m, i_i * kmer_num)
						anchor_com_a_pre[0] >>= (tempMove_new + 2);
						if(buf_t >> 2)
							anchor_com_a_pre[0] |= ((buf_t >> 2) << tempMove_r_new3);

						buf_t_tmp = (buf_t&0X3);

						if(pre_flag)
						{
							if(buf_t_tmp < 2)  	//(buf_t_tmp == 0) || (buf_t_tmp == 1) 0->1, end, forward Y
							{
								memcpy(write_fy[num] + f_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
								f_num++;
							}
							else  				//2->3, start + end, x
							{
								memcpy(write_x[num] + x_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
								x_num++;
							}
						}
						else
						{
							if(buf_t_tmp == 2)  	//start, reverse Y
							{
								memcpy(write_ry[num] + r_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
								r_num++;
							}
							else if(buf_t_tmp == 1)  	//0->1, end, forward Y
							{
								memcpy(write_fy[num] + f_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
								f_num++;
							}
							else if(buf_t_tmp == 3)  	//2->3, start + end, x
							{
								memcpy(write_x[num] + x_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
								x_num++;
							}
						}
					}
				}
				else
				{
					buf_t = bufMK2_t[pre_i];
					buf_t &= 0Xf;

					getkmerno_s(anchor_com_a_pre, bufMK2_m, pre_i * kmer_num)
					anchor_com_a_pre[0] >>= (tempMove_new + 2);
					if(buf_t >> 2)
						anchor_com_a_pre[0] |= ((buf_t >> 2) << tempMove_r_new3);

					buf_t_tmp = (buf_t&0X3);

					if(buf_t_tmp == 2)  	//start, reverse Y
					{
						memcpy(write_ry[num] + r_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
						r_num++;
					}
					else if(buf_t_tmp == 1)  	//0->1, end, forward Y
					{
						memcpy(write_fy[num] + f_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
						f_num++;
					}
					else if(buf_t_tmp == 3)  	//2->3, start + end, x
					{
						memcpy(write_x[num] + x_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
						x_num++;
					}
				}
				pre_i = i;
				pres_i = 0;
			}

			memcpy(anchor_com_pre, anchor_com, kmer_num_tmp);

			pres[pres_i] = (bufMK2_t[i] >> 4);
			pres_i++;
		}

#ifdef	BOUNDARY
		if(num != totalNum - 1)
		{
#endif
			if(i - pre_i > 1)
			{
				for(i_i = pre_i; i_i < i; i_i++)   //0,2
				{
					buf_t = bufMK2_t[i_i];
					buf_t_s = (buf_t >> 4);
					buf_t &= 0Xf;

					pre_flag = 0;
					for(pre_k = 0; pre_k < pres_i; pre_k++)
						if((pres[pre_k] & buf_t_s) && (pre_k != (i_i - pre_i)))
							pre_flag = 1;

					getkmerno_s(anchor_com_a_pre, bufMK2_m, i_i * kmer_num)
					anchor_com_a_pre[0] >>= (tempMove_new + 2);
					if(buf_t >> 2)
						anchor_com_a_pre[0] |= ((buf_t >> 2) << tempMove_r_new3);

					buf_t_tmp = (buf_t&0X3);

					if(pre_flag)
					{
						if(buf_t_tmp < 2)  	//(buf_t_tmp == 0) || (buf_t_tmp == 1) 0->1, end, forward Y
						{
							memcpy(write_fy[num] + f_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
							f_num++;
						}
						else  				//2->3, start + end, x
						{
							memcpy(write_x[num] + x_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
							x_num++;
						}
					}
					else
					{
						if(buf_t_tmp == 2)  	//start, reverse Y
						{
							memcpy(write_ry[num] + r_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
							r_num++;
						}
						else if(buf_t_tmp == 1)  	//0->1, end, forward Y
						{
							memcpy(write_fy[num] + f_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
							f_num++;
						}
						else if(buf_t_tmp == 3)  	//2->3, start + end, x
						{
							memcpy(write_x[num] + x_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
							x_num++;
						}
					}
				}
			}
			else
			{
				buf_t = bufMK2_t[pre_i];
				buf_t &= 0Xf;

				getkmerno_s(anchor_com_a_pre, bufMK2_m, pre_i * kmer_num)
				anchor_com_a_pre[0] >>= (tempMove_new + 2);
				if(buf_t >> 2)
					anchor_com_a_pre[0] |= ((buf_t >> 2) << tempMove_r_new3);

				buf_t_tmp = (buf_t&0X3);

				if(buf_t_tmp == 2)  	//start, reverse Y
				{
					memcpy(write_ry[num] + r_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
					r_num++;
				}
				else if(buf_t_tmp == 1)  	//0->1, end, forward Y
				{
					memcpy(write_fy[num] + f_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
					f_num++;
				}
				else if(buf_t_tmp == 3)  	//2->3, start + end, x
				{
					memcpy(write_x[num] + x_num * kmer_num, anchor_com_a_pre, kmer_num_tmp);
					x_num++;
				}
			}

#ifdef	BOUNDARY
		}
		else	last_block_num = pre_i;
#endif
	}

	w_f_p[num] = f_num;
	w_r_p[num] = r_num;
	w_x_p[num] = x_num;

	//fprintf(stderr, "num: %u last_block_num : %u %u %u %u %u\n", num, last_block_num, start, end, i, pre_i);
	//fflush(stderr);

	return (void *)NULL;
}

void multiMergeTagcheck_addkmer(uint64_t i, uint64_t pre_i, FILE* fp_fy, FILE* fp_ry, FILE* fp_x)
{
	uint8_t pres[4];
	uint8_t pres_i = 0;
	uint8_t pre_k = 0;
	uint8_t pre_flag = 0;
	uint8_t buf_t_s = 0;
	uint8_t buf_t_tmp = 0;

	uint64_t i_i = 0;
	uint64_t buf_t = 0;

	for(i_i = pre_i; i_i < i; i_i++)
	{
		pres[pres_i] = (bufMK2_t[i_i] >> 4);
		pres_i++;
	}


	if(kmer_num == 1)
	{
		uint64_t anchor_a_pre = 0;

		if(i - pre_i > 1)
		{
			for(i_i = pre_i; i_i < i; i_i++)   //0,2
			{
				buf_t = bufMK2_t[i_i];
				buf_t_s = (buf_t >> 4);
				buf_t &= 0Xf;

				pre_flag = 0;
				for(pre_k = 0; pre_k < pres_i; pre_k++)
					if((pres[pre_k] & buf_t_s) && (pre_k != (i_i - pre_i)))
						pre_flag = 1;

				anchor_a_pre = bufMK2_m[i_i];
				anchor_a_pre >>= 2;

				if(buf_t >> 2)
					anchor_a_pre |= ((buf_t >> 2) << 62);

				buf_t_tmp = (buf_t&0X3);

				if(pre_flag)
				{
					if(buf_t_tmp < 2)  	//(buf_t_tmp == 0) || (buf_t_tmp == 1) 0->1, end, forward Y
					{
						fwrite(&anchor_a_pre, 8, 1, fp_fy);
						fy_n++;
					}
					else  				//2->3, start + end, x
					{
						fwrite(&anchor_a_pre, 8, 1, fp_x);
						x_n++;
					}
				}
				else
				{
					if(buf_t_tmp == 2)  	//start, reverse Y
					{
						fwrite(&anchor_a_pre, 8, 1, fp_ry);
						ry_n++;
					}
					else if(buf_t_tmp == 1)  	//0->1, end, forward Y
					{
						fwrite(&anchor_a_pre, 8, 1, fp_fy);
						fy_n++;
					}
					else if(buf_t_tmp == 3)  	//2->3, start + end, x
					{
						fwrite(&anchor_a_pre, 8, 1, fp_x);
						x_n++;
					}
				}
			}
		}
		else
		{
			buf_t = bufMK2_t[pre_i];
			buf_t &= 0Xf;

			anchor_a_pre = bufMK2_m[pre_i];
			anchor_a_pre >>= 2;
			if(buf_t >> 2)
				anchor_a_pre |= ((buf_t >> 2) << 62);

			buf_t_tmp = (buf_t&0X3);
			if(buf_t_tmp == 2)  	//start, reverse Y
			{
				fwrite(&anchor_a_pre, 8, 1, fp_ry);
				ry_n++;
			}
			else if(buf_t_tmp == 1)  	//0->1, end, forward Y
			{
				fwrite(&anchor_a_pre, 8, 1, fp_fy);
				fy_n++;
			}
			else if(buf_t_tmp == 3)  	//2->3, start + end, x
			{
				fwrite(&anchor_a_pre, 8, 1, fp_x);
				x_n++;
			}
		}
	}
	else
	{
		uint8_t kmer_num_tmp = (kmer_num << 3);
		uint64_t anchor_com_a_pre[kmer_num];

		if(i - pre_i > 1)
		{
			for(i_i = pre_i; i_i < i; i_i++)   //0,2
			{
				buf_t = bufMK2_t[i_i];
				buf_t_s = (buf_t >> 4);
				buf_t &= 0Xf;

				pre_flag = 0;
				for(pre_k = 0; pre_k < pres_i; pre_k++)
					if((pres[pre_k] & buf_t_s) && (pre_k != (i_i - pre_i)))
						pre_flag = 1;

				getkmerno_s(anchor_com_a_pre, bufMK2_m, i_i * kmer_num)
				anchor_com_a_pre[0] >>= (tempMove_new + 2);
				if(buf_t >> 2)
					anchor_com_a_pre[0] |= ((buf_t >> 2) << tempMove_r_new3);

				buf_t_tmp = (buf_t&0X3);

				if(pre_flag)
				{
					if(buf_t_tmp < 2)  	//(buf_t_tmp == 0) || (buf_t_tmp == 1) 0->1, end, forward Y
					{
						fwrite(anchor_com_a_pre, kmer_num_tmp, 1, fp_fy);
						fy_n++;
					}
					else  				//2->3, start + end, x
					{
						fwrite(anchor_com_a_pre, kmer_num_tmp, 1, fp_x);
						x_n++;
					}
				}
				else
				{
					if(buf_t_tmp == 2)  	//start, reverse Y
					{
						fwrite(anchor_com_a_pre, kmer_num_tmp, 1, fp_ry);
						ry_n++;
					}
					else if(buf_t_tmp == 1)  	//0->1, end, forward Y
					{
						fwrite(anchor_com_a_pre, kmer_num_tmp, 1, fp_fy);
						fy_n++;
					}
					else if(buf_t_tmp == 3)  	//2->3, start + end, x
					{
						fwrite(anchor_com_a_pre, kmer_num_tmp, 1, fp_x);
						x_n++;
					}
				}
			}
		}
		else
		{
			buf_t = bufMK2_t[pre_i];
			buf_t &= 0Xf;

			getkmerno_s(anchor_com_a_pre, bufMK2_m, pre_i * kmer_num)
			anchor_com_a_pre[0] >>= (tempMove_new + 2);
			if(buf_t >> 2)
				anchor_com_a_pre[0] |= ((buf_t >> 2) << tempMove_r_new3);

			buf_t_tmp = (buf_t&0X3);

			if(buf_t_tmp == 2)  	//start, reverse Y
			{
				fwrite(anchor_com_a_pre, kmer_num_tmp, 1, fp_ry);
				ry_n++;
			}
			else if(buf_t_tmp == 1)  	//0->1, end, forward Y
			{
				fwrite(anchor_com_a_pre, kmer_num_tmp, 1, fp_fy);
				fy_n++;
			}
			else if(buf_t_tmp == 3)  	//2->3, start + end, x
			{
				fwrite(anchor_com_a_pre, kmer_num_tmp, 1, fp_x);
				x_n++;
			}
		}
	}
}


void getbranchkmer(uint64_t THREAD_NUM, char* tmp_route)
{
	uint8_t key_len = 0;

	//there may sth error here
	KMER_LENGTH_PlusOne = KMER_LENGTH + 1;
	KMER_LENGTH_PlusTwo = KMER_LENGTH + 2;

	if(KMER_LENGTH_PlusTwo < 32)	OneUnit = 1;
	else	OneUnit = 0;

	kmer_num = (KMER_LENGTH_PlusTwo >> 5);
	kmer_rem = (KMER_LENGTH_PlusTwo & 0X1f);
	if(kmer_rem)	kmer_num++;

	kmer_num_new = kmer_num;

	uint64_t key_tmp_re = (((KMER_LENGTH << 1) + 1) << 2);
	kmer_num_br = (key_tmp_re >> 6);
	if(key_tmp_re & 0X3f)	kmer_num_br++;

	uint64_t tempMove_r_new = ((KMER_LENGTH_PlusTwo & 0X3) << 1);		//0 6 4 2
	tempMove_r_new2 = (tempMove_r_new == 0 ? 6:(tempMove_r_new - 2));	//6 4 2 0
	tempMove_r_new3 = tempMove_r_new2 + 56;								//62 60 58 56
	tempMove_new = (tempMove_r_new == 0 ? 0:(8 - tempMove_r_new));		//0 2 4 6
	tempMove_new2 = (tempMove_r_new3 - (BUCKET_LENGTH << 1));

	temphash = 3;
	br_off = (tempMove_new >> 1);										//0 1 2 3
	br_offm = br_off + 1;												//1 2 3 4

	uint64_t key_bits = 0;
	key_bits = (KMER_LENGTH_PlusTwo << 1);
	key_len = key_bits / 8 + !! (key_bits % 8);

	key_tmp_re = (key_len & 0X7);
	keyMove = 64 - (key_tmp_re << 3);
	if(keyMove == 64)	keyMove	= 0;

	leftmove32 = 0;
	if(key_tmp_re)
	{
		KMask = (~(((uint64_t)1 << ((33 - (key_tmp_re << 2)) << 1)) - 1));
		KMasks = (~(((uint64_t)1 << ((34 - (key_tmp_re << 2)) << 1)) - 1));
	}
	else
	{
		KMask = ~((uint64_t)3);
		KMasks = ~((uint64_t)15);

		if(KMER_LENGTH_PlusTwo > 31)	kmer_num_new++;

		if(KMER_LENGTH_PlusTwo == 32)
		{
			KMask <<= 2;
			KMasks <<= 2;
			leftmove32 = 2;
		}
	}

	if(flag_filter)	kmer_num_filter = (kmer_num_new == kmer_num ? kmer_num + 1:kmer_num_new);
	else kmer_num_filter = kmer_num_new;

	if(OneUnit)	kmer_num_filtermerge = kmer_num << 4;
	else	kmer_num_filtermerge = kmer_num_new << 3;

	if(OneUnit)	kmer_num_fm = kmer_num << 1;
	else	kmer_num_fm = kmer_num_new;

	first_bp_move = tempMove_r_new3 - 2;

	one_out = 62 - (KMER_LENGTH_PlusTwo << 1);
	one_suf = one_out + 2;
	one_sufp = one_suf + 2;

	if(flag_filter)
	{
		firstls[0] = 0;
		firstls[1] = ((uint64_t)1 << 62);
		firstls[2] = ((uint64_t)2 << 62);
		firstls[3] = ((uint64_t)3 << 62);

		firstls[4] = 0;
		firstls[5] = ((uint64_t)1 << 60);
		firstls[6] = ((uint64_t)2 << 60);
		firstls[7] = ((uint64_t)3 << 60);

		firstls[8] = 0;
		firstls[9] = ((uint64_t)1 << keyMove);
		firstls[10] = ((uint64_t)2 << keyMove);
		firstls[11] = ((uint64_t)3 << keyMove);

		type_one = ((uint64_t)1 << 62);
		type_two = ((uint64_t)2 << 62);
	}

	long start1 = 0;
	long end1 = 0;
	float cost_time = 0;

	uint64_t bufferSize = 0;
	bufferSize = (uint64_t )((double)((uint64_t)memoryKmer<<20) / (double)(((KMER_LENGTH << 2) + 1) * (kmer_num << 4)));

	uint64_t bufferSizeM = bufferSize * THREAD_NUM;
	uint64_t bufferSize_tmp = 0;
	uint64_t eliminate=~(ELIMINATE);
	uint64_t start = 0, end = 0;
	uint64_t char_tmp = 0;
	uint64_t well_tmp = 4;

	uint16_t w_i = 0;
	uint64_t w_in = 0;

	int i = 0;
	int j = 0, num = 0;
	int L_NUM = 1;//6 3
	int temp = 0;
	char cNum[6];
	char string_tmp[MAX_ROUTE];

	uint8_t judge_tmp = 0;
	int64_t heapTail = 0;
	uint64_t pBufM = 0;
	uint64_t f1_tmp64 = 0;
	uint64_t s1_tmp64 = 0;
	uint64_t f2_tmp64 = 0;
	uint64_t s2_tmp64 = 0;
	uint64_t fy_size = 0;
	uint64_t x_size = 0;
	uint64_t totalKmerNum = 0;
	uint64_t kmer_num_br_char = kmer_num_br << 3;
	uint64_t* write_tmp_array = calloc(kmer_num_br, 8);

	clock_t start_c = 0, end_c = 0;
	time_t start_t = 0, end_t = 0;

#ifndef	LAST_HASH_DIS
	uint64_t* last_hashKmer_write = NULL;
#endif

	wta_n = (((KMER_LENGTH + 1) & 0Xf) == 0 ? ((KMER_LENGTH + 1) >> 4):(((KMER_LENGTH + 1) >> 4) + 1));

	fprintf(stderr, "Begin generating branching kmers wta_n: %"PRId64"\n", wta_n);

	FILE* fp_first = NULL;
	FILE* fp_second = NULL;
	FILE* fpKmer = NULL;

#ifdef	BEF_BRANCH
	//ry_size = ry_n;
	fy_size = fy_n;
	x_size = x_n;

	if((fy_size + x_size) < bufferSize)	bufferSize = (fy_size + x_size);

	strcpy(string_tmp, "cat ");
	strcat(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerFY ");
	strcat(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerX > ");
	strcat(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerFYN");
	system(string_tmp);

	strcpy(string_tmp, "cat ");
	strcat(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerRY ");
	strcat(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerX > ");
	strcat(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerRYN");
	system(string_tmp);

	strcpy(string_tmp, "rm -f ");
	strcat(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerFY");
	system(string_tmp);

	strcpy(string_tmp, "rm -f ");
	strcat(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerRY");
	system(string_tmp);

	strcpy(string_tmp, "rm -f ");
	strcat(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerX");
	system(string_tmp);
#endif

	strcpy(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerFYN");

	//fprintf(stderr, "%s\n", string_tmp);
	fp_first = fopen(string_tmp,"rb");
	if(fp_first == NULL) fprintf(stderr,"Fail to open kmerFYN\n"),exit(1);

	strcpy(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerRYN");
	fp_second = fopen(string_tmp,"rb");
	if(fp_second == NULL) fprintf(stderr,"Fail to open kmerRYN\n"),exit(1);


	l_segCount = (uint64_t* )calloc(THREAD_NUM+1,sizeof(uint64_t));
	l_segCount[THREAD_NUM] = BUCKET_CAPACITY;

	fprintf(stderr, "Begin branching using %"PRId64" threads %"PRId64"\n", THREAD_NUM, bufferSize);

	gettimeofday(&t_start1, NULL);
	start1 = ((long)t_start1.tv_sec)*1000+(long)t_start1.tv_usec/1000;

	memset(write_tmp_array, 0, kmer_num_br << 3);

	uint8_t kmer_num_tmp = (kmer_num << 3);

	first_buffer = (uint64_t* )calloc(bufferSize+1, kmer_num_tmp);
	second_buffer = (uint64_t* )calloc(bufferSize+1, kmer_num_tmp);

	uint64_t hash_c = 0;
	uint64_t hash_dis = 0;
	uint64_t locateNum = 0;
	uint64_t locate = 0;
	uint8_t fp_flag = 0;

	thread_data* tt = NULL;

	pthread_t* myThread = NULL;

	last_hashKmer_dis = (uint64_t** )calloc(THREAD_NUM, sizeof(uint64_t* ));

	last_countKmer_dis = (uint64_t** )calloc(THREAD_NUM,sizeof(uint64_t* ));

	while(!(feof(fp_first) || feof(fp_second)))   //
	{
		read_num = fread(first_buffer, kmer_num_tmp, bufferSize, fp_first);
		read_num = fread(second_buffer, kmer_num_tmp, read_num, fp_second);

		if(!read_num)	break;

		myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
		tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));
		for(i=0; i<THREAD_NUM; i++)
		{
			tt[i].num = i;
			tt[i].thread = THREAD_NUM;

			int check=pthread_create(&myThread[i], NULL, last_multiCount, tt + i);
			if(check)
			{
				fprintf(stderr,"threadNum:%d, Error - pthread_create() return code: %d\n",i,check);
				exit(EXIT_FAILURE);
			}
		}
		for(i=0; i<THREAD_NUM; i++)
		{
			pthread_join(myThread[i], NULL);
		}
		free(myThread);
		free(tt);

		fprintf(stderr, "Block %d kmers counting and distrbuting finish\n", L_NUM);

		last_countKmer = (uint64_t* )calloc(BUCKET_CAPACITY+1,sizeof(uint64_t));

		last_countKmer[0] = 0;
		for(i=1; i<BUCKET_CAPACITY+1; i++)
		{
			for(j=0; j<THREAD_NUM; j++)
			{
				last_countKmer[i] += (last_countKmer_dis[j][i] - last_countKmer_dis[j][i - 1]);
			}
			last_countKmer[i] += last_countKmer[i - 1];
		}
		totalKmerNum = last_countKmer[BUCKET_CAPACITY];

		last_hashKmer = (uint64_t* )calloc(totalKmerNum, kmer_num_br << 3);

		for(i=0; i<BUCKET_CAPACITY; i++)
		{
			for(j=0, hash_c = 0; j<THREAD_NUM; j++)
			{
				hash_dis = last_countKmer_dis[j][i+1] - last_countKmer_dis[j][i];

				if(hash_dis)
					memcpy(last_hashKmer + (last_countKmer[i] + hash_c) * kmer_num_br, last_hashKmer_dis[j] + last_countKmer_dis[j][i] * kmer_num_br, hash_dis * kmer_num_br_char);

				hash_c += hash_dis;
			}
		}

		fprintf(stderr, "\ntotalKmerNum: %"PRId64"\n", totalKmerNum);

		l_segCount[0]=0;
		for(i=1; i<THREAD_NUM; i++)
		{
			locateNum = i * (totalKmerNum/THREAD_NUM);
			locate = BinarySearch(locateNum,last_countKmer,BUCKET_CAPACITY-1);
			l_segCount[i] = locate;
		}

#ifdef	LAST_HASH_DIS
		last_hashKmer_write = (uint64_t* )malloc((totalKmerNum * wta_n) << 3);
#endif
		myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
		tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));
		for(i=0; i<THREAD_NUM; i++)
		{
			tt[i].thread = THREAD_NUM;
			tt[i].num = i;
			int check=pthread_create(&myThread[i], NULL, last_multiThreadSort, tt + i);

			if(check)
			{
				fprintf(stderr,"threadNum:%d, Error - pthread_create() return code: %d\n",i,check);
				exit(EXIT_FAILURE);
			}
		}

		for(i=0; i<THREAD_NUM; i++)
		{
			pthread_join( myThread[i], NULL);
		}
		free(myThread);
		free(tt);

		fprintf(stderr, "Block %d kmers sorting finish\n", L_NUM);//: %u, totalKmerNum

		strcpy(string_tmp, tmp_route);
		strcat(string_tmp, "/last_kmer.");

		sprintf(cNum,"%d", L_NUM);
		strcat(string_tmp, cNum);
		fpKmer = fopen(string_tmp,"wb");

#ifndef	LAST_HASH_DIS
		last_hashKmer_write = (uint64_t* )malloc((totalKmerNum * wta_n) << 3);
		for(i = 0; i < totalKmerNum; i++)
			memcpy(last_hashKmer_write + (i * wta_n), last_hashKmer + (i * kmer_num_br), wta_n << 3);
#endif

		fwrite(last_hashKmer_write, wta_n << 3, totalKmerNum, fpKmer);

		if(fpKmer)	fclose(fpKmer);
		if(last_hashKmer_write)	free(last_hashKmer_write);
		if(last_countKmer)	free(last_countKmer);
		if(last_hashKmer)	free(last_hashKmer);

		L_NUM++;
	}

	if(first_buffer)	free(first_buffer);
	if(second_buffer)	free(second_buffer);

	if(last_hashKmer_dis)	free(last_hashKmer_dis);
	if(last_countKmer_dis)	free(last_countKmer_dis);

	strcpy(string_tmp, "rm -f ");
	strcat(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerFYN");
	system(string_tmp);

	strcpy(string_tmp, "rm -f ");
	strcat(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerRYN");
	system(string_tmp);

	gettimeofday(&t_end1, NULL);
	end1 = ((long)t_end1.tv_sec)*1000+(long)t_end1.tv_usec/1000;

	cost_time = end1 - start1;
	fprintf(stderr, "Generating branching sorted kmers time is %.2lf s\n", cost_time/1000);

	fprintf(stderr, "Begin merge and output BWT seq to file using %"PRId64" threads\n", THREAD_NUM);//, L_NUM

	gettimeofday(&t_start1, NULL);
	start1 = ((long)t_start1.tv_sec)*1000+(long)t_start1.tv_usec/1000;

	int k = 0;

	uint64_t* pBuf = NULL;
	uint64_t* readNum_b = NULL;
	uint64_t** bufK2_b = NULL;
	int* heap = NULL;

	FILE *fpK2_b[L_NUM];
	FILE* fpK2_b_t = NULL;
#ifdef	GFA_COM
	FILE* fpK2_e_t = NULL;
#endif

	heap = (int *)calloc(L_NUM,sizeof(int));
	pBuf = (uint64_t* )calloc(L_NUM,sizeof(uint64_t));

	readNum_b = (uint64_t* )calloc(L_NUM,sizeof(uint64_t));

	bufferSize = (uint64_t )((double)((uint64_t)(memoryKmer - 256) << 20) / (double)(((wta_n << 3) * L_NUM) + kmer_num_tmp + 2)) - 1;//8 * t * wta_n + kmer_num_tmp * t + t + 256M
	fprintf(stderr, "Branch: %"PRId64" %"PRId64"\n", ((uint64_t)memoryKmer<<20), ((wta_n << 3) + kmer_num_tmp + 1));//2 * 8 + 8 + 1 = 25


	bufferSizeM = bufferSize * THREAD_NUM;

	bufK2_b = (uint64_t** )calloc(L_NUM, sizeof(uint64_t* ));
	for(num=0; num<L_NUM; num++)
		bufK2_b[num]=(uint64_t* )calloc((bufferSize+1) * wta_n,sizeof(uint64_t));

	uint64_t* bufk = (uint64_t* )calloc(bufferSize + 1, kmer_num_tmp);
	uint8_t* bufk_t = (uint8_t* )calloc(bufferSize + 1, 1);

#ifdef	GFA_COM
	uint8_t* bufk_e_t = (uint8_t* )calloc(bufferSize + 1, 1);
#endif

	strcpy(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerMergewhole");
	fpK2_b[0]=fopen(string_tmp,"rb");
	if(fpK2_b[0]==NULL) fprintf(stderr,"Fail to open last kmerMergewhole\n"),exit(1);

	strcpy(string_tmp, tmp_route);
	strcat(string_tmp, "/kmerMergewhole_t");
	fpK2_b_t = fopen(string_tmp,"rb");
	if(fpK2_b_t == NULL) fprintf(stderr,"Fail to open last kmerMergewhole_t\n"),exit(1);

#ifdef	GFA_COM
	strcpy(string_tmp, tmp_route);
	strcat(string_tmp, "/edges_t");
	fpK2_e_t = fopen(string_tmp,"rb");
	if(fpK2_e_t == NULL) fprintf(stderr,"Fail to open edges_t\n"),exit(1);
#endif

	uint64_t* startp = NULL;

	if(kmer_num == 1)
	{
		br_off = 0;
		br_offm = 0;
	}

	readNum_b[0] = fread(bufk, kmer_num_tmp, bufferSize, fpK2_b[0]);
	fread(bufk_t, 1, bufferSize, fpK2_b_t);

#ifdef	GFA_COM
	fread(bufk_e_t, 1, bufferSize, fpK2_e_t);
#endif

	int check = 0;
	int twofourcnt = 0;
#ifdef	TWOBITTOTOFOURBIT

	if(readNum_b[0] > 0)
	{
		uint64_t segment=readNum_b[0]/THREAD_NUM;

		splitMK2[0]=0; //lower bound (can reach)
		for(i=1; i<THREAD_NUM; i++)
		{
			splitMK2[i]=i*segment;
		}
		splitMK2[THREAD_NUM]=readNum_b[0];

		fprintf(stderr, "original sorted 2-bit k-mers transforming to 4-bit k-mers start for block %d\n", twofourcnt);
		twofourcnt++;

		myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
		thread_branch* tt = (thread_branch* )calloc(THREAD_NUM, sizeof(thread_branch ));
		for(i=0; i<THREAD_NUM; i++)
		{
			tt[i].wta_n = wta_n;
			tt[i].num = i;
			tt[i].readp = bufk + splitMK2[i] * kmer_num;
			tt[i].readp_t = bufk_t + splitMK2[i];
			tt[i].writep = bufK2_b[0] + splitMK2[i] * wta_n;

			check=pthread_create(&myThread[i], NULL, twobittofourbit, tt + i);

			if(check)
			{
				fprintf(stderr, "Thread Num: %d, Error: pthread create return code: %d\n", i, check);
				exit(EXIT_FAILURE);
			}
		}

		for(i=0; i<THREAD_NUM; i++)
		{
			pthread_join(myThread[i], NULL);
		}
		free(myThread);
		free(tt);
	}


#else
	for(i = 0; i < readNum_b[0]; i++)
	{
		startp = bufk + i * kmer_num;
		end = (bufk_t[i] & 0Xf);

		if(((end & 0X3) == 2) || ((end & 0X3) == 3))	write_tmp_array[0] = 4;
		else	write_tmp_array[0] = (end >> 2);

		if(kmer_num != 1)
			startp[0] >>= (tempMove_new + 2);

		for(j = 0, k = 1; j < KMER_LENGTH; j++, k++)
		{
			char_tmp = ((startp[(j + br_offm) >> 5] >> ((31 - ((j + br_offm) & 0X1f)) << 1)) & 0X3);

			write_tmp_array[k >> 4] <<= 4;
			if(char_tmp)	write_tmp_array[k >> 4] |= char_tmp;

		}

		if(k & 0Xf)	write_tmp_array[k >> 4] <<= ((16 - (k & 0Xf)) << 2);

		memcpy(bufK2_b[0] + i * wta_n, write_tmp_array, wta_n << 3);
	}

#endif

	pBuf[0]=0;
	for(num=1; num<L_NUM; num++)   //1
	{
		sprintf(cNum,"%u",num);
		strcpy(string_tmp, tmp_route);
		strcat(string_tmp, "/last_kmer.");
		strcat(string_tmp, cNum);

		fpK2_b[num]=fopen(string_tmp,"rb");
		if(fpK2_b[num]==NULL) fprintf(stderr,"Fail to open last_kmer.%u\n",num),exit(1);

		//readNum[num]=fread(&bufK2[num][1],sizeof(uint64_t)*2,bufferSize,fpK2[num]);
		readNum_b[num]=fread(bufK2_b[num], wta_n << 3, bufferSize, fpK2_b[num]);// sizeof(kmer_s )

		pBuf[num]=0;
	}

	FILE* fp_bwt_seq = fopen(out_route, "wb");

	if(fp_bwt_seq == NULL)	fprintf(stderr, "Fail to open bwt seq file\n"), exit(1);
	fseek(fp_bwt_seq, 8, SEEK_SET);

	uint64_t* bwt_seq_buffer = (uint64_t* )calloc(MAX_MWA_L, 8);
	uint64_t seq_c_n = 0;
	uint64_t seq_n = 0;
	uint64_t max_mwa_l_char = MAX_MWA_L << 4;

#ifdef	GFA_COM
	strcpy(string_tmp, tmp_route);
	strcat(string_tmp, "/edges.seq");
	FILE* fp_edges_seq = fopen(string_tmp, "wb");

	if(fp_edges_seq == NULL)	fprintf(stderr, "Fail to open bwt seq file\n"), exit(1);
	fseek(fp_edges_seq, 8, SEEK_SET);
	uint64_t* edges_seq_buffer = (uint64_t* )calloc(MAX_MWA_L, 8);
	uint64_t seq_e_n = 0;
	uint64_t seq_n_e = 0;
	uint64_t max_edges_char = MAX_MWA_L << 3;
#endif


	for(num=0; num<L_NUM; num++)
	{
		heapTail=num;
		for(i=heapTail; i>0; i=(i-1)>>1)
		{
			//if(L_cmpMK2(bufK2_b[num] + wta_n, bufK2_b[heap[(i-1)>>1]] + wta_n))
			if(L_cmpMK2(bufK2_b[num], bufK2_b[heap[(i-1)>>1]]))
			{
				heap[i]=heap[(i-1)>>1];
			}
			else break;
		}
		heap[i]=num;
	}

	uint64_t ubwt_cnt_tol[5];
	memset(ubwt_cnt_tol, 0, 40);

	int cdi1=0,cdi2=0;
	int heap_tmp = 0;

	while(heapTail>=0)
	{
		pBufM=0;

		while(pBufM<bufferSizeM&&heapTail>=0)   //get bufMK2
		{
			num=heap[0];

			char_tmp = ((bufK2_b[num][pBuf[num] * wta_n] >> 60) & 0Xf);
			bwt_seq_buffer[seq_c_n >> 4] <<= 4;//not need set 0
			bwt_seq_buffer[seq_c_n >> 4] |= char_tmp;

			++seq_c_n;
			if(seq_c_n == max_mwa_l_char)
			{
				fwrite(bwt_seq_buffer, 8, MAX_MWA_L, fp_bwt_seq);
				seq_c_n = 0;
				seq_n++;
			}

#ifdef	GFA_COM
			if(num == 0)
				char_tmp = bufk_e_t[pBuf[0]];
			else	char_tmp = 0;

			edges_seq_buffer[seq_e_n >> 3] <<= 8;//not need set 0
			edges_seq_buffer[seq_e_n >> 3] |= char_tmp;

			++seq_e_n;
			if(seq_e_n == max_edges_char)
			{
				fwrite(edges_seq_buffer, 8, MAX_MWA_L, fp_edges_seq);
				seq_e_n = 0;
				seq_n_e++;
			}
#endif

			pBufM++;
			++pBuf[num];

			if(pBuf[num]==readNum_b[num])
			{
				if(num == 0)
				{
					readNum_b[0] = fread(bufk, kmer_num_tmp, bufferSize, fpK2_b[0]);
					fread(bufk_t, 1, bufferSize, fpK2_b_t);
#ifdef	GFA_COM
					fread(bufk_e_t, 1, bufferSize, fpK2_e_t);
#endif

#ifdef	TWOBITTOTOFOURBIT
					if(readNum_b[0] > 0)
					{
						uint64_t segment=readNum_b[0]/THREAD_NUM;

						splitMK2[0]=0; //lower bound (can reach)
						for(i=1; i<THREAD_NUM; i++)
						{
							splitMK2[i]=i*segment;
						}
						splitMK2[THREAD_NUM]=readNum_b[0];

						fprintf(stderr, "original sorted 2-bit k-mers transforming to 4-bit k-mers start for block %d\n", twofourcnt);
						twofourcnt++;

						myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
						thread_branch* tt = (thread_branch* )calloc(THREAD_NUM, sizeof(thread_branch ));
						for(i=0; i<THREAD_NUM; i++)
						{
							tt[i].wta_n = wta_n;
							tt[i].num = i;
							tt[i].readp = bufk + splitMK2[i] * kmer_num;
							tt[i].readp_t = bufk_t + splitMK2[i];
							tt[i].writep = bufK2_b[0] + splitMK2[i] * wta_n;

							check=pthread_create(&myThread[i], NULL, twobittofourbit, tt + i);

							if(check)
							{
								fprintf(stderr, "Thread Num: %d, Error: pthread create return code: %d\n", i, check);
								exit(EXIT_FAILURE);
							}
						}

						for(i=0; i<THREAD_NUM; i++)
						{
							pthread_join(myThread[i], NULL);
						}
						free(myThread);
						free(tt);

					}
#else
					for(i = 0; i < readNum_b[0]; i++)
					{
						startp = bufk + (i * kmer_num);
						end = (bufk_t[i] & 0Xf);

						if(((end & 0X3) == 2) || ((end & 0X3) == 3))	write_tmp_array[0] = 4;
						else	write_tmp_array[0] = (end >> 2);

						//
						//if(!OneUnit)
						if(kmer_num != 1)
							startp[0] >>= (tempMove_new + 2);

						for(j = 0, k = 1; j < KMER_LENGTH; j++, k++)
						{
							char_tmp = ((startp[(j + br_offm) >> 5] >> ((31 - ((j + br_offm) & 0X1f)) << 1)) & 0X3);

							write_tmp_array[k >> 4] <<= 4;
							if(char_tmp)	write_tmp_array[k >> 4] |= char_tmp;

						}

						if(k & 0Xf)	write_tmp_array[k >> 4] <<= ((16 - (k & 0Xf)) << 2);

						memcpy(bufK2_b[0] + i * wta_n, write_tmp_array, wta_n << 3);
					}
#endif
				}
				else	readNum_b[num]=fread(bufK2_b[num], wta_n << 3, bufferSize, fpK2_b[num]);

				pBuf[num]=0;

				if(readNum_b[num]==0)
				{
					fclose(fpK2_b[num]);

					num=heap[heapTail];
					heapTail--;
				}
			}

			for(i=0; (i<<1)+1<=heapTail; i=temp)
			{
				cdi1=(i<<1)+2,cdi2=cdi1-1;

				temp=(cdi1<=heapTail&&L_cmpMK2(bufK2_b[heap[cdi1]] + pBuf[heap[cdi1]] * wta_n, bufK2_b[heap[cdi2]] + pBuf[heap[cdi2]] * wta_n))?cdi1:cdi2;

				heap_tmp = heap[temp];
				if(heap_tmp == num)	heap[i] = num;
				else
				{
					if(L_cmpMK2(bufK2_b[heap_tmp] + pBuf[heap_tmp] * wta_n ,bufK2_b[num] + pBuf[num] * wta_n))
					{
						heap[i]=heap_tmp;
					}
					else break;
				}

			}
			heap[i]=num;
		}
	}

	if(seq_c_n)
	{
		if(seq_c_n & 0Xf)
			fwrite(bwt_seq_buffer, 8, ((seq_c_n >> 4) + 1), fp_bwt_seq);
		else	fwrite(bwt_seq_buffer, 8, (seq_c_n >> 4), fp_bwt_seq);
	}

	seq_n = (seq_n * max_mwa_l_char) + seq_c_n;
	fseek(fp_bwt_seq, 0, SEEK_SET);
	fwrite(&seq_n, 8, 1, fp_bwt_seq);


#ifdef	GFA_COM
	if(seq_e_n)
	{
		if(seq_e_n & 0X7)
			fwrite(edges_seq_buffer, 8, ((seq_e_n >> 3) + 1), fp_edges_seq);
		else	fwrite(edges_seq_buffer, 8, (seq_e_n >> 3), fp_edges_seq);
	}

	seq_n_e = (seq_n_e * max_edges_char) + seq_e_n;
	fseek(fp_edges_seq, 0, SEEK_SET);
	fwrite(&seq_n_e, 8, 1, fp_edges_seq);


#endif

	gettimeofday(&t_end1, NULL);
	end1 = ((long)t_end1.tv_sec)*1000+(long)t_end1.tv_usec/1000;

	cost_time = end1 - start1;
	fprintf(stderr, "Last merging and outputing BWT seq time is %.2lf s %"PRId64" %"PRId64"\n", cost_time/1000, seq_n, seq_c_n);

	gettimeofday(&t_start1, NULL);
	start1 = ((long)t_start1.tv_sec)*1000+(long)t_start1.tv_usec/1000;

	if(bwt_seq_buffer)	free(bwt_seq_buffer);
	if(heap)	free(heap);
	if(pBuf)	free(pBuf);

	for(num=0; num<L_NUM; num++)
		if(bufK2_b[num])	free(bufK2_b[num]);
	if(bufK2_b)	free(bufK2_b);

	if(fp_bwt_seq)	fclose(fp_bwt_seq);
	if(bufk)	free(bufk);
	if(bufk_t)	free(bufk_t);

#ifdef	GFA_COM
	if(bufk_e_t)	free(bufk_e_t);
	if(edges_seq_buffer)	free(edges_seq_buffer);
#endif

	//if(fpK2_b)	fclose(fpK2_b);
	if(fpK2_b_t)	fclose(fpK2_b_t);
#ifdef	GFA_COM
	if(fpK2_e_t)	fclose(fpK2_e_t);
	if(fp_edges_seq)	fclose(fp_edges_seq);
#endif


	for(num=1; num<L_NUM; num++)
	{
		sprintf(cNum,"%u",num);
		strcpy(string_tmp, "rm -f ");
		strcat(string_tmp, tmp_route);
		strcat(string_tmp, "/last_kmer.");
		strcat(string_tmp, cNum);
		system(string_tmp);
	}
}


void *last_multiCount(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	uint64_t num = argT->num;
	uint64_t THREAD_NUM = argT->thread;

	uint8_t b_i = 0;
	uint16_t w_i = 0;
	int j = 0, k = 0;
	uint64_t segment = (read_num / THREAD_NUM);
	uint64_t start_p = (num * segment);
	uint64_t end_p = 0;

	uint64_t i,i_in,tempSeq;
	uint64_t well_tmp = 4;
	uint64_t char_tmp = 0;
	uint64_t w_cnt = 0;

	if(num < THREAD_NUM - 1)
	{
		end_p = ((num + 1) * segment);
	}
	else
	{
		end_p = read_num;
	}

	uint64_t tmp = BUCKET_CAPACITY;

	last_countKmer_dis[num] = (uint64_t* )calloc(BUCKET_CAPACITY+1,sizeof(uint64_t ));
	uint64_t* write_tmp_array = calloc(kmer_num_br, 8);

	uint64_t kmer_num_br_tmp = kmer_num_br << 3;
	uint64_t* write_tmp = (uint64_t* )calloc((end_p - start_p) * KMER_LENGTH, kmer_num_br_tmp);

	if(write_tmp == NULL)	fprintf(stderr, "Fail to allocate memory for write\n"), exit(1);

	if(kmer_num == 1)
	{
		uint64_t start = 0, end = 0;

		memset(write_tmp_array, 0, kmer_num_br_tmp);

		for(i_in = start_p; i_in < end_p; i_in++)
		{
			start = first_buffer[i_in];
			end = second_buffer[i_in];
			start <<= 2;
			end <<= 2;

			for(b_i = 0; b_i < KMER_LENGTH; b_i++)   //add prefix node
			{
				for(i = b_i, w_i = 0; i < KMER_LENGTH; i++)
				{
					char_tmp = ((start >> ((31 - i) << 1)) & 0X3);

					write_tmp_array[w_i >> 4] <<= 4;
					if(char_tmp)	write_tmp_array[w_i >> 4] |= char_tmp;

					w_i++;
				}

				write_tmp_array[w_i >> 4] <<= 4;
				write_tmp_array[w_i >> 4] |= well_tmp;
				w_i++;

				for(j = 0; j < KMER_LENGTH; j++)   //b_i
				{
					char_tmp = ((end >> ((31 - j) << 1)) & 0X3);

					write_tmp_array[w_i >> 4] <<= 4;
					if(char_tmp)	write_tmp_array[w_i >> 4] |= char_tmp;

					w_i++;

				}

				if(w_i & 0Xf)	write_tmp_array[w_i >> 4] <<= ((16 - (w_i & 0Xf)) << 2);

				memcpy(write_tmp + w_cnt * kmer_num_br, write_tmp_array, kmer_num_br_tmp);//kmer_num_br << 3

				w_cnt++;

				last_countKmer_dis[num][(write_tmp_array[0]>>l_tempMove)&l_maskBKT]++;
			}
		}
	}
	else
	{
		uint64_t* start = NULL;
		uint64_t* end = NULL;

		for(i_in = start_p; i_in < end_p; i_in++)
		{
			start = first_buffer + (i_in * kmer_num);
			end = second_buffer + (i_in * kmer_num);

			for(b_i = 0; b_i < KMER_LENGTH; b_i++)   //add prefix node
			{
				for(i = b_i, w_i = 0; i < KMER_LENGTH; i++)
				{
					char_tmp = ((start[(i + br_offm) >> 5] >> ((31 - ((i + br_offm) & 0X1f)) << 1)) & 0X3);

					write_tmp_array[w_i >> 4] <<= 4;
					if(char_tmp)	write_tmp_array[w_i >> 4] |= char_tmp;

					w_i++;
				}

				write_tmp_array[w_i >> 4] <<= 4;
				write_tmp_array[w_i >> 4] |= well_tmp;

				w_i++;

				for(j = 0; j < KMER_LENGTH; j++)   //b_i
				{
					char_tmp = ((end[(j + br_offm) >> 5] >> ((31 - ((j + br_offm) & 0X1f)) << 1)) & 0X3);

					write_tmp_array[w_i >> 4] <<= 4;
					if(char_tmp)	write_tmp_array[w_i >> 4] |= char_tmp;

					w_i++;
				}

				if(w_i & 0Xf)	write_tmp_array[w_i >> 4] <<= ((16 - (w_i & 0Xf)) << 2);

				memcpy(write_tmp + w_cnt * kmer_num_br, write_tmp_array, kmer_num_br << 3);
				w_cnt++;

				last_countKmer_dis[num][(write_tmp_array[0]>>l_tempMove)&l_maskBKT]++;
			}
		}
	}

	if(write_tmp_array)	free(write_tmp_array);

	for(i=1; i<BUCKET_CAPACITY; i++)
		last_countKmer_dis[num][i] = last_countKmer_dis[num][i-1] + last_countKmer_dis[num][i];

	last_countKmer_dis[num][BUCKET_CAPACITY] = last_countKmer_dis[num][BUCKET_CAPACITY-1];

	last_hashKmer_dis[num] = (uint64_t* )calloc(last_countKmer_dis[num][BUCKET_CAPACITY], kmer_num_br_tmp);//kmer_num_br << 3

	for(i = 0; i < w_cnt; i++)
	{
		tempSeq = ((write_tmp[i * kmer_num_br]>>l_tempMove)&l_maskBKT);

		last_countKmer_dis[num][tempSeq]--;

		memcpy(last_hashKmer_dis[num] + last_countKmer_dis[num][tempSeq] * kmer_num_br, write_tmp + i * kmer_num_br, kmer_num_br_tmp);//kmer_num_br << 3

	}

	if(write_tmp)	free(write_tmp);

	return (void*)NULL;
}


void *last_multiThreadSort(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	uint64_t num=argT->num;
	uint64_t i, low, up;
	uint64_t kmer_num_br_char = kmer_num_br << 3;

	low = l_segCount[num];
	up = l_segCount[num+1];

	if(last_countKmer_dis[num])	free(last_countKmer_dis[num]);
	if(last_hashKmer_dis[num])	free(last_hashKmer_dis[num]);

	for(i = low; i < up; i++)
		qsort(last_hashKmer + last_countKmer[i] * kmer_num_br, last_countKmer[i+1] - last_countKmer[i], kmer_num_br_char, L_cmpKmers);

#ifdef	LAST_HASH_DIS
	uint64_t wta_n_char = wta_n << 3, j = 0;
	for(i = low; i < up; i++)
		for(j = last_countKmer[i]; j < last_countKmer[i+1]; j++)
			memcpy(last_hashKmer_write + (j * wta_n), last_hashKmer + (j * kmer_num_br), wta_n_char);
#endif

	return (void *)NULL;
}

void *twobittofourbit(void *arg)
{
	thread_branch* argT=(thread_branch* )arg;
	int num = argT->num;
	uint64_t wta_n = argT->wta_n;

	uint64_t start=splitMK2[num], end=splitMK2[num+1];
	uint64_t bufMK1_len=end-start;

	if(bufMK1_len==0) return (void *)NULL;

	uint64_t* write_tmp_array = calloc(kmer_num_br, 8);

	memset(write_tmp_array, 0, kmer_num_br << 3);

	uint64_t i = 0;
	uint64_t char_tmp = 0;
	int j = 0;
	int k = 0;
	uint64_t* startp = NULL;
	uint64_t endp = 0;
	uint64_t* readp = NULL;
	uint8_t* readp_t = NULL;
	uint64_t* writep = NULL;

	readp = argT->readp;
	readp_t = argT->readp_t;
	writep = argT->writep;

	for(i = 0; i < bufMK1_len; i++)
	{
		startp = readp + i * kmer_num;
		endp = (readp_t[i] & 0Xf);

		if(((endp & 0X3) == 2) || ((endp & 0X3) == 3))	write_tmp_array[0] = 4;
		else	write_tmp_array[0] = (endp >> 2);

		if(kmer_num != 1)
			startp[0] >>= (tempMove_new + 2);

		for(j = 0, k = 1; j < KMER_LENGTH; j++, k++)
		{
			char_tmp = ((startp[(j + br_offm) >> 5] >> ((31 - ((j + br_offm) & 0X1f)) << 1)) & 0X3);

			write_tmp_array[k >> 4] <<= 4;
			if(char_tmp)	write_tmp_array[k >> 4] |= char_tmp;

		}

		if(k & 0Xf)	write_tmp_array[k >> 4] <<= ((16 - (k & 0Xf)) << 2);

		memcpy(writep + i * wta_n, write_tmp_array, wta_n << 3);
	}

	if(write_tmp_array)	free(write_tmp_array);

	return (void *)NULL;
}


int L_cmpKmer(const void *a, const void *b)
{
	kmer_s* vav = (kmer_s* )a;
	kmer_s* vbv = (kmer_s* )b;

	uint64_t va = 0;
	uint64_t vb = 0;
	uint8_t i = 0;

	va = vav->seq[0];
	vb = vbv->seq[0];

	va = va << 4;
	vb = vb << 4;

	if(va < vb)	return -1;
	else if(va > vb)	return 1;
	else
	{
		for(i = 1; i < MAX_KMER_SEQ_BIT; i++)
		{
			va = vav->seq[i];
			vb = vbv->seq[i];

			if(va < vb)	return -1;
			else if(va > vb)	return 1;
		}
	}

	return 0;
}

int L_cmpKmers(const void *a, const void *b)
{

	uint64_t* vav = (uint64_t* )a;
	uint64_t* vbv = (uint64_t* )b;

	uint64_t va = 0;
	uint64_t vb = 0;
	uint8_t i = 0;

	va = (*vav);
	vb = (*vbv);

	va <<= 4;
	vb <<= 4;

	if(va < vb)	return -1;
	else if(va > vb)	return 1;
	else
	{
		for(i = 1; i < kmer_num_br; i++)   //MAX_KMER_SEQ_BIT
		{
			va = vav[i];
			vb = vbv[i];

			if(va < vb)	return -1;
			else if(va > vb)	return 1;
		}
	}


	return 0;
}

inline int cmpMK2(uint64_t a, uint64_t b)
{
	a=a<<2, b=b<<2;
	if(a<b) return 1;
	else return 0;
}

inline int cmpMK2_s(uint64_t* va, uint64_t* vb)
{
	uint64_t a = (*va);
	uint64_t b = (*vb);

	if(OneUnit)
	{
#ifdef	FIRST_LAST_KMER
		a=a<<4, b=b<<4;
#else
		a=a<<2, b=b<<2;
#endif

		if(a<b) return 1;
		else return 0;
	}
	else
	{
		a <<= (tempMove_new + 2);
		b <<= (tempMove_new + 2);

		if(a < b)	return 1;
		else if(a > b)	return 0;
		else
		{
			if(kmer_num > 1)
			{
				uint8_t i = 0;

				for(i = 1; i < kmer_num - 1; i++)
				{
					a = (*(va + i));
					b = (*(vb + i));

					if(a<b) return 1;
					else if(a>b) return 0;
				}
				if(kmer_num == kmer_num_new)
				{
					a = (va[kmer_num - 1] >> 2) << 2;
					b = (vb[kmer_num - 1] >> 2) << 2;

					if(a<b) return 1;
					else	return 0;
				}
				else
				{
					if(va[i] < vb[i]) return 1;
					else	return 0;
				}
			}
			else	return 0;
		}
	}
}

inline int cmpMK2_sm(uint64_t* va, uint64_t* vb)
{
	uint64_t a = (*va);
	uint64_t b = (*vb);

	if(kmer_num == 1)
	{
		a=a<<2, b=b<<2;
		if(a<b) return 1;
		else return 0;
	}
	else
	{
		a <<= 2;
		b <<= 2;

		if(a < b)	return 1;
		else if(a > b)	return 0;
		else
		{
			uint8_t i = 0;

			for(i = 1; i < kmer_num; i++)
			{
				a = (*(va + i));
				b = (*(vb + i));
				if(a<b) return 1;
				else if(a>b) return 0;
			}
			return 0;
		}
	}
}

inline int cmpMK2_s_filter(uint64_t* va, uint64_t* vb)
{
	uint64_t a = (*va);
	uint64_t b = (*vb);

	//if(kmer_num == 1)
	if(OneUnit)
	{
		a=a<<2, b=b<<2;
		if(a<b) return 1;
		else return 0;
	}
	else
	{
		a <<= 2;
		b <<= 2;

		if(a < b)	return 1;
		else if(a > b)	return 0;
		else
		{
			if(kmer_num > 1)
			{
				uint8_t i = 0;

				for(i = 1; i < kmer_num - 1; i++)
				{
					a = (*(va + i));
					b = (*(vb + i));

					if(a<b) return 1;
					else if(a>b) return 0;
				}
				if(kmer_num == kmer_num_new)
				{
					a = (va[kmer_num - 1] >> 8) << 8;
					b = (vb[kmer_num - 1] >> 8) << 8;

					if(a<b) return 1;
					else	return 0;
				}
				else
				{
					if(va[i] < vb[i]) return 1;
					else	return 0;
				}
			}
			else	return 0;
		}
	}
}

inline int L_cmpMK2(uint64_t* a, uint64_t* b)
{
	uint64_t va = (*a);
	uint64_t vb = (*b);
	uint8_t i = 0;

	va = va << 4;
	vb = vb << 4;

	if(va < vb) return 1;
	else if(va == vb)
	{
		for(i=1; i<wta_n; i++)
		{
			va = (*(a+i));
			vb = (*(b+i));
			if(va < vb)	return 1;
			else if(va > vb)	return 0;
		}
		return 0;
	}
	else return 0;
}

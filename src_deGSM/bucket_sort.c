#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
//#include "collect#$.h"
//#include "INandOut.h"
#include "bucket_sort.h"
//#include "geUnipath.h"


void *multiThreadSort(void *arg);
void *multiThreadSort_new(void *arg);

int cmpKmer(const void * , const void * );
void *multiDistri(void * );
void *multiCount(void * );
void *multiDistri_lfkmers(void * );
void *multiCount_lfkmers(void * );
void *multiThreadSort_lfkmers(void *);
int cmpKmer_lfkmers(const void *, const void *);
int cmpMK2_s_lfkmers(uint64_t* , uint64_t* );
void *multiTagcheck_lfkmers(void *);
void multiTagcheck_lfkmers_addkmer(uint64_t , uint64_t , FILE* );

void bucket_sort_lfkmers(uint64_t, uint64_t, char*, char* );
void lfkmer_fa(char* );
void lfkmer_fq(char*,  char* , char* , char , uint8_t );

uint64_t *countKmer = NULL;

uint64_t trans[256];
uint64_t totalKmerNum;
uint64_t *segCount=NULL;
uint8_t* readBuf=NULL;//
uint64_t readNum;
uint64_t tempMove=0;
uint64_t maskBKT;
uint8_t jelly_f = 0;
uint64_t key_len = 0;
uint64_t pair_len = 0;
uint64_t readNum_new = 0;
uint64_t bufferSize_new = 0;

uint64_t* splitIN = NULL;
char* readBuf_new = NULL;
uint64_t** readBuf_p = NULL;
uint64_t* readBuf_p_num = NULL;

uint64_t** countKmer_dis = NULL;
uint64_t** hashKmer_dis = NULL;

uint64_t* hashKmer = NULL;

uint64_t* splitMK2_lf = NULL;
uint64_t* bufMK2_lf = NULL;
uint64_t w_lf_n[32];

//bit operation
uint8_t charToDna5[] =
{
	/*   0 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  16 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  32 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  48 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*  64 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 3, 0,
	/*    A     C           G                    N */
	/*  80 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 127, 0, 0, 0, 0, 0,//'Z' 6 127
	/*             T */
	/*  96 */ 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 3, 0,
	/*    a     c           g                    n */
	/* 112 */ 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/*             t */
	/* 128 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 144 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 160 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 176 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 192 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 208 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 224 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
	/* 240 */ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
};

uint8_t weirdn = 22;
char weirds[23] = {'N', 'n', 'R', 'r', 'M', 'm', 'S', 's', 'Y', 'y', 'K', 'k', 'W', 'w', 'B', 'b', 'D', 'd', 'H', 'h', 'V', 'v'};

uint8_t weirdn_fq = 22;
char weirds_fq[23] = {'N', 'n', 'R', 'r', 'M', 'm', 'S', 's', 'Y', 'y', 'K', 'k', 'W', 'w', 'B', 'b', 'D', 'd', 'H', 'h', 'V', 'v'};

#ifdef	FIRST_LAST_KMER

void FLkmeradd_F(char* input_kmer)
{
	uint8_t char_tmp = 0;
	uint16_t k = 0;
	uint8_t write_kmer[pair_len];
	int16_t fl_i = 0;

	memset(write_kmer, 0, key_len);
	for(fl_i = KMER_LENGTH, k = 0; fl_i >= 0; fl_i--, k++)
	{
		char_tmp = charToDna5[(uint8_t)input_kmer[fl_i]];
		if(char_tmp)	write_kmer[k >> 2] |= (char_tmp << ((k & 0X3) << 1));
	}

	write_kmer[pair_len - 1] = 1;

	fwrite(write_kmer, pair_len, 1, fp_lfk_l);

}
void FLkmeradd_BFQ(char* input_kmer)
{
	uint8_t char_tmp = 0;
	uint8_t k = 0;
	uint8_t write_kmer[pair_len];
	int16_t fl_i = 0;

	memset(write_kmer, 0, key_len);
	for(fl_i = KMER_LENGTH - 1, k = 1; fl_i >= 0; fl_i--, k++)
	{
		char_tmp = charToDna5[(uint8_t)input_kmer[fl_i]];
		if(char_tmp)	write_kmer[k >> 2] |= (char_tmp << ((k & 0X3) << 1));
	}

	write_kmer[pair_len - 1] = 3;

	fwrite(write_kmer, pair_len, 1, fp_lfk_b);

	bfq_n++;

}
void FLkmeradd_FFQ(char* input_kmer)
{
	int16_t fl_i = 0;
	char* pch_tmp = NULL;

	input_kmer[KMER_LENGTH_PlusOne] = '\0';
	for(fl_i = 0; fl_i < weirdn_fq; fl_i++)
	{
		if((pch_tmp = strchr(input_kmer, weirds_fq[fl_i])))
			break;
	}

	if(pch_tmp == NULL)
	{
		uint8_t char_tmp = 0;
		uint8_t k = 0;
		uint8_t write_kmer[pair_len];


		memset(write_kmer, 0, key_len);
		for(fl_i = KMER_LENGTH, k = 0; fl_i >= 0; fl_i--, k++)
		{
			char_tmp = charToDna5[(uint8_t)input_kmer[fl_i]];
			if(char_tmp)	write_kmer[k >> 2] |= (char_tmp << ((k & 0X3) << 1));

		}

		write_kmer[pair_len - 1] = 1;

		fwrite(write_kmer, pair_len, 1, fp_lfk_l);
	}
}

void FLkmeradd_L(char* input_kmer)
{
	uint8_t char_tmp = 0;
	uint16_t k = 0;
	uint8_t write_kmer[pair_len];
	int16_t fl_i = 0;

	memset(write_kmer, 0, key_len);
	for(fl_i = KMER_LENGTH, k = 1; fl_i >= 0; fl_i--, k++)
	{
		char_tmp = charToDna5[(uint8_t)input_kmer[fl_i]];
		if(char_tmp)	write_kmer[k >> 2] |= (char_tmp << ((k & 0X3) << 1));
	}

	write_kmer[pair_len - 1] = 2;

	fwrite(write_kmer, pair_len, 1, fp_lfk_r);
}

void FLkmeradd_LFQ(char* input_kmer)
{
	int16_t fl_i = 0;
	char* pch_tmp = NULL;
	input_kmer[KMER_LENGTH_PlusOne] = '\0';
	for(fl_i = 0; fl_i < weirdn_fq; fl_i++)
	{
		if((pch_tmp = strchr(input_kmer, weirds_fq[fl_i])))
			break;
	}

	if(pch_tmp == NULL)
	{
		uint8_t char_tmp = 0;
		uint8_t k = 0;
		uint8_t write_kmer[pair_len];

		memset(write_kmer, 0, key_len);
		for(fl_i = KMER_LENGTH, k = 1; fl_i >= 0; fl_i--, k++)
		{
			char_tmp = charToDna5[(uint8_t)input_kmer[fl_i]];
			if(char_tmp)	write_kmer[k >> 2] |= (char_tmp << ((k & 0X3) << 1));
		}

		write_kmer[pair_len - 1] = 2;

		fwrite(write_kmer, pair_len, 1, fp_lfk_r);
	}
}

void FLkmeradd_LFQOFF(char* input_kmer, uint16_t offset)
{
	uint8_t char_tmp = 0;
	uint8_t k = 0;
	uint8_t write_kmer[pair_len];
	int16_t fl_i = 0;

	memset(write_kmer, 0, key_len);
	for(fl_i = KMER_LENGTH + offset, k = 1; fl_i >= offset; fl_i--, k++)
	{
		char_tmp = charToDna5[(uint8_t)input_kmer[fl_i]];
		if(char_tmp)	write_kmer[k >> 2] |= (char_tmp << ((k & 0X3) << 1));
	}

	write_kmer[pair_len - 1] = 2;

	fwrite(write_kmer, pair_len, 1, fp_lfk_r);
}

uint8_t forward_kmeradd(char* fa_inputs, char* pch)
{
	int16_t num_tmp = 0;
	int16_t num_tmp_a = 0;
	int16_t num_tmp_b = 0;
	int16_t fl_i = 0;
	char input_kmer[KMER_LENGTH_PlusTwo + 2];
	char fa_inputs_tmp[MAX_FA_LINE + 1];
	char* pch_tmp = NULL;
	char* fa_in_p = NULL;

	num_tmp = strlen(fa_inputs) - 1 - (pch - fa_inputs);

	if(num_tmp < KMER_LENGTH_PlusOne)
	{
		num_tmp_a = ((KMER_LENGTH_PlusOne - num_tmp) / fa_line_len);
		num_tmp_b = ((KMER_LENGTH_PlusOne - num_tmp) % fa_line_len);

		strncpy(input_kmer, pch, num_tmp);

		for(fl_i = 0; fl_i < num_tmp_a; fl_i++)
		{
			fa_in_p = fgets(fa_inputs_tmp, MAX_FA_LINE, fp_fa);

			if(!fa_in_p)	return 0;
			else if((fa_inputs_tmp[0] == '>') || (strlen(fa_inputs_tmp) - 1 < fa_line_len))
			{
				fseek(fp_fa, -(((fa_line_len + 1) * fl_i) + strlen(fa_inputs_tmp)), SEEK_CUR);
				return 0;
			}

			strncpy(input_kmer + num_tmp, fa_inputs_tmp, fa_line_len);
			num_tmp += fa_line_len;
		}

		if(num_tmp_b)
		{
			fa_in_p = fgets(fa_inputs_tmp, MAX_FA_LINE, fp_fa);

			if(!fa_in_p)	return 0;
			else if((fa_inputs_tmp[0] == '>') || (strlen(fa_inputs_tmp) - 1 < num_tmp_b))
			{
				fseek(fp_fa, -(((fa_line_len + 1) * fl_i) + strlen(fa_inputs_tmp)), SEEK_CUR);
				return 0;
			}

			strncpy(input_kmer + num_tmp, fa_inputs_tmp, num_tmp_b);
		}

		fseek(fp_fa, -(((fa_line_len + 1) * num_tmp_a) + (num_tmp_b == 0 ? 0:strlen(fa_inputs_tmp))), SEEK_CUR);
	}
	else
	{
		strncpy(input_kmer, pch, KMER_LENGTH_PlusOne);
	}

	input_kmer[KMER_LENGTH_PlusOne] = '\0';

	pch_tmp = NULL;
	for(fl_i = 0; fl_i < weirdn; fl_i++)
	{
		if((pch_tmp = strchr(input_kmer, weirds[fl_i])))
			break;
	}

	if(pch_tmp == NULL)
		FLkmeradd_F(input_kmer);
	else if(flag_filter)
	{
		char ter_end = input_kmer[KMER_LENGTH];
		char ter_start = input_kmer[0];
		if((ter_start == 'N') || (ter_end == 'N'))
		{
			char input_kmerb[KMER_LENGTH + 1];
			if(ter_start == 'N')
			{
				strncpy(input_kmerb, input_kmer + 1, KMER_LENGTH);
			}
			else if(ter_end == 'N')
			{
				strncpy(input_kmerb, input_kmer, KMER_LENGTH);
			}
			input_kmerb[KMER_LENGTH] = '\0';

			for(fl_i = 0; fl_i < weirdn; fl_i++)
			{
				if((pch_tmp = strchr(input_kmerb, weirds[fl_i])))
					break;
			}
			if(pch_tmp == NULL)	FLkmeradd_BFQ(input_kmerb);
		}
	}

	return 1;
}

uint8_t reverse_kmerdadd(char* fa_inputs, char* pch, int16_t fap_off)
{
	int16_t num_tmp = 0;
	int16_t num_tmp_a = 0;
	int16_t num_tmp_b = 0;
	int16_t fl_i = 0;

	char input_kmer[KMER_LENGTH_PlusTwo + 2];
	char fa_inputs_tmp[MAX_FA_LINE + 1];

	char* pch_tmp = NULL;
	char* pch_tmp_tmp = NULL;

	num_tmp = pch - fa_inputs;

	if(num_tmp < KMER_LENGTH_PlusOne)
	{
		num_tmp_a = ((KMER_LENGTH_PlusOne - num_tmp) / fa_line_len);
		num_tmp_b = ((KMER_LENGTH_PlusOne - num_tmp) % fa_line_len);

		fl_i = (num_tmp_a + (num_tmp_b == 0 ? 0:1));

		if(fl_i >= line_c)	return 0;

		fseek(fp_fa, -fap_off, SEEK_CUR);
		fseek(fp_fa, -((fa_line_len + 1) * fl_i), SEEK_CUR);//jump to start of a previous line

		if(num_tmp_b)
		{
			fgets(fa_inputs_tmp, MAX_FA_LINE, fp_fa);
			strncpy(input_kmer, fa_inputs_tmp + fa_line_len - num_tmp_b, num_tmp_b);
		}

		for(fl_i = 0; fl_i < num_tmp_a; fl_i++)
		{
			fgets(fa_inputs_tmp, MAX_FA_LINE, fp_fa);
			strncpy(input_kmer + num_tmp_b, fa_inputs_tmp, fa_line_len);
			num_tmp_b += fa_line_len;
		}

		strncpy(input_kmer + num_tmp_b, fa_inputs, num_tmp);

		fgets(fa_inputs_tmp, MAX_FA_LINE, fp_fa);
	}
	else
	{
		strncpy(input_kmer, pch - KMER_LENGTH_PlusOne, KMER_LENGTH_PlusOne);
	}

	input_kmer[KMER_LENGTH_PlusOne] = '\0';

	pch_tmp = NULL;
	for(fl_i = 0; fl_i < weirdn; fl_i++)
	{
		if((pch_tmp = strchr(input_kmer, weirds[fl_i])))
			break;
	}

	if(pch_tmp == NULL)
		FLkmeradd_L(input_kmer);
	else if(flag_filter)
	{
		char ter_end = input_kmer[KMER_LENGTH];
		char ter_start = input_kmer[0];
		if((ter_start == 'N') || (ter_end == 'N'))
		{
			char input_kmerb[KMER_LENGTH + 1];
			if(ter_start == 'N')
			{
				strncpy(input_kmerb, input_kmer + 1, KMER_LENGTH);
			}
			else if(ter_end == 'N')
			{
				strncpy(input_kmerb, input_kmer, KMER_LENGTH);
			}
			input_kmerb[KMER_LENGTH] = '\0';

			for(fl_i = 0; fl_i < weirdn; fl_i++)
			{
				if((pch_tmp = strchr(input_kmerb, weirds[fl_i])))
					break;
			}
			if(pch_tmp == NULL)	FLkmeradd_BFQ(input_kmerb);
		}
	}

	return 1;
}

#endif


void lfkmer_fa(char* fa_inputs)
{
	char* pch = NULL;
	char* pch_tmp = NULL;
	char* pch_tmp_tmp = NULL;
	char char_tmp_f, char_tmp_l;
	char char_tmp, char_tmp_le, char_tmp_ri;

	uint8_t first_flag = 1;
	uint8_t more_flag = 0;
	uint8_t allN_flag = 0;
	uint16_t min_pch = 0;
	uint8_t go_flag = 0;
	uint8_t k = 0;
	uint16_t n_i = 0;
	int16_t num_tmp = 0;
	int16_t num_tmp_more = 0;

	fgets(fa_inputs, MAX_FA_LINE, fp_fa);
	fgets(fa_inputs, MAX_FA_LINE, fp_fa);

	fa_line_len = strlen(fa_inputs) - 1;

	fseek(fp_fa, 0, SEEK_SET);

	line_c = 0;
	while(!feof(fp_fa))
	{
		pch = fgets(fa_inputs, MAX_FA_LINE, fp_fa);
		if(pch == NULL)	break;

		num_tmp_more = strlen(fa_inputs);

		char_tmp_f = fa_inputs[0];
		char_tmp_l = fa_inputs[num_tmp_more - 2];

		if(line_c == 0)	fa_line_len = num_tmp_more - 1;

		line_c++;

		if(char_tmp_f == '>')
		{
			if((!first_flag) && (!allN_flag))
			{
				line_c--;
				num_tmp = num_tmp_more + fa_line_len + 3;//
				fseek(fp_fa, -num_tmp, SEEK_CUR);
				fgets(fa_inputs, MAX_FA_LINE, fp_fa);
				fgets(fa_inputs, MAX_FA_LINE, fp_fa);
				pch = fa_inputs + strlen(fa_inputs) - 1;

				reverse_kmerdadd(fa_inputs, pch, (pch - fa_inputs) + 1);

				fseek(fp_fa, num_tmp_more, SEEK_CUR);
			}
			first_flag = 0;
			allN_flag = 0;
			more_flag = 1;
			line_c = 0;
			continue;
		}

		if((char_tmp_f == 'N') && (char_tmp_l == 'N'))
		{
			for(n_i = 0; n_i < num_tmp_more - 1; n_i++)
			{
				char_tmp = fa_inputs[n_i];

				if((char_tmp == 'A') || (char_tmp == 'C') || (char_tmp == 'G') || (char_tmp == 'T') || (char_tmp == 'a') || (char_tmp == 'c') || (char_tmp == 'g') || (char_tmp == 't'))//may need to be modified
				{
					char_tmp_le = fa_inputs[n_i - 1];
					if((!((char_tmp_le == 'A') || (char_tmp_le == 'C') || (char_tmp_le == 'G') || (char_tmp_le == 'T') || (char_tmp_le == 'a') || (char_tmp_le == 'c') || (char_tmp_le == 'g') || (char_tmp_le == 't'))) && (num_tmp_more - 1 - n_i > KMER_LENGTH))
						forward_kmeradd(fa_inputs, fa_inputs + n_i);
					char_tmp_ri = fa_inputs[n_i + 1];
					if((!((char_tmp_ri == 'A') || (char_tmp_ri == 'C') || (char_tmp_ri == 'G') || (char_tmp_ri == 'T') || (char_tmp_ri == 'a') || (char_tmp_ri == 'c') || (char_tmp_ri == 'g') || (char_tmp_ri == 't'))) && (n_i > KMER_LENGTH))
						reverse_kmerdadd(fa_inputs, fa_inputs + n_i + 1, num_tmp_more);
				}
			}

			if((!allN_flag) && (!more_flag))
			{
				num_tmp = num_tmp_more + fa_line_len + 3;
				fseek(fp_fa, -num_tmp, SEEK_CUR);
				fgets(fa_inputs, MAX_FA_LINE, fp_fa);
				fgets(fa_inputs, MAX_FA_LINE, fp_fa);
				pch = fa_inputs + strlen(fa_inputs) - 1;

				reverse_kmerdadd(fa_inputs, pch, (pch - fa_inputs) + 1);

				fseek(fp_fa, num_tmp_more, SEEK_CUR);
			}
			allN_flag = 1;
			more_flag = 0;
			continue;
		}

		if((more_flag) || (allN_flag))
		{
			forward_kmeradd(fa_inputs, fa_inputs);
			more_flag = 0;
			allN_flag = 0;
		}

		pch = NULL;
		min_pch = MAX_FA_LINE;
		for(k = 0; k < weirdn; k++)
		{
			pch_tmp = strchr(fa_inputs, weirds[k]);
			if(pch_tmp)
			{
				if(pch_tmp - fa_inputs < min_pch)
				{
					min_pch = pch_tmp - fa_inputs;
					pch = pch_tmp;
				}
			}
		}

		if(pch)
		{
			while(pch)
			{
				reverse_kmerdadd(fa_inputs, pch, num_tmp_more);

				++pch;

				while(1)
				{
					go_flag = 0;
					for(k = 0; k < weirdn; k++)
						if((*pch) == weirds[k])
						{
							go_flag = 1;
							break;
						}

					if(go_flag)	pch++;
					else	break;
				}

				forward_kmeradd(fa_inputs, pch);

				if(pch - fa_inputs == num_tmp_more - 1)
					break;

				pch_tmp = NULL;
				min_pch = MAX_FA_LINE;
				for(k = 0; k < weirdn; k++)
				{
					pch_tmp_tmp = strchr(pch, weirds[k]);
					if(pch_tmp_tmp)
					{
						if(pch_tmp_tmp - pch < min_pch)
						{
							min_pch = pch_tmp_tmp - pch;
							pch_tmp = pch_tmp_tmp;
						}
					}
				}

				pch = pch_tmp;
			}
			allN_flag = 0;
			more_flag = 0;
		}
	}

	num_tmp_more = strlen(fa_inputs);
	pch = fa_inputs + num_tmp_more - 1;

	reverse_kmerdadd(fa_inputs, pch, num_tmp_more);

}

void lfkmer_fq(char* fa_inputs, char* fq_inputs, char* fq_quality, char quality_char, uint8_t quality_flag)
{
	uint8_t go_flag = 0;
	uint8_t w_flag = 0;
	uint8_t k = 0;

	int16_t num_tmp_more = 0;
	int min_pch = 0;
	uint64_t key_bits = 0;
	uint8_t line_tmp = 0;

	char* pch = NULL;
	char* pch_tmp = NULL;
	char* pch_tmp_tmp = NULL;
	char input_kmer[KMER_LENGTH_PlusTwo + 2];

	char* pch_pre = NULL;
	int16_t fl_i = 0;

	while(!feof(fp_fa))
	{
		fgets(fq_inputs, MAX_FQ_LINE, fp_fa);

		line_tmp = (key_bits & 0X3);
		if(line_tmp == 1)
		{
			strcpy(fa_inputs, fq_inputs);
			num_tmp_more = strlen(fa_inputs) - 1;
		}
		else if(line_tmp == 3)
		{
			if(quality_flag)
			{
				strcpy(fq_quality, fq_inputs);
				for(fl_i = 0; fl_i < num_tmp_more; fl_i++)
					if(fq_quality[fl_i] < quality_char)
						fa_inputs[fl_i] = 'N';
			}

			if(num_tmp_more < KMER_LENGTH_PlusOne)
			{
				if((num_tmp_more == KMER_LENGTH) && flag_filter)
				{
					strncpy(input_kmer, fa_inputs, KMER_LENGTH);

					input_kmer[KMER_LENGTH] = '\0';
					for(fl_i = 0; fl_i < weirdn_fq; fl_i++)
					{
						if((pch_tmp = strchr(input_kmer, weirds_fq[fl_i])))
							break;
					}
					if(pch_tmp == NULL)	FLkmeradd_BFQ(input_kmer);
				}

				key_bits++;
				continue;
			}

			//if((pch = strchr(fa_inputs, 'n')) || (pch = strchr(fa_inputs, 'N')))
			pch = NULL;
			min_pch = MAX_FQ_LINE;
			for(k = 0; k < weirdn_fq; k++)
			{
				pch_tmp = strchr(fa_inputs, weirds_fq[k]);
				if(pch_tmp)
				{
					if(pch_tmp - fa_inputs < min_pch)
					{
						min_pch = pch_tmp - fa_inputs;
						pch = pch_tmp;
					}
				}
			}

			if(flag_filter)	pch_pre = fa_inputs;

			if(pch)
			{
				strncpy(input_kmer, fa_inputs, KMER_LENGTH_PlusOne);
				FLkmeradd_FFQ(input_kmer);

				w_flag = 1;
				while(pch)
				{
					if(pch - fa_inputs > KMER_LENGTH)
					{
						strncpy(input_kmer, pch - KMER_LENGTH_PlusOne, KMER_LENGTH_PlusOne);
						FLkmeradd_LFQ(input_kmer);
					}

					if((pch - pch_pre == KMER_LENGTH) && flag_filter)
					{
						strncpy(input_kmer, pch - KMER_LENGTH, KMER_LENGTH);
						input_kmer[KMER_LENGTH] = '\0';
						FLkmeradd_BFQ(input_kmer);
					}

					++pch;

					while(1)
					{
						go_flag = 0;
						for(k = 0; k < weirdn_fq; k++)
							if((*pch) == weirds_fq[k])
							{
								go_flag = 1;
								break;
							}

						if(go_flag)	pch++;
						else	break;
					}

					if(flag_filter)	pch_pre = pch;

					if(pch - fa_inputs < num_tmp_more - KMER_LENGTH)
					{
						strncpy(input_kmer, pch, KMER_LENGTH_PlusOne);
						FLkmeradd_FFQ(input_kmer);

						pch_tmp = NULL;
						min_pch = MAX_FQ_LINE;
						for(k = 0; k < weirdn_fq; k++)
						{
							pch_tmp_tmp = strchr(pch, weirds_fq[k]);
							if(pch_tmp_tmp)
							{
								if(pch_tmp_tmp - pch < min_pch)
								{
									min_pch = pch_tmp_tmp - pch;
									pch_tmp = pch_tmp_tmp;
								}
							}
						}

						pch = pch_tmp;
					}
					else if((pch - fa_inputs == num_tmp_more - KMER_LENGTH) && flag_filter)
					{
						strncpy(input_kmer, pch, KMER_LENGTH);
						input_kmer[KMER_LENGTH] = '\0';
						for(fl_i = 0; fl_i < weirdn_fq; fl_i++)
						{
							if((pch_tmp = strchr(input_kmer, weirds_fq[fl_i])))
								break;
						}
						if(pch_tmp == NULL)	FLkmeradd_BFQ(input_kmer);

						w_flag = 0;
						break;
					}
					else
					{
						w_flag = 0;
						break;
					}
				}

				if(w_flag)
				{
					strncpy(input_kmer, fa_inputs + num_tmp_more - KMER_LENGTH_PlusOne, KMER_LENGTH_PlusOne);
					FLkmeradd_LFQ(input_kmer);
				}
			}
			else
			{
				FLkmeradd_F(fa_inputs);
				FLkmeradd_LFQOFF(fa_inputs, num_tmp_more - KMER_LENGTH_PlusOne);
			}
		}
		key_bits++;
	}
}

void *bucket_sort(void* arg)
{
	void **argT=(void **)arg;
	char *bin=(char *)argT[0];
	uint64_t THREAD_NUM=(uint64_t )argT[1];
	char* source = (char *)argT[2];
	char* tmp_route = (char *)argT[3];

	trans['A']=trans['a']=0;
	trans['C']=trans['c']=1;
	trans['N']=trans['G']=trans['g']=2;
	trans['T']=trans['t']=3;
	trans['#']=4;// patch for special characters
	trans['$']=5;

	uint64_t i = 0;
	int j = 0;
	char segKmerName[300];
	char cNum[6];

	//////////////////////////////////////kmercounting///////////////////////////////////////////
	uint64_t bufferSize=BUFFERSIZE<<1;//1<<22

	tempMove=((KMER_LENGTH_PlusTwo - (BUCKET_LENGTH+1))<<1);

	maskBKT=(1<<(BUCKET_LENGTH<<1))-1;

	pthread_t* myThread = NULL;

	uint64_t key_bits = 0;
	uint64_t val_len = 0;
	uint64_t key_ct = 0;
	uint64_t k = 0;
	uint64_t header_size = 0;
	uint64_t segStart=0, segment = 0;
	uint64_t loopNum=0;
	uint64_t seq = 0;
	uint64_t locateNum = 0;
	uint64_t locate = 0;
	uint64_t i_ini = 0;

	FILE* fpKmer_w = NULL;
	strcpy(segKmerName, tmp_route);
	strcat(segKmerName, "/kmerInfo");


#ifdef	KEEP_TMP
	FILE* fpKmer = fopen(segKmerName,"rb");//kmerPath

	char* readBuf_head = (char *)calloc(64,sizeof(char));

	fread(readBuf_head,sizeof(char),64 ,fpKmer);

	//back
	/**/
	char *DATABASE_FILE_TYPE = "JFLISTDN";
	if (strncmp(readBuf_head, DATABASE_FILE_TYPE, strlen(DATABASE_FILE_TYPE)))
	{
		jelly_f = 1;
		fseek(fpKmer,0,SEEK_SET);
	}

	if(jelly_f)//version 2.x
	{
		int buff_n = 0;
		fread(readBuf_head, 1, 9, fpKmer);
		buff_n = atoi(readBuf_head);

		key_bits = (KMER_LENGTH_PlusTwo << 1);
		key_len = key_bits / 8 + !! (key_bits % 8);
		val_len = 4;
		header_size = buff_n + 9;
	}
	else  //version 1.x
	{
		memcpy(&key_bits, readBuf_head + 8, 8);
		memcpy(&val_len, readBuf_head + 16, 8);
		memcpy(&key_ct, readBuf_head + 48, 8);

		//k = key_bits / 2;
		key_len = key_bits / 8 + !! (key_bits % 8);

		header_size=72+2*(4+8*key_bits);
	}
	free(readBuf_head);

	fseek(fpKmer,header_size,SEEK_SET);

#endif

	pair_len=key_len+val_len;

	KMER_LENGTH_PlusOne = KMER_LENGTH + 1;

	if(KMER_LENGTH_PlusTwo < 32)	OneUnit = 1;
	else	OneUnit = 0;

	kmer_num = (KMER_LENGTH_PlusTwo >> 5);
	kmer_rem = (KMER_LENGTH_PlusTwo & 0X1f);
	if(kmer_rem)	kmer_num++;

	kmer_num_new = kmer_num;


#ifdef	JELLY_LFKMER

#ifdef	FIRST_LAST_KMER_EXTRACT

	strcpy(segKmerName, tmp_route);
	strcat(segKmerName, "/lfInfo_l");
	fp_lfk_l = fopen(segKmerName, "wb");
	strcpy(segKmerName, tmp_route);
	strcat(segKmerName, "/lfInfo_r");
	fp_lfk_r = fopen(segKmerName, "wb");
	strcpy(segKmerName, tmp_route);
	strcat(segKmerName, "/lfInfo_b");
	fp_lfk_b = fopen(segKmerName, "wb");
	bfq_n = 0;

	fprintf(stderr, "Begin loading start and end kmers for every chr\n");

	if(quality_flag)
		fprintf(stderr, "Do quality filtering for characters in input file\n");
	i = MAX_FQ_LINE;
	fprintf(stderr, "The maximal length of line in input file should not be more than %"PRId64"\n", i);

	char fa_inputs[MAX_FQ_LINE];
	char fq_inputs[MAX_FQ_LINE];
	char fq_quality[MAX_FQ_LINE];

	struct stat buf;
	if(lstat(source, &buf) < 0)
	{
		fprintf(stderr, "lstat error for %s\n", source);
		exit(1);
	}

	if(S_ISDIR(buf.st_mode))
	{
		DIR *directory_pointer;

		fprintf(stderr, "%s is DIR\n", source);

		if((directory_pointer = opendir(source))==NULL)
			fprintf(stderr, "Error opening %s\n", source);
		else
		{
			fprintf(stderr, "%s is DIR\n", source);

			struct dirent *entry;

			char file_ref[MAX_FILE_NMAE];

			int filen = 0;
			while((entry = readdir(directory_pointer)))
			{
				//begin with each chr file
				strcpy(file_ref, source);
				strcat(file_ref, "/");
				strcat(file_ref, entry->d_name);

				//do with first_last kmers here
				if(strstr(entry->d_name, ".fq") || strstr(entry->d_name, ".fastq"))
				{
					fp_fa = fopen(file_ref, "r");
					if(fp_fa == NULL)
					{
						fprintf(stderr, "File error of openning ref file %s\n", file_ref);
						exit(1);
					}

					lfkmer_fq(fq_inputs, fa_inputs, fq_quality, quality_char, quality_flag);//,

					filen++;
				}
				else if(strstr(entry->d_name, ".fna") || strstr(entry->d_name, ".fa") || strstr(entry->d_name, ".fasta"))
				{
					fp_fa = fopen(file_ref, "r");
					if(fp_fa == NULL)
					{
						fprintf(stderr, "File error of openning ref file %s\n", file_ref);
						exit(1);
					}

					lfkmer_fa(fa_inputs);//,
					filen++;
				}
				else
				{
					fprintf(stderr, "Wrong input file type\n");
					continue;
				}

				fclose(fp_fa);
			}
			if(filen == 0)
			{
				fprintf(stderr, "No available files under %s\n", source);
				exit(1);
			}
		}
	}
	else if(S_ISREG(buf.st_mode))
	{

		fp_fa = fopen(source, "r");
		if(fp_fa == NULL)
		{
			fprintf(stderr, "File error of openning ref file %s\n", source);
			exit(1);
		}

		if(strstr(source, ".fq") || strstr(source, ".fastq"))
		{
			fprintf(stderr, "%s is FQ FILE\n", source);
			lfkmer_fq(fq_inputs, fa_inputs, fq_quality, quality_char, quality_flag);//,
		}
		else if(strstr(source, ".fna") || strstr(source, ".fa") || strstr(source, ".fasta"))
		{
			fprintf(stderr, "%s is FA FILE\n", source);
			lfkmer_fa(fa_inputs);
		}
		else
		{
			fprintf(stderr, "Wrong input file type\n");
			exit(1);
		}
		fclose(fp_fa);
	}
	else
	{
		fprintf(stderr, "File error for %s\n", source);
		exit(1);
	}

	if(fp_lfk_l)	fclose(fp_lfk_l);
	if(fp_lfk_r)	fclose(fp_lfk_r);
	if(fp_lfk_b)	fclose(fp_lfk_b);

#endif

#endif

	thread_data* tt = NULL;

	leftmove32 = 0;

	uint64_t key_tmp_re = (((KMER_LENGTH << 1) + 1) << 2);
	kmer_num_br = (key_tmp_re >> 6);
	if(key_tmp_re & 0X3f)	kmer_num_br++;

	uint64_t tempMove_r_new = ((KMER_LENGTH_PlusTwo & 0X3) << 1);
	tempMove_r_new2 = (tempMove_r_new == 0 ? 6:(tempMove_r_new - 2));
	tempMove_r_new3 = tempMove_r_new2 + 56;
	tempMove_new = (tempMove_r_new == 0 ? 0:(8 - tempMove_r_new));
	tempMove_new2 = (tempMove_r_new3 - (BUCKET_LENGTH << 1));

	temphash = 3;
	br_off = (tempMove_new >> 1);										//0 1 2 3
	br_offm = br_off + 1;												//1 2 3 4

	key_tmp_re = (key_len & 0X7);
	keyMove = 64 - (key_tmp_re << 3);
	if(keyMove == 64)	keyMove	= 0;

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

	uint8_t kmer_num_char = kmer_num_filter << 3;

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

	//get unique k-mer lf-kmers and add them into kmerInfo

	uint64_t mem_tol = 0;
	double mem_tmp = 0;

	mem_tmp = ((double )(kmer_num_new << 4)) / ((double )pair_len) + 1;

	mem_tol = (uint64_t )((double)(((uint64_t)memoryKmer<<20) - (((uint64_t )1 << 19) * (THREAD_NUM + 1))) / mem_tmp);

	mem_tol = ((mem_tol / pair_len) * pair_len);


	fprintf(stderr, "key_bits=%lu, key_len=%lu, val_len=%lu, pair_len=%lu\n",key_bits,key_len,val_len,pair_len);//8 4
	fprintf(stderr, "header_size: %lu tempMove: %"PRId64" mem_tol: %"PRId64"\n", header_size, tempMove, mem_tol);
	fprintf(stderr, "flag_filter: %u filter_min: %u filter_max: %u %"PRId64" %"PRId64" %lf\n\n", flag_filter, filter_min, filter_max, ((uint64_t)memoryKmer<<20), (((uint64_t )1 << 19) * (THREAD_NUM + 1)), mem_tmp);

#ifdef	JELLY_LFKMER
	lf_tage = 1;
	strcpy(segKmerName, tmp_route);
	strcat(segKmerName, "/lfInfo_l");
	bucket_sort_lfkmers(mem_tol, THREAD_NUM, segKmerName, tmp_route);

	lf_tage = 2;
	strcpy(segKmerName, tmp_route);
	strcat(segKmerName, "/lfInfo_r");
	bucket_sort_lfkmers(mem_tol, THREAD_NUM, segKmerName, tmp_route);

#ifdef	ADD_TMP
	bfq_n = 1;
#endif
	if(bfq_n && flag_filter)
	{
		lf_tage = 3;
		strcpy(segKmerName, tmp_route);
		strcat(segKmerName, "/lfInfo_b");
		bucket_sort_lfkmers(mem_tol, THREAD_NUM, segKmerName, tmp_route);
	}
#endif

	if(flag_filter)	mem_tmp = ((double )(kmer_num_filter << 4)) / ((double )pair_len) + 1;
	else	mem_tmp = ((double )(kmer_num_new << 4)) / ((double )pair_len) + 1;

	mem_tol = (uint64_t )((double)(((uint64_t)memoryKmer<<20) - (((uint64_t )1 << 19) * (THREAD_NUM + 1))) / mem_tmp);

	mem_tol = ((mem_tol / pair_len) * pair_len);

#ifdef	MULTI_COUNT
	uint64_t hash_c = 0;
	uint64_t hash_dis = 0;
	uint64_t i_in = 0;

	hashKmer_dis = (uint64_t** )calloc(THREAD_NUM, sizeof(uint64_t* ));

	countKmer_dis = (uint64_t** )calloc(THREAD_NUM, sizeof(uint64_t* ));
#endif

	hashCount=0;

	readBuf = (uint8_t *)calloc(mem_tol + 1, 1);
	if(readBuf == NULL)	fprintf(stderr, "Fail to allocate memory for readBuf\n"), exit(1);
	else fprintf(stderr, "Load %"PRId64" bytes into memory every time\n", mem_tol);

	segCount=(uint64_t *)calloc(THREAD_NUM+1,sizeof(uint64_t));
	segCount[THREAD_NUM]=BUCKET_CAPACITY;

	time_t start_t = 0, end_t = 0;
	clock_t start_c = 0, end_c = 0;

	while(!feof(fpKmer))//loopFlag
	{
		if((readNum = fread(readBuf, 1, mem_tol, fpKmer))>0)
		{
			countKmer=(uint64_t *)calloc(BUCKET_CAPACITY+1,sizeof(uint64_t));

			myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
			tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));
			for(i=0; i<THREAD_NUM; i++)
			{
				tt[i].num = i;
				tt[i].thread = THREAD_NUM;

				int check=pthread_create(&myThread[i], NULL, multiCount, tt + i);
				if(check)
				{
					fprintf(stderr,"threadNum:%lu, Error count - pthread_create() return code: %d\n",i,check);
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
		else
		{
			break;
		}

		myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
		tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));
		for(i=0; i<THREAD_NUM; i++)
		{
			tt[i].thread = THREAD_NUM;
			tt[i].num = i;

			int check=pthread_create(&myThread[i], NULL, multiDistri, tt + i);

			if(check)
			{
				fprintf(stderr,"threadNum:%lu, Error distri - pthread_create() return code: %d\n",i,check);
				exit(EXIT_FAILURE);
			}
		}
		for(i=0; i<THREAD_NUM; i++)
		{
			pthread_join( myThread[i], NULL);
		}
		free(myThread);
		free(tt);

		countKmer[0] = 0;

		for(i=1; i<BUCKET_CAPACITY+1; i++)
		{
			for(i_in=0; i_in<THREAD_NUM; i_in++)
			{
				countKmer[i] += (countKmer_dis[i_in][i] - countKmer_dis[i_in][i - 1]);
			}
			countKmer[i] += countKmer[i - 1];
		}
		totalKmerNum = countKmer[BUCKET_CAPACITY];
		hashKmer = (uint64_t* )calloc(totalKmerNum, kmer_num_char);

		for(i=0; i<BUCKET_CAPACITY; i++)
		{
			for(i_in=0, hash_c = 0; i_in<THREAD_NUM; i_in++)
			{
				hash_dis = countKmer_dis[i_in][i+1] - countKmer_dis[i_in][i];

				if(hash_dis)
				{
					memcpy(hashKmer + (countKmer[i] + hash_c) * kmer_num_filter, hashKmer_dis[i_in] + countKmer_dis[i_in][i] * kmer_num_filter, hash_dis * kmer_num_char);//(hash_dis * kmer_num) << 3
				}

				hash_c += hash_dis;
			}
		}

		segCount[0]=0;
		for(i=1; i<THREAD_NUM; i++)
		{
			locateNum = i*(totalKmerNum/THREAD_NUM);
			locate = BinarySearch(locateNum,countKmer,BUCKET_CAPACITY-1);
			segCount[i] = locate;
		}

		myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
		tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));
		for(i=0; i<THREAD_NUM; i++)
		{
			tt[i].thread = THREAD_NUM;
			tt[i].num = i;
			int check=pthread_create(&myThread[i], NULL, multiThreadSort, tt + i);

			if(check)
			{
				fprintf(stderr,"threadNum:%lu, Error sort - pthread_create() return code: %d\n",i,check);
				exit(EXIT_FAILURE);
			}
		}

		for(i=0; i<THREAD_NUM; i++)
		{
			pthread_join( myThread[i], NULL);
		}
		free(myThread);
		free(tt);

		strcpy(segKmerName, tmp_route);
		strcat(segKmerName, "/kmerInfo.");
		sprintf(cNum,"%u", hashCount);
		strcat(segKmerName, cNum);

		fpKmer_w = fopen(segKmerName,"wb");

#ifndef	TIME_TEST
		fwrite(hashKmer, kmer_num_char, totalKmerNum,fpKmer_w);
#endif
		fclose(fpKmer_w);
		free(hashKmer);
		free(countKmer);

		fprintf(stderr, "Block %d counting and sorting finish\n", hashCount + 1);

		hashCount++;
	}

	if(hashKmer_dis)	free(hashKmer_dis);
	if(countKmer_dis)	free(countKmer_dis);

	if(readBuf)	free(readBuf);
	if(segCount)	free(segCount);

	if(fpKmer)	fclose(fpKmer);
	
	/**/
	strcpy(segKmerName, "rm -f ");
	strcat(segKmerName, tmp_route);
	strcat(segKmerName, "/kmerInfo");
	system(segKmerName);
	

	return (void *)NULL;
}

void bucket_sort_lfkmers(uint64_t mem_tol, uint64_t THREAD_NUM, char* lfinfo, char* tmp_route)
{
	uint8_t i_in = 0;
	uint8_t kmer_num_char = kmer_num_new << 3;

	int i = 0;
	uint64_t hash_c = 0;
	uint64_t hash_dis = 0;
	unsigned int hashC = 0;
	uint64_t totalKmerNum = 0;
	uint64_t locateNum = 0;
	uint64_t locate = 0;
	char cNum[10];
	char segKmerName[50];

	//pthread_t myThread[THREAD_NUM];
	pthread_t* myThread = NULL;

	thread_data* tt = NULL;
	FILE* fpKmer_w = NULL;

	FILE* fpKmer = fopen(lfinfo,"rb");
	if(fpKmer==NULL) fprintf(stderr,"Fail to open %s\n", lfinfo),exit(1);

	hashKmer_dis = (uint64_t** )calloc(THREAD_NUM + 1, sizeof(uint64_t* ));
	countKmer_dis = (uint64_t** )calloc(THREAD_NUM + 1, sizeof(uint64_t* ));

	readBuf = (uint8_t *)calloc(mem_tol + 1, 1);
	if(readBuf == NULL)	fprintf(stderr, "Fail to allocate memory for readBuf\n"), exit(1);

	segCount = (uint64_t* )calloc(THREAD_NUM+1,sizeof(uint64_t));
	segCount[THREAD_NUM]=BUCKET_CAPACITY;

	while(!feof(fpKmer))//loopFlag
	{
		if((readNum = fread(readBuf, 1, mem_tol, fpKmer))>0)
		{
			countKmer=(uint64_t *)calloc(BUCKET_CAPACITY+1,sizeof(uint64_t));

			myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
			tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));
			for(i=0; i<THREAD_NUM; i++)
			{
				tt[i].num = i;
				tt[i].thread = THREAD_NUM;

				int check=pthread_create(&myThread[i], NULL, multiCount_lfkmers, tt + i);
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
		}
		else
		{
			break;
		}

		fprintf(stderr, "LF-kmers Block %d reading kmers and counting finish\n", hashC + 1);

		myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
		tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));
		for(i=0; i<THREAD_NUM; i++)
		{
			tt[i].thread = THREAD_NUM;
			tt[i].num = i;

			int check=pthread_create(&myThread[i], NULL, multiDistri_lfkmers, tt + i);

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

		fprintf(stderr, "LF-kmers Block %d distributing finish\n", hashC + 1);

		countKmer[0] = 0;

		for(i=1; i<BUCKET_CAPACITY+1; i++)
		{
			for(i_in=0; i_in<THREAD_NUM; i_in++)
			{
				countKmer[i] += (countKmer_dis[i_in][i] - countKmer_dis[i_in][i - 1]);
			}
			countKmer[i] += countKmer[i - 1];
		}
		totalKmerNum = countKmer[BUCKET_CAPACITY];
		hashKmer = (uint64_t* )calloc(totalKmerNum, kmer_num_char);

		for(i=0; i<BUCKET_CAPACITY; i++)
		{
			for(i_in=0, hash_c = 0; i_in<THREAD_NUM; i_in++)
			{
				hash_dis = countKmer_dis[i_in][i+1] - countKmer_dis[i_in][i];

				if(hash_dis)
					memcpy(hashKmer + (countKmer[i] + hash_c) * kmer_num_new, hashKmer_dis[i_in] + countKmer_dis[i_in][i] * kmer_num_new, hash_dis * kmer_num_char);//(hash_dis * kmer_num) << 3

				hash_c += hash_dis;
			}
		}

		segCount[0]=0;
		for(i=1; i<THREAD_NUM; i++)
		{
			locateNum = i*(totalKmerNum/THREAD_NUM);
			locate = BinarySearch(locateNum,countKmer,BUCKET_CAPACITY-1);
			segCount[i] = locate;
		}

		myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
		tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));
		for(i=0; i<THREAD_NUM; i++)
		{
			tt[i].thread = THREAD_NUM;
			tt[i].num = i;
			int check=pthread_create(&myThread[i], NULL, multiThreadSort_lfkmers, tt + i);

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

		fprintf(stderr, "LF-kmers Block %d sorting finish\n", hashC + 1);

		strcpy(segKmerName, tmp_route);
		strcat(segKmerName, "/lfInfo.");
		sprintf(cNum,"%u", hashC);
		strcat(segKmerName, cNum);

		fpKmer_w = fopen(segKmerName,"wb");

		fwrite(hashKmer, kmer_num_char, totalKmerNum, fpKmer_w);

		fclose(fpKmer_w);
		free(hashKmer);
		free(countKmer);

		hashC++;
	}

		
	

	
	if(hashKmer_dis)	free(hashKmer_dis);
	if(countKmer_dis)	free(countKmer_dis);
	if(readBuf)	free(readBuf);
	if(fpKmer)	fclose(fpKmer);
	if(segCount) free(segCount);

#ifdef	DEL_LFINFO
	strcpy(segKmerName, "rm -f ");
	strcat(segKmerName, lfinfo);
	system(segKmerName);
#endif

	//multi-tag on lf-kmers
	uint8_t kmer_num_tmp = (kmer_num << 3);
	uint8_t kmer_num_tmp_new = (kmer_num_new << 3);

	uint32_t num = 0;
	int block_num = 0;
	int heapTail = 0;
	int temp = 0;
	int cdi1 = 0;
	int cdi2 = 0;

	int *heap=NULL;

	uint64_t pBufM = 0;
	uint64_t pBufM_pre = 0;
	uint64_t bufferSizeM = 0;
	uint64_t segment = 0;
	uint64_t tempK = 0;
	uint64_t bufferSize = 0;

	uint64_t eliminate=~(ELIMINATE);
	uint64_t eliminate_mi = (eliminate << 2);

	splitMK2_lf = (uint64_t* )calloc(THREAD_NUM+1,sizeof(uint64_t));

	uint64_t readNum[hashC];
	uint64_t pBuf[hashC];
	uint64_t* bufK2[hashC];
	FILE *fpK2[hashC];

	uint64_t tmp_comf[kmer_num];
	uint64_t tmp_coms[kmer_num];

	bufferSize = (uint64_t )((double)((uint64_t)memoryKmer<<20) / (double)(THREAD_NUM * (kmer_num_tmp_new + pair_len) + hashC * kmer_num_tmp_new));

	bufferSizeM = bufferSize * THREAD_NUM;

	bufMK2_lf = (uint64_t* )calloc(bufferSizeM, kmer_num_tmp_new);
	heap = (int *)calloc(hashC, sizeof(int));

	write_lf_buff = (uint8_t** )calloc(THREAD_NUM, sizeof(uint8_t* ));
	for(i=0; i<THREAD_NUM; i++)
		if((write_lf_buff[i] = (uint8_t* )calloc(bufferSize + THREAD_NUM, pair_len)) == NULL)
			fprintf(stderr, "Fail to allocate memory for lf_kmers write buff\n"), exit(1);

	strcpy(segKmerName, tmp_route);
	strcat(segKmerName, "/kmerInfo");
	FILE* fp_allinfo = fopen(segKmerName,"ab+");
	if(fp_allinfo==NULL) fprintf(stderr,"Fail to open kmerInfo to add\n"),exit(1);

	if(hashC == 1)
	{
		strcpy(segKmerName, tmp_route);
		strcat(segKmerName, "/lfInfo.0");

		fpK2[0]=fopen(segKmerName,"rb");
		if(fpK2[0]==NULL) fprintf(stderr,"Fail to open lfInfo.0\n"),exit(1);
	}
	else
	{
		for(num=0; num<hashC; num++)
		{
			sprintf(cNum,"%u",num);

			strcpy(segKmerName, tmp_route);
			strcat(segKmerName, "/lfInfo.");
			strcat(segKmerName, cNum);

			bufK2[num]=(uint64_t* )calloc(bufferSize+1, kmer_num_tmp_new);

			fpK2[num]=fopen(segKmerName,"rb");
			if(fpK2[num]==NULL) fprintf(stderr,"Fail to open lfInfo.%u\n",num),exit(1);

			readNum[num]=fread(bufK2[num], kmer_num_tmp_new, bufferSize, fpK2[num]);

			pBuf[num]=0;
		}
		for(num=0; num<hashC; num++)
		{
			heapTail=num;
			for(i=heapTail; i>0; i=(i-1)>>1)
			{
				if(cmpMK2_s_lfkmers(bufK2[num], bufK2[heap[(i-1)>>1]]))
				{
					heap[i]=heap[(i-1)>>1];
				}
				else break;
			}
			heap[i]=num;
		}
	}

	fprintf(stderr, "hashC: %u %"PRId64"\n", hashC, bufferSizeM);

	int bound_off = 0;

	while(heapTail>=0)
	{
		pBufM=bound_off;

		if(hashC == 1)
		{
			pBufM = fread(bufMK2_lf + pBufM * kmer_num_new, kmer_num_tmp_new, bufferSizeM - pBufM, fpK2[0]);
			if(pBufM == 0)	break;

			fprintf(stderr, "hashC==1 %"PRId64"\n", pBufM);
		}
		else
		{
			while((pBufM < bufferSizeM) && (heapTail >= 0))//get bufMK2_lf
			{
				num=heap[0];

				memcpy(bufMK2_lf + pBufM * kmer_num_new, bufK2[num] + pBuf[num] * kmer_num_new, kmer_num_tmp_new);

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

					temp = (((cdi1<=heapTail) && cmpMK2_s_lfkmers(bufK2[heap[cdi1]] + pBuf[heap[cdi1]] * kmer_num_new, bufK2[heap[cdi2]] + pBuf[heap[cdi2]] * kmer_num_new)) ? cdi1:cdi2);

					if(cmpMK2_s_lfkmers(bufK2[heap[temp]] + pBuf[heap[temp]] * kmer_num_new, bufK2[num] + pBuf[num] * kmer_num_new))
					{
						heap[i]=heap[temp];
					}
					else break;
				}
				heap[i]=num;

				if(pBufM>=bufferSizeM-15)
				{
					if(OneUnit)
					{
						if((bufMK2_lf[pBufM-1]<<2)!=(bufK2[heap[0]][pBuf[heap[0]]]<<2))
							break;
					}
					else
					{
						copykmer(tmp_comf, bufMK2_lf, (pBufM-1) * kmer_num_new)
						copykmer(tmp_coms, bufK2[heap[0]], pBuf[heap[0]] * kmer_num_new)
						if(memcmp(tmp_comf, tmp_coms, kmer_num_tmp))	break;
					}
				}
			}
		}

		if(pBufM == 0)	break;

		segment=pBufM/THREAD_NUM;//4194289

		uint64_t tempK = 0;
		splitMK2_lf[0]=0; //lower bound (can reach)
		for(i=1; i<THREAD_NUM; i++)
		{
			splitMK2_lf[i]=i*segment;

			if(pBufM < THREAD_NUM)	continue;

				if(OneUnit)
				{
					tempK = (bufMK2_lf[splitMK2_lf[i]]<<2);

					while((bufMK2_lf[splitMK2_lf[i]+1]<<2)==tempK)
					{
						splitMK2_lf[i]++;
					}
				}
				else
				{
					copykmer(tmp_comf, bufMK2_lf, splitMK2_lf[i] * kmer_num_new)

					while(1)
					{
						copykmer(tmp_coms, bufMK2_lf, (splitMK2_lf[i] + 1) * kmer_num_new)
						if(memcmp(tmp_comf, tmp_coms, kmer_num_tmp))	break;
						else	splitMK2_lf[i]++;
					}
				}

			splitMK2_lf[i]++;
		}
		splitMK2_lf[THREAD_NUM]=pBufM;

		fprintf(stderr, "LF-kmers multi-Tagcheck start:-------------------------------------------- %d\n", block_num);

		myThread = (pthread_t* )calloc(THREAD_NUM, sizeof(pthread_t ));
		thread_data* tt = (thread_data* )calloc(THREAD_NUM, sizeof(thread_data ));
		for(i=0; i<THREAD_NUM; i++)
		{
			w_lf_n[i] = 0;

			tt[i].thread = THREAD_NUM;
			tt[i].num = i;
			int check=pthread_create(&myThread[i], NULL, multiTagcheck_lfkmers, tt + i);

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
			fwrite(write_lf_buff[i], pair_len, w_lf_n[i], fp_allinfo);

		//get num of last k-1 mers and copy these kmers to new bufMK2_lf
		if(heapTail >= 0)
		{
			bound_off = pBufM - last_block_num_lfkmers;
			fprintf(stderr, "LF-kmers BOUNDARY: %d %"PRId64" %"PRId64"\n", bound_off, pBufM, last_block_num_lfkmers);
			memcpy(bufMK2_lf, bufMK2_lf + (pBufM - bound_off) * kmer_num_new, bound_off * kmer_num_tmp_new);
		}

		pBufM_pre = pBufM;

		fprintf(stderr, "LF-kmers Merging block %d and multiTagcheck finish\n", block_num + 1);

		block_num++;
	}

	multiTagcheck_lfkmers_addkmer(pBufM_pre, last_block_num_lfkmers, fp_allinfo);

	fprintf(stderr, "LF-kmers MultiTag-check different k-1 mers finish\n");

	if(bufMK2_lf)	free(bufMK2_lf);
	if(heap)	free(heap);
	if(splitMK2_lf)	free(splitMK2_lf);

	for(i=0; i<THREAD_NUM; i++)
		if(write_lf_buff[i])	free(write_lf_buff[i]);
	if(write_lf_buff)	free(write_lf_buff);

	if(fp_allinfo)	fclose(fp_allinfo);

#ifdef	COM_DEL
	for(i=0; i<hashC; i++)
	{
		strcpy(segKmerName, "rm -f ");
		strcat(segKmerName, tmp_route);
		strcat(segKmerName, "/lfInfo.");
		sprintf(cNum,"%d", i);
		strcat(segKmerName, cNum);
		system(segKmerName);
	}
#endif
}

void multiTagcheck_lfkmers_addkmer(uint64_t i, uint64_t anchor_km_pre_i, FILE* fp_allinfo)
{
	uint32_t kmer_n = i - anchor_km_pre_i;
	int k = 0;
	uint8_t write_kmer[pair_len];

	if(OneUnit)
	{
		uint8_t kmerlen_minus = (32 - KMER_LENGTH_PlusTwo) << 1;

		uint64_t anchor_km_pre = bufMK2_lf[anchor_km_pre_i];
		uint64_t anchor_tmp = (anchor_km_pre << 2) >> kmerlen_minus;
		for(k = 0; k < key_len; k++)
			write_kmer[k] = ((anchor_tmp >> (k << 3)) & 0Xff);

		if(kmer_n > 0Xffffff)	kmer_n = 0Xffffff;

		write_kmer[key_len] = (kmer_n & 0Xff);
		write_kmer[key_len + 1] = ((kmer_n >> 8) & 0Xff);
		write_kmer[key_len + 2] = ((kmer_n >> 16) & 0Xff);

		write_kmer[pair_len - 1] = lf_tage;

		fwrite(write_kmer, pair_len, 1, fp_allinfo);

	}
	else
	{
		uint8_t kmer_num_tmp = (kmer_num << 3);
		uint64_t start_get = (anchor_km_pre_i * kmer_num_new);
		uint64_t anchor_kms_pre[kmer_num];

		copykmer(anchor_kms_pre, bufMK2_lf, start_get)

		for(k = 0; k < key_len; k++)
			write_kmer[key_len - 1 - k] = ((anchor_kms_pre[k >> 3] >> (56 - ((k & 0X7) << 3))) & 0Xff);

		if(kmer_n > 0Xffffff)	kmer_n = 0Xffffff;

		write_kmer[key_len] = (kmer_n & 0Xff);
		write_kmer[key_len + 1] = ((kmer_n >> 8) & 0Xff);
		write_kmer[key_len + 2] = ((kmer_n >> 16) & 0Xff);

		write_kmer[pair_len - 1] = lf_tage;

		fwrite(write_kmer, pair_len, 1, fp_allinfo);

	}

}

int cmpMK2_s_lfkmers(uint64_t* va, uint64_t* vb)
{
	uint64_t a = (*va);
	uint64_t b = (*vb);

	if(OneUnit)
	{
		a=a<<2, b=b<<2;

		if(a<b) return 1;
		else return 0;
	}
	else
	{
		a <<= tempMove_new;
		b <<= tempMove_new;

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
					//else if(a>b) return 0;
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

void *multiTagcheck_lfkmers(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	int num = argT->num;
	int totalNum = argT->thread;

	uint64_t i = 0;
	uint64_t anchor_km_pre_i = 0;
	uint64_t w_lf_t = 0;

	uint64_t start=splitMK2_lf[num], end=splitMK2_lf[num+1];
	uint64_t bufMK1_len=end-start;

	if(bufMK1_len==0) return (void *)NULL;

	uint8_t write_kmer[pair_len];

	int k = 0;
	uint32_t kmer_n = 0;

	anchor_km_pre_i = start;

	if(OneUnit)
	{
		uint64_t anchor_km = 0;
		uint64_t anchor_km_pre = 0;
		uint64_t anchor_tmp = 0;

		uint8_t kmerlen_minus = (32 - KMER_LENGTH_PlusTwo) << 1;
		for(i = start; i < end; i++)
		{
			anchor_km = bufMK2_lf[i];

			if((anchor_km != anchor_km_pre) && (i > start))
			{
				anchor_tmp = (anchor_km_pre << 2) >> kmerlen_minus;

				for(k = 0; k < key_len; k++)
					write_kmer[k] = ((anchor_tmp >> (k << 3)) & 0Xff);

				if(kmer_n > 0Xffffff)	kmer_n = 0Xffffff;

				write_kmer[key_len] = (kmer_n & 0Xff);
				write_kmer[key_len + 1] = ((kmer_n >> 8) & 0Xff);
				write_kmer[key_len + 2] = ((kmer_n >> 16) & 0Xff);

				write_kmer[pair_len - 1] = lf_tage;

				memcpy(write_lf_buff[num] + w_lf_t * pair_len, write_kmer, pair_len);
				w_lf_t++;

				kmer_n = 0;
				anchor_km_pre_i = i;
			}

			kmer_n++;
			anchor_km_pre = anchor_km;
		}

		if(num != totalNum - 1)
		{
			anchor_km_pre = bufMK2_lf[anchor_km_pre_i];
			anchor_tmp = (anchor_km_pre << 2) >> kmerlen_minus;
			for(k = 0; k < key_len; k++)
				write_kmer[k] = ((anchor_tmp >> (k << 3)) & 0Xff);

			if(kmer_n > 0Xffffff)	kmer_n = 0Xffffff;

			write_kmer[key_len] = (kmer_n & 0Xff);
			write_kmer[key_len + 1] = ((kmer_n >> 8) & 0Xff);
			write_kmer[key_len + 2] = ((kmer_n >> 16) & 0Xff);

			write_kmer[pair_len - 1] = lf_tage;

			memcpy(write_lf_buff[num] + w_lf_t * pair_len, write_kmer, pair_len);
			w_lf_t++;
		}
		else	last_block_num_lfkmers = anchor_km_pre_i;
	}
	else
	{
		uint64_t start_get = 0;
		uint64_t anchor_kms[kmer_num];
		uint64_t anchor_kms_pre[kmer_num];

		uint8_t kmer_num_tmp = (kmer_num << 3);

		for(i = start; i < end; i++)
		{
			start_get = (i * kmer_num_new);

			copykmer(anchor_kms, bufMK2_lf, start_get)

			if(memcmp(anchor_kms, anchor_kms_pre, kmer_num_tmp) && (i > start))
			{
				for(k = 0; k < key_len; k++)
					write_kmer[key_len - 1 - k] = ((anchor_kms_pre[k >> 3] >> (56 - ((k & 0X7) << 3))) & 0Xff);

				if(kmer_n > 0Xffffff)	kmer_n = 0Xffffff;

				write_kmer[key_len] = (kmer_n & 0Xff);
				write_kmer[key_len + 1] = ((kmer_n >> 8) & 0Xff);
				write_kmer[key_len + 2] = ((kmer_n >> 16) & 0Xff);

				write_kmer[pair_len - 1] = lf_tage;

				memcpy(write_lf_buff[num] + w_lf_t * pair_len, write_kmer, pair_len);
				w_lf_t++;

				kmer_n = 0;
				anchor_km_pre_i = i;
			}

			kmer_n++;
			memcpy(anchor_kms_pre, anchor_kms, kmer_num_tmp);
		}

		if(num != totalNum - 1)
		{
			start_get = (anchor_km_pre_i * kmer_num_new);

			copykmer(anchor_kms_pre, bufMK2_lf, start_get)

			for(k = 0; k < key_len; k++)
				write_kmer[key_len - 1 - k] = ((anchor_kms_pre[k >> 3] >> (56 - ((k & 0X7) << 3))) & 0Xff);

			if(kmer_n > 0Xffffff)	kmer_n = 0Xffffff;

			write_kmer[key_len] = (kmer_n & 0Xff);
			write_kmer[key_len + 1] = ((kmer_n >> 8) & 0Xff);
			write_kmer[key_len + 2] = ((kmer_n >> 16) & 0Xff);

			write_kmer[pair_len - 1] = lf_tage;

			memcpy(write_lf_buff[num] + w_lf_t * pair_len, write_kmer, pair_len);
			w_lf_t++;
		}
		else	last_block_num_lfkmers = anchor_km_pre_i;
	}

	w_lf_n[num] = w_lf_t;

	return (void *)NULL;
}

void *multiThreadSort(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	uint64_t num=argT->num;
	uint64_t THREAD_NUM=argT->thread;

	uint64_t i,low,up;
	//uint64_t kmer_num_tmp = (kmer_num_new << 3);
	uint64_t kmer_num_tmp = (kmer_num_filter << 3);

	low=segCount[num];
	up=segCount[num+1];

#ifdef	MULTI_COUNT
	if(countKmer_dis[num])	free(countKmer_dis[num]);
	if(hashKmer_dis[num])	free(hashKmer_dis[num]);
#endif

	for(i=low; i<up; i++)
	{
		qsort(hashKmer + countKmer[i] * kmer_num_filter, countKmer[i+1] - countKmer[i], kmer_num_tmp, cmpKmer);
	}
	return (void *)NULL;
}

void *multiThreadSort_lfkmers(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	uint64_t num=argT->num;
	uint64_t THREAD_NUM=argT->thread;

	uint64_t i,low,up;
	uint64_t kmer_num_tmp = (kmer_num_new << 3);

	low=segCount[num];
	up=segCount[num+1];

	if(countKmer_dis[num])	free(countKmer_dis[num]);
	if(hashKmer_dis[num])	free(hashKmer_dis[num]);

	for(i=low; i<up; i++)
	{
		qsort(hashKmer + countKmer[i] * kmer_num_new, countKmer[i+1] - countKmer[i], kmer_num_tmp, cmpKmer_lfkmers);
	}
	return (void *)NULL;
}

int cmpKmer(const void *a, const void *b)
{
	if(OneUnit)
	{
		uint64_t* vav = (uint64_t* )a;
		uint64_t* vbv = (uint64_t* )b;
		uint64_t va = *vav;
		uint64_t vb = *vbv;

		//printf("cmp: %u\n", va);
		//decode(va);
#ifdef	FIRST_LAST_KMER
		va <<= 4;
		vb <<= 4;
#else
		va <<= 2;
		vb <<= 2;
#endif

		//printf("cmp: %u\n", va);
		//decode(va);

		//printf("va=%lu; vb=%lu\n",va,vb );
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

		va_t <<= (tempMove_new + 2);
		vb_t <<= (tempMove_new + 2);

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
					va_t = (va[kmer_num - 1] >> 2) << 2;
					vb_t = (vb[kmer_num - 1] >> 2) << 2;
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


int cmpKmer_lfkmers(const void *a, const void *b)
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

		va_t <<= tempMove_new;
		vb_t <<= tempMove_new;

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
					va_t = (va[kmer_num - 1] >> 2) << 2;
					vb_t = (vb[kmer_num - 1] >> 2) << 2;
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

void *multiDistri(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	uint64_t num=argT->num;
	uint64_t THREAD_NUM=argT->thread;

	uint64_t segment=(readNum/pair_len)/THREAD_NUM;
	uint64_t start=(num*segment)*pair_len,i;
	uint64_t end;
	if(num<THREAD_NUM-1)
	{
		end=((num+1)*segment)*pair_len;
	}
	else
	{
		end=readNum;
	}

	int8_t j = 0;
	uint8_t k = 0;
	uint8_t kmerlen_minus = 0;
#ifdef	FIRST_LAST_KMER
	kmerlen_minus = (31 - KMER_LENGTH_PlusTwo) << 1;
#else
	kmerlen_minus = (32 - KMER_LENGTH_PlusTwo) << 1;
#endif

	uint64_t seq=0;
	uint64_t tempSeq = 0;

	uint64_t* seq_a = (uint64_t* )calloc(kmer_num_filter, 8);
	uint8_t kmer_num_tmp = kmer_num_filter << 3;
	uint64_t filter_num = 0;
	uint64_t start_get = 0;

	for(i=start; i<end; i+=pair_len)
	{
		if(jelly_f)
		{
			if(OneUnit)
			{
				seq = 0;
				for(j = key_len - 1; j>=0; j--)
				{
					seq <<= 8;
					seq |= readBuf[i + j];
				}
				tempSeq=(seq>>tempMove)&maskBKT;

				if(readBuf[i + pair_len - 1])
					seq |= ((uint64_t )readBuf[i + pair_len - 1] << (KMER_LENGTH_PlusTwo << 1));

				if(flag_filter)
				{
					filter_num = 0;
					filter_num |= readBuf[i + key_len];
					filter_num |= ((uint32_t )readBuf[i + key_len + 1] << 8);
					filter_num |= ((uint32_t )readBuf[i + key_len + 2] << 16);
				}
			}
			else
			{
				if(flag_filter)
					seq_a[kmer_num_filter - 1] = 0;

				seq_a[kmer_num_new - 1] = 0;
				for(j = key_len - 1, k=0; j>=0; j--, k++)
				{
					seq_a[k >> 3] <<= 8;
					seq_a[k >> 3] |= readBuf[i + j];
				}
				seq_a[kmer_num - 1] <<= keyMove;
				seq_a[kmer_num_new - 1] |= readBuf[i + pair_len - 1];

				tempSeq=(seq_a[0]>>tempMove_new2)&maskBKT;

				if(flag_filter)
				{
					filter_num = 0;
					filter_num |= readBuf[i + key_len];
					filter_num |= ((uint32_t )readBuf[i + key_len + 1] << 8);
					filter_num |= ((uint32_t )readBuf[i + key_len + 2] << 16);

					seq_a[kmer_num_filter - 1] |= (filter_num << 32);
				}
			}
		}
		else
		{
			memcpy(&seq,&readBuf[i],key_len);
			tempSeq=(seq>>tempMove)&maskBKT;
		}

		countKmer_dis[num][tempSeq]--;

		if(OneUnit)
		{
			if(flag_filter)
			{
				start_get = countKmer_dis[num][tempSeq] * kmer_num_filter;
				hashKmer_dis[num][start_get]=(seq << kmerlen_minus);
				hashKmer_dis[num][start_get + 1]=filter_num;
			}
			else	hashKmer_dis[num][countKmer_dis[num][tempSeq]]=(seq << kmerlen_minus);
		}
		else
		{
			memcpy(hashKmer_dis[num] + countKmer_dis[num][tempSeq] * kmer_num_filter, seq_a, kmer_num_tmp);
		}
	}

	free(seq_a);

	return (void *)NULL;
}



void *multiDistri_lfkmers(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	uint64_t num=argT->num;
	uint64_t THREAD_NUM=argT->thread;
	uint64_t segment=(readNum/pair_len)/THREAD_NUM;
	uint64_t start=(num*segment)*pair_len,i;
	uint64_t end;
	if(num<THREAD_NUM-1)
	{
		end=((num+1)*segment)*pair_len;
	}
	else
	{
		end=readNum;
	}

	int8_t j = 0;
	uint8_t k = 0;
	uint8_t kmerlen_minus = 0;

	kmerlen_minus = (31 - KMER_LENGTH_PlusTwo) << 1;

	uint64_t seq=0;
	uint64_t tempSeq = 0;
	uint64_t* seq_a = (uint64_t* )calloc(kmer_num_new, 8);
	uint8_t kmer_num_tmp = kmer_num_new << 3;

	for(i=start; i<end; i+=pair_len)
	{
		if(jelly_f)
		{
			if(OneUnit)
			{
				seq = 0;
				for(j = key_len - 1; j>=0; j--)
				{
					seq <<= 8;
					seq |= readBuf[i + j];
				}
				tempSeq=(seq>>(tempMove + 2))&maskBKT;

				if(readBuf[i + pair_len - 1])
					seq |= ((uint64_t )readBuf[i + pair_len - 1] << (KMER_LENGTH_PlusTwo << 1));

			}
			else
			{
				seq_a[kmer_num_new - 1] = 0;
				for(j = key_len - 1, k=0; j>=0; j--, k++)
				{
					seq_a[k >> 3] <<= 8;
					seq_a[k >> 3] |= readBuf[i + j];

				}
				seq_a[kmer_num - 1] <<= keyMove;
				seq_a[kmer_num_new - 1] |= readBuf[i + pair_len - 1];

				tempSeq=(seq_a[0]>>(tempMove_new2 + 2))&maskBKT;
			}
		}
		else
		{
			memcpy(&seq,&readBuf[i],key_len);
			tempSeq=(seq>>tempMove)&maskBKT;
		}

		countKmer_dis[num][tempSeq]--;

		if(OneUnit)
			hashKmer_dis[num][countKmer_dis[num][tempSeq]]=(seq << kmerlen_minus);
		else	memcpy(hashKmer_dis[num] + countKmer_dis[num][tempSeq] * kmer_num_new, seq_a, kmer_num_tmp);
	}

	free(seq_a);

	return (void *)NULL;
}


void *multiCount(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	uint64_t num = argT->num;
	uint64_t THREAD_NUM = argT->thread;

	uint64_t segment=(readNum/pair_len)/THREAD_NUM;
	uint64_t start=(num*segment)*pair_len;
	uint64_t end;
	uint64_t i,seq,tempSeq;
	uint8_t kmer_num_tmp = kmer_num_filter << 3;

	int j = 0;

	if(num<THREAD_NUM-1)
	{
		end=((num+1)*segment)*pair_len;
	}
	else
	{
		end=readNum;
	}

	countKmer_dis[num] = (uint64_t* )calloc(BUCKET_CAPACITY+1,sizeof(uint64_t ));

	if(jelly_f)
	{
		for(i=start; i<end; i+=pair_len)
		{
			seq = 0;

			for(j = key_len - 1; j> key_len - temphash - 1; j--)
			{
				seq <<= 8;
				seq |= readBuf[i + j];
			}

			countKmer_dis[num][(seq>>tempMove_r_new2)&maskBKT]++;
		}
	}
	else
	{
		for(i=start; i<end; i+=pair_len)
		{
			memcpy(&seq, &readBuf[i], key_len);
			tempSeq=(seq>>tempMove)&maskBKT;
			countKmer_dis[num][tempSeq]++;
		}
	}

	for(i=1; i<BUCKET_CAPACITY; i++)
		countKmer_dis[num][i]=countKmer_dis[num][i-1]+countKmer_dis[num][i];

	countKmer_dis[num][BUCKET_CAPACITY] = countKmer_dis[num][BUCKET_CAPACITY-1];

	hashKmer_dis[num] = (uint64_t* )calloc(countKmer_dis[num][BUCKET_CAPACITY], kmer_num_tmp);// * kmer_num, sizeof(uint64_t)

	return (void*)NULL;
}


void *multiCount_lfkmers(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	uint64_t num = argT->num;
	uint64_t THREAD_NUM = argT->thread;

	uint64_t segment=(readNum/pair_len)/THREAD_NUM;
	uint64_t start=(num*segment)*pair_len;
	uint64_t end;
	uint64_t i,seq,tempSeq;
	uint8_t kmer_num_tmp = kmer_num_new << 3;

	int j = 0;

	if(num<THREAD_NUM-1)
	{
		end=((num+1)*segment)*pair_len;
	}
	else
	{
		end=readNum;
	}

	countKmer_dis[num] = (uint64_t* )calloc(BUCKET_CAPACITY+1,sizeof(uint64_t ));

	if(jelly_f)
	{
		for(i=start; i<end; i+=pair_len)
		{
			seq = 0;

			for(j = key_len - 1; j> key_len - temphash - 1; j--)
			{
				seq <<= 8;
				seq |= readBuf[i + j];
			}

			countKmer_dis[num][(seq>>(tempMove_r_new2 + 2))&maskBKT]++;
		}
	}
	else
	{
		for(i=start; i<end; i+=pair_len)
		{
			memcpy(&seq, &readBuf[i], key_len);
			tempSeq=(seq>>tempMove)&maskBKT;
			countKmer_dis[num][tempSeq]++;
		}
	}

	for(i=1; i<BUCKET_CAPACITY; i++)
		countKmer_dis[num][i]=countKmer_dis[num][i-1]+countKmer_dis[num][i];

	countKmer_dis[num][BUCKET_CAPACITY] = countKmer_dis[num][BUCKET_CAPACITY-1];

	hashKmer_dis[num] = (uint64_t* )calloc(countKmer_dis[num][BUCKET_CAPACITY], kmer_num_tmp);// * kmer_num, sizeof(uint64_t)

	return (void*)NULL;
}


char *getPath(char *DIR, char *name)
{
	char *path=(char *)calloc(PATH_LEN,sizeof(char));
	strcpy(path,DIR);
	strcat(path,name);
	return path;
}

uint64_t BinarySearch(uint64_t mk, uint64_t* target, int64_t up)//Search the # and $ point
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

void *multiThreadSort_new(void *arg)
{
	thread_data* argT=(thread_data* )arg;
	uint64_t num=argT->num;
	uint64_t THREAD_NUM=argT->thread;

	uint64_t i, start, end, len, r_p = 0, adjMove;
	adjMove = ((32-KMER_LENGTH_PlusTwo) << 1);

	start = splitIN[num];
	end = splitIN[num + 1];

	for(i = start; i < end; i++)
	{
		memcpy(readBuf_p[num] + r_p, readBuf_new + i * pair_len, key_len);
		readBuf_p[num][r_p] <<= adjMove;
		r_p++;
	}

	qsort(readBuf_p[num], r_p, sizeof(uint64_t ), cmpKmer);
	readBuf_p_num[num] = r_p;

	return (void *)NULL;
}


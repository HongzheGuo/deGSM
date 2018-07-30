#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include "ubwt.h"
#include "utils.h"

unsigned char nst_nt4_table[256] = {
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

void ubwt_init(ubwt_t *ubwt, ubwt_count_t ubwt_l)
{
    ubwt->ubwt_l = ubwt_l;
    ubwt_count_t n_occ = (ubwt_l + _OCC_INV - 1) / _OCC_INV; // do NOT store last one C[5]
    ubwt->ubwt_size = (ubwt_l + _BWT_INV - 1) / _BWT_INV + n_occ * _OCC_C;
    ubwt->ubwt = (ubwt_int_t*)_err_calloc(ubwt->ubwt_size, sizeof(ubwt_int_t));
    ubwt->ubwt_map = NULL;
    ubwt->ubwt_unit = 0, ubwt->ubwt_i = 0, ubwt->ubwt_k = 0;
    int i; for (i = 0; i < _OCC_C; ++i) ubwt->C[i] = 0;
}

void ubwt_free(ubwt_t *ubwt) { 
    if (ubwt->ubwt) free(ubwt->ubwt);
    if (ubwt->ubwt_map) free(ubwt->ubwt_map);
    if (ubwt) free(ubwt); 
}

// calculate C[] & OCC, based on bwt_str
void ubwt_cal_occ(uint8_t bwt_nt, ubwt_t *ubwt)
{
    bwt_nt = bwt_nt > nt_N ? nt_N : bwt_nt;
    //stdout_printf("BWT: %c\n", "ACGTN"[bwt_nt]);
    int i;
    if (ubwt->ubwt_i % _OCC_INV == 0) {
        for (i = 0; i < _OCC_C; ++i) ubwt->ubwt[ubwt->ubwt_k++] = ubwt->C[i];
    }
    ubwt->ubwt_unit |= (ubwt_int_t)(bwt_nt) << ((~ubwt->ubwt_i&_BWT_INV_M)<<_BWT_NT_B);
    ++ubwt->C[bwt_nt];
    ++ubwt->ubwt_i;
    
    if (ubwt->ubwt_i % _BWT_INV == 0 || ubwt->ubwt_i == ubwt->ubwt_l) {
        ubwt->ubwt[ubwt->ubwt_k++] = ubwt->ubwt_unit;
        ubwt->ubwt_unit = 0;
    }
}

int ubwt_cal(ubwt_t *ubwt, uint8_t *ubwt_bstr, ubwt_count_t ubwt_l)
{
    ubwt_count_t i, uni_c = 0;
    for (i = 0; i < ubwt_l; ++i) {
        if (ubwt_bstr[i] >= nt_N) ++uni_c;
        ubwt_cal_occ(ubwt_bstr[i], ubwt);
    }
    return uni_c;
}

void ubwt_gen_bit_table16(ubwt_t *ubwt)
{
    int i;
    ubwt->bit_table16[0] = 0;
    for (i = 0; i != 65536; ++i) 
        ubwt->bit_table16[i] = (i&1) + ubwt->bit_table16[i/2];
}

void ubwt_update(ubwt_t *ubwt)
{
    int i; 
    // update C[]
    for (i = 1; i < _OCC_C; ++i) ubwt->C[i] += ubwt->C[i-1];
    for (i = _OCC_C-1; i > 0; --i) ubwt->C[i] = ubwt->C[i-1];
    ubwt->C[0] = 0;
    // gen cnt table
    ubwt_gen_bit_table16(ubwt);
}

#define ubwt_bwt(ubwt, bwt_i) ((ubwt)->ubwt[_OCC_C + ubwt_bwt_occ_a((bwt_i)>>_OCC_INV_B) + (((bwt_i)&(_OCC_INV_M))>>_BWT_INV_B)])
#define ubwt_bwt_nt(ubwt, bwt_i) ((ubwt_bwt(ubwt, bwt_i)>>((~(bwt_i)&_BWT_INV_M)<<_BWT_NT_B))&_BWT_NT_M)
#define ubwt_occ_a(ubwt, bwt_i) ((ubwt)->ubwt + ubwt_bwt_occ_a((bwt_i)>>_OCC_INV_B))
#define __occ_cnt4(table, b) (table[(b)&0x1111] + table[(b)>>16&0x1111] + table[(b)>>32&0x1111] + table[(b)>>48&0x1111])

// count numbre of '1'
static inline int __occ_aux(ubwt_int_t b, const ubwt_t *ubwt, uint8_t c)
{
    b = ((c&4)?b:~b)>>2 & ((c&2)?b:~b)>>1 & ((c&1)?b:~b) & 0x1111111111111111ull;
    return __occ_cnt4(ubwt->bit_table16, b);
}

static inline ubwt_count_t ubwt_occ(const ubwt_t *ubwt, ubwt_count_t bwt_i, uint8_t nt)
{
    ubwt_int_t *p, *end;
    ubwt_count_t occ = (p = ubwt_occ_a(ubwt, bwt_i))[nt];
    p += _OCC_C;

    // OCC up to bwt_i/_BWT_INV
    end = p + ((bwt_i & _OCC_INV_M) >> _BWT_INV_B);
    for (; p < end; ++p) occ += __occ_aux(*p, ubwt, nt);

    occ += __occ_aux((*p >> ((~bwt_i&_BWT_INV_M)<<_BWT_NT_B))>>(1<<_BWT_NT_B), ubwt, nt);
    if (nt == 0) occ -= ((~bwt_i&_BWT_INV_M)+1);
    return occ;
}

// requiring k <= l
void ubwt_2occ(const ubwt_t *ubwt, ubwt_count_t k, ubwt_count_t l, uint8_t nt, ubwt_count_t *ok, ubwt_count_t *ol)
{
    if (k >> _OCC_INV_B != l >> _OCC_INV_B) {
        *ok = ubwt_occ(ubwt, k, nt);
        *ol = ubwt_occ(ubwt, l, nt);
    } else {
        ubwt_count_t occ_k, occ_l, *p0, *p, *end_k, *end_l;
        occ_k = (p = ubwt_occ_a(ubwt, k))[nt];
        p += _OCC_C, p0 = p;
        // calculate occ for k
        end_k = p + ((k & _OCC_INV_M) >> _BWT_INV_B);
        for (; p < end_k; ++p) occ_k += __occ_aux(*p, ubwt, nt);
        occ_l = occ_k;
        occ_k += __occ_aux((*p >> ((~k&_BWT_INV_M)<<_BWT_NT_B))>>(1<<_BWT_NT_B), ubwt, nt);
        // calculate occ for l
        end_l = p0 + ((l & _OCC_INV_M) >> _BWT_INV_B);
        for (; p < end_l; ++p) occ_l += __occ_aux(*p, ubwt, nt);
        occ_l += __occ_aux((*p >> ((~l&_BWT_INV_M)<<_BWT_NT_B))>>(1<<_BWT_NT_B), ubwt, nt);

        if (nt == 0) occ_k -= ((~k&_BWT_INV_M)+1), occ_l -= ((~l&_BWT_INV_M)+1);
        *ok = occ_k, *ol = occ_l;
    }
}

#define ubwt_set_intv(d, nt, k, l) {(k) = (d)->C[(nt)], (l) = (d)->C[(nt)+1]-1;}

ubwt_count_t ubwt_exact_match(const ubwt_t *ubwt, int qlen, const uint8_t *query, ubwt_count_t *bwt_k, ubwt_count_t *bwt_l)
{
    if (qlen < 1) return 0;
    ubwt_count_t k, l, occ_k, occ_l;

    if (query[qlen-1] >= nt_N) return 0;
    ubwt_set_intv(ubwt, query[qlen-1], k, l);

    int i;
    for (i = qlen-2; i >= 0; --i) {
        uint8_t c = query[i];
        if (c >= nt_N) return 0; // 'N' is not in bwt
        ubwt_2occ(ubwt, k, l+1, c, &occ_k, &occ_l);
        if (occ_k < occ_l) {
            k = ubwt->C[c] + occ_k;
            l = ubwt->C[c] + occ_l - 1;
        } else return 0;
    }
    *bwt_k = k, *bwt_l = l;
    return l-k+1;
}
ubwt_count_t ubwt_uid(ubwt_t *ubwt, ubwt_count_t k)
{
    ubwt_count_t occ_k;
    uint8_t nt;

    while (1) {
        nt = ubwt_bwt_nt(ubwt, k);
        occ_k = ubwt_occ(ubwt, k, nt);
        if (nt >= nt_N) break;
        k = ubwt->C[nt] + occ_k;
    }
    return ubwt->ubwt_map[occ_k];
}

ubwt_count_t ubwt_cal_off(ubwt_t *ubwt, ubwt_count_t k, ubwt_count_t *off)
{
    ubwt_count_t u_off=0, occ_k;
    uint8_t nt;

    while (1) {
        nt = ubwt_bwt_nt(ubwt, k);
        occ_k = ubwt_occ(ubwt, k, nt);
        if (nt >= nt_N) break;
        ++u_off;
        k = ubwt->C[nt] + occ_k;
    }
    *off = u_off;
    return ubwt->ubwt_map[occ_k];
}

// return unipath length, excluding 'N' symbol
ubwt_count_t ubwt_uni_len(ubwt_t *ubwt, ubwt_count_t uid) {
    ubwt_count_t uni_len = 0, k, occ_k;
    uint8_t nt;
    k = ubwt->C[nt_N] + uid;
    while (ubwt_bwt_nt(ubwt, k) < nt_N) {
        ++uni_len;
        nt  = ubwt_bwt_nt(ubwt, k);
        occ_k = ubwt_occ(ubwt, k, nt);
        k = ubwt->C[nt] + occ_k;
    }
    return uni_len;
}

void ubwt_out_unipath1(ubwt_t *ubwt, ubwt_count_t uid, FILE *out)
{
    char *unipath = (char*)_err_calloc(1000, sizeof(char));
    uint8_t nt;
    int64_t i; ubwt_count_t k, occ_k, uni_i = 0, uni_len = 1000;
    
    k = ubwt->C[nt_N] + uid;
    while (1) {
        if (uni_i == uni_len) {
            uni_len <<= 1;
            unipath = (char*)_err_realloc(unipath, uni_len * sizeof(char));
        }
        if((nt = ubwt_bwt_nt(ubwt, k)) >= nt_N) break;
        unipath[uni_i++] = "ACGTN"[nt];
        occ_k = ubwt_occ(ubwt, k, nt);
        k = ubwt->C[nt] + occ_k;
    }
    // store in struct, output for every N unipaths XXX
    fprintf(out, ">%lld_%lld\n", (long long)uid+1, (long long)uni_i);
    for (i = uni_i-1; i >= 0; --i) fprintf(out, "%c", unipath[i]);
    fprintf(out, "%c\n", unipath[0]);
    free(unipath);
}

void ubwt_gen_map1(ubwt_t *ubwt, ubwt_count_t uid) {
    ubwt_count_t k, occ_k; uint8_t nt;
    k = ubwt->C[nt_N]+uid;
    while (1) {
        ;
        if ((nt = ubwt_bwt_nt(ubwt, k)) >= nt_N) break;
        occ_k = ubwt_occ(ubwt, k, nt);
        k = ubwt->C[nt] + occ_k;
    }
    occ_k = ubwt_occ(ubwt, k, nt_N);
    ubwt->ubwt_map[occ_k] = uid; // occ_k+1 => i+1
}
void ubwt_gen_unipath1(ubwt_t *ubwt, ubwt_count_t uid, ubwt_count_t uni_s, char **out_unipath)
{
    char *unipath = (char*)_err_calloc(1000, sizeof(char));
    int64_t i; ubwt_count_t j, uni_i = 0, uni_len = 1000; uint8_t nt;
    ubwt_count_t k, occ_k;

    k = ubwt->C[nt_N] + uid;
    while (1) {
        if (uni_i == uni_len) {
            uni_len <<= 1;
            unipath = (char*)_err_realloc(unipath, uni_len * sizeof(char));
        }
        if ((nt = ubwt_bwt_nt(ubwt, k)) >= nt_N) break; 
        unipath[uni_i++] = "ACGTN"[nt];
        occ_k = ubwt_occ(ubwt, k, nt);
        k = ubwt->C[nt] + occ_k;
    }
    out_unipath[uid-uni_s] = (char*)_err_malloc(uni_i+1);
    j = 0;
    for (i = uni_i-1; i >= 0; --i)
        out_unipath[uid-uni_s][j++] = unipath[i];
    out_unipath[uid-uni_s][j] = '\0';
    free(unipath);
}

// TODO edge is 4-bits
int ubwt_get_out_edge(uint8_t edge, uint8_t edge_seq[4]) {
    int i, edge_n = 0; uint8_t seq[4] = {0x1, 0x2, 0x4, 0x8};
    for (i = 0; i < 4; ++i) {
        if (edge & seq[i]) edge_seq[edge_n++] = i;
    }
    return edge_n;
}

void ubwt_gen_gfa1(ubwt_t *ubwt, uint8_t *edge, ubwt_count_t uid, ubwt_count_t uni_s, ubwt_count_t klen, char **out_unipath, int *out_n, ubwt_count_t **out_uid)
{
    char *unipath = (char*)_err_calloc(1000, sizeof(char));
    int64_t i; ubwt_count_t j, uni_i = 0, uni_len = 1000; 
    ubwt_count_t k, occ_k;
    uint8_t nt, *query_seq = (uint8_t*)_err_calloc(klen, sizeof(uint8_t));
    ubwt_count_t out_n1=1000, out_uid1[4], query_k, query_l; uint8_t out_edge_seq[4];

    k = ubwt->C[nt_N] + uid;
    while (1) {
        if (uni_i == uni_len) {
            uni_len <<= 1;
            unipath = (char*)_err_realloc(unipath, uni_len * sizeof(char));
        }

        nt = ubwt_bwt_nt(ubwt, k);
        if (uni_i <= klen-2) { // last k-mer of current unipath
            query_seq[klen-2 - uni_i] = nt;
        } else if (uni_i == klen) { // construct k-mer on out-edge
            out_n1 = ubwt_get_out_edge(edge[k], out_edge_seq);
            for (j = 0; j < out_n1; ++j) {
                query_seq[klen-1] = out_edge_seq[j];
                if (ubwt_exact_match(ubwt, klen, query_seq, &query_k, &query_l) != 1)
                    err_fatal(__func__, "Error: no exact match of k-mer. (k=%d)", (int)klen);
                out_uid1[j] = ubwt_uid(ubwt, query_k);
            }
        }
        if (nt >= nt_N) break; 

        unipath[uni_i++] = "ACGTN"[nt];
        occ_k = ubwt_occ(ubwt, k, nt);
        k = ubwt->C[nt] + occ_k;
    }
    // dump unipath sequence
    out_unipath[uid-uni_s] = (char*)_err_malloc(uni_i+1);
    j = 0;
    for (i = uni_i-1; i >= 0; --i)
        out_unipath[uid-uni_s][j++] = unipath[i];
    out_unipath[uid-uni_s][j] = '\0';
    // dump link
    out_n[uid-uni_s] = out_n1;
    for (j = 0; j < out_n1; ++j) out_uid[uid-uni_s][j] = out_uid1[j];

    free(unipath); free(query_seq);
}

uint64_t THREAD_I;
pthread_rwlock_t RWLOCK;

static void *ubwt_thread_gen_map(void *aux)
{
    ubwt_thread_aux_t *a = (ubwt_thread_aux_t*)aux;

    ubwt_count_t i;
    while (1) {
        pthread_rwlock_wrlock(&RWLOCK);
        i = THREAD_I++;
        pthread_rwlock_unlock(&RWLOCK);
        if (i >= a->uni_c) break;

        ubwt_gen_map1(a->ubwt, i); // XXX
    }
    return 0;
}

static void *ubwt_thread_gen_unipath(void *aux)
{
    ubwt_thread_aux_t *a = (ubwt_thread_aux_t*)aux;

    ubwt_count_t i;
    while (1) {
        pthread_rwlock_wrlock(&RWLOCK);
        i = THREAD_I++;
        pthread_rwlock_unlock(&RWLOCK);
        if (i >= a->uni_c) break;

        ubwt_gen_unipath1(a->ubwt, i, a->uni_s, a->unipath); // XXX
    }
    return 0;
}

static void *ubwt_thread_gen_gfa(void *aux)
{
    ubwt_thread_aux_t *a = (ubwt_thread_aux_t*)aux;

    ubwt_count_t i;
    while (1) {
        pthread_rwlock_wrlock(&RWLOCK);
        i = THREAD_I++;
        pthread_rwlock_unlock(&RWLOCK);
        if (i >= a->uni_c) break;

        ubwt_gen_gfa1(a->ubwt, a->edge, i, a->uni_s, a->klen, a->unipath, a->out_n, a->out_uid); // XXX
    }
    return 0;
}

void ubwt_gen_map(ubwt_t *ubwt, ubwt_count_t uni_c, int t, int chunk_size)
{
    ubwt->ubwt_map = (ubwt_count_t*)_err_malloc(uni_c * sizeof(ubwt_count_t));
    ubwt_count_t i, k, occ_k; uint8_t nt;
    if (t > 1) {
        int i, chunk_n = chunk_size;
        ubwt_count_t uni_s = 0, remain_uni = uni_c;
        ubwt_thread_aux_t *aux = (ubwt_thread_aux_t*)_err_malloc(t * sizeof(ubwt_thread_aux_t));
        for (i = 0; i < t; ++i) {
            aux[i].tid = i;
            aux[i].ubwt = ubwt;
        }
        while (remain_uni > 0) {
            if (remain_uni < (unsigned long)chunk_n) chunk_n = (int)remain_uni;

            for (i = 0; i < t; ++i) {
                aux[i].uni_s = uni_s; aux[i].uni_c = uni_s + chunk_n;
            }

            THREAD_I = uni_s;
            pthread_t *tid = (pthread_t*)_err_malloc(t * sizeof(pthread_t)); pthread_attr_t attr;
            pthread_attr_init(&attr); pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            for (i = 0; i < t; ++i) {
                pthread_create(&tid[i], &attr, ubwt_thread_gen_map, aux+i);
            }
            for (i = 0; i < t; ++i) pthread_join(tid[i], 0);
            free(tid);

            uni_s += chunk_n;              
            remain_uni -= chunk_n;
        }
        free(aux);
    } else {
        for (i = 0; i < uni_c; ++i) {
            k = ubwt->C[nt_N]+i;
            while (1) {
                if ((nt = ubwt_bwt_nt(ubwt, k)) >= nt_N) break;
                occ_k = ubwt_occ(ubwt, k, nt);
                k = ubwt->C[nt] + occ_k;
            }
            occ_k = ubwt_occ(ubwt, k, nt_N);
            ubwt->ubwt_map[occ_k] = i; // occ_k+1 => i+1
            //printf("%d %d\n", occ_k+1, i+1);
        }
    }
}

// fasta and binary format output
void ubwt_gen_unipath(ubwt_t *ubwt, ubwt_count_t uni_c, FILE *out, int b_out, int t, int chunk_size)
{
    char f[2][10] = {"fasta", "binary"};
    err_func_format_printf(__func__, "Generating unipath in %s format ...\n", f[b_out]);
    ubwt_count_t max_len=0;
    if (t > 1) {
        int i, chunk_n = chunk_size;
        ubwt_count_t uni_s = 0, remain_uni = uni_c;
        char **unipath = (char**)_err_malloc(chunk_n * sizeof(char*));
        ubwt_thread_aux_t *aux = (ubwt_thread_aux_t*)_err_malloc(t * sizeof(ubwt_thread_aux_t));
        for (i = 0; i < t; ++i) {
            aux[i].tid = i;
            aux[i].ubwt = ubwt;
            aux[i].unipath = unipath;
        }
        while (remain_uni > 0) {
            if (remain_uni < (unsigned long)chunk_n) chunk_n = (int)remain_uni;

            for (i = 0; i < t; ++i) {
                aux[i].uni_s = uni_s; aux[i].uni_c = uni_s + chunk_n;
            }
            
            THREAD_I = uni_s;
            pthread_t *tid = (pthread_t*)_err_malloc(t * sizeof(pthread_t)); pthread_attr_t attr;
            pthread_attr_init(&attr); pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            for (i = 0; i < t; ++i) {
                pthread_create(&tid[i], &attr, ubwt_thread_gen_unipath, aux+i);
            }
            for (i = 0; i < t; ++i) pthread_join(tid[i], 0);
            free(tid);

            // output chunk_n unipaths
            for (i = 0; i < (int)chunk_n; ++i) {
                ubwt_count_t len = strlen(unipath[i]);
                if (len > max_len) max_len = len;
                if (b_out) {
                } else fprintf(out, ">%lld_%lld\n%s\n", (long long)uni_s+i+1, (long long)strlen(unipath[i]), unipath[i]);
            }
            for(i = 0; i < chunk_n; ++i) free(unipath[i]);
            uni_s += chunk_n;              
            remain_uni -= chunk_n;
        }
        free(aux); free(unipath);
    } else {
        char *unipath = (char*)_err_calloc(1000, sizeof(char));
        ubwt_count_t uni_i, i, j, k, occ_k, uni_len = 1000; uint8_t nt;

        for (i = 0; i < uni_c; ++i) {
            uni_i = 0;
            k = ubwt->C[nt_N]+i; // i-th unipath
            while (1) {
                if (uni_i == uni_len) {
                    uni_len <<= 1;
                    unipath = (char*)_err_realloc(unipath, uni_len * sizeof(char));
                }
                if ((nt = ubwt_bwt_nt(ubwt,k)) >= nt_N) break; 
                unipath[uni_i++] = "ACGTN"[nt];
                occ_k = ubwt_occ(ubwt, k, nt);
                k = ubwt->C[nt] + occ_k;
                //unipath[uni_i++] = "ACGTN"[ubwt_bstr[k]];
                //occ_k = ubwt_occ(ubwt, k, ubwt_bstr[k]);
                //k = ubwt->C[ubwt_bstr[k]] + occ_k;
                //if (ubwt_bstr[k] >= nt_N) break; 
            }
            if (b_out) {
            } else {
            fprintf(out, ">%lld_%lld\n", (long long)i+1, (long long)uni_i);
            for (j = uni_i-1; j != 0; --j)
                fprintf(out, "%c", unipath[j]);
            fprintf(out, "%c\n", unipath[0]);
            }
            if (uni_i > max_len) max_len = uni_i;
        }
        free(unipath);
    }
    err_func_format_printf(__func__, "Generating unipath in %s format done.\n", f[b_out]);
    // maximum unipath length
    fprintf(stderr, "MAX: %lld\n", (long long)max_len);
}

// GFA format output
void ubwt_gen_gfa(ubwt_t *ubwt, uint8_t *edge, ubwt_count_t uni_c, ubwt_count_t klen, FILE *out, int t, int chunk_size) {
    err_func_format_printf(__func__, "Generating unipath in GFA format ...\n");
    ubwt_count_t max_len=0;
    if (t > 1) {
        int i, j, chunk_n = chunk_size;
        ubwt_count_t uni_s = 0, remain_uni = uni_c;
        char **unipath = (char**)_err_malloc(chunk_n * sizeof(char*));
        int *out_n = (int*)_err_malloc(chunk_n * sizeof(int));
        ubwt_count_t **out_uid = (ubwt_count_t**)_err_malloc(chunk_n * sizeof(ubwt_count_t*));
        for (i = 0; i < chunk_n; ++i) out_uid[i] = (ubwt_count_t*)_err_malloc(4 * sizeof(ubwt_count_t));
        ubwt_thread_aux_t *aux = (ubwt_thread_aux_t*)_err_malloc(t * sizeof(ubwt_thread_aux_t));

        for (i = 0; i < t; ++i) {
            aux[i].tid = i;
            aux[i].klen = klen;
            aux[i].ubwt = ubwt; aux[i].edge = edge;
            aux[i].unipath = unipath;
            aux[i].out_n = out_n; aux[i].out_uid = out_uid;
        }
        while (remain_uni > 0) {
            if (remain_uni < (unsigned long)chunk_n) chunk_n = (int)remain_uni;

            for (i = 0; i < t; ++i) {
                aux[i].uni_s = uni_s; aux[i].uni_c = uni_s + chunk_n;
            }
            
            THREAD_I = uni_s;
            pthread_t *tid = (pthread_t*)_err_malloc(t * sizeof(pthread_t)); pthread_attr_t attr;
            pthread_attr_init(&attr); pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
            for (i = 0; i < t; ++i) {
                pthread_create(&tid[i], &attr, ubwt_thread_gen_gfa, aux+i);
            }
            for (i = 0; i < t; ++i) pthread_join(tid[i], 0);
            free(tid);

            // output chunk_n unipaths and links
            for (i = 0; i < (int)chunk_n; ++i) {
                ubwt_count_t len = strlen(unipath[i]);
                if (len > max_len) max_len = len;
                fprintf(out, "S\t%lld\t%s\tLN:i:%lld\n", (long long)uni_s+i+1, unipath[i], (long long)len);
                for (j = 0; j < out_n[i]; ++j) {
                    fprintf(out, "L\t%lld\t+\t%lld\t+\t%lldM\n", (long long)uni_s+i+1, (long long)out_uid[i][j]+1, (long long)klen-1);
                }
            }
            for(i = 0; i < chunk_n; ++i) free(unipath[i]);
            uni_s += chunk_n;              
            remain_uni -= chunk_n;
        }
        free(aux); free(unipath); free(out_n); 
        for (i = 0; i < chunk_size; ++i) free(out_uid[i]); free(out_uid);
    } else {
        char *unipath = (char*)_err_calloc(1000, sizeof(char));
        ubwt_count_t uni_i, i, j, k, occ_k, uni_len = 1000;

        uint8_t nt, *query_seq = (uint8_t*)_err_calloc(klen, sizeof(uint8_t)); ubwt_count_t query_k, query_l, out_edge_n, out_uid[4]; uint8_t out_edge_seq[4];
        for (i = 0; i < uni_c; ++i) {
            uni_i = 0;
            k = ubwt->C[nt_N]+i; // i-th unipath
            out_edge_n = 0;
            while (1) {
                if (uni_i == uni_len) {
                    uni_len <<= 1;
                    unipath = (char*)_err_realloc(unipath, uni_len * sizeof(char));
                }
                nt = ubwt_bwt_nt(ubwt, k);
                if (uni_i <= klen-2) { // last k-mer of current unipath
                    query_seq[klen-2 - uni_i] = nt;
                } else if (uni_i == klen) { // construct k-mer on out-edge
                    out_edge_n = ubwt_get_out_edge(edge[k], out_edge_seq);
                    for (j = 0; j < out_edge_n; ++j) {
                        query_seq[klen-1] = out_edge_seq[j];
                        if (ubwt_exact_match(ubwt, klen, query_seq, &query_k, &query_l) != 1)
                            err_fatal(__func__, "Error: no exact match of k-mer. (k=%d)", (int)klen);
                        out_uid[j] = ubwt_uid(ubwt, query_k);
                    }
                }
                if (nt >= nt_N) break; 
                unipath[uni_i++] = "ACGTN"[nt];
                occ_k = ubwt_occ(ubwt, k, nt);
                k = ubwt->C[nt] + occ_k;
            }
            fprintf(out, "S\t%lld\t", (long long)i+1); // Segment
            for (j = uni_i-1; j != 0; --j)
                fprintf(out, "%c", unipath[j]);
            fprintf(out, "%c\tLN:i:%lld\n", unipath[0], (long long)uni_i);
            for (j = 0; j < out_edge_n; ++j) {
                fprintf(out, "L\t%lld\t+\t%lld\t+\t%dM\n", (long long)i+1, (long long)out_uid[j]+1, (int)klen-1); // Link
            }
            if (uni_i > max_len) max_len = uni_i;
        }
        free(unipath); free(query_seq);
    }
    err_func_format_printf(__func__, "Generating unipath in GFA format done.\n");
    // maximum unipath length
    fprintf(stderr, "MAX: %lld\n", (long long)max_len);
}

uint8_t *ubwt_read_seq(char *fn, ubwt_count_t *seq_l)
{
    uint8_t *bseq = (uint8_t*)_err_malloc(100 * sizeof(uint8_t));
    char ch; ubwt_count_t i = 0, m = 100;
    FILE *fp = xopen(fn, "r");
    while ((ch = fgetc(fp)) != EOF) {
        if (isspace(ch) || ch == '\n') continue;
        if (i == m) {
            m <<= 1;
            bseq = (uint8_t*)_err_realloc(bseq, m * sizeof(uint8_t));
        }
        bseq[i++] = nst_nt4_table[(int)ch];
    }
    *seq_l = i;
    err_fclose(fp);
    return bseq;
}

uint8_t *ubwt_read_bwt_str(char *fn, int input_b, ubwt_count_t *ubwt_l)
{
    err_func_format_printf(__func__, "Reading ubwt string from %s ...\n", fn);
    uint64_t bwt_int;
    uint8_t *ubwt_bstr; ubwt_count_t bwt_i, i; int j;
    if (input_b) { // binary file, 4-bit per bp, first 64-bit: length
        FILE *fp = xopen(fn, "rb");
        err_fread_noeof(ubwt_l, sizeof(ubwt_count_t), 1, fp);
        ubwt_bstr = (uint8_t*)_err_malloc(*ubwt_l * sizeof(uint8_t));
        bwt_i = 0;
        for (i = 0; i < *ubwt_l / 16; ++i) {
            err_fread_noeof(&bwt_int, sizeof(uint64_t), 1, fp);
            for (j = 15; j >= 0; --j) {
                ubwt_bstr[bwt_i++] = (bwt_int >> (4*j)) & 0x7;
            }
        }
        // last one bwt_int
        if (*ubwt_l != (*ubwt_l/16)*16) {
            err_fread_noeof(&bwt_int, sizeof(uint64_t), 1, fp);
            j = *ubwt_l % 16;
            for (i = 0; i < *ubwt_l % 16; ++i) {
                ubwt_bstr[bwt_i+j-i-1] = bwt_int & 0x7;
                bwt_int >>= 4;
            }
        }
        err_fclose(fp);
    } else {       // plain text
        // read bwt-str
        ubwt_bstr = ubwt_read_seq(fn, ubwt_l);
    }
    err_func_format_printf(__func__, "Reading ubwt string done.\n");
    return ubwt_bstr;
}

uint8_t *ubwt_read_edge(char *fn, ubwt_count_t *edge_l)
{
    err_func_format_printf(__func__, "Reading ubwt edge sequence from %s ...\n", fn);
    uint64_t edge_int;
    uint8_t *edge; ubwt_count_t edge_i, i; int j;
    // binary file, 8-bit per bp, first 64-bit: length
    // TGCA   TGCA
    // ---- | ----
    //  in     out
    FILE *fp = xopen(fn, "rb");
    err_fread_noeof(edge_l, sizeof(ubwt_count_t), 1, fp);
    edge = (uint8_t*)_err_malloc(*edge_l * sizeof(uint8_t));
    edge_i = 0;
    for (i = 0; i < *edge_l / 8; ++i) {
        err_fread_noeof(&edge_int, sizeof(uint64_t), 1, fp);
        for (j = 7; j >= 0; --j) {
            edge[edge_i++] = (edge_int >> (8*j)) & 0xff;
        }
    }
    // last one bwt_int
    if (*edge_l != (*edge_l/8)*8) {
        err_fread_noeof(&edge_int, sizeof(uint64_t), 1, fp);
        j = *edge_l % 8;
        for (i = 0; i < *edge_l % 8; ++i) {
            edge[edge_i+j-i-1] = edge_int & 0xff;
            edge_int >>= 8;
        }
    }
    err_fclose(fp);
    err_func_format_printf(__func__, "Reading ubwt edge sequence done.\n");
    return edge;
}

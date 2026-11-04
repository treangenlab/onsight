/* pore_model.c: the pore model and reference signal; see pore_model.h. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pore_model.h"
#include "utils.h"

static char *next_field(char **stringp, const char *delim)
{
    char* start = *stringp;
    char* p;

    p = (start != NULL) ? strpbrk(start, delim) : NULL;

    if (p == NULL) {
        *stringp = NULL;
    } else {
        *p = '\0';
        *stringp = p + 1;
    }

    return start;
}

static uint32_t rev_complement(uint32_t x, const short k) {
    uint32_t y = 0;
    for (short i = 0; i < k; i++) {
        y = (y<<2) | ((x&3)^3);
        x >>= 2;
    }
    return y;
}

static int compare_value_index_pairs(const void* a, const void* b) {
    float diff = ((ri_porei_t*)a)->pore_val - ((ri_porei_t*)b)->pore_val;
    return (diff < 0)?-1:(diff > 0);
}

static ri_porei_t* create_sorted_pairs(const ri_pore_t* pore) {
    double sum = 0, sum2 = 0, mean, std_dev;

    // First pass: Calculate mean and standard deviation
    for (uint32_t i = 0; i < pore->n_pore_vals; i++) {
        sum += pore->pore_vals[i];
        sum2 += pore->pore_vals[i] * pore->pore_vals[i];
    }
    mean = sum / pore->n_pore_vals;
    std_dev = sqrt(sum2 / pore->n_pore_vals - mean * mean);

    // Allocate memory for pairs
    ri_porei_t* pairs = (ri_porei_t*)malloc(pore->n_pore_vals * sizeof(ri_porei_t));

    // Second pass: Normalize values and create pairs
    for (uint32_t i = 0; i < pore->n_pore_vals; i++) {
        pairs[i].pore_val = (pore->pore_vals[i] - mean) / std_dev; // Normalization
        pairs[i].ind = i;
        pairs[i].rev_ind = rev_complement(i, pore->k);
    }

    // Sort pairs based on normalized pore_val
    qsort(pairs, pore->n_pore_vals, sizeof(ri_porei_t), compare_value_index_pairs);
    return pairs;
}

void ri_pore_load(const char* fpore, const short k, const short lev_col, ri_pore_t* pore){
	FILE* fp = fopen(fpore, "r");
	if(fp == NULL){
		fprintf(stderr, "Error: cannot open file %s\n", fpore);
		return;
	}

	pore->pore_vals = (float*)malloc(sizeof(float) * (1U<<(2*k)));
	char line[1024];
	char* token;
	int i = 0;
    double mean = 0, std_dev = 0, sum = 0, sum2 = 0;

	while(fgets(line, sizeof(line), fp) != NULL){
		if(!strncmp(line, "kmer", 4)) continue;
    	char* rest = line;
    	int j = 0;
		while((token = next_field(&rest, "\t")) != NULL){
			if(j++ == lev_col){
				float value;
				if (sscanf(token, "%f", &value) == 1) {
					pore->pore_vals[i] = value;
                    sum += value;
		            sum2 += value*value;
				} else {
					fprintf(stderr, "Error: cannot convert '%s' to float\n", token);
					free(pore->pore_vals); pore->pore_vals = NULL;
					fclose(fp);
					return;
				}
				break;
			}
		}
		i++;
	}
	fclose(fp);

    mean = sum/i;
	std_dev = sqrt(sum2/i - (mean)*(mean));

    for(int j = 0; j < i; ++j) pore->pore_vals[j] = (pore->pore_vals[j]-mean)/std_dev;

    pore->n_pore_vals = 1U<<(2*k);
    pore->k = k;
    pore->pore_inds = create_sorted_pairs(pore);
}


void ri_reference_signal(const char *str,
				   int len,
				   const ri_pore_t* pore,
				   const int k,
				   const int strand, 
				   uint32_t* s_len,
				   float* s_values){

	int i, j, l, pos, n = 0;
	uint64_t mask = (1ULL<<2*k) - 1, kmer = 0;

	for (i = l = j = n = 0; i < len; ++i) {
		if(strand) pos = len-i-1;
		else pos = i;
		int c = seq_nt4_table[(uint8_t)str[pos]];
		if (c < 4) { // not an ambiguous base
			if(!strand) kmer = (kmer << 2 | c) & mask;    // forward k-mer
			// else kmer = (kmer >> 2) | (3ULL^c) << shift1; // reverse k-mer
			else kmer = ((kmer << 2) | (3ULL^c)) & mask; // reverse k-mer
		}

		if(i+1 < k) continue;

		s_values[j++] = pore->pore_vals[kmer];;
	}

	*s_len = j;
}


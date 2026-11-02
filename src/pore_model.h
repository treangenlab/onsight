/* pore_model.h: the pore model and the reference signal of a sequence (index build).
 *
 * The pore model gives the expected current of each k-mer in the pore (6 bases for R9.4.1, 9 for R10.4.1), normalized
 * over all k-mers. The reference signal of a sequence is the expected level of each of its k-mers in order, on either
 * strand.
 *
 * In: the pore model file (-p, levels in column --level_column); a signature's sequence. Out: its expected event
 * levels. State: the loaded model (ri_pore_t), kept in the index.
 */
#ifndef PORE_MODEL_H
#define PORE_MODEL_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ri_porei_s{float pore_val; unsigned int ind; unsigned int rev_ind;} ri_porei_t;
typedef struct ri_pore_s{ri_porei_t* pore_inds; float* pore_vals; unsigned int n_pore_vals; short k; float max_val, min_val;} ri_pore_t;

/* Loads the level of every k-mer from column lev_col of fpore; pore->pore_vals is NULL on error. */
void ri_pore_load(const char* fpore, const short k, const short lev_col, ri_pore_t* pore);

/**
 * Converts the sequence into its expected event values
 *
 * @param str		sequence to convert to its expected event values
 * @param len		length of the $str
 * @param pore		a pore model including expected event values for each possible k-mer
 * @param pore_kmer	k-mer size of each event
 * @param strand	directin of strand. 1 is forward, 0 is reverse direction.
 * @param s_len		length of $s_values
 * @param s_values	expected event values of each k-mer in $str
 */
void ri_reference_signal(const char *str,
				   int len,
				   const ri_pore_t* pore,
				   const int pore_kmer,
				   const int strand,
				   uint32_t* s_len,
				   float* s_values);

#ifdef __cplusplus
}
#endif
#endif // PORE_MODEL_H

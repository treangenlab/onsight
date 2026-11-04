#ifndef SIGNAL_INPUT_H
#define SIGNAL_INPUT_H

/* signal_input.h: each read's raw signal from FAST5, POD5 and SLOW5/BLOW5 files.
 *
 * In: signal files, or folders of them. Out: per read, its name and raw signal in pA, samples outside 30-200 pA
 * dropped. State: the open file and the next read in it.
 */

#include "utils.h"
#ifndef NHDF5RH
#include "hdf5_tools.hpp"
#endif
#ifndef NPOD5RH
#include "pod5_format/c_api.h"
#endif
#ifndef NSLOW5RH
#include <slow5/slow5.h>
#include <slow5/slow5_mt.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct ri_signal_s{
	uint32_t rid, l_sig; //read id and length of the signal values
	char *name; //name of the read

	uint64_t offset; // offset in ri_idx_t::S

	// float dig, ran, offset; //digitalisation, range, offset
	float* sig; //signal values of a read
#ifdef ONS_PROF
	uint64_t prof_t_read; // ns to read and convert it
#endif
} ri_signal_t;

typedef struct { size_t n, m; ri_signal_t **a; } ri_signal_v;

typedef struct ri_signal_file_s {
	// gzFile fp;
	// kseq_t *ks;
	// ri_signal_t s;
	int num_read; //Number of reads
	int cur_read; //Number of processed reads by RawHash (shows the id of the next read to process)
	
	//HDF5-related
	char** raw_path; //List of paths to raw values
	char** ch_path; //List of paths to channel
	#ifndef NHDF5RH
	hdf5_tools::File* fp; //FAST5 file pointer
	#endif
	
	//POD5-related
	unsigned long int pod5_batch_count;
	unsigned long int pod5_batch;
	unsigned long int pod5_row_count;
	unsigned long int pod5_row;
	#ifndef NPOD5RH
	Pod5ReadRecordBatch_t* batch;
	Pod5FileReader_t* pp; //POD5 file pointer
	#endif

	#ifndef NSLOW5RH
	slow5_mt_t* slow5_mt;
	slow5_batch_t* slow5_batch;
	// slow5_rec_t** slow5_rec;
	slow5_file_t* sp; //SLOW5 file pointer
	#endif
} ri_signal_file_t;

/* Closes the signal file and frees its buffers. */
void ri_signal_file_close(ri_signal_file_t *fp);

/* Opens a FAST5, POD5, SLOW5 or BLOW5 file; NULL if it cannot be opened. */
ri_signal_file_t *ri_signal_file_open(const char *fn, int io_n_threads);

/* The next read's name and signal (s->name, s->sig, s->l_sig); samples outside 30-200 pA are dropped. */
void ri_read_signal(ri_signal_file_t* fp, ri_signal_t* s);

/* A itself if it is a signal file, else every .fast5, .pod5, .slow5 and .blow5 file under folder A, recursively,
 * appended to fnames. */
void ri_find_signal_files(const char *A, ri_char_v *fnames);

#ifdef __cplusplus
}
#endif
#endif //SIGNAL_INPUT_H
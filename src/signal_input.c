/* signal_input.c: each read's raw signal; see signal_input.h. */
#include "signal_input.h"
#include "rh_kvec.h"
#include <math.h>
#ifndef NHDF5RH
#include <cstring>
#include <cassert>
#endif
#include <errno.h>
#include <assert.h>
#include <sys/stat.h>
#include <dirent.h>

#ifndef NHDF5RH
static inline ri_signal_file_t *signal_file_open_fast5(const char *fn)
{
	ri_signal_file_t *fp;

	hdf5_tools::File* fast5_file = new hdf5_tools::File();
	fast5_file->open(std::string(fn));
	if (!fast5_file->is_open()) return 0;

	fp = (ri_signal_file_t*)calloc(1, sizeof(ri_signal_file_t));
	#ifndef NPOD5RH
	fp->pp = NULL;
	#endif

	#ifndef NSLOW5RH
	fp->sp = NULL;
	#endif

	fp->fp = fast5_file;

	bool is_single = false;
	std::vector<std::string> fast5_file_groups = fast5_file->list_group("/");
	fp->num_read = fast5_file_groups.size();
	fp->ch_path = (char**)calloc(fp->num_read, sizeof(char*));
	fp->raw_path = (char**)calloc(fp->num_read, sizeof(char*));

	for (std::string &group : fast5_file_groups) {
		if (group == "Raw") {
			is_single = true;
			break;
		}
	}

	std::string raw_path;
	std::string ch_path;
	int i = 0;

	if (is_single) {
		ch_path = "/UniqueGlobalKey/channel_id";
		for (std::string &read : fast5_file->list_group("/Raw/Reads")) {
			raw_path = "/Raw/Reads/" + read;
			if(i == fp->num_read){
				fprintf(stderr, "ERROR: More reads than previously predicted (%d). Stopped reading the reads here.\n", fp->num_read);
				break;
			}
			fp->ch_path[i] = strdup(ch_path.c_str());
			fp->raw_path[i++] = strdup(raw_path.c_str());
		}
	} else {
		for (std::string &read : fast5_file_groups) {
			raw_path = "/" + read + "/Raw";
			ch_path = "/" + read + "/channel_id";
			fp->ch_path[i] = strdup(ch_path.c_str());
			fp->raw_path[i++] = strdup(raw_path.c_str());
		}
	}

	fp->num_read = i;
	fp->cur_read = 0;
	return fp;
}
#endif

#ifndef NPOD5RH
static inline ri_signal_file_t *signal_file_open_pod5(const char *fn){
	pod5_init();

	ri_signal_file_t *fp;

	// Read the file rather than map it: the pages of a mapped file count toward memory as they are read, which on a
	// POD5 file of tens of GB adds tens of GB.
	Pod5ReaderOptions_t opts = {1}; // force_disable_file_mapping
	Pod5FileReader_t* pod5_file = pod5_open_file_options(fn, &opts);

	if (!pod5_file) {
		fprintf(stderr, "ERROR: Failed to open file (%s) %s\n", fn, pod5_get_error_string());
		return 0;
    }

	long unsigned int batch_count = 0;
    if (pod5_get_read_batch_count(&batch_count, pod5_file) != POD5_OK) {
        fprintf(stderr, "Failed to query batch count: %s\n", pod5_get_error_string());
		pod5_close_and_free_reader(pod5_file);
		return 0;
    }

	fp = (ri_signal_file_t*)calloc(1, sizeof(ri_signal_file_t));
	fp->pp = pod5_file;

	#ifndef NHDF5RH
	fp->fp = NULL;
	#endif

	#ifndef NSLOW5RH
	fp->sp = NULL;
	#endif

	fp->num_read = batch_count; //num_read is the number of batches for pod5
	fp->cur_read = 0; //cur_read is cur_batch for pod5

	Pod5ReadRecordBatch_t* batch = NULL;
	if(pod5_get_read_batch(&batch, pod5_file, fp->cur_read) != POD5_OK){
		fprintf(stderr, "Failed to get batch: %s\n", pod5_get_error_string());
		pod5_close_and_free_reader(pod5_file);
		free(fp);
		return 0;
	}

	long unsigned int batch_row_count = 0;
	if(pod5_get_read_batch_row_count(&batch_row_count, batch) != POD5_OK) {
		fprintf(stderr, "Failed to get batch row count\n");
		pod5_close_and_free_reader(pod5_file);
		pod5_free_read_batch(batch);
		free(fp);
		return 0;
	}

	fp->pod5_row_count = batch_row_count;
	fp->pod5_row = 0;
	
	fp->batch = batch;
	fp->num_read = batch_count;

	return fp;
}
#endif

#ifndef NSLOW5RH
static inline ri_signal_file_t *signal_file_open_slow5(const char *fn, int io_n_threads){
	
	ri_signal_file_t *fp;

	//open the SLOW5 file
	slow5_file_t *sp = slow5_open(fn, "r");
	if(sp==NULL){
       fprintf(stderr, "ERROR: Failed to open file %s\n", fn);
	   return 0;
    }

	fp = (ri_signal_file_t*)calloc(1, sizeof(ri_signal_file_t));

	#ifndef NHDF5RH
	fp->fp = NULL;
	#endif

	#ifndef NPOD5RH
	fp->pp = NULL;
	#endif
	
	fp->sp = sp;
	fp->slow5_mt = slow5_init_mt(io_n_threads, fp->sp);
	fp->slow5_batch = slow5_init_batch(4096);   /* reads are taken 4096 at a time */
	fp->num_read = slow5_get_next_batch(fp->slow5_mt, fp->slow5_batch, 4096);
	fp->cur_read = 0;

	if(fp->num_read <= 0){
		fp->num_read = 0;
		fprintf(stderr, "ERROR: Failed to read the first batch of reads from %s\n", fn);
		slow5_close(fp->sp);
		free(fp);
		return 0;
	}

	return fp;
}
#endif

ri_signal_file_t *signal_file_open(const char *fn, int io_n_threads){
	if (strstr(fn, ".fast5")) {
		#ifndef NHDF5RH
		return signal_file_open_fast5(fn);
		#endif
	} else if (strstr(fn, ".pod5") || strstr(fn, ".pod")) {
		#ifndef NPOD5RH
		return signal_file_open_pod5(fn);
		#endif
	} else if (strstr(fn, ".slow5") || strstr(fn, ".blow5")) {
		#ifndef NSLOW5RH
		return signal_file_open_slow5(fn, io_n_threads);
		#endif
	}

	return 0;
}

void ri_signal_file_close(ri_signal_file_t *fp)
{
	if(!fp) return;

	#ifndef NHDF5RH
	if(fp->fp){
		fp->fp->close();
		for(int i = 0; i < fp->num_read; ++i){
			if(fp->ch_path[i])free(fp->ch_path[i]);
			if(fp->raw_path[i])free(fp->raw_path[i]);
		}
		free(fp->ch_path);
		free(fp->raw_path);
		delete fp->fp;
	}
	#endif
	#ifndef NPOD5RH
	if(fp->pp){
		pod5_close_and_free_reader(fp->pp);
	}
	#endif

	#ifndef NSLOW5RH
	if(fp->sp){
		slow5_free_batch(fp->slow5_batch);
    	slow5_free_mt(fp->slow5_mt);
		slow5_close(fp->sp);
	}
	#endif

	free(fp);
}

ri_signal_file_t *ri_signal_file_open(const char *fn, int io_n_threads)
{
	ri_signal_file_t *fp = signal_file_open(fn, io_n_threads);
	if (fp == 0) fprintf(stderr, "ERROR: failed to open file '%s': %s\n", fn, strerror(errno));
	return fp;
}

//Check if the input const char* A is a directory
int is_dir(const char *A)
{
	struct stat st;
	if (stat(A, &st) == -1) return 0;
	return S_ISDIR(st.st_mode);
}

void ri_find_signal_files(const char *A, ri_char_v *fnames)
{
	if (!is_dir(A)) {
		if (strstr(A, ".fast5") || strstr(A, ".pod5") || strstr(A, ".pod") || strstr(A, ".slow5") || strstr(A, ".blow5")) {
			char** cur_fname;
			rh_kv_pushp(char*, 0, *fnames, &cur_fname);
			(*cur_fname) = strdup(A);
		}
		return;
	}

	DIR *dir;
	struct dirent *ent;
	if ((dir = opendir(A)) != NULL) {
		while ((ent = readdir(dir)) != NULL) {
			char *tmp = (char*)malloc(strlen(A) + strlen(ent->d_name) + 2);
			sprintf(tmp, "%s/%s", A, ent->d_name);
			if (is_dir(tmp)) {
				if (strcmp(ent->d_name, ".") && strcmp(ent->d_name, ".."))
					ri_find_signal_files(tmp, fnames);
			} else {
				if (strstr(ent->d_name, ".fast5") || strstr(ent->d_name, ".pod5") || strstr(ent->d_name, ".pod") || strstr(ent->d_name, ".slow5") || strstr(ent->d_name, ".blow5")) {
					char** cur_fname;
					rh_kv_pushp(char*, 0, *fnames, &cur_fname);
					(*cur_fname) = strdup(tmp);
				}
			}
			free(tmp);
		}
		closedir(dir);
	}
}

#ifndef NHDF5RH
static inline void read_signal_fast5(ri_signal_file_t* fp, ri_signal_t* s){

	if(fp->cur_read >= fp->num_read) return;
	
	s->name = 0;
	for (auto a : fp->fp->get_attr_map(fp->raw_path[fp->cur_read])) {
		if (a.first == "read_id") {
			s->name = strdup(a.second.c_str());
		}
	}

	assert(s->name);

	float dig = 0, ran = 0, offset = 0;

	for (auto a : fp->fp->get_attr_map(fp->ch_path[fp->cur_read])) {
		if (a.first == "digitisation") {
			dig = atof(a.second.c_str());
		} else if (a.first == "range") {
			ran = atof(a.second.c_str());
		} else if (a.first == "offset") {
			offset = atof(a.second.c_str());
		}
	}

	std::string sig_path = std::string(fp->raw_path[fp->cur_read]) + "/Signal";
	std::vector<int16_t> sig;
	try {
		fp->fp->read(sig_path, sig);
	} catch (hdf5_tools::Exception &e) {
		if (H5Zfilter_avail(32020) <= 0) // 32020: ONT's VBZ filter
			fprintf(stderr, "[ERROR] cannot read the signal of read %s: VBZ plugin not found; run 'conda activate onsight' first, or set HDF5_PLUGIN_PATH to the plugin's folder\n", s->name);
		else
			fprintf(stderr, "[ERROR] cannot read the signal of read %s: %s\n", s->name, e.what());
		exit(1);
	}
	// convert to pA
	uint32_t l_sig = 0;
	float scale = ran/dig;
	float pa = 0;
	for (size_t i = 0; i < sig.size(); i++) {
		pa = (sig[i]+offset)*scale;
		if (pa > 30.0f && pa < 200.0f) {
			sig[l_sig++] = pa;
		}
	}

	s->sig = (float*)calloc(l_sig, sizeof(float));
	s->l_sig = l_sig;
	std::copy(sig.begin(), sig.begin() + l_sig, s->sig);
	fp->cur_read++;
}
#endif

#ifndef NPOD5RH
static inline void read_signal_pod5(ri_signal_file_t* fp, ri_signal_t* s){

	if(fp->cur_read >= fp->num_read) return;
	if(fp->pod5_row >= fp->pod5_row_count) return;

	pod5_init();

    // Open the file ready for walking:
    Pod5FileReader_t* pod5_file = fp->pp;

    if (!pod5_file) {
        fprintf(stderr, "Failed to open file %s\n", pod5_get_error_string());
		return;
    }

	uint16_t read_table_version = 0;
	ReadBatchRowInfo_t read_data;
	if (pod5_get_read_batch_row_info_data(fp->batch, fp->pod5_row, READ_BATCH_ROW_INFO_VERSION, &read_data, &read_table_version) != POD5_OK) {
		fprintf(stderr, "Failed to get read %lu\n", fp->pod5_row);
		return;
	}

	//Retrieve global information for the run
	RunInfoDictData_t* run_info_data;
	if (pod5_get_run_info(fp->batch, read_data.run_info, &run_info_data) != POD5_OK) {
		fprintf(stderr, "Failed to get Run Info %lu %s\n", fp->pod5_row, pod5_get_error_string());
		pod5_release_run_info(run_info_data);
		return;
	}

	int16_t *sig = (int16_t*)malloc(read_data.num_samples * sizeof(int16_t));
	float *sigF = (float*)malloc(read_data.num_samples * sizeof(float));
	
	if (pod5_get_read_complete_signal(pod5_file, fp->batch, fp->pod5_row, read_data.num_samples, sig) != POD5_OK) {
		fprintf(stderr, "Failed to get read %lu signal: %s\n", fp->pod5_row, pod5_get_error_string());
		pod5_release_run_info(run_info_data);
		if(sig) free(sig);
		return;
	}

	char read_id_tmp[37];
	pod5_format_read_id(read_data.read_id, read_id_tmp);

	uint32_t l_sig = 0;
	float pa = 0.0f;
	for(uint64_t i = 0; i < read_data.num_samples; i++){
		pa = (sig[i]+read_data.calibration_offset)*read_data.calibration_scale;
		if (pa > 30 && pa < 200) {
			sigF[l_sig++] = pa;
		}
	}

	free(sig);
	s->sig = (float*)calloc(l_sig, sizeof(float));
	s->l_sig = l_sig;
	memcpy(s->sig, sigF, l_sig*sizeof(float));
	s->name = strdup(read_id_tmp);
	free(sigF);

	pod5_release_run_info(run_info_data);

	fp->pod5_row++;
	if(fp->pod5_row >= fp->pod5_row_count){

		if (pod5_free_read_batch(fp->batch) != POD5_OK) {
            fprintf(stderr, "Failed to release batch\n");
			// pod5_close_and_free_reader(fp->pp);
			return;
        }

		fp->cur_read++;
		if(fp->cur_read < fp->num_read){
			Pod5ReadRecordBatch_t* batch = NULL;
			if(pod5_get_read_batch(&batch, fp->pp, fp->cur_read) != POD5_OK){
				fprintf(stderr, "Failed to get batch: %s\n", pod5_get_error_string());
				return;
			}

			long unsigned int batch_row_count = 0;
			if(pod5_get_read_batch_row_count(&batch_row_count, batch) != POD5_OK) {
				fprintf(stderr, "Failed to get batch row count\n");
				return;
			}

			fp->pod5_row_count = batch_row_count;
			fp->pod5_row = 0;
			
			fp->batch = batch;
		}
	}
}
#endif

#ifndef NSLOW5RH
static inline void read_signal_slow5(ri_signal_file_t* fp, ri_signal_t* s){
	
	if(fp->cur_read >= fp->num_read) return;

	slow5_file_t *sp = fp->sp;
	slow5_rec_t *rec = NULL;

	rec = (fp->slow5_batch->slow5_rec)[fp->cur_read];

	s->name = strdup(rec->read_id);
	float *sigF = (float*)malloc(rec->len_raw_signal * sizeof(float));
	uint32_t l_sig = 0;
	float pa = 0.0f;
	float scale = rec->range/rec->digitisation;
	
	for(int i = 0; i < rec->len_raw_signal; ++i){
		pa = (rec->raw_signal[i]+rec->offset)*scale;
		if (pa > 30.0f && pa < 200.0f) {
			sigF[l_sig++] = pa;
		}
	}

	fp->cur_read++;

	if((fp->cur_read >= fp->num_read) && !(fp->num_read < 4096)){

		fp->num_read = slow5_get_next_batch(fp->slow5_mt, fp->slow5_batch, 4096);
		fp->cur_read = 0;

		if(fp->num_read <= 0){
			fp->num_read = 0;
		}
	}


	s->sig = (float*)calloc(l_sig, sizeof(float));
	s->l_sig = l_sig;
	memcpy(s->sig, sigF, l_sig*sizeof(float));
	free(sigF);
	// slow5_rec_free(rec);
}
#endif

void ri_read_signal(ri_signal_file_t* fp, ri_signal_t* s){

	assert(fp->cur_read < fp->num_read);

	#ifndef NHDF5RH
	if(fp->fp) read_signal_fast5(fp, s);
	#endif
	#ifndef NPOD5RH
	if(fp->pp) read_signal_pod5(fp, s);
	#endif
	#ifndef NSLOW5RH
	if(fp->sp) read_signal_slow5(fp, s);
	#endif
}
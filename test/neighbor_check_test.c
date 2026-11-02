/* Test of neighbor_check.c against the Python check's results.
 *
 *   neighbor_check_test check <regions.lcr> <pore model> <margin> <reads dir> <out.tsv> [threads]
 *     reads dir from neighbor_regions.py reads (reads.tsv, events.f32). Per read: nc_keep without and with costs, and
 *     every region's cost, compared with the Python costs (written to 4 decimals) and kept_all. With threads > 1 the
 *     decisions are made again on that many threads sharing one db, and must not change.
 *   neighbor_check_test model <pore model> <z.f32>
 *     the model's z-scores against encoding_lib.load_model's, bit for bit.
 *   neighbor_check_test cuts <regions.lcr> <pore model> <margin> <cuts.tsv>
 *     region cuts against neighbor_check.near_neighbor_region's (neighbor_regions.py cuts).
 *
 * gcc -std=c99 -O3 -march=x86-64-v3 -ffp-contract=off -pthread test/neighbor_check_test.c -lm -o neighbor_check_test
 * (includes src/neighbor_check.c to reach its internal steps) */
#include "../src/neighbor_check.c"

#include <pthread.h>
#include <time.h>

#define TOL 0.0000501

typedef struct {
	char *cell, *index, *read, *py_own, *py_near, *py_near_genome;
	uint32_t rid, n_levels;
	uint64_t offset;
	nc_band_t band;
	int band_hits, py_events, py_band_hits, py_strand, py_kept, fast;
	double py_sig;
} row_t;

static double now(void)
{
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	return t.tv_sec + 1e-9 * t.tv_nsec;
}

static char *slurp(const char *path, size_t *n)
{
	FILE *fp = fopen(path, "rb");
	char *buf = NULL;
	long size;
	if (!fp) return NULL;
	fseek(fp, 0, SEEK_END);
	size = ftell(fp);
	fseek(fp, 0, SEEK_SET);
	buf = (char*)malloc(size + 1);
	if (buf && fread(buf, 1, size, fp) != (size_t)size) {
		free(buf);
		buf = NULL;
	}
	fclose(fp);
	if (buf) buf[size] = 0;
	*n = size;
	return buf;
}

static int split(char *line, char **f, int max)
{
	int n = 0;
	f[n++] = line;
	for (; *line && n < max; ++line)
		if (*line == '\t') {
			*line = 0;
			f[n++] = line + 1;
		}
	return n;
}

static void memory(const char *when)
{
	char line[256];
	FILE *fp = fopen("/proc/self/status", "r");
	printf("memory %s:", when);
	while (fp && fgets(line, sizeof line, fp))
		if (!strncmp(line, "VmHWM:", 6) || !strncmp(line, "RssAnon:", 8) || !strncmp(line, "RssFile:", 8)) {
			line[strcspn(line, "\n")] = 0;
			printf(" %s;", line);
		}
	printf("\n");
	if (fp) fclose(fp);
}

/* cost of genome name in a Python list "g:c;g:c" ("-" when empty); -1 if absent */
static double py_cost(const char *list, const char *name)
{
	size_t len = strlen(name);
	const char *p = list;
	while (p && *p && *p != '-') {
		if (strncmp(p, name, len) == 0 && p[len] == ':') return atof(p + len + 1);
		p = strchr(p, ';');
		if (p) ++p;
	}
	return -1.0;
}

static int py_count(const char *list)
{
	int n = 1;
	if (!*list || *list == '-') return 0;
	for (; *list; ++list) n += *list == ';';
	return n;
}

static int by_value(const void *a, const void *b)
{
	double x = *(const double*)a, y = *(const double*)b;
	return (x > y) - (x < y);
}

static void times(const char *what, double *t, int n)
{
	double sum = 0;
	int i;
	for (i = 0; i < n; ++i) sum += t[i];
	qsort(t, n, sizeof(double), by_value);
	printf("time per call, %s: median %.2f ms, mean %.2f ms, 90%% %.2f ms, max %.0f ms, total %.1f s\n", what,
	       1e3 * t[n / 2], 1e3 * sum / n, 1e3 * t[(int)(0.9 * n)], 1e3 * t[n - 1], sum);
}

/* every compared region's cost, the signature first (genome -1), with no early stop */
static int all_costs(const nc_db_t *db, uint32_t rid, const float *levels, uint32_t n_levels, const nc_band_t *band,
		     void *buf, double *cost, int32_t *genome, int *own, uint32_t *events)
{
	const nc_sig_t *s = &db->sigs[rid];
	const nc_reg_t *g;
	double *prev = (double*)buf, *cur = prev + db->max_cut + 1;
	float *y = (float*)(cur + db->max_cut + 1);
	uint32_t a, n, m, len, r;
	int64_t b0, b1, t0, t1;
	int k = 0;
	if (rid >= db->n_sigs || s->len == 0 || !band_span(db, s, band, n_levels, &a, &n, &b0, &b1)) return -1;
	*events = n;
	t0 = b0 > db->pad ? b0 - db->pad : 0;
	t1 = b1 + db->pad < s->len ? b1 + db->pad : s->len;
	m = cut_levels(db, s->base + (uint64_t)t0, (uint32_t)(t1 - t0), band->strand != 0, y);
	cost[k] = subseq_dtw(levels + a, n, y, m, HUGE_VAL, prev, cur);
	genome[k] = -1;
	own[k++] = 1;
	for (r = 0; r < s->n_reg; ++r) {
		g = &db->regs[s->rd + r];
		len = region_cut(db, g, b0, b1, &t0);
		if (len == 0) continue;
		m = cut_levels(db, g->base + (uint64_t)t0, len, ((g->flags & NC_MINUS) != 0) != (band->strand != 0), y);
		cost[k] = subseq_dtw(levels + a, n, y, m, HUGE_VAL, prev, cur);
		genome[k] = (int32_t)(g->flags & NC_GENOME);
		own[k++] = (g->flags & NC_OWN) != 0;
	}
	return k;
}

typedef struct {
	const nc_db_t *db;
	row_t *rows;
	const float *ev;
	int n, start, step, changed;
} job_t;

static void *decide(void *arg)
{
	job_t *j = (job_t*)arg;
	void *buf = malloc(nc_buf_size(j->db));
	int i;
	for (i = j->start; i < j->n; i += j->step) {
		row_t *r = &j->rows[i];
		if (r->band.strand >= 0 && nc_keep(j->db, r->rid, j->ev + r->offset, r->n_levels, &r->band, buf, NULL) != r->fast)
			++j->changed;
	}
	free(buf);
	return NULL;
}

static int run_check(int argc, char **argv)
{
	const char *dir = argv[5];
	int threads = argc > 7 ? atoi(argv[7]) : 1, n = 0, cap = 0, i, j, k, checked = 0, no_check = 0;
	int agree = 0, fast_ok = 0, sig_ok = 0, sets_ok = 0, costs_ok = 0, near_ok = 0, span_ok = 0;
	int kk = 0, kd = 0, dk = 0, dd = 0;
	char path[4096], *tsv, *line, *next, *f[32];
	size_t n_tsv, n_ev;
	double t0, *t_fast, *t_exact, *cost;
	int32_t *genome;
	int *own;
	row_t *rows = NULL;
	float *ev;
	void *buf;
	nc_db_t *db;
	FILE *fo;

	t0 = now();
	db = nc_load(argv[2], argv[3], atoi(argv[4]));
	if (!db) return 1;
	printf("loaded %s (k %d, margin %d) in %.2f s; scratch per thread %.0f kB\n", argv[2], db->k, db->pad, now() - t0,
	       nc_buf_size(db) / 1e3);
	memory("after load");
	snprintf(path, sizeof path, "%s/reads.tsv", dir);
	tsv = slurp(path, &n_tsv);
	snprintf(path, sizeof path, "%s/events.f32", dir);
	ev = (float*)slurp(path, &n_ev);
	if (!tsv || !ev) {
		fprintf(stderr, "cannot read %s\n", dir);
		return 1;
	}
	for (line = strchr(tsv, '\n') + 1; line && *line; line = next) {
		row_t *r;
		next = strchr(line, '\n');
		if (next) *next++ = 0;
		if (split(line, f, 32) < 22) continue;
		if (n == cap) {
			cap = cap ? 2 * cap : 1024;
			rows = (row_t*)realloc(rows, cap * sizeof(row_t));
		}
		r = &rows[n++];
		r->cell = f[0]; r->index = f[1]; r->read = f[2];
		r->rid = (uint32_t)atol(f[3]);
		r->n_levels = (uint32_t)atol(f[5]);
		r->offset = (uint64_t)atoll(f[6]);
		r->band.q0 = (uint32_t)atol(f[7]); r->band.q1 = (uint32_t)atol(f[8]);
		r->band.r0 = (uint32_t)atol(f[9]); r->band.r1 = (uint32_t)atol(f[10]);
		r->band.strand = atoi(f[11]);
		r->band_hits = atoi(f[12]);
		r->py_events = atoi(f[13]); r->py_band_hits = atoi(f[14]); r->py_strand = atoi(f[15]);
		r->py_sig = atof(f[16]);
		r->py_own = f[17]; r->py_near = f[18]; r->py_near_genome = f[20];
		r->py_kept = atoi(f[21]);
	}
	printf("reads %d\n", n);
	buf = malloc(nc_buf_size(db));
	t_fast = (double*)malloc(n * sizeof(double));
	t_exact = (double*)malloc(n * sizeof(double));
	cost = (double*)malloc(65536 * sizeof(double));
	genome = (int32_t*)malloc(65536 * sizeof(int32_t));
	own = (int*)malloc(65536 * sizeof(int));
	fo = fopen(argv[6], "w");
	fprintf(fo, "cell\tindex\tread\trid\tevents\tpy_events\tc_own\tpy_own\tc_near\tc_near_genome\tpy_near_genome\t"
		    "c_kept\tpy_kept\tc_kept_fast\tregions\tpy_regions\tcosts_off\tmax_diff\tms_fast\tms_exact\n");
	for (i = 0; i < n; ++i) {
		row_t *r = &rows[i];
		nc_costs_t c;
		uint32_t events = 0;
		int kept, m, off = 0, sets = 1;
		double t, diff, max_diff, py_own = r->py_sig, p;
		const char *ng;
		if (r->band.strand < 0) continue;	/* no band: no Python row either */
		t = now();
		r->fast = nc_keep(db, r->rid, ev + r->offset, r->n_levels, &r->band, buf, NULL);
		t_fast[checked] = now() - t;
		t = now();
		kept = nc_keep(db, r->rid, ev + r->offset, r->n_levels, &r->band, buf, &c);
		t_exact[checked++] = now() - t;
		m = all_costs(db, r->rid, ev + r->offset, r->n_levels, &r->band, buf, cost, genome, own, &events);
		if (kept < 0 || m <= 0) {
			++no_check;
			fprintf(fo, "%s\t%s\t%s\t%u\tno check\n", r->cell, r->index, r->read, r->rid);
			continue;
		}
		max_diff = fabs(cost[0] - r->py_sig);
		for (k = 1; k < m; ++k) {
			p = py_cost(own[k] ? r->py_own : r->py_near, nc_genome_name(db, genome[k]));
			if (p < 0) {
				sets = 0;
				continue;
			}
			diff = fabs(cost[k] - p);
			if (diff > max_diff) max_diff = diff;
			off += diff > TOL;
			if (own[k] && p < py_own) py_own = p;
		}
		if (m != 1 + py_count(r->py_own) + py_count(r->py_near)) sets = 0;
		ng = nc_genome_name(db, c.near_genome);
		agree += kept == r->py_kept;
		fast_ok += r->fast == kept;
		sig_ok += fabs(cost[0] - r->py_sig) <= TOL;
		sets_ok += sets;
		costs_ok += off == 0 && max_diff <= TOL;
		near_ok += strcmp(ng, r->py_near_genome) == 0;
		span_ok += r->band_hits == r->py_band_hits && r->band.strand == r->py_strand && (int)events == r->py_events;
		kk += kept && r->py_kept;
		kd += kept && !r->py_kept;
		dk += !kept && r->py_kept;
		dd += !kept && !r->py_kept;
		fprintf(fo, "%s\t%s\t%s\t%u\t%u\t%d\t%.6f\t%.4f\t%.6f\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%.6f\t%.3f\t%.3f\n",
			r->cell, r->index, r->read, r->rid, events, r->py_events, c.own, py_own, c.near, ng,
			r->py_near_genome, kept, r->py_kept, r->fast, m, 1 + py_count(r->py_own) + py_count(r->py_near), off,
			max_diff, 1e3 * t_fast[checked - 1], 1e3 * t_exact[checked - 1]);
	}
	fclose(fo);
	printf("checked %d, no check possible %d\n", checked, no_check);
	printf("band and read span as in the Python results: %d\n", span_ok);
	printf("keep decision as Python: %d of %d (both kept %d, C kept/Python dropped %d, C dropped/Python kept %d, "
	       "both dropped %d)\n", agree, checked - no_check, kk, kd, dk, dd);
	printf("decision-only call as the call with costs: %d\n", fast_ok);
	printf("signature cost within 4 decimals: %d; same regions compared: %d; every cost within 4 decimals: %d; "
	       "same best near neighbor: %d\n", sig_ok, sets_ok, costs_ok, near_ok);
	times("decision only (out NULL)", t_fast, checked);
	times("with costs", t_exact, checked);
	memory("after the check");
	if (threads > 1) {
		pthread_t *tid = (pthread_t*)malloc(threads * sizeof(pthread_t));
		job_t *jobs = (job_t*)malloc(threads * sizeof(job_t));
		int changed = 0;
		t0 = now();
		for (j = 0; j < threads; ++j) {
			jobs[j].db = db; jobs[j].rows = rows; jobs[j].ev = ev; jobs[j].n = n;
			jobs[j].start = j; jobs[j].step = threads; jobs[j].changed = 0;
			pthread_create(&tid[j], NULL, decide, &jobs[j]);
		}
		for (j = 0; j < threads; ++j) {
			pthread_join(tid[j], NULL);
			changed += jobs[j].changed;
		}
		printf("%d threads on one db: decisions changed %d; wall %.1f s\n", threads, changed, now() - t0);
		free(tid);
		free(jobs);
	}
	nc_free(db);
	free(rows); free(tsv); free(ev); free(buf); free(t_fast); free(t_exact); free(cost); free(genome); free(own);
	return 0;
}

static int run_model(char **argv)
{
	size_t n_py, i;
	int k, diff = 0;
	float *z = load_model(argv[2], &k), *py = (float*)slurp(argv[3], &n_py);
	if (!z || !py || n_py != ((size_t)4 << 2 * k)) {
		fprintf(stderr, "cannot compare %s with %s\n", argv[2], argv[3]);
		return 1;
	}
	for (i = 0; i < n_py / 4; ++i) diff += memcmp(&z[i], &py[i], 4) != 0;
	printf("%s: k %d, %zu levels, %d differ from load_model's\n", argv[2], k, n_py / 4, diff);
	free(z);
	free(py);
	return diff != 0;
}

/* cuts.tsv: rid, genome, b0, b1, strand, the cut's bases ("-": None). The region's cut is decoded and compared. */
static int run_cuts(char **argv)
{
	static const char *acgtn = "ACGTN";
	nc_db_t *db = nc_load(argv[2], argv[3], atoi(argv[4]));
	char *tsv, *line, *next, *f[8], *s;
	size_t n_tsv;
	uint64_t p, j, lo, hi, from;
	uint32_t rid, len, r, i, b;
	int64_t b0, b1, t0 = 0;
	int strand, n = 0, same = 0, rc, found, with_n = 0;
	tsv = slurp(argv[5], &n_tsv);
	if (!db || !tsv) return 1;
	s = (char*)malloc(db->max_cut + 1);
	for (line = tsv; line && *line; line = next) {
		next = strchr(line, '\n');
		if (next) *next++ = 0;
		if (split(line, f, 8) < 6) continue;
		rid = (uint32_t)atol(f[0]);
		b0 = atoll(f[2]);
		b1 = atoll(f[3]);
		strand = atoi(f[4]);
		for (i = 0; f[5][i]; ++i)	/* letters other than ACGT read as N here */
			if (!strchr("ACGT-", f[5][i])) f[5][i] = 'N';
		for (r = 0, found = 0; r < db->sigs[rid].n_reg && !found; ++r) {
			const nc_reg_t *g = &db->regs[db->sigs[rid].rd + r];
			if (strcmp(nc_genome_name(db, (int32_t)(g->flags & NC_GENOME)), f[1]) != 0) continue;
			found = 1;
			len = region_cut(db, g, b0, b1, &t0);
			rc = ((g->flags & NC_MINUS) != 0) != (strand != 0);
			from = g->base + (uint64_t)t0;
			for (lo = 0, hi = db->n_runs; lo < hi;) {	/* runs that overlap the cut, as cut_levels finds them */
				j = (lo + hi) / 2;
				if ((uint64_t)db->runs[j].start + db->runs[j].len <= from) lo = j + 1;
				else hi = j;
			}
			for (hi = lo; hi < db->n_runs && db->runs[hi].start < from + len; ++hi) {}
			with_n += len > 0 && hi > lo;
			for (i = 0; i < len; ++i) {
				p = rc ? from + len - 1 - i : from + i;
				b = (db->bases[p >> 2] >> ((p & 3) << 1)) & 3;
				for (j = lo; j < hi; ++j)
					if (p >= db->runs[j].start && p - db->runs[j].start < db->runs[j].len) b = 4;
				s[i] = acgtn[b < 4 && rc ? 3 - b : b];
			}
			s[len] = 0;
			/* the Python check uses a cut only if it has more than k + 10 bases */
			if (len == 0) found += f[5][0] == '-' || strlen(f[5]) <= (size_t)db->k + 10;
			else found += strcmp(s, f[5]) == 0;
			if (found == 1 && n - same < 5)
				printf("differs: rid %u %s b0 %lld b1 %lld strand %d: C %u bases, Python %zu\n", rid, f[1],
				       (long long)b0, (long long)b1, strand, len, strlen(f[5]));
		}
		n += found > 0;
		same += found == 2;
	}
	printf("cuts %d (with bases other than ACGT %d), same as near_neighbor_region %d\n", n, with_n, same);
	nc_free(db);
	free(tsv);
	free(s);
	return same != n;
}

int main(int argc, char **argv)
{
	if (argc >= 7 && strcmp(argv[1], "check") == 0) return run_check(argc, argv);
	if (argc == 4 && strcmp(argv[1], "model") == 0) return run_model(argv);
	if (argc == 6 && strcmp(argv[1], "cuts") == 0) return run_cuts(argv);
	fprintf(stderr, "usage: neighborcheck_test check <regions.lcr> <pore model> <margin> <reads dir> <out.tsv> [threads]\n"
			"       neighborcheck_test model <pore model> <z.f32>\n"
			"       neighborcheck_test cuts <regions.lcr> <pore model> <margin> <cuts.tsv>\n");
	return 1;
}

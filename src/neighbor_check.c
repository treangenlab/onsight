/* neighbor_check.c: the near-neighbor check; see neighbor_check.h. A port of neighbor_check.py's check() (subseq_dtw,
 * near_neighbor_region and the kept_all rule) on the regions written by panel-src/near_neighbors/neighbor_regions.py. Costs are summed in double as in the Python (numba) code, so the decisions match it. */

#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif
#include "neighbor_check.h"

#include <fcntl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

/* region file records (neighbor_regions.py). LCR1 holds base positions in 32 bits; LCR2 (panels past 2^32 stored bases)
 * appends the high 32 bits to each record (sig and region: one more uint32; N run: start_hi after len), so the first
 * fields read the same in both. */
typedef struct { uint32_t base, len, rd, n_reg; } nc_sig_t;
typedef struct { uint32_t base; int32_t off; uint16_t len, qs, qe, flags; } nc_reg_t;

#define NC_GENOME 0x3fffu
#define NC_MINUS  0x4000u	/* hit on the genome's '-' strand */
#define NC_OWN    0x8000u	/* genome of the signature's own species */

struct nc_db_s {
	int k, pad;
	uint32_t n_sigs, n_genomes, max_cut;
	uint64_t n_runs;
	float *z;			/* expected level of each k-mer */
	uint32_t *name_at;
	int wide;			/* LCR2: 64-bit base positions */
	size_t sig_size, reg_size, run_size;
	const uint8_t *sigs, *regs;
	const uint8_t *runs;		/* bases other than A, C, G, T, sorted */
	const char *names;
	const uint8_t *bases;		/* 2 bits per base */
	void *map;
	size_t map_size;
};

static inline const nc_sig_t *sig_at(const nc_db_t *db, uint64_t i)
{
	return (const nc_sig_t*)(db->sigs + i * db->sig_size);
}

static inline const nc_reg_t *reg_at(const nc_db_t *db, uint64_t i)
{
	return (const nc_reg_t*)(db->regs + i * db->reg_size);
}

/* a sig's or region's base position: the high 32 bits follow the LCR1 fields in LCR2 */
static inline uint64_t base_of(const nc_db_t *db, const void *rec, uint32_t lo)
{
	return db->wide ? (uint64_t)((const uint32_t*)rec)[4] << 32 | lo : lo;
}

static inline uint64_t run_start(const nc_db_t *db, uint64_t j)
{
	const uint32_t *r = (const uint32_t*)(db->runs + j * db->run_size);
	return db->wide ? (uint64_t)r[2] << 32 | r[0] : r[0];
}

static inline uint32_t run_len(const nc_db_t *db, uint64_t j)
{
	return ((const uint32_t*)(db->runs + j * db->run_size))[1];
}

/* numpy's pairwise sum, so the z-scores equal encoding_lib.load_model's bit for bit */
static double pairwise_sum(const double *a, size_t n)
{
	size_t i, n2;
	double r[8], s;
	if (n < 8) {
		for (i = 0, s = -0.0; i < n; ++i) s += a[i];
		return s;
	}
	if (n <= 128) {
		for (i = 0; i < 8; ++i) r[i] = a[i];
		for (i = 8; i < n - n % 8; i += 8) {
			r[0] += a[i];     r[1] += a[i + 1]; r[2] += a[i + 2]; r[3] += a[i + 3];
			r[4] += a[i + 4]; r[5] += a[i + 5]; r[6] += a[i + 6]; r[7] += a[i + 7];
		}
		s = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
		for (; i < n; ++i) s += a[i];
		return s;
	}
	n2 = n / 2;
	n2 -= n2 % 8;
	return pairwise_sum(a, n2) + pairwise_sum(a + n2, n - n2);
}

/* Pore model read as encoding_lib.load_model reads it: "kmer level" lines in ACGT order, an optional "kmer" header;
 * z = (level - mean) / SD over all k-mers, levels parsed to float, sums in double. */
static float *load_model(const char *path, int *k)
{
	static const char *acgt = "ACGT";
	FILE *fp = fopen(path, "r");
	char line[1024], *e;
	const char *p, *q;
	double *lev = NULL, *sq, mean, sd;
	float *z = NULL;
	size_t n = 0, cap = 0, i, code;
	int ok = 1;
	if (!fp) {
		fprintf(stderr, "[nc_load] cannot open %s\n", path);
		return NULL;
	}
	while (ok && fgets(line, sizeof line, fp)) {
		for (p = line; *p == ' ' || *p == '\t'; ++p) {}
		if (*p == '\n' || *p == '\r' || !*p || (strncmp(p, "kmer", 4) == 0 && (p[4] == ' ' || p[4] == '\t')))
			continue;
		for (code = 0; *p && strchr(" \t\r\n", *p) == NULL; ++p) {
			q = strchr(acgt, *p);
			ok = q != NULL;
			if (!ok) break;
			code = code << 2 | (size_t)(q - acgt);
		}
		ok = ok && code == n;	/* the k-mer must be the line's index */
		if (ok && n == cap) {
			cap = cap ? cap << 1 : 4096;
			lev = (double*)realloc(lev, cap * sizeof(double));
			ok = lev != NULL;
		}
		if (ok) {
			lev[n] = (double)(float)strtod(p, &e);
			ok = e != p;
			++n;
		}
	}
	fclose(fp);
	for (*k = 0; *k < 15 && ((size_t)1 << 2 * *k) < n; ++*k) {}
	if (!ok || n == 0 || ((size_t)1 << 2 * *k) != n) {
		fprintf(stderr, "[nc_load] %s: not a pore model of 4^k levels in ACGT order\n", path);
		free(lev);
		return NULL;
	}
	sq = (double*)malloc(n * sizeof(double));
	z = (float*)malloc(n * sizeof(float));
	if (sq && z) {
		for (i = 0; i < n; ++i) sq[i] = lev[i] * lev[i];
		mean = pairwise_sum(lev, n) / (double)n;
		sd = sqrt(pairwise_sum(sq, n) / (double)n - mean * mean);
		for (i = 0; i < n; ++i) z[i] = (float)((lev[i] - mean) / sd);
	} else {
		free(z);
		z = NULL;
	}
	free(sq);
	free(lev);
	return z;
}

nc_db_t *nc_load(const char *path, const char *pore_model, int pad)
{
	const uint8_t *h;
	uint32_t max_pad, i, j;
	uint64_t n_regs, n_bases, name_bytes, size;
	struct stat st;
	nc_db_t *db;
	int fd = open(path, O_RDONLY);
	void *map;
	if (fd < 0 || fstat(fd, &st) != 0 || st.st_size < 48) {
		fprintf(stderr, "[nc_load] cannot read %s\n", path);
		if (fd >= 0) close(fd);
		return NULL;
	}
	map = mmap(NULL, (size_t)st.st_size, PROT_READ, MAP_SHARED, fd, 0);
	close(fd);
	if (map == MAP_FAILED) {
		fprintf(stderr, "[nc_load] cannot map %s\n", path);
		return NULL;
	}
	db = (nc_db_t*)calloc(1, sizeof(nc_db_t));
	h = (const uint8_t*)map;
	memcpy(&max_pad, h + 4, 4);
	memcpy(&db->n_sigs, h + 8, 4);
	memcpy(&db->n_genomes, h + 12, 4);
	memcpy(&n_regs, h + 16, 8);
	memcpy(&n_bases, h + 24, 8);
	memcpy(&db->n_runs, h + 32, 8);
	memcpy(&name_bytes, h + 40, 8);
	db->wide = memcmp(h, "LCR2", 4) == 0;
	db->sig_size = db->wide ? 20 : 16;
	db->reg_size = db->wide ? 20 : 16;
	db->run_size = db->wide ? 12 : 8;
	size = 48 + db->sig_size * (uint64_t)db->n_sigs + db->reg_size * n_regs + db->run_size * db->n_runs + name_bytes
	       + (n_bases + 3) / 4;
	db->map = map;
	db->map_size = (size_t)st.st_size;
	if ((memcmp(h, "LCR1", 4) != 0 && !db->wide) || size != (uint64_t)st.st_size) {
		fprintf(stderr, "[nc_load] %s is not a near-neighbor region file\n", path);
		nc_free(db);
		return NULL;
	}
	if (pad < 0 || (uint32_t)pad > max_pad) {
		fprintf(stderr, "[nc_load] margin %d: %s has windows for a margin of at most %u\n", pad, path, max_pad);
		nc_free(db);
		return NULL;
	}
	db->pad = pad;
	db->sigs = h + 48;
	db->regs = db->sigs + db->sig_size * db->n_sigs;
	db->runs = db->regs + db->reg_size * n_regs;
	db->names = (const char*)(db->runs + db->run_size * db->n_runs);
	db->bases = (const uint8_t*)db->names + name_bytes;
	db->name_at = (uint32_t*)malloc(db->n_genomes * sizeof(uint32_t) + 1);
	for (i = 0, j = 0; i < db->n_genomes && j < name_bytes; ++i) {
		db->name_at[i] = j;
		j += (uint32_t)strlen(db->names + j) + 1;
	}
	for (i = 0; i < db->n_sigs; ++i)
		if (sig_at(db, i)->len > db->max_cut) db->max_cut = sig_at(db, i)->len;
	db->max_cut += 2 * pad;		/* a cut is the band's span plus the margins, so no longer than this */
	db->z = load_model(pore_model, &db->k);
	if (!db->z) {
		nc_free(db);
		return NULL;
	}
	return db;
}

void nc_free(nc_db_t *db)
{
	if (!db) return;
	if (db->map) munmap(db->map, db->map_size);
	free(db->z);
	free(db->name_at);
	free(db);
}

size_t nc_buf_size(const nc_db_t *db)
{
	return 2 * ((size_t)db->max_cut + 1) * sizeof(double) + (size_t)db->max_cut * sizeof(float);
}

const char *nc_genome_name(const nc_db_t *db, int32_t genome)
{
	return genome >= 0 && (uint32_t)genome < db->n_genomes ? db->names + db->name_at[genome] : "-";
}

uint32_t nc_sig_len(const nc_db_t *db, uint32_t rid)
{
	return rid < db->n_sigs ? sig_at(db, rid)->len : 0;
}

/* Expected levels of bases [from, from + n), of their reverse complement if rc. As encoding_lib.kmer_codes, a base
 * other than A, C, G, T leaves the k-mer as it was. */
static uint32_t cut_levels(const nc_db_t *db, uint64_t from, uint32_t n, int rc, float *y)
{
	uint64_t lo = 0, hi = db->n_runs, mid, p, j;
	uint32_t i, m = 0, kmer = 0, mask = (1u << 2 * db->k) - 1, b;
	while (lo < hi) {	/* runs that overlap the cut: [lo, hi) */
		mid = (lo + hi) / 2;
		if (run_start(db, mid) + run_len(db, mid) <= from) lo = mid + 1;
		else hi = mid;
	}
	for (hi = lo; hi < db->n_runs && run_start(db, hi) < from + n; ++hi) {}
	for (i = 0; i < n; ++i) {
		p = rc ? from + n - 1 - i : from + i;
		b = (db->bases[p >> 2] >> ((p & 3) << 1)) & 3;
		for (j = lo; j < hi; ++j)
			if (p >= run_start(db, j) && p - run_start(db, j) < run_len(db, j)) b = 4;
		if (b < 4) kmer = (kmer << 2 | (rc ? 3 - b : b)) & mask;
		if (i + 1 >= (uint32_t)db->k) y[m++] = db->z[kmer];
	}
	return m;
}

/* neighbor_check.subseq_dtw: cost of the best alignment of all of x to any stretch of y (squared difference), over n.
 * A row's minimum is never below the previous row's, so once it reaches stop * n the final cost cannot fall below
 * stop and the alignment is given up (HUGE_VAL). */
static double subseq_dtw(const float *x, uint32_t n, const float *y, uint32_t m, double stop, double *prev,
			 double *cur)
{
	uint32_t i, j;
	double d, best, left, xi, low = HUGE_VAL, *t;
	for (j = 0; j <= m; ++j) prev[j] = 0.0;
	for (i = 0; i < n; ++i) {
		xi = x[i];
		left = cur[0] = 1e30;
		low = HUGE_VAL;
		for (j = 1; j <= m; ++j) {
			d = xi - (double)y[j - 1];
			d *= d;
			best = prev[j - 1];
			if (prev[j] < best) best = prev[j];
			if (left < best) best = left;
			left = d + best;
			cur[j] = left;
			if (left < low) low = left;
		}
		if (low / n >= stop) return HUGE_VAL;
		t = prev;
		prev = cur;
		cur = t;
	}
	return low / n;
}

/* The read span ev[max(q0 - pad, 0):q1 + pad] and the signature bases [b0, b1) the band covers, on the forward
 * strand. 0 if either is empty. */
static int band_span(const nc_db_t *db, const nc_sig_t *s, const nc_band_t *band, uint32_t n_levels, uint32_t *a,
		     uint32_t *n, int64_t *b0, int64_t *b1)
{
	int64_t lo = (int64_t)band->q0 - db->pad, hi = (int64_t)band->q1 + db->pad;
	if (lo < 0) lo = 0;
	if (hi > n_levels) hi = n_levels;
	if (band->strand == 0) {
		*b0 = band->r0;
		*b1 = (int64_t)band->r1 + db->k;
	} else {
		*b0 = (int64_t)s->len - ((int64_t)band->r1 + db->k);
		*b1 = (int64_t)s->len - band->r0;
	}
	*a = (uint32_t)lo;
	*n = hi > lo ? (uint32_t)(hi - lo) : 0;
	return *n > 0 && *b0 >= 0 && *b1 - *b0 >= db->k && *b1 <= s->len;
}

/* neighbor_check.near_neighbor_region on the stored window: the genome bases homologous to signature bases [b0, b1), plus the
 * margin. Returns the cut's length; 0 if the hit misses [b0, b1) or the cut has k + 10 bases or fewer, which the
 * Python check skips. */
static uint32_t region_cut(const nc_db_t *db, const nc_reg_t *g, int64_t b0, int64_t b1, int64_t *t0)
{
	int64_t t1;
	if (b1 <= g->qs || b0 >= g->qe) return 0;
	if (g->flags & NC_MINUS) {
		*t0 = (int64_t)g->off - b1 - db->pad;
		t1 = (int64_t)g->off - b0 + db->pad;
	} else {
		*t0 = (int64_t)g->off + b0 - db->pad;
		t1 = (int64_t)g->off + b1 + db->pad;
	}
	if (*t0 < 0) *t0 = 0;
	if (t1 > g->len) t1 = g->len;
	return t1 - *t0 > db->k + 10 ? (uint32_t)(t1 - *t0) : 0;
}

int nc_keep(const nc_db_t *db, uint32_t rid, const float *levels, uint32_t n_levels, const nc_band_t *band,
	    void *buf, nc_costs_t *out)
{
	const nc_sig_t *s;
	const nc_reg_t *g;
	int64_t b0, b1, t0, t1;
	uint32_t a, n, m, len, r;
	double *prev = (double*)buf, *cur = prev + db->max_cut + 1, own, near = HUGE_VAL, cost, stop;
	float *y = (float*)(cur + db->max_cut + 1);
	int32_t near_genome = -1;
	int own_side;

	if (rid >= db->n_sigs || sig_at(db, rid)->len == 0) return -1;
	s = sig_at(db, rid);
	if (!band_span(db, s, band, n_levels, &a, &n, &b0, &b1)) return -1;
	t0 = b0 > db->pad ? b0 - db->pad : 0;
	t1 = b1 + db->pad < s->len ? b1 + db->pad : s->len;
	m = cut_levels(db, base_of(db, s, s->base) + (uint64_t)t0, (uint32_t)(t1 - t0), band->strand != 0, y);
	own = subseq_dtw(levels + a, n, y, m, HUGE_VAL, prev, cur);

	/* own strains first: their best cost is the cost the near neighbors must stay above */
	for (own_side = 1; own_side >= 0; --own_side) {
		for (r = 0; r < s->n_reg; ++r) {
			g = reg_at(db, (uint64_t)s->rd + r);
			if (((g->flags & NC_OWN) != 0) != own_side) continue;
			len = region_cut(db, g, b0, b1, &t0);
			if (len == 0) continue;
			m = cut_levels(db, base_of(db, g, g->base) + (uint64_t)t0, len,
				       ((g->flags & NC_MINUS) != 0) != (band->strand != 0), y);
			/* a region is followed only while it can still change the result */
			if (own_side) stop = own;
			else stop = out ? near : nextafter(own, HUGE_VAL);
			cost = subseq_dtw(levels + a, n, y, m, stop, prev, cur);
			if (own_side) {
				if (cost < own) own = cost;
			} else if (cost < near) {
				near = cost;
				near_genome = (int32_t)(g->flags & NC_GENOME);
				if (!out && near <= own) break;
			}
		}
	}
	if (out) {
		out->own = own;
		out->near = near;
		out->near_genome = near_genome;
	}
	return own < near;	/* HUGE_VAL near: no near neighbor overlaps the band, kept */
}

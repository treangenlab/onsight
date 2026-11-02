# onsight - Root Makefile
#
# 'make' builds bin/onsight through src/Makefile.
#
# Options (passed through to src/Makefile):
#   ENABLE_POD5=OFF    Disable POD5 (.pod5)
#   ENABLE_HDF5=OFF    Disable HDF5 (.fast5)
#   ENABLE_SLOW5=OFF   Disable SLOW5 (.slow5 / .blow5)
#   DEBUG=1            Debug build (-O2 -fsanitize=address -g)
#   PROFILE=1          Profiling build (-g -fno-omit-frame-pointer -DPROFILERH=1)
#
# Hardware profiling (doc/HW_PROFILING.md):
#   make prof          bin/onsight-prof (per-stage timers and work counts, -DONS_PROF)
#   make attr          bin/onsight-attr (as prof, hot helpers out of line for perf)
#   asan=1 / tsan=1    Sanitizers

.PHONY: all onsight subset prof attr clean help install uninstall

define build_done
	@echo "Build complete: bin/onsight"
	@echo "Run it as ./bin/onsight, or 'make install' to put it on PATH."
endef

# Install location. Inside a conda env, default to that env so the command is
# on PATH whenever it is active; otherwise the usual /usr/local. Override with
# `make install PREFIX=/somewhere`.
PREFIX ?= $(if $(CONDA_PREFIX),$(CONDA_PREFIX),/usr/local)
BINDIR ?= $(PREFIX)/bin

all: onsight

onsight:
	@mkdir -p bin
	+$(MAKE) -C src
	@cp src/onsight bin/
	$(build_done)

# Built after the normal binary, so the external libraries are already there.
prof: onsight
	+$(MAKE) -C src subset ONS_PROF=1
	@cp src/onsight-prof bin/
	@echo "Build complete: bin/onsight-prof"

attr: onsight
	+$(MAKE) -C src subset ONS_PROF=1 ONS_ATTRIB=1
	@cp src/onsight-attr bin/
	@echo "Build complete: bin/onsight-attr"

subset:
	@mkdir -p bin
	+$(MAKE) -C src subset
	@cp src/onsight bin/
	$(build_done)

# --- Install the command onto PATH, and the pore models where --pore finds them ---
SHAREDIR ?= $(PREFIX)/share/onsight
PORE_MODELS = extern/local_kmer_models/r94_means_only.txt extern/local_kmer_models/uncalled_r1041_model_only_means.txt

install: all
	@mkdir -p $(BINDIR) $(SHAREDIR)
	@install -m 755 bin/onsight $(BINDIR)/onsight
	@install -m 644 $(PORE_MODELS) $(SHAREDIR)/
	@echo "Installed $(BINDIR)/onsight, pore models in $(SHAREDIR)"

uninstall:
	@rm -f $(BINDIR)/onsight
	@rm -f $(addprefix $(SHAREDIR)/,$(notdir $(PORE_MODELS)))
	@echo "Removed $(BINDIR)/onsight, the pore models from $(SHAREDIR)"

clean:
	rm -rf bin
	+$(MAKE) -C src clean

help:
	@echo "onsight build"
	@echo ""
	@echo "  make                  Build bin/onsight (FAST5, POD5, SLOW5/BLOW5)"
	@echo "  make subset           Build without recompiling external dependencies"
	@echo "  make ENABLE_POD5=OFF  Build without POD5 (likewise ENABLE_HDF5, ENABLE_SLOW5)"
	@echo "  make DEBUG=1          Debug build with AddressSanitizer"
	@echo "  make PROFILE=1        Build with profiling support"
	@echo "  make prof             Build bin/onsight-prof (hardware profiling, doc/HW_PROFILING.md)"
	@echo "  make attr             Build bin/onsight-attr (prof with hot helpers out of line, for perf)"
	@echo "  make install          Put onsight on PATH"
	@echo "  make clean            Remove all build artifacts"
	@echo "  asan=1 / tsan=1       AddressSanitizer / ThreadSanitizer"

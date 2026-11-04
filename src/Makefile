# src/Makefile - Standalone build for onsight
#
# All signal formats are enabled by default.
#
# Options:
#   ENABLE_POD5=OFF    Disable POD5 (.pod5)
#   ENABLE_HDF5=OFF    Disable HDF5 (.fast5)
#   ENABLE_SLOW5=OFF   Disable SLOW5 (.slow5 / .blow5)
#   DEBUG=1            Debug build (-O2 -fsanitize=address -g)
#   FP_CONTRACT=off    ARM: no fused multiply-add, as x86 builds (see Architecture-specific flags)
#   ONS_PROF=1         Hardware-characterization build (see onsprof.h); objects go
#                      to obj-prof/ and the binary is onsight-prof
#   ONS_ATTRIB=1       With ONS_PROF: hot helpers out of line for perf attribution;
#                      obj-attr/, binary onsight-attr
#   asan=1 / tsan=1    Sanitizers
#
# Note: a pure-C build requires disabling BOTH POD5 and HDF5 (see NEED_CXX below).
#
# When CONDA_PREFIX is set, HDF5 and zstd are taken from the conda environment
# instead of extern/, so those submodules never need initializing.

# =====================================================================
# Feature toggles - public interface is ENABLE_*=ON|OFF, default ON
# =====================================================================
ENABLE_POD5  ?= ON
ENABLE_HDF5  ?= ON
ENABLE_SLOW5 ?= ON

# The internal NOxxx flags tested below; `override`, so ENABLE_* decides even if NOxxx is on the command line.
override NOPOD5  := $(if $(filter ON on On,$(ENABLE_POD5)),0,1)
override NOHDF5  := $(if $(filter ON on On,$(ENABLE_HDF5)),0,1)
override NOSLOW5 := $(if $(filter ON on On,$(ENABLE_SLOW5)),0,1)

ifeq ($(NOPOD5)$(NOHDF5)$(NOSLOW5),111)
    $(error At least one signal format must be enabled \
      (ENABLE_POD5 / ENABLE_HDF5 / ENABLE_SLOW5) - the binary could read nothing)
endif

# =====================================================================
# System detection
# =====================================================================
SYSTEM_NAME      ?= $(shell uname -s)
SYSTEM_PROCESSOR ?= $(shell uname -m)

# =====================================================================
# Directories
# =====================================================================
SRCDIR  := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
WORKDIR := $(abspath $(SRCDIR)/../extern)

# =====================================================================
# Compiler selection
# When POD5 or HDF5 is enabled, ALL .c files must be compiled as C++:
#   - signal_input.c contains C++ code (std::string, new, etc.) in #ifndef NHDF5RH
#   - POD5 static libs (libarrow.a) require C++ linker
# When neither is enabled, pure C compilation suffices.
# =====================================================================
NEED_CXX := 0
ifneq ($(NOPOD5),1)
    NEED_CXX := 1
endif
ifneq ($(NOHDF5),1)
    NEED_CXX := 1
endif

ifeq ($(NEED_CXX),1)
    COMPILER    = $(CXX)
    LANG_STD    = -std=c++11
    LANG_MODE   = -x c++
    COMPAT_WARN =
else
    COMPILER    = $(CC)
    LANG_STD    = -std=c99
    LANG_MODE   =
    COMPAT_WARN = -Wc++-compat
endif

# =====================================================================
# Base compiler flags
# =====================================================================
BASE_DEFS = -DHAVE_KALLOC

ifdef DEBUG
    OPT_FLAGS = -O2 -fsanitize=address -g
else
    OPT_FLAGS = -O3
endif


ifdef ONS_PROF
    OPT_FLAGS += -g -fno-omit-frame-pointer
    BASE_DEFS += -DONS_PROF=1
endif
ifdef ONS_ATTRIB
    BASE_DEFS += -DONS_ATTRIB=1
endif

# =====================================================================
# Architecture-specific flags
# =====================================================================
ARCH_FLAGS =
ARCH_DEFS  =

# --- x86 / x86_64 ---
ifneq ($(filter x86_64 i686 i386,$(SYSTEM_PROCESSOR)),)
    # x86-64-v3 (AVX2, BMI, FMA: Intel Haswell and AMD Zen or newer); -march=native could emit AVX-512, which
    # other x86 machines cannot run. No fused multiply-add (FP_CONTRACT below), so every x86 host gives the same floats.
    ARCH_FLAGS += -march=x86-64-v3 -mtune=generic
    FP_CONTRACT ?= off
    # Detect immintrin.h availability
    HAVE_IMMINTRIN := $(shell echo '\#include <immintrin.h>' | \
        $(if $(filter 1,$(NEED_CXX)),$(CXX),$(CC)) -x c -fsyntax-only - 2>/dev/null && echo 1 || echo 0)
    ifeq ($(HAVE_IMMINTRIN),1)
        ARCH_DEFS += -DHAVE_IMMINTRIN_H
    endif
endif

# --- AArch64 (arm64) ---
ifneq ($(filter aarch64 arm64 ARM64,$(SYSTEM_PROCESSOR)),)
    ARCH_FLAGS += -fsigned-char
    ARCH_DEFS  += -D_FILE_OFFSET_BITS=64
    # -march=native may not work on all Apple Silicon toolchains
    MARCH_OK := $(shell $(if $(filter 1,$(NEED_CXX)),$(CXX),$(CC)) -march=native -x c -c /dev/null -o /dev/null 2>/dev/null && echo 1 || echo 0)
    ifeq ($(MARCH_OK),1)
        ARCH_FLAGS += -march=native
    endif
endif

# --- ARM 32-bit ---
ifneq ($(filter arm armv7l armv6l,$(SYSTEM_PROCESSOR)),)
    ARCH_FLAGS += -fsigned-char -mfpu=neon
    ARCH_DEFS  += -D_FILE_OFFSET_BITS=64
    MARCH_OK := $(shell $(CC) -march=native -x c -c /dev/null -o /dev/null 2>/dev/null && echo 1 || echo 0)
    ifeq ($(MARCH_OK),1)
        ARCH_FLAGS += -march=native
    endif
endif

# Fused multiply-add. x86 builds turn it off (off above); ARM builds keep the compiler's default, which fuses, as in
# the paper's Jetson runs. FP_CONTRACT=off gives ARM the x86 setting (an x86 build with fast gives the Jetson's calls).
ifneq ($(FP_CONTRACT),)
    ARCH_FLAGS += -ffp-contract=$(FP_CONTRACT)
endif

# =====================================================================
# Sanitizer flags (standalone, outside DEBUG)
# =====================================================================
ifneq ($(asan),)
    OPT_FLAGS     += -fsanitize=address
    EXTRA_LDFLAGS += -fsanitize=address
endif

ifneq ($(tsan),)
    OPT_FLAGS     += -fsanitize=thread
    EXTRA_LDFLAGS += -fsanitize=thread
endif

# =====================================================================
# Combined compiler flags
# =====================================================================
ALL_FLAGS = $(LANG_STD) -Wall $(COMPAT_WARN) $(OPT_FLAGS) $(ARCH_FLAGS) \
            $(BASE_DEFS) $(ARCH_DEFS) -pthread

# =====================================================================
# POD5 configuration (v0.3.36)
# =====================================================================
POD5_VERSION = 0.3.36
POD5_REPO    = https://github.com/nanoporetech/pod5-file-format

ifneq ($(NOPOD5),1)
    # Platform-specific URL and library directory
    ifeq ($(SYSTEM_NAME),Linux)
        ifneq ($(filter aarch64 arm64 ARM64,$(SYSTEM_PROCESSOR)),)
            POD5_ARCHIVE = lib_pod5-$(POD5_VERSION)-linux-arm64.tar.gz
        else
            POD5_ARCHIVE = lib_pod5-$(POD5_VERSION)-linux-x64.tar.gz
        endif
        POD5_LIB_SUBDIR      = lib64
        POD5_EXTRA_LIBS_FILE = libjemalloc_pic.a
    else ifeq ($(SYSTEM_NAME),Darwin)
        ifneq ($(filter arm64 ARM64 aarch64,$(SYSTEM_PROCESSOR)),)
            POD5_ARCHIVE = lib_pod5-$(POD5_VERSION)-osx-14.0-arm64.tar.gz
        else
            $(error POD5 $(POD5_VERSION) has no pre-built binaries for macOS x86_64. Use ENABLE_POD5=OFF.)
        endif
        POD5_LIB_SUBDIR      = lib
        POD5_EXTRA_LIBS_FILE =
    endif

    POD5_URL = $(POD5_REPO)/releases/download/$(POD5_VERSION)/$(POD5_ARCHIVE)
    POD5_DIR = $(WORKDIR)/pod5-$(POD5_VERSION)-$(SYSTEM_NAME)

    # Allow user to override include/lib paths (skip download if set)
    ifdef POD5_INCLUDE_DIR
        POD5_DOWNLOAD = 0
    endif
    POD5_INCLUDE_DIR ?= $(POD5_DIR)/include
    POD5_LIB_DIR     ?= $(POD5_DIR)/$(POD5_LIB_SUBDIR)

    POD5_LIBRARIES = $(POD5_LIB_DIR)/libpod5_format.a \
                     $(POD5_LIB_DIR)/libarrow.a
    ifneq ($(POD5_EXTRA_LIBS_FILE),)
        POD5_LIBRARIES += $(POD5_LIB_DIR)/$(POD5_EXTRA_LIBS_FILE)
    endif

    INCLUDES += -I$(POD5_INCLUDE_DIR)
    LIBS     += $(POD5_LIBRARIES)
else
    ALL_FLAGS += -DNPOD5RH=1
endif

# =====================================================================
# HDF5 configuration
# =====================================================================
HDF5_DIR       ?= $(WORKDIR)/hdf5
HDF5_BUILD_DIR ?= $(HDF5_DIR)/build

# Prefer conda's HDF5 when the active environment provides it. Setting
# HDF5_INCLUDE_DIR implies HDF5_COMPILE=0 below, so the ~118 MB autotools build
# of extern/hdf5 is skipped entirely and that submodule never needs initializing.
ifneq ($(CONDA_PREFIX),)
    ifneq ($(wildcard $(CONDA_PREFIX)/include/H5pubconf.h),)
        HDF5_INCLUDE_DIR ?= $(CONDA_PREFIX)/include
        HDF5_LIB_DIR     ?= $(CONDA_PREFIX)/lib
    endif
endif

ifdef HDF5_INCLUDE_DIR
    HDF5_COMPILE = 0
endif
HDF5_INCLUDE_DIR ?= $(HDF5_BUILD_DIR)/include
HDF5_LIB_DIR     ?= $(HDF5_BUILD_DIR)/lib
HDF5_LIB         ?= hdf5

ifneq ($(NOHDF5),1)
    INCLUDES += -I$(HDF5_INCLUDE_DIR)
    # The static archive when one exists (the extern/ build), otherwise the shared library (conda ships only .so)
    ifneq ($(wildcard $(HDF5_LIB_DIR)/lib$(HDF5_LIB).a),)
        LIBS += $(HDF5_LIB_DIR)/lib$(HDF5_LIB).a
    else
        LIBS += -L$(HDF5_LIB_DIR) -l$(HDF5_LIB) -Wl,-rpath,$(HDF5_LIB_DIR)
    endif
else
    ALL_FLAGS += -DNHDF5RH=1
endif

# =====================================================================
# SLOW5 configuration
# =====================================================================
SLOW5_DIR ?= $(WORKDIR)/slow5lib

ifdef SLOW5_INCLUDE_DIR
    SLOW5_COMPILE = 0
endif
SLOW5_INCLUDE_DIR ?= $(SLOW5_DIR)/include
SLOW5_LIB_DIR     ?= $(SLOW5_DIR)/lib

ifneq ($(NOSLOW5),1)
    INCLUDES += -I$(SLOW5_INCLUDE_DIR)
    LIBS     += $(SLOW5_LIB_DIR)/libslow5.a
else
    ALL_FLAGS += -DNSLOW5RH=1
endif

# =====================================================================
# zstd configuration (needed by both POD5 and SLOW5)
# =====================================================================
# POD5's libpod5_format.a has 8 undefined ZSTD_* symbols (static archives do not
# bundle their dependencies), so zstd is required whenever POD5 or SLOW5 is on.
#
# Prefer conda's zstd. This must move together with HDF5 above: linking a shared
# libhdf5 while zstd stayed static from extern/ would put two zstd versions in one process.
ZSTD_NEEDED = 0
ifneq ($(NOPOD5),1)
    ZSTD_NEEDED = 1
endif
ifneq ($(NOSLOW5),1)
    ZSTD_NEEDED = 1
endif

ifeq ($(ZSTD_NEEDED),1)
    ifneq ($(CONDA_PREFIX),)
        ifneq ($(wildcard $(CONDA_PREFIX)/include/zstd.h),)
            ZSTD_FROM_CONDA = 1
        endif
    endif
    ifeq ($(ZSTD_FROM_CONDA),1)
        INCLUDES += -I$(CONDA_PREFIX)/include
        LIBS     += -L$(CONDA_PREFIX)/lib -lzstd -Wl,-rpath,$(CONDA_PREFIX)/lib
    else
        LIBS     += $(WORKDIR)/zstd/lib/libzstd.a
    endif
endif

# =====================================================================
# System libraries (must come after static libs for correct link order)
# =====================================================================
LIBS += -lm -lz $(EXTRA_LDFLAGS)
ifeq ($(SYSTEM_NAME),Linux)
    LIBS += -ldl
endif

# =====================================================================
# C++ runtime libraries (must come LAST for correct link ordering)
# POD5 v0.3.36 static archives (libarrow.a, libpod5_format.a) reference
# C++ stdlib symbols that must be resolved after those archives.
# =====================================================================
ifeq ($(NEED_CXX),1)
    ifeq ($(SYSTEM_NAME),Linux)
        # POD5 v0.3.36 uses std::filesystem; -lstdc++fs provides these
        # symbols on GCC < 9 and is a harmless stub on GCC 9+.
        LIBS += -lstdc++fs -lstdc++
    endif
endif

# =====================================================================
# Source files and objects
# =====================================================================
# In pipeline order after the vendored code (third_party/: klib, minimap2's FASTA reader, threads and memory pool)
SRCS = third_party/kthread.c third_party/kalloc.c third_party/bseq.c \
       options.c utils.c pore_model.c index.c \
       signal_input.c event_detection.c seeds.c seed_lookup.c seed_hits.c \
       diagonal_bands.c slope_range.c calls.c neighbor_check.c streaming.c \
       sig_summary.c pipeline.c main.c

INCLUDES += -Ithird_party

ifdef ONS_PROF
SRCS += onsprof.c
ifdef ONS_ATTRIB
PROF_OBJDIR = obj-attr
PROG = onsight-attr
else
PROF_OBJDIR = obj-prof
PROG = onsight-prof
endif
OBJS = $(addprefix $(PROF_OBJDIR)/,$(SRCS:.c=.o))
else
OBJS = $(SRCS:.c=.o)
PROG = onsight
endif

# =====================================================================
# Determine which external dependencies need building
# =====================================================================
EXT_TARGETS =

NEED_ZSTD = 0
ifneq ($(NOPOD5),1)
    NEED_ZSTD = 1
endif
ifneq ($(NOSLOW5),1)
    NEED_ZSTD = 1
endif
# When conda supplies zstd we link it directly and never build extern/zstd,
# so that submodule does not need initializing either.
ifeq ($(ZSTD_FROM_CONDA),1)
    NEED_ZSTD = 0
endif
ifeq ($(NEED_ZSTD),1)
    EXT_TARGETS += zstd
endif

ifneq ($(NOHDF5),1)
    EXT_TARGETS += hdf5
endif

ifneq ($(NOSLOW5),1)
    EXT_TARGETS += slow5
endif

ifneq ($(NOPOD5),1)
    EXT_TARGETS += pod5
endif

# Validation targets (run after deps are built)
CHECK_TARGETS =
ifneq ($(NOPOD5),1)
    CHECK_TARGETS += check_pod5
endif
ifneq ($(NOHDF5),1)
    CHECK_TARGETS += check_hdf5
endif
ifneq ($(NOSLOW5),1)
    CHECK_TARGETS += check_slow5
endif

# =====================================================================
# Top-level targets
# =====================================================================
.PHONY: all subset clean help zstd hdf5 slow5 pod5 \
        check_hdf5 check_slow5 check_pod5

all: $(EXT_TARGETS) $(CHECK_TARGETS) $(PROG)

subset: $(CHECK_TARGETS) $(PROG)

# With make -j: the checks after the downloads and builds, the objects after the checks
$(CHECK_TARGETS): | $(EXT_TARGETS)
$(OBJS): | $(CHECK_TARGETS)


# =====================================================================
# Link
# =====================================================================
$(PROG): $(OBJS)
	$(COMPILER) $(ALL_FLAGS) $(OBJS) -o $@ $(LIBS)

# =====================================================================
# Compile: all .c files use the same compiler (CC or CXX)
# =====================================================================
%.o: %.c
	$(COMPILER) -c $(LANG_MODE) $(ALL_FLAGS) $(INCLUDES) -MMD -MP $< -o $@

$(PROF_OBJDIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(COMPILER) -c $(LANG_MODE) $(ALL_FLAGS) $(INCLUDES) -MMD -MP $< -o $@

# =====================================================================
# External dependency: zstd
# =====================================================================
zstd:
ifeq ($(NEED_ZSTD),1)
	@echo "=== Building zstd ==="
	+$(MAKE) -C $(WORKDIR)/zstd lib-release
endif

# =====================================================================
# External dependency: HDF5
# =====================================================================
hdf5:
ifneq ($(NOHDF5),1)
ifneq ($(HDF5_COMPILE),0)
	@echo "=== Building HDF5 ==="
	cd $(HDF5_DIR) && \
		mkdir -p build && \
		./configure --enable-threadsafe --disable-hl --prefix="$(HDF5_BUILD_DIR)" && \
		$(MAKE) && \
		$(MAKE) install
endif
endif

# =====================================================================
# External dependency: slow5lib (depends on zstd)
# =====================================================================
slow5: zstd
ifneq ($(NOSLOW5),1)
ifneq ($(SLOW5_COMPILE),0)
	@echo "=== Building slow5lib ==="
	+$(MAKE) -C $(SLOW5_DIR) \
		slow5_mt=1 \
		zstd_local=$(WORKDIR)/zstd/lib \
		lib/libslow5.a
endif
endif

# =====================================================================
# External dependency: POD5 (download pre-built, depends on zstd)
# =====================================================================
pod5: zstd
ifneq ($(NOPOD5),1)
ifneq ($(POD5_DOWNLOAD),0)
	@echo "=== Downloading POD5 v$(POD5_VERSION) ==="
	@mkdir -p $(POD5_DIR)
	@if [ ! -f "$(POD5_INCLUDE_DIR)/pod5_format/c_api.h" ]; then \
		echo "Fetching $(POD5_URL)"; \
		curl -fSL $(POD5_URL) | tar -xz -C $(POD5_DIR); \
	else \
		echo "POD5 already downloaded, skipping."; \
	fi
endif
endif

# =====================================================================
# Validation checks
# =====================================================================
check_pod5:
ifneq ($(NOPOD5),1)
	@[ -f "$(POD5_INCLUDE_DIR)/pod5_format/c_api.h" ] || \
		{ echo "ERROR: POD5 headers not found at $(POD5_INCLUDE_DIR)" >&2; exit 1; }
	@[ -f "$(POD5_LIB_DIR)/libpod5_format.a" ] || \
		{ echo "ERROR: POD5 library not found at $(POD5_LIB_DIR)" >&2; exit 1; }
endif

check_hdf5:
ifneq ($(NOHDF5),1)
	@[ -f "$(HDF5_INCLUDE_DIR)/H5pubconf.h" ] || \
		{ echo "ERROR: HDF5 headers not found at $(HDF5_INCLUDE_DIR)" >&2; exit 1; }
	@[ -f "$(HDF5_LIB_DIR)/lib$(HDF5_LIB).a" ] || \
		[ -f "$(HDF5_LIB_DIR)/lib$(HDF5_LIB).so" ] || \
		{ echo "ERROR: HDF5 library not found at $(HDF5_LIB_DIR)" >&2; exit 1; }
endif

check_slow5:
ifneq ($(NOSLOW5),1)
	@[ -f "$(SLOW5_INCLUDE_DIR)/slow5/slow5.h" ] || \
		{ echo "ERROR: SLOW5 headers not found at $(SLOW5_INCLUDE_DIR)" >&2; exit 1; }
	@[ -f "$(SLOW5_LIB_DIR)/libslow5.a" ] || \
		{ echo "ERROR: SLOW5 library not found at $(SLOW5_LIB_DIR)" >&2; exit 1; }
endif

# =====================================================================
# Clean
# =====================================================================
clean:
	rm -f *.o *.d third_party/*.o third_party/*.d version.h onsight onsight-prof onsight-attr *~
	rm -rf obj-prof obj-attr

# =====================================================================
# Help
# =====================================================================
help:
	@echo ""
	@echo "src/Makefile options:"
	@echo "  DEBUG=1       Debug build (-O2 -fsanitize=address -g)"
	@echo "  FP_CONTRACT=off   ARM: no fused multiply-add, as x86 builds"
	@echo "  ENABLE_POD5=OFF   Disable POD5 (likewise ENABLE_HDF5, ENABLE_SLOW5)"
	@echo "  asan=1        Enable AddressSanitizer"
	@echo "  tsan=1        Enable ThreadSanitizer"
	@echo ""
	@echo "  make subset   Skip recompiling external dependencies"

# =====================================================================
# --version: the number in ../VERSION
# =====================================================================
version.h: ../VERSION
	@echo '#define ONSIGHT_VERSION "v$(shell cat ../VERSION)"' > $@
main.o: version.h
$(PROF_OBJDIR)/main.o: version.h

# =====================================================================
# Header dependencies, written by the compiler (-MMD) next to each object
# =====================================================================
-include $(OBJS:.o=.d)

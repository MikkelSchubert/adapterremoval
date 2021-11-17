###############################################################################
# Makefile options: Edit / comment / uncomment to change build behavior

# Installation destinations
PREFIX := /usr/local

# Default compilation flags
CXXFLAGS := ${CXXFLAGS} -std=c++11 -O3

## Optional features; comment out or set to value other than 'yes' to disable

# Enable SSE, SSE2, and AVX2 instructions for a significant performance gain (YMMV)
VECTORIZE := yes

# Use Intelligent Storage Acceleration Library (ISA-L) for gzip decompression
LIBISAL := yes

# Use libdeflate for block compression
LIBDEFLATE := yes

# Hide individual commands during build; only shows summaries instead.
QUIET_BUILD := yes

# Use of colored output during build
COLOR_BUILD := yes

# Debug build; adds warnings and debugging symbols
DEBUG_BUILD := no

# Include coverage instrumentation in build
COVERAGE := no


###############################################################################
# Makefile internals. Normally you do not need to touch these.
INSTALLEXE = install -m 0755
INSTALLDAT = install -m 0644
INSTALLDOC = install -m 0644
MKDIR      = install -d  # act as mkdir -p

# Libraries required by AdapterRemoval
LIBRARIES := -pthread -lz

# Build directory; modified depending on build options
BDIR     := build/main

ifeq ($(strip ${QUIET_BUILD}),yes)
QUIET := @
endif

ifeq ($(strip ${COLOR_BUILD}),yes)
COLOR_YELLOW := "\033[0;33m"
COLOR_GREEN := "\033[0;32m"
COLOR_CYAN := "\033[0;36m"
COLOR_END := "\033[0m"
endif


ifeq ($(strip ${VECTORIZE}),yes)
$(info Building AdapterRemoval with SSE/SSE2/AVX2 extensions: yes)
CXXFLAGS := $(CXXFLAGS) -msse -msse2 -mavx2
else
$(info Building AdapterRemoval with SSE/SSE2/AVX2 extensions: no)
endif


ifeq ($(strip ${LIBISAL}),yes)
$(info Building AdapterRemoval with isa-l: yes)
CXXFLAGS := $(CXXFLAGS) -DUSE_LIBISAL
LIBRARIES := $(LIBRARIES) -lisal
else
$(info Building AdapterRemoval with isa-l: no)
endif


ifeq ($(strip ${LIBDEFLATE}),yes)
$(info Building AdapterRemoval with libdeflate: yes)
CXXFLAGS := $(CXXFLAGS) -DUSE_LIBDEFLATE
LIBRARIES := $(LIBRARIES) -ldeflate
else
$(info Building AdapterRemoval with libdeflate: no)
endif

ifeq ($(strip ${COVERAGE}), yes)
$(info Building AdapterRemoval with coverage instrumentation: yes)
CXXFLAGS := ${CXXFLAGS} --coverage
DEBUG_BUILD := yes
else
$(info Building AdapterRemoval with coverage instrumentation: no)
endif

ifeq ($(strip ${DEBUG_BUILD}), yes)
$(info Building AdapterRemoval with debug information: yes)
CXXFLAGS := ${CXXFLAGS} -g -pedantic -Wall -Wextra -Wcast-align -Wcast-qual \
	-Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self \
	-Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wsign-promo \
	-Wstrict-overflow=2 -Wswitch-default -Wundef -Weffc++ -Wdeprecated
else
$(info Building AdapterRemoval with debug information: no)
endif

PROG     := AdapterRemoval
LIBOBJS  := $(BDIR)/adapterset.o \
            $(BDIR)/alignment.o \
            $(BDIR)/alignment_tables.o \
            $(BDIR)/argparse.o \
            $(BDIR)/barcode_table.o \
            $(BDIR)/debug.o \
            $(BDIR)/demultiplexing.o \
            $(BDIR)/fastq.o \
            $(BDIR)/fastq_enc.o \
            $(BDIR)/fastq_io.o \
            $(BDIR)/json.o \
            $(BDIR)/linereader.o \
            $(BDIR)/linereader_joined.o \
            $(BDIR)/main_adapter_id.o \
            $(BDIR)/main_adapter_rm.o \
            $(BDIR)/main_demultiplex.o \
            $(BDIR)/main_fastq_ro.o \
            $(BDIR)/managed_writer.o \
            $(BDIR)/reports.o \
            $(BDIR)/scheduler.o \
            $(BDIR)/statistics.o \
            $(BDIR)/strutils.o \
            $(BDIR)/threads.o \
            $(BDIR)/timer.o \
            $(BDIR)/trimming.o \
            $(BDIR)/userconfig.o \
            $(BDIR)/utilities.o
OBJS     := ${LIBOBJS} $(BDIR)/main.o
DFILES   := $(OBJS:.o=.deps)


.PHONY: all install clean test clean_tests static regression docs

all: build/$(PROG)

everything: all static test regression docs

# Clean
clean: clean_tests clean_docs
	@echo $(COLOR_GREEN)"Cleaning ..."$(COLOR_END)
	$(QUIET) rm -f build/$(PROG)
	$(QUIET) rm -rvf build/regression
	$(QUIET) rm -rvf $(BDIR)

# Install
install: build/$(PROG)
	@echo $(COLOR_GREEN)"Installing AdapterRemoval .."$(COLOR_END)
	@echo $(COLOR_GREEN)"  .. binary into ${PREFIX}/bin/"$(COLOR_END)
	$(QUIET) $(MKDIR) ${PREFIX}/bin/
	$(QUIET) $(INSTALLEXE) build/$(PROG) ${PREFIX}/bin/

	@echo $(COLOR_GREEN)"  .. man-page into ${PREFIX}/share/man/man1/"$(COLOR_END)
	$(QUIET) $(MKDIR) ${PREFIX}/share/man/man1/
	$(QUIET) $(INSTALLDOC) $(PROG).1 ${PREFIX}/share/man/man1/

	@echo $(COLOR_GREEN)"  .. README into ${PREFIX}/share/adapterremoval/"$(COLOR_END)
	$(QUIET) $(MKDIR) ${PREFIX}/share/adapterremoval/
	$(QUIET) $(INSTALLDOC) README.md ${PREFIX}/share/adapterremoval/

	@echo $(COLOR_GREEN)"  .. examples into ${PREFIX}/share/adapterremoval/examples/"$(COLOR_END)
	$(QUIET) $(MKDIR) ${PREFIX}/share/adapterremoval/examples/
	$(QUIET) $(INSTALLDAT) examples/*.* ${PREFIX}/share/adapterremoval/examples/

# Object files
$(BDIR)/%.o: src/%.cpp
	@echo $(COLOR_CYAN)"Building $@ from $<"$(COLOR_END)
	$(QUIET) $(MKDIR) $(BDIR)
	$(QUIET) $(CXX) $(CXXFLAGS) -pthread -c -o $@ $<
	$(QUIET) $(CXX) $(CXXFLAGS) -w -MM -MT $@ -MF $(@:.o=.deps) $<

# Executable
build/$(PROG): $(OBJS)
	@echo $(COLOR_GREEN)"Linking executable $@"$(COLOR_END)
	$(QUIET) $(CXX) $(CXXFLAGS) ${LDFLAGS} $^ ${LIBRARIES} -o $@

# Automatic header depencencies
-include $(DFILES)


###############################################################################
# Unit testing
TEST_DIR := build/tests
TEST_OBJS := $(TEST_DIR)/main_test.o \
             $(TEST_DIR)/debug.o \
             $(TEST_DIR)/alignment.o \
             $(TEST_DIR)/alignment_tables.o \
             $(TEST_DIR)/alignment_test.o \
             $(TEST_DIR)/argparse.o \
             $(TEST_DIR)/argparse_test.o \
             $(TEST_DIR)/barcodes_test.o \
             $(TEST_DIR)/barcode_table.o \
             $(TEST_DIR)/fastq.o \
             $(TEST_DIR)/fastq_test.o \
             $(TEST_DIR)/fastq_enc.o \
             $(TEST_DIR)/json.o \
             $(TEST_DIR)/json_test.o \
             $(TEST_DIR)/strutils.o \
             $(TEST_DIR)/strutils_test.o
TEST_DEPS := $(TEST_OBJS:.o=.deps)

TEST_CXXFLAGS := -Isrc -DAR_TEST_BUILD -g

test: $(TEST_DIR)/main
	@echo $(COLOR_GREEN)"Running unit tests"$(COLOR_END)
	$(QUIET) $(TEST_DIR)/main

clean_tests:
	@echo $(COLOR_GREEN)"Cleaning tests ..."$(COLOR_END)
	$(QUIET) rm -rvf $(TEST_DIR)

$(TEST_DIR)/main: $(TEST_OBJS)
	@echo $(COLOR_GREEN)"Linking executable $@"$(COLOR_END)
	$(QUIET) $(CXX) $(CXXFLAGS) ${LIBRARIES} $^ -o $@

$(TEST_DIR)/%.o: tests/unit/%.cpp
	@echo $(COLOR_CYAN)"Building $@ from $<"$(COLOR_END)
	$(QUIET) $(MKDIR) $(TEST_DIR)
	$(QUIET) $(CXX) $(CXXFLAGS) $(TEST_CXXFLAGS) -c -o $@ $<
	$(QUIET) $(CXX) $(CXXFLAGS) $(TEST_CXXFLAGS) -w -MM -MT $@ -MF $(@:.o=.deps) $<

$(TEST_DIR)/%.o: src/%.cpp
	@echo $(COLOR_CYAN)"Building $@ from $<"$(COLOR_END)
	$(QUIET) $(MKDIR) $(TEST_DIR)
	$(QUIET) $(CXX) $(CXXFLAGS) $(TEST_CXXFLAGS) -c -o $@ $<
	$(QUIET) $(CXX) $(CXXFLAGS) $(TEST_CXXFLAGS) -w -MM -MT $@ -MF $(@:.o=.deps) $<


###############################################################################
# Validation
VALIDATION_BDIR=./build/regression
VALIDATION_SDIR=./tests/regression

regression: build/$(PROG)
	@echo $(COLOR_GREEN)"Running regression tests"$(COLOR_END)
	@$(MKDIR) $(VALIDATION_BDIR)
	@$(VALIDATION_SDIR)/run $(VALIDATION_BDIR) $(VALIDATION_SDIR)


# Automatic header dependencies for tests
-include $(TEST_DEPS)


###############################################################################
# Documentation
SPHINXOPTS  = -n -q
SPHINXBUILD = sphinx-build

docs:
	$(QUIET) @$(SPHINXBUILD) -M html docs build/docs $(SPHINXOPTS)
	$(QUIET) @$(SPHINXBUILD) -M man docs build/docs $(SPHINXOPTS)
	$(QUIET) cp -v "build/docs/man/AdapterRemoval.1" .

clean_docs:
	@echo $(COLOR_GREEN)"Cleaning documentation ..."$(COLOR_END)
	$(QUIET) rm -rvf build/docs

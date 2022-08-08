###############################################################################
# Makefile options: Edit / comment / uncomment to change build behavior

# Installation destinations
PREFIX := /usr/local

# Default compilation flags
CXXFLAGS := ${CXXFLAGS} -std=c++14 -O3

## Optional features; comment out or set to value other than 'yes' to disable

# Generate statically linked binary
STATIC := no

# Use Intelligent Storage Acceleration Library (ISA-L) for gzip decompression
LIBISAL := yes

# Use libdeflate for block compression
LIBDEFLATE := yes

# Show individual commands during build; otherwise shows summaries instead.
VERBOSE := no

# Use of colored output during build
COLOR_BUILD := yes

# Debug build; adds warnings and debugging symbols
DEBUG_BUILD := no

# Include coverage instrumentation in build
COVERAGE := no


###############################################################################
# Makefile internals. Normally you do not need to touch these.

EXEC_MAIN   := adapterremoval3
EXEC_TEST   := unit_tests

BUILD_DIR   := build
OBJS_DIR    := $(BUILD_DIR)/objects
EXECUTABLE  := $(BUILD_DIR)/$(EXEC_MAIN)
TEST_RUNNER := $(BUILD_DIR)/$(EXEC_TEST)

REGRESSION_TESTS := tests/regression
REGRESSION_DIR   := $(BUILD_DIR)/regression

INSTALLEXE = install -m 0755
INSTALLDAT = install -m 0644
INSTALLDOC = install -m 0644
MKDIR      = install -d  # act as mkdir -p

# Libraries required by AdapterRemoval
LIBRARIES := -pthread -lz

ifeq ($(strip ${VERBOSE}),no)
QUIET := @
endif

ifeq ($(strip ${COLOR_BUILD}),yes)
ifneq ($(strip $(MAKE_TERMOUT)),)
COLOR_YELLOW := "\033[0;33m"
COLOR_GREEN := "\033[0;32m"
COLOR_CYAN := "\033[0;36m"
COLOR_END := "\033[0m"
endif
endif

ifeq ($(strip ${STATIC}),yes)
$(info Building static AdapterRemoval binary: yes)
# https://gcc.gnu.org/bugzilla/show_bug.cgi?id=58909
CXXFLAGS := $(CXXFLAGS) -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
OBJS_DIR := $(BUILD_DIR)/static
EXECUTABLE := $(BUILD_DIR)/$(EXEC_MAIN).static
else
$(info Building static AdapterRemoval binary: no)
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
CXXFLAGS := ${CXXFLAGS} -g -DDEBUG \
	-pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy \
	-Wdisabled-optimization -Wformat=2 -Winit-self -Wold-style-cast \
	-Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-overflow=2 \
	-Wswitch-default -Wundef -Weffc++ -Wdeprecated
else
$(info Building AdapterRemoval with debug information: no)
endif

################################################################################

# Build objects shared between unit tests and executable
CORE_OBJS := \
	$(OBJS_DIR)/alignment_tables.o \
	$(OBJS_DIR)/alignment.o \
	$(OBJS_DIR)/argparse.o \
	$(OBJS_DIR)/barcode_table.o \
	$(OBJS_DIR)/debug.o \
	$(OBJS_DIR)/fastq_enc.o \
	$(OBJS_DIR)/fastq.o \
	$(OBJS_DIR)/json.o \
	$(OBJS_DIR)/logging.o \
	$(OBJS_DIR)/strutils.o \
	$(OBJS_DIR)/utilities.o

# Build objects used only by the executable
EXEC_OBJS := \
	$(OBJS_DIR)/adapterset.o \
	$(OBJS_DIR)/demultiplexing.o \
	$(OBJS_DIR)/fastq_io.o \
	$(OBJS_DIR)/linereader_joined.o \
	$(OBJS_DIR)/linereader.o \
	$(OBJS_DIR)/main_adapter_id.o \
	$(OBJS_DIR)/main_adapter_rm.o \
	$(OBJS_DIR)/main_demultiplex.o \
	$(OBJS_DIR)/main_fastq_ro.o \
	$(OBJS_DIR)/main.o \
	$(OBJS_DIR)/managed_writer.o \
	$(OBJS_DIR)/reports_html.o \
	$(OBJS_DIR)/reports_json.o \
	$(OBJS_DIR)/reports_terminal.o \
	$(OBJS_DIR)/reports_template_html.o \
	$(OBJS_DIR)/scheduler.o \
	$(OBJS_DIR)/statistics.o \
	$(OBJS_DIR)/threads.o \
	$(OBJS_DIR)/timer.o \
	$(OBJS_DIR)/progress.o \
	$(OBJS_DIR)/trimming.o \
	$(OBJS_DIR)/userconfig.o

# Build objects used only by the test runner
TEST_OBJS := \
	$(OBJS_DIR)/alignment_test.o \
	$(OBJS_DIR)/argparse_test.o \
	$(OBJS_DIR)/barcodes_test.o \
	$(OBJS_DIR)/debug_test.o \
	$(OBJS_DIR)/counts_test.o \
	$(OBJS_DIR)/fastq_test.o \
	$(OBJS_DIR)/json_test.o \
	$(OBJS_DIR)/logging_test.o \
	$(OBJS_DIR)/main_test.o \
	$(OBJS_DIR)/strutils_test.o \
	$(OBJS_DIR)/utilities_test.o

################################################################################

.PHONY: all clean docs everything examples install regression test

all: $(EXECUTABLE)

clean:
	@echo $(COLOR_GREEN)"Cleaning"$(COLOR_END)
	$(QUIET) rm -rf $(BUILD_DIR)
	$(QUIET) make -C examples clean

everything: all test regression docs examples

examples: $(EXECUTABLE)
	@echo $(COLOR_GREEN)"Running examples"$(COLOR_END)
	$(QUIET) make -C examples EXE=$(PWD)/$(EXECUTABLE)

install: $(EXECUTABLE)
	@echo $(COLOR_GREEN)"Installing AdapterRemoval .."$(COLOR_END)
	@echo $(COLOR_GREEN)"  .. binary into ${PREFIX}/bin/"$(COLOR_END)
	$(QUIET) $(MKDIR) ${PREFIX}/bin/
	$(QUIET) $(INSTALLEXE) $(EXECUTABLE) ${PREFIX}/bin/

	@echo $(COLOR_GREEN)"  .. man-page into ${PREFIX}/share/man/man1/"$(COLOR_END)
	$(QUIET) $(MKDIR) ${PREFIX}/share/man/man1/
	$(QUIET) $(INSTALLDOC) $(EXEC_MAIN).1 ${PREFIX}/share/man/man1/

	@echo $(COLOR_GREEN)"  .. README into ${PREFIX}/share/adapterremoval/"$(COLOR_END)
	$(QUIET) $(MKDIR) ${PREFIX}/share/adapterremoval/
	$(QUIET) $(INSTALLDOC) README.md ${PREFIX}/share/adapterremoval/

	@echo $(COLOR_GREEN)"  .. examples into ${PREFIX}/share/adapterremoval/examples/"$(COLOR_END)
	$(QUIET) $(MKDIR) ${PREFIX}/share/adapterremoval/examples/
	$(QUIET) $(INSTALLDAT) examples/*.* ${PREFIX}/share/adapterremoval/examples/

# HTML report templates
.INTERMEDIATE: src/reports_template.intermediate

src/reports_template_html.hpp src/reports_template_html.cpp: src/reports_template.intermediate ;
src/reports_template.intermediate: src/reports_template.html src/reports_template_html.py
	@echo $(COLOR_CYAN)"Building HTML templates from $<"$(COLOR_END)
	$(QUIET) python3 src/reports_template_html.py $< src/reports_template_html

regression: $(EXECUTABLE)
	@echo $(COLOR_GREEN)"Running regression tests"$(COLOR_END)
	$(QUIET) $(MKDIR) $(REGRESSION_DIR)
	$(QUIET) python3 $(REGRESSION_TESTS).py $(REGRESSION_DIR) $(REGRESSION_TESTS) \
		--executable $(EXECUTABLE)

test: $(TEST_RUNNER)
	@echo $(COLOR_GREEN)"Running unit tests"$(COLOR_END)
	$(QUIET) $(TEST_RUNNER)

$(EXECUTABLE): $(CORE_OBJS) $(EXEC_OBJS)
	@echo $(COLOR_GREEN)"Linking executable $@"$(COLOR_END)
	$(QUIET) $(CXX) $(CXXFLAGS) ${LDFLAGS} $^ ${LIBRARIES} -o $@

$(TEST_RUNNER): $(CORE_OBJS) $(TEST_OBJS)
	@echo $(COLOR_GREEN)"Linking executable $@"$(COLOR_END)
	$(QUIET) $(CXX) $(CXXFLAGS) ${LIBRARIES} $^ -o $@

# Main object files 
$(OBJS_DIR)/%.o: src/%.cpp
	@echo $(COLOR_CYAN)"Building $@ from $<"$(COLOR_END)
	$(QUIET) $(MKDIR) $(OBJS_DIR)
	$(QUIET) $(CXX) $(CXXFLAGS) -pthread -c -o $@ $<
	$(QUIET) $(CXX) $(CXXFLAGS) -pthread -w -MM -MT $@ -MF $(@:.o=.d) $<

# Unit test object files
$(OBJS_DIR)/%.o: tests/unit/%.cpp
	@echo $(COLOR_CYAN)"Building $@ from $<"$(COLOR_END)
	$(QUIET) $(MKDIR) $(OBJS_DIR)
	$(QUIET) $(CXX) $(CXXFLAGS) -Isrc -pthread -c -o $@ $<
	$(QUIET) $(CXX) $(CXXFLAGS) -Isrc -pthread -w -MM -MT $@ -MF $(@:.o=.d) $<


###############################################################################
# Documentation

SPHINXOPTS  = -n -q
SPHINXBUILD = sphinx-build

docs:
	$(QUIET) @$(SPHINXBUILD) -M html docs build/docs $(SPHINXOPTS)
	$(QUIET) @$(SPHINXBUILD) -M man docs build/docs $(SPHINXOPTS)
	$(QUIET) cp -v "build/docs/man/AdapterRemoval.1" .

###############################################################################

# Automatic tracking of include-file dependencies
-include $(CORE_OBJS:.o=.d)
-include $(EXEC_OBJS:.o=.d)
-include $(TEST_OBJS:.o=.d)

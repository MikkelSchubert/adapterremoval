###############################################################################
# Makefile options: Edit here or specify on command-line (e.g. make STATIC=yes)

# Installation destinations
PREFIX := /usr/local

## Optional features; comment out or set to value other than 'yes' to disable

# Generate statically linked binary
STATIC := no

# Show individual commands during build; otherwise shows summaries instead.
VERBOSE := no

# Use of colored output during build (yes/no/auto)
COLOR := auto

# Debug build; adds warnings and debugging symbols
DEBUG := no

# Include coverage instrumentation in build
COVERAGE := no

# Enable address and undefined behavior sanitation
SANITIZE := no


###############################################################################
# Makefile internals. Normally you do not need to touch these.

EXEC_MAIN   := adapterremoval3
EXEC_TEST   := unit_tests

BUILD_DIR   := build
OBJS_DIR    := $(BUILD_DIR)/objects
COV_DIR     := $(BUILD_DIR)/coverage
DOCS_DIR    := $(BUILD_DIR)/docs

MAN_PAGE    := $(DOCS_DIR)/man/$(EXEC_MAIN).1

REGRESSION_TESTS := tests/regression
REGRESSION_DIR   := $(BUILD_DIR)/regression

INSTALLEXE = install -m 0755
INSTALLDAT = install -m 0644
INSTALLDOC = install -m 0644
MKDIR      = install -d  # act as mkdir -p

# Default compilation flags
CXXFLAGS := ${CXXFLAGS} -std=c++17 -O3
LDLIBS := -pthread -lisal -ldeflate ${LDLIBS}
LDFLAGS := ${LDFLAGS}

ifeq ($(strip ${VERBOSE}),no)
QUIET := @
endif

ifeq ($(strip ${COLOR}),auto)
ifeq (${NO_COLOR},)
ifneq ($(strip $(MAKE_TERMOUT)),)
COLOR := yes
endif
endif
endif

ifeq ($(strip ${COLOR}),yes)
COLOR_YELLOW := "\033[0;33m"
COLOR_GREEN := "\033[0;32m"
COLOR_CYAN := "\033[0;36m"
COLOR_END := "\033[0m"
else
COLOR := no
endif

BUILD_NAME := release
BUILD_NAME_PREFIX :=
BUILD_NAME_POSTFIX :=

ifeq ($(strip ${STATIC}),yes)
$(info Building static AdapterRemoval binary: yes)
# https://gcc.gnu.org/bugzilla/show_bug.cgi?id=58909
LDFLAGS := $(LDFLAGS) -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -static
OBJS_DIR := $(BUILD_DIR)/static
BUILD_NAME_PREFIX := static-
else
$(info Building static AdapterRemoval binary: no)
endif

ifeq ($(strip ${COVERAGE}), yes)
$(info Building AdapterRemoval with coverage instrumentation: yes)
CXXFLAGS := ${CXXFLAGS} --coverage
DEBUG := yes
BUILD_NAME := coverage
else
$(info Building AdapterRemoval with coverage instrumentation: no)
endif

ifeq ($(strip ${DEBUG}), yes)
$(info Building AdapterRemoval with debug information: yes)
CXXFLAGS := ${CXXFLAGS} -g -DDEBUG \
	-pedantic -Wall -Wextra -Wcast-align -Wcast-qual -Wctor-dtor-privacy \
	-Wdisabled-optimization -Wformat=2 -Winit-self -Wold-style-cast \
	-Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-overflow=2 \
	-Wswitch-default -Wswitch-enum -Wundef -Weffc++ -Wdeprecated

ifneq ($(strip ${COVERAGE}), yes)
BUILD_NAME := debug
endif
else
$(info Building AdapterRemoval with debug information: no)
endif

ifeq ($(strip ${SANITIZE}), yes)
$(info Building AdapterRemoval with sanitizers: yes)
CXXFLAGS := ${CXXFLAGS} -fno-sanitize-recover=all -fsanitize=undefined,address
ifeq ($(findstring clang,${CXX}),)
LDLIBS := $(LDLIBS) -lasan
endif
BUILD_NAME_POSTFIX := ${BUILD_NAME_POSTFIX}-sanitize
else
$(info Building AdapterRemoval with sanitizers: no)
endif

# Use a different folder for each build configuration
BUILD_NAME  := $(BUILD_NAME_PREFIX)$(BUILD_NAME)$(BUILD_NAME_POSTFIX)
EXECUTABLE  := $(BUILD_DIR)/$(BUILD_NAME)/$(EXEC_MAIN)
TEST_RUNNER := $(BUILD_DIR)/$(BUILD_NAME)/$(EXEC_TEST)

OBJS_DIR := $(BUILD_DIR)/$(BUILD_NAME)

$(info Build configuration: $(BUILD_NAME))

################################################################################

# Build objects shared between unit tests and executable
CORE_OBJS := \
	$(OBJS_DIR)/alignment.o \
	$(OBJS_DIR)/argparse.o \
	$(OBJS_DIR)/barcode_table.o \
	$(OBJS_DIR)/debug.o \
	$(OBJS_DIR)/errors.o \
	$(OBJS_DIR)/fastq_enc.o \
	$(OBJS_DIR)/fastq.o \
	$(OBJS_DIR)/json.o \
	$(OBJS_DIR)/linereader.o \
	$(OBJS_DIR)/logging.o \
	$(OBJS_DIR)/managed_io.o \
	$(OBJS_DIR)/mathutils.o \
	$(OBJS_DIR)/simd_avx2.o \
	$(OBJS_DIR)/simd_avx512bw.o \
	$(OBJS_DIR)/simd_sse2.o \
	$(OBJS_DIR)/simd_std.o \
	$(OBJS_DIR)/simd.o \
	$(OBJS_DIR)/strutils.o \
	$(OBJS_DIR)/utilities.o

# Build objects used only by the executable
EXEC_OBJS := \
	$(OBJS_DIR)/adapter_id.o \
	$(OBJS_DIR)/adapterset.o \
	$(OBJS_DIR)/benchmarking.o \
	$(OBJS_DIR)/demultiplexing.o \
	$(OBJS_DIR)/fastq_io.o \
	$(OBJS_DIR)/linereader_joined.o \
	$(OBJS_DIR)/main_adapter_id.o \
	$(OBJS_DIR)/main_adapter_rm.o \
	$(OBJS_DIR)/main_benchmark.o \
	$(OBJS_DIR)/main_demultiplex.o \
	$(OBJS_DIR)/main.o \
	$(OBJS_DIR)/reports_html.o \
	$(OBJS_DIR)/reports_json.o \
	$(OBJS_DIR)/reports_terminal.o \
	$(OBJS_DIR)/reports_template_html.o \
	$(OBJS_DIR)/scheduler.o \
	$(OBJS_DIR)/statistics.o \
	$(OBJS_DIR)/timer.o \
	$(OBJS_DIR)/progress.o \
	$(OBJS_DIR)/trimming.o \
	$(OBJS_DIR)/userconfig.o

# Build objects used only by the test runner
TEST_OBJS := \
	$(OBJS_DIR)/alignment_test.o \
	$(OBJS_DIR)/argparse_test.o \
	$(OBJS_DIR)/barcodes_test.o \
	$(OBJS_DIR)/buffer_test.o \
	$(OBJS_DIR)/counts_test.o \
	$(OBJS_DIR)/debug_test.o \
	$(OBJS_DIR)/fastq_test.o \
	$(OBJS_DIR)/json_test.o \
	$(OBJS_DIR)/linereader_test.o \
	$(OBJS_DIR)/logging_test.o \
	$(OBJS_DIR)/main_test.o \
	$(OBJS_DIR)/mathutils_test.o \
	$(OBJS_DIR)/strutils_test.o \
	$(OBJS_DIR)/utilities_test.o

################################################################################

.PHONY: all clean coverage docs man everything examples install regression test

all: $(EXECUTABLE)

clean:
	@echo $(COLOR_GREEN)"Cleaning"$(COLOR_END)
	$(QUIET) $(RM) -r $(BUILD_DIR)
	$(QUIET) $(MAKE) -C examples clean

everything: all test regression docs examples

examples: $(EXECUTABLE)
	@echo $(COLOR_GREEN)"Running examples"$(COLOR_END)
	$(QUIET) $(MAKE) -j1 -C examples EXE=$(CURDIR)/$(EXECUTABLE)

install: $(EXECUTABLE) $(MAN_PAGE)
	@echo $(COLOR_GREEN)"Installing AdapterRemoval .."$(COLOR_END)
	@echo $(COLOR_GREEN)"  .. binary into ${DESTDIR}${PREFIX}/bin/"$(COLOR_END)
	$(QUIET) $(MKDIR) ${DESTDIR}${PREFIX}/bin/
	$(QUIET) $(INSTALLEXE) $(EXECUTABLE) ${DESTDIR}${PREFIX}/bin/

	@echo $(COLOR_GREEN)"  .. man-page into ${DESTDIR}${PREFIX}/share/man/man1/"$(COLOR_END)
	$(QUIET) $(MKDIR) ${DESTDIR}${PREFIX}/share/man/man1/
	$(QUIET) $(INSTALLDOC) "$(MAN_PAGE)" ${DESTDIR}${PREFIX}/share/man/man1/

	@echo $(COLOR_GREEN)"  .. README into ${DESTDIR}${PREFIX}/share/adapterremoval3/"$(COLOR_END)
	$(QUIET) $(MKDIR) ${DESTDIR}${PREFIX}/share/adapterremoval3/
	$(QUIET) $(INSTALLDOC) README.md ${DESTDIR}${PREFIX}/share/adapterremoval3/

	@echo $(COLOR_GREEN)"  .. examples into ${DESTDIR}${PREFIX}/share/adapterremoval3/examples/"$(COLOR_END)
	$(QUIET) $(MKDIR) ${DESTDIR}${PREFIX}/share/adapterremoval3/examples/
	$(QUIET) $(INSTALLDAT) examples/Makefile examples/*.* ${DESTDIR}${PREFIX}/share/adapterremoval3/examples/

# HTML report templates
.INTERMEDIATE: src/reports_template.intermediate

src/reports_template_html.hpp src/reports_template_html.cpp: src/reports_template.intermediate ;
src/reports_template.intermediate: src/reports_template.html scripts/html_template_to_cpp.py
	@echo $(COLOR_CYAN)"Building HTML templates from $<"$(COLOR_END)
	$(QUIET) python3 scripts/html_template_to_cpp.py $< src/reports_template_html

regression: $(EXECUTABLE)
	@echo $(COLOR_GREEN)"Running regression tests"$(COLOR_END)
	$(QUIET) $(MKDIR) $(REGRESSION_DIR)
	$(QUIET) python3 scripts/regression_test_runner.py $(REGRESSION_DIR) $(REGRESSION_TESTS) \
		--json-schema "schema.json" \
		--executable $(EXECUTABLE)

test: $(TEST_RUNNER)
	@echo $(COLOR_GREEN)"Running unit tests"$(COLOR_END)
	$(QUIET) $(TEST_RUNNER) --invisibles --use-colour $(COLOR)

coverage:
ifeq ($(strip ${COVERAGE}), yes)
ifneq ($(shell which gcovr), )
	$(MAKE) test
	@echo $(COLOR_GREEN)"Running coverage analysis"$(COLOR_END)
	$(QUIET) mkdir -p "$(COV_DIR)"
	$(QUIET) gcovr \
		--txt $(COV_DIR)/index.txt \
		--xml $(COV_DIR)/coverage.xml \
		--exclude-lines-by-pattern ".*\\sREQUIRE[_A-Z]*\\(.*" \
		--sort-percentage \
		--exclude src/robin_hood.hpp \
		--exclude tests/unit/catch.hpp  \
		--delete
	$(QUIET) cat "$(COV_DIR)/index.txt"
else
	$(error Cannot analyse coverage: gcovr not found)
endif
else
	$(MAKE) coverage COVERAGE=yes
endif

$(EXECUTABLE): $(CORE_OBJS) $(EXEC_OBJS)
	@echo $(COLOR_GREEN)"Linking executable $@"$(COLOR_END)
	$(QUIET) $(CXX) $(CXXFLAGS) ${LDFLAGS} $^ ${LDLIBS} -o $@

$(TEST_RUNNER): $(CORE_OBJS) $(TEST_OBJS)
	@echo $(COLOR_GREEN)"Linking executable $@"$(COLOR_END)
	$(QUIET) $(CXX) $(CXXFLAGS) ${LDFLAGS} $^ ${LDLIBS} -o $@

# Main object files
$(OBJS_DIR)/%.o: src/%.cpp
	@echo $(COLOR_CYAN)"Building $@ from $<"$(COLOR_END)
	$(QUIET) $(MKDIR) $(OBJS_DIR)
	$(QUIET) $(CXX) $(CXXFLAGS) -pthread -c -MMD -MQ $@ -MF $(@:.o=.d) -o $@ $<

# Objects built with support for specific CPU instructions
$(OBJS_DIR)/simd_avx2.o $(OBJS_DIR)/simd_avx512bw.o $(OBJS_DIR)/simd_sse2.o : $(OBJS_DIR)/simd_%.o: src/simd_%.cpp
	@echo $(COLOR_CYAN)"Building $@ from $< (-m$*)"$(COLOR_END)
	$(QUIET) $(MKDIR) $(OBJS_DIR)
	$(QUIET) $(CXX) $(CXXFLAGS) -pthread -c -MMD -MQ $@ -MF $(@:.o=.d) -o $@ $< -m$*

# Unit test object files
$(OBJS_DIR)/%.o: tests/unit/%.cpp
	@echo $(COLOR_CYAN)"Building $@ from $<"$(COLOR_END)
	$(QUIET) $(MKDIR) $(OBJS_DIR)
	$(QUIET) $(CXX) $(CXXFLAGS) -pthread -c -MMD -MQ $@ -MF $(@:.o=.d) -o $@ $< -Isrc


###############################################################################
# Documentation

SPHINXOPTS  = -n -q
SPHINXBUILD = sphinx-build

docs: man
	@echo $(COLOR_CYAN)"Building manual"$(COLOR_END)
	$(QUIET) @$(SPHINXBUILD) -M html docs $(BUILD_DIR)/docs $(SPHINXOPTS)

man $(MAN_PAGE):
	@echo $(COLOR_CYAN)"Building man-page"$(COLOR_END)
	$(QUIET) @$(SPHINXBUILD) -M man docs $(BUILD_DIR)/docs $(SPHINXOPTS)

###############################################################################

# Automatic tracking of include-file dependencies
-include $(CORE_OBJS:.o=.d)
-include $(EXEC_OBJS:.o=.d)
-include $(TEST_OBJS:.o=.d)

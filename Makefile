###############################################################################
# Makefile options: Edit / comment / uncomment to change build behavior
#

# Installation destinations
PREFIX := /usr/local

# Default compilation flags
CXXFLAGS := -O3 -ansi -pedantic -Wall -Wextra
CXXFLAGS := ${CXXFLAGS} -Wcast-align -Wcast-qual -Wctor-dtor-privacy \
	-Wdisabled-optimization -Wformat=2 -Winit-self -Wold-style-cast \
	-Woverloaded-virtual -Wredundant-decls -Wsign-promo -Wstrict-overflow=5 \
	-Wswitch-default -Wundef -Weffc++

## Optional features; comment out or set to value other than 'yes' to disable

# Enable reading writing of gzip compressed files using libz.
ENABLE_GZIP_SUPPORT := yes

# Enable reading writing of bzip2 compressed files using libbz2.
ENABLE_BZIP2_SUPPORT := yes

# Enable multi-threading support using pthreads.
ENABLE_PTHREAD_SUPPORT := yes

# Hide individual commands during build; only shows summaries instead.
ENABLE_QUIET_BUILD := yes

# Use of colored output during build
ENABLE_COLOR_BUILD := yes


###############################################################################
# Makefile internals. Normally you do not need to touch these.

# Libraries required by AdapterRemoval
LIBRARIES :=

# Build directory; modified depending on build options
BDIR     := build/main


ifeq ($(strip ${ENABLE_QUIET_BUILD}),yes)
QUIET := @
endif

ifeq ($(strip ${ENABLE_COLOR_BUILD}),yes)
COLOR_YELLOW := "\033[0;33m"
COLOR_GREEN := "\033[0;32m"
COLOR_CYAN := "\033[0;36m"
COLOR_END := "\033[0m"
endif

ifeq ($(strip ${ENABLE_GZIP_SUPPORT}),yes)
$(info Building AdapterRemoval with gzip support: yes)
CXXFLAGS := ${CXXFLAGS} -DAR_GZIP_SUPPORT
LIBRARIES := ${LIBRARIES} -lz
BDIR := ${BDIR}_gz
else
$(info Building AdapterRemoval with gzip support: no)
endif

ifeq ($(strip ${ENABLE_BZIP2_SUPPORT}),yes)
$(info Building AdapterRemoval with bzip2 support: yes)
CXXFLAGS := ${CXXFLAGS} -DAR_BZIP2_SUPPORT
LIBRARIES := ${LIBRARIES} -lbz2
BDIR := ${BDIR}_bz2
else
$(info Building AdapterRemoval with bzip2 support: no)
endif

ifeq ($(strip ${ENABLE_PTHREAD_SUPPORT}),yes)
$(info Building AdapterRemoval with pthreads support: yes)
CXXFLAGS := ${CXXFLAGS} -DAR_PTHREAD_SUPPORT
LIBRARIES := ${LIBRARIES} -lpthread
BDIR := ${BDIR}_pthreads
else
$(info Building AdapterRemoval with pthreads support: no)
endif


PROG     := AdapterRemoval
LIBNAME  := libadapterremoval
LIBOBJS  := $(BDIR)/adapterset.o \
            $(BDIR)/alignment.o \
            $(BDIR)/argparse.o \
            $(BDIR)/debug.o \
            $(BDIR)/demultiplex.o \
            $(BDIR)/fastq.o \
            $(BDIR)/fastq_enc.o \
            $(BDIR)/fastq_io.o \
            $(BDIR)/linereader.o \
            $(BDIR)/main_adapter_id.o \
            $(BDIR)/main_adapter_rm.o \
            $(BDIR)/scheduler.o \
            $(BDIR)/strutils.o \
            $(BDIR)/threads.o \
            $(BDIR)/timer.o \
            $(BDIR)/userconfig.o
OBJS     := ${LIBOBJS} $(BDIR)/main.o
DFILES   := $(OBJS:.o=.deps)


.PHONY: all install clean test clean_tests static

all: build/$(PROG) build/$(PROG).1

# Clean
clean: clean_tests
	@echo $(COLOR_GREEN)"Cleaning ..."$(COLOR_END)
	$(QUIET) rm -f build/$(PROG) build/$(PROG).1
	$(QUIET) rm -rvf $(BDIR)

# Install
install: build/$(PROG) build/$(PROG).1
	@echo $(COLOR_GREEN)"Installing AdapterRemoval .."$(COLOR_END)
	@echo $(COLOR_GREEN)"  .. binary into ${PREFIX}/bin/"$(COLOR_END)
	$(QUIET) mkdir -p ${PREFIX}/bin/
	$(QUIET) mv -f build/$(PROG) ${PREFIX}/bin/
	$(QUIET) chmod a+x ${PREFIX}/bin/$(PROG)

	@echo $(COLOR_GREEN)"  .. man-page into ${PREFIX}/share/man/man1/"$(COLOR_END)
	$(QUIET) mkdir -p ${PREFIX}/share/man/man1/
	$(QUIET) mv -f build/$(PROG).1 ${PREFIX}/share/man/man1/
	$(QUIET) chmod a+r ${PREFIX}/share/man/man1/$(PROG).1

	@echo $(COLOR_GREEN)"  .. README into ${PREFIX}/share/adapterremoval/"$(COLOR_END)
	$(QUIET) mkdir -p ${PREFIX}/share/adapterremoval/
	$(QUIET) cp -a README.md ${PREFIX}/share/adapterremoval/
	$(QUIET) chmod a+r ${PREFIX}/share/adapterremoval/README.md

	@echo $(COLOR_GREEN)"  .. examples into ${PREFIX}/share/adapterremoval/examples/"$(COLOR_END)
	$(QUIET) mkdir -p ${PREFIX}/share/adapterremoval/examples/
	$(QUIET) cp -a examples/*.{txt,fq} ${PREFIX}/share/adapterremoval/examples/
	$(QUIET) chmod a+r ${PREFIX}/share/adapterremoval/examples/*.{txt,fq}

static: build/$(LIBNAME).a

# Object files
$(BDIR)/%.o: src/%.cc
	@echo $(COLOR_CYAN)"Building '$@' from '$<'"$(COLOR_END)
	$(QUIET) mkdir -p $(BDIR)
	$(QUIET) $(CXX) $(CXXFLAGS) -c -o $@ $<
	$(QUIET) $(CXX) $(CXXFLAGS) -w -MM -MT $@ -MF $(@:.o=.deps) $<

# Executable
build/$(PROG): $(OBJS)
	@echo $(COLOR_GREEN)"Linking executable '$@'"$(COLOR_END)
	$(QUIET) $(CXX) $(CXXFLAGS) $^ ${LIBRARIES} -o $@

# Static library
build/$(LIBNAME).a: $(LIBOBJS)
	@echo $(COLOR_GREEN)"Linking static library '$@'"$(COLOR_END)
	$(AR) rcs build/$(LIBNAME).a $(LIBOBJS)

build/%.1: %.pod
	@echo $(COLOR_GREEN)"Constructing man-page '$@' from '$<'"$(COLOR_END)
	$(QUIET) mkdir -p $(BDIR)
	$(QUIET) pod2man $< > $@

# Automatic header depencencies
-include $(DFILES)


#
# Unit testing
#
TEST_DIR := build/tests
TEST_OBJS := $(TEST_DIR)/alignment.o \
             $(TEST_DIR)/alignment_test.o \
             $(TEST_DIR)/argparse.o \
             $(TEST_DIR)/argparse_test.o \
             $(TEST_DIR)/debug.o \
             $(TEST_DIR)/fastq.o \
             $(TEST_DIR)/fastq_enc.o \
             $(TEST_DIR)/fastq_test.o \
             $(TEST_DIR)/strutils.o \
             $(TEST_DIR)/strutils_test.o
TEST_DEPS := $(TEST_OBJS:.o=.deps)

GTEST_DIR := googletest-release-1.7.0
GTEST_OBJS := $(TEST_DIR)/gtest-all.o $(TEST_DIR)/gtest_main.o
GTEST_LIB :=$(TEST_DIR)/libgtest.a

TEST_CXXFLAGS := -isystem $(GTEST_DIR)/include -I$(GTEST_DIR) -Isrc -DAR_TEST_BUILD
GTEST_CXXFLAGS := $(TEST_CXXFLAGS)

test: $(TEST_DIR)/main
	@echo $(COLOR_GREEN)"Running tests"$(COLOR_END)
	$(QUIET) $< --gtest_print_time=0 --gtest_shuffle

clean_tests:
	@echo $(COLOR_GREEN)"Cleaning tests ..."$(COLOR_END)
	$(QUIET) rm -rvf $(TEST_DIR)

$(TEST_DIR)/main: $(GTEST_LIB) $(TEST_OBJS)
	@echo $(COLOR_GREEN)"Linking executable $@\033[0m"$(COLOR_END)
	$(QUIET) $(CXX) $(CXXFLAGS) -pthread $^ -o $@

$(TEST_DIR)/libgtest.a: $(GTEST_OBJS)
	@echo $(COLOR_GREEN)"Linking GTest library '$@'"$(COLOR_END)
	$(QUIET) ar -rv $@ $^

$(TEST_DIR)/%.o: tests/%.cc
	@echo $(COLOR_CYAN)"Building $@ from $<"$(COLOR_END)
	mkdir -p $(TEST_DIR)
	$(QUIET) $(CXX) $(CXXFLAGS) $(TEST_CXXFLAGS) -c -o $@ $<
	$(QUIET) $(CXX) $(CXXFLAGS) $(TEST_CXXFLAGS) -w -MM -MT $@ -MF $(@:.o=.deps) $<

$(TEST_DIR)/%.o: src/%.cc
	@echo $(COLOR_CYAN)"Building $@ from $<"$(COLOR_END)
	mkdir -p $(TEST_DIR)
	$(QUIET) $(CXX) $(CXXFLAGS) $(TEST_CXXFLAGS) -c -o $@ $<
	$(QUIET) $(CXX) $(CXXFLAGS) $(TEST_CXXFLAGS) -w -MM -MT $@ -MF $(@:.o=.deps) $<

$(TEST_DIR)/gtest%.o: $(GTEST_DIR)/src/gtest%.cc
	@echo $(COLOR_CYAN)"Building '$@' from '$<'"$(COLOR_END)
	$(QUIET) mkdir -p $(TEST_DIR)
	$(QUIET) $(CXX) $(GTEST_CXXFLAGS) -pthread -c $< -o $@

.PRECIOUS: $(GTEST_DIR)/src/gtest%.cc
$(GTEST_DIR)/src/gtest%.cc: googletest-release-1.7.0.zip
	$(QUIET) if ! test -e "$@"; \
	then \
		echo $(COLOR_CYAN)"Unpacking Google Test library"$(COLOR_END); \
		unzip -qo googletest-release-1.7.0.zip; \
	fi

googletest-release-1.7.0.zip:
ifneq ("$(shell which wget)", "")
	@echo $(COLOR_CYAN)"Fetching Google Test library using wget"$(COLOR_END)
	$(QUIET) wget -q https://github.com/google/googletest/archive/release-1.7.0.zip -O googletest-release-1.7.0.zip
else ifneq ("$(shell which curl)", "")
	@echo $(COLOR_CYAN)"Fetching Google Test library using curl"$(COLOR_END)
	$(QUIET) curl -L https://github.com/google/googletest/archive/release-1.7.0.zip -o googletest-release-1.7.0.zip
else
	@echo $(COLOR_YELLOW)"To run tests, first download and unpack GoogleTest 1.7.0 in this folder:"$(COLOR_END)
	@echo $(COLOR_YELLOW)"  $$ wget https://github.com/google/googletest/archive/release-1.7.0.zip -O googletest-release-1.7.0.zip"$(COLOR_END)
	@echo $(COLOR_YELLOW)"  $$ unzip googletest-release-1.7.0.zip"$(COLOR_END)
	@exit 1
endif

# Automatic header dependencies for tests
-include $(TEST_DEPS)

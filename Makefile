###############################################################################
# Makefile options: Edit / comment / uncomment to change build behavior
#

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
OBJS     := $(BDIR)/main.o \
			$(BDIR)/main_adapter_id.o \
			$(BDIR)/main_adapter_rm.o \
			$(BDIR)/argparse.o \
            $(BDIR)/alignment.o \
            $(BDIR)/fastq.o \
            $(BDIR)/fastq_io.o \
            $(BDIR)/fastq_enc.o \
            $(BDIR)/userconfig.o \
            $(BDIR)/timer.o \
            $(BDIR)/linereader.o \
            $(BDIR)/scheduler.o \
            $(BDIR)/threads.o \
            $(BDIR)/strutils.o

DFILES   := $(OBJS:.o=.deps)


.PHONY: all install clean test clean_tests

all: build/$(PROG) build/$(PROG).1

# Clean
clean: clean_tests
	@echo $(COLOR_GREEN)"Cleaning ..."$(COLOR_END)
	$(QUIET) rm -f build/$(PROG) build/$(PROG).1
	$(QUIET) rm -rvf $(BDIR)

# Install
install:
	@echo $(COLOR_GREEN)"Installing ..."$(COLOR_END)
	$(QUIET) mv -f build/$(PROG) /usr/local/bin/
	$(QUIET) mv -f build/$(PROG).1 /usr/share/man/man1/

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
TEST_OBJS := $(TEST_DIR)/fastq_test.o $(BDIR)/fastq.o $(BDIR)/fastq_enc.o \
	$(TEST_DIR)/alignment_test.o $(BDIR)/alignment.o \
	$(TEST_DIR)/argparse_test.o $(BDIR)/argparse.o \
	$(TEST_DIR)/strutils_test.o $(BDIR)/strutils.o
TEST_DEPS := $(TEST_OBJS:.o=.deps)

GTEST_DIR := gtest-1.7.0
GTEST_OBJS := $(TEST_DIR)/gtest-all.o $(TEST_DIR)/gtest_main.o
GTEST_LIB :=$(TEST_DIR)/libgtest.a

TEST_CXXFLAGS := -isystem $(GTEST_DIR)/include -I$(GTEST_DIR) -Isrc
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

$(TEST_DIR)/gtest%.o: $(GTEST_DIR)/src/gtest%.cc
	@echo $(COLOR_CYAN)"Building '$@' from '$<'"$(COLOR_END)
	$(QUIET) mkdir -p $(TEST_DIR)
	$(QUIET) $(CXX) $(GTEST_CXXFLAGS) -pthread -c $< -o $@

$(GTEST_DIR)/src/gtest%.cc:
	@echo $(COLOR_YELLOW)"To run tests, first download and unpack GoogleTest 1.7.0 in this folder."$(COLOR_END)
	@exit 1

# Automatic header dependencies for tests
-include $(TEST_DEPS)

# Put objects into subdirs
AUTOMAKE_OPTIONS = subdir-objects

# Default extension to use when no bla_SOURCES exists.
AM_DEFAULT_SOURCE_EXT = .cpp

# Macros (like boost) are in the m4 folder
ACLOCAL_AMFLAGS = -I m4

# Common flags for preprocessor and linker
AM_CPPFLAGS = \
    -I$(top_srcdir)/include \
    -DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN \
    $(BOOST_CPPFLAGS) \
    $(ARB_CPPFLAGS) \
    $(CODE_COVERAGE_CPPFLAGS) \
    $(TBB_CPPFLAGS) \
    $(TBB_MALLOC_CPPFLAGS)

AM_LDFLAGS = \
    $(BOOST_LDFLAGS) \
    $(ARB_LDFLAGS) \
    $(TBB_LDFLAGS) \
    $(TBB_MALLOC_LDFLAGS)

AM_CXXFLAGS = \
    $(CODE_COVERAGE_CXXFLAGS)

AM_CFLAGS = \
    $(CODE_COVERAGE_CFLAGS)


# Files to copy into distribution packages
EXTRA_BINDIST = \
    LICENSE \
    CHANGELOG \
    README \
    arb_LICENSE.txt \
    arb_PUBLIB_LICENSE.txt

EXTRA_DIST = \
    $(EXTRA_BINDIST) \
    Makefile.am \
    configure.ac \
    autogen.sh \
    tools \
    test_data/ltp_reduced.arb.xz \
    include \
    tests

CLEANFILES = \
    CHANGELOG \
    README

### TARGETS ###

# calling pkglibdir our ARBHOME (.../lib/sina), we install 
# internal libs and binaries to .../bin and .../lib subdirs

bin_PROGRAMS = src/sina
lib_LTLIBRARIES = src/libsina.la

if ENABLE_FAT_TAR 
arblibdir = $(libdir)/arb/lib

bin_PROGRAMS += arb_pt_server arb_message
arblib_DATA = \
    arb_tcp.dat

arb_message_SOURCES=
arb_pt_server_SOURCES=
endif

### TARGET CONFIGS ###

src_libsina_la_SOURCES = \
    src/align.cpp \
    src/align.h \
    src/aligned_base.cpp \
    src/aligned_base.h \
    src/alignment_stats.cpp \
    src/alignment_stats.h \
    src/buffer.h \
    src/cache.h \
    src/cseq.cpp \
    src/cseq.h \
    src/cseq_impl.h \
    src/cseq_comparator.cpp \
    src/cseq_comparator.h \
    src/famfinder.h \
    src/famfinder.cpp \
    src/graph.h \
    src/helpers.h \
    src/idset.h \
    src/mesh.h \
    src/mesh_debug.h \
    src/mseq.cpp \
    src/mseq.h \
    src/progress.h \
    src/pseq.cpp \
    src/pseq.h \
    src/query_arb.cpp \
    src/query_arb.h \
    src/query_pt.cpp \
    src/query_pt.h \
    src/rw_fasta.cpp \
    src/rw_fasta.h \
    src/rw_arb.cpp \
    src/rw_arb.h \
    src/rw_csv.cpp \
    src/rw_csv.h \
    src/log.cpp \
    src/log.h \
    src/scoring_schemes.h \
    src/search.h \
    src/search.cpp \
    src/search_filter.cpp \
    src/search_filter.h \
    src/tempfile.h \
    src/timer.h \
    src/tray.cpp \
    src/tray.h \
    src/kmer.h \
    src/kmer_search.cpp \
    src/kmer_search.h

src_libsina_la_LDFLAGS = \
    $(ARB_LDFLAGS) \
    $(BOOST_LDFLAGS)

src_libsina_la_LIBADD = \
    $(Z_LIB) \
    $(ARB_LIBS) \
    $(BOOST_PROGRAM_OPTIONS_LIB) \
    $(BOOST_THREAD_LIB) \
    $(BOOST_IOSTREAMS_LIB) \
    $(BOOST_SYSTEM_LIB) \
    $(BOOST_FILESYSTEM_LIB) \
    $(TBB_LIB) \
    $(TBB_MALLOC_LIB)

sina_libs = \
    src/libsina.la \
    $(ARB_PROBE_LIBS) \
    $(ARB_HELIX_LIBS) \
    $(ARB_LIBS) \
    $(CODE_COVERAGE_LIBS)

src_sina_LDADD = $(sina_libs)

arb_pt_server$(EXEEXT): $(ARBHOME)/bin/arb_pt_server$(EXEEXT)
	cp $< $@
arb_message$(EXEEXT): $(ARBHOME)/bin/arb_message$(EXEXT)
	cp $< $@
arb_tcp.dat: $(ARBHOME)/lib/arb_tcp.dat
	cp $< $@


@CODE_COVERAGE_RULES@
CODE_COVERAGE_DIRECTORY = $(top_builddir)/src
CODE_COVERAGE_OUTPUT_DIRECTORY = $(top_builddir)/coverage
CODE_COVERAGE_OUTPUT_FILE =  $(top_builddir)/coverage.info
header_paths = $(abspath $(patsubst -I%,%,$(filter -I%,$(AM_CPPFLAGS) $(CPPFLAGS))))
CODE_COVERAGE_IGNORE_PATTERN = \
  $(patsubst %,%/\*,$(header_paths) /usr) \
  \*/c++/\* \
  \*unit_tests\* \
  \*include\*

### DOCUMENTATION ###

DOC_SOURCE = \
    doc/source/conf.py \
    doc/source/index.rst \
    doc/source/install.rst \
    doc/source/commandline.rst \
    doc/source/CHANGELOG.rst \
    doc/source/README.rst \
    README.rst

doc_DATA := $(filter %.txt,$(DOC_SOURCE:doc/source/%.rst=doc/text/%.txt))

EXTRA_DIST += $(DOC_SOURCE)

man1_MANS = doc/man/sina.1


$(doc_DATA): doc/build.stamp
doc/build.stamp: $(DOC_SOURCE)
	$(SPHINX_BUILD) -M text $$(dirname $<) doc
	touch $@
CLEANFILES += doc/build.stamp

$(man1_MANS): $(DOC_SOURCE)
	$(SPHINX_BUILD) -M man $$(dirname $<) doc

%: doc/text/%.txt
	cp $< $@



CLEANFILES += sina-manual.log
LATEXTOMAN=$(top_srcdir)/tools/latex2man

export TEXINPUTS=$(top_srcdir):$(top_srcdir)/tools/:

CLEANFILES += $(doc_DATA) $(check_DATA) $(man1_MANS)

clean-local-docs:
	rm -rf doc/doctrees


### PACKAGE ###


dist-hook:
	find $(distdir) \( -name \*~ -or -name \*.log -or -name \*.arb \) -delete


BINDISTFILES =
bindistdir = $(PACKAGE)-$(VERSION)-$(OSNAME)

bindistdir: all
	rm -rf $(bindistdir)
	mkdir $(bindistdir)
	$(MAKE) DESTDIR=$(abspath $(bindistdir)) install
	mv -f $(bindistdir)/$(prefix)/* $(bindistdir)
	cd $(bindistdir); for f in bin/sina share/doc/sina*; do \
	  test -f $$f && $(LN_S) $$f; done || true
	for f in $(EXTRA_BINDIST); do \
          cp -fpR $(top_srcdir)/$$f $(bindistdir)/$$f 2>/dev/null || \
          cp -fpR $(top_builddir)/$$f $(bindistdir)/$$f; done; \
	find $(bindistdir) -depth -type d -empty -exec rmdir {} \;
	binaries=$$(find $(bindistdir) -type f -exec file --mime {} \; | \
                $(GREP) -E "application/x-.*(mach|executable|sharedlib)" | cut -d: -f1); \
	export LDPATHS="$(subst -L,,$(filter -L%, $(BOOST_LDFLAGS)))"; \
        for binary in $$binaries; do $(top_srcdir)/tools/fix_libpaths.sh "$$binary"; done; \
	find $(bindistdir) \( -name \*.la -or -name \*.a \) -delete

tar_allow_lns = $${TAR-tar} cof - "$$tardir"
bindist-gzip: bindistdir
	tardir=$(bindistdir) && $(tar_allow_lns) | GZIP=$(GZIP_ENV) gzip -c >$(bindistdir).tar.gz
	rm -rf $(bindistdir)
bindist-bzip2: bindistdir
	tardir=$(bindistdir) && $(tar_allow_lns) | BZIP2=$${BZIP2--9} bzip2 -c >$(bindistdir).tar.bz2
	rm -rf $(bindistdir)


.PHONY: bindistdir bindist-gzip bindist-bzip2


### Unit Tests ###

test_libs = \
    $(sina_libs) \
    $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
    $(BOOST_SYSTEM_LIB)

LOG_COMPILER = $(top_srcdir)/tools/test_driver.sh

TEST_LOG_DRIVER = env AM_TAP_AWK='$(AWK)' $(SHELL) \
                  $(top_srcdir)/build-aux/tap-driver.sh --comments

check_PROGRAMS = \
    src/unit_tests/cseq_test \
    src/unit_tests/cseq_comparator_test \
    src/unit_tests/idset_test \
    src/unit_tests/kmer_test \
    src/unit_tests/kmer_search_test \
    src/unit_tests/famfinder_test \
    src/unit_tests/progress_test \
    src/unit_tests/rw_csv_test \
    src/unit_tests/query_arb_test

src_unit_tests_cseq_test_LDADD = $(test_libs)
src_unit_tests_cseq_test_LDFLAGS = $(ARB_LDFLAGS)
src_unit_tests_cseq_comparator_test_LDADD = $(test_libs)
src_unit_tests_cseq_comparator_test_LDFLAGS = $(ARB_LDFLAGS)
src_unit_tests_idset_test_LDADD = $(test_libs)
src_unit_tests_idset_test_LDFLAGS = $(ARB_LDFLAGS)
src_unit_tests_kmer_test_LDADD = $(test_libs)
src_unit_tests_kmer_test_LDFLAGS = $(ARB_LDFLAGS)
src_unit_tests_kmer_search_test_LDADD = $(test_libs)
src_unit_tests_kmer_search_test_LDFLAGS = $(ARB_LDFLAGS)
src_unit_tests_famfinder_test_LDADD = $(test_libs)
src_unit_tests_famfinder_test_LDFLAGS = $(ARB_LDFLAGS)
src_unit_tests_rw_csv_test_LDADD = $(test_libs)
src_unit_tests_rw_csv_test_LDFLAGS = $(ARB_LDFLAGS)
src_unit_tests_query_arb_test_LDADD = $(test_libs)
src_unit_tests_query_arb_test_LDFLAGS = $(ARB_LDFLAGS)
src_unit_tests_progress_test_LDADD = $(BOOST_UNIT_TEST_FRAMEWORK_LIB)

test_scripts = \
    tests/readwrite.test \
    tests/align.test \
    tests/search.test \
    tests/infocommands.test \
    tests/accuracy_kmer.test \
    tests/accuracy.test


check_DATA = test_data/ltp_reduced.arb

TESTS = $(check_PROGRAMS) $(test_scripts)

## add `make check-filtered P=<pattern>`
P ?= /
# call contains, MATCH, WORDS - filters by substr
contains = $(foreach v,$(2),$(if $(findstring $(1),$(v)),$(v),))
# call strip_term - removes leading ^ and trailing $
strip_term = $(patsubst ^%,%,$(patsubst %$$,%,$(1)))
# call add_term - adds leading ^ and trailing $
add_term = $(patsubst %,^%$$,$(1))
# call contains_term MATCH WORDS - filters by substr with ^ and $ matching start/end
contains_term = $(call strip_term,$(call contains,$(1),$(call add_term,$(2))))
filtered_tests = $(call contains_term,$(P),$(TEST_LOGS))
check-filtered:
	$(MAKE) $(AM_MAKEFLAGS) check TEST_LOGS="$(filtered_tests)"
.PHONY: check-filtered

## pass ARGS to tests
TEST_ARGS := $(check_DATA) $(TEST_EXTRA_ARGS)
export TEST_ARGS
$(test_scripts): $(bin_PROGRAMS) $(check_DATA) tests/test_helper.sh

## unzip source test data into build folder
test_data/%: test_data/%.xz
	mkdir -p $$(dirname $@)
	unxz -c $< > $@

clean-local: clean-local-check clean-local-docs

clean-local-check:
	rm -f test.*.xml
	rm -f test_data/ltp_reduced.arb.*
	rm -f test_data/ltp_reduced.sidx

.PHONY: clean-local-check clean-local-docs



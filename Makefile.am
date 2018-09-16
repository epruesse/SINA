# Put objects into subdirs
AUTOMAKE_OPTIONS = subdir-objects

# Default extension to use when no bla_SOURCES exists.
AM_DEFAULT_SOURCE_EXT = .cpp

# Macros (like boost) are in the m4 folder
ACLOCAL_AMFLAGS = -I m4

# Common flags for preprocessor and linker
AM_CPPFLAGS = \
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
    README.md \
    arb_LICENSE.txt \
    arb_PUBLIB_LICENSE.txt

EXTRA_DIST = \
    $(EXTRA_BINDIST) \
    LICENSE \
    CHANGELOG \
    README.md \
    arb_LICENSE.txt \
    arb_PUBLIB_LICENSE.txt \
    Makefile.am \
    configure.ac \
    autogen.sh \
    tools/latex2man \
    tools/latex2man.sty \
    tools/latex2man.cfg \
    tools/test_driver.sh

CLEANFILES =


### TARGETS ###

# calling pkglibdir our ARBHOME (.../lib/sina), we install 
# internal libs and binaries to .../bin and .../lib subdirs

bin_PROGRAMS = \
	src/sina \
	src/sina_probematch \
	src/kmer_test
lib_LTLIBRARIES = src/libsina.la

if ENABLE_FAT_TAR 
arbbindir = $(libdir)/sina/arb/bin
arblibdir = $(libdir)/sina/arb/lib

arbbin_PROGRAMS = arb_pt_server arb_message
arblib_DATA = \
    arb_tcp.dat \
    $(ARBHOME)/lib/libCORE$(LIBEXT) \
    $(ARBHOME)/lib/libARBDB$(LIBEXT) 

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
    src/copy_alignment.h \
    src/cseq.cpp \
    src/cseq.h \
    src/cseq_impl.h \
    src/cseq_comparator.cpp \
    src/cseq_comparator.h \
    src/famfinder.h \
    src/famfinder.cpp \
    src/graph.h \
    src/mesh.h \
    src/mesh_debug.h \
    src/mseq.cpp \
    src/mseq.h \
    src/null_filter.h \
    src/object_pipe.h \
    src/pipe.h \
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
    src/log.cpp \
    src/log.h \
    src/scoring_schemes.h \
    src/search.h \
    src/search_filter.cpp \
    src/search_filter.h \
    src/timer.h \
    src/tray.cpp \
    src/tray.h \
    src/kmer.h \
    src/kmer_search.cpp \
    src/kmer_search.h \
    include/pstream.h

src_libsina_la_LDFLAGS = \
    $(ARB_LDFLAGS) \
    $(BOOST_LDFLAGS)

src_libsina_la_LIBADD = \
    $(Z_LIB) \
    $(ARB_LIBS) \
    $(BOOST_PROGRAM_OPTIONS_LIB) \
    $(BOOST_THREAD_LIB) \
    $(BOOST_SERIALIZATION_LIB) \
    $(BOOST_SYSTEM_LIB) \
    $(TBB_LIB) \
    $(TBB_MALLOC_LIB)

sina_libs = \
    src/libsina.la \
    $(ARB_PROBE_LIBS) \
    $(ARB_HELIX_LIBS)

src_sina_probematch_LDADD = $(sina_libs)
src_kmer_test_LDADD = $(sina_libs)
src_sina_LDADD = $(sina_libs)

arb_pt_server$(EXEEXT): $(ARBHOME)/bin/arb_pt_server$(EXEEXT)
	cp $< $@
arb_message$(EXEEXT): $(ARBHOME)/bin/arb_message$(EXEXT)
	cp $< $@
libARBDB$(LIBEXT): $(ARBHOME)/lib/libARBDB$(LIBEXT) 
	cp $< $@
libCORE$(LIBEXT): $(ARBHOME)/lib/libCORE$(LIBEXT)
	cp $< $@
arb_tcp.dat: $(ARBHOME)/lib/arb_tcp.dat
	cp $< $@


### DOCUMENTATION ###


if ENABLE_DOCS
doc_DATA = doc/sina-manual.pdf doc/sina-manual.html
man1_MANS = doc/sina.man
endif
dist_doc_DATA = doc/sina.tex


CLEANFILES += sina-manual.log
LATEXTOMAN=$(top_srcdir)/tools/latex2man

export TEXINPUTS=$(top_srcdir):$(top_srcdir)/tools/:

CLEANFILES += $(doc_DATA) $(man1_MANS) doc/sina-manual.log  doc/sina-manual.toc doc/sina-manual.aux doc/sina-manual.out


%.txt: sina_%.txt
	-@mkdir `dirname $@`
	groff -T utf8 $< > $@

%.man: %.tex
	-@mkdir `dirname $@`
	$(SED) 's/{}//g; s/``/"/g; s/'\'\''/"/g;' $< | ${LATEXTOMAN} - $@

%-manual.tex: %.tex
	-@mkdir `dirname $@`
	${LATEXTOMAN} -L $< $@

%.pdf: %.tex
	-@mkdir `dirname $@`
	pdflatex -output-directory `dirname $@` $< 
	pdflatex -output-directory `dirname $@` $< 

%-manual.html: %.tex
	-@mkdir `dirname $@`
	$(SED) 's/{}//g; s/``/"/g; s/'\'\''/"/g;' $< | ${LATEXTOMAN} -H - $@


### PACKAGE ###

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
          cp -fpR $(top_srcdir)/$$f $(bindistdir)/$$f; done
	find $(bindistdir) -depth -type d -empty -exec rmdir {} \;
	binaries=$$(find $(bindistdir) -type f -exec file --mime {} \; | \
                grep -E "application/x-(mach|executable)" | cut -d: -f1); \
	export LDPATHS="$(subst -L,,$(filter -L%, $(BOOST_LDFLAGS)))"; \
        for binary in $$binaries; do $(top_srcdir)/tools/fix_libpaths.sh "$$binary"; done; \
	find $(bindistdir) -name \*.la -delete -o -name \*.a -delete

bindist-gzip: bindistdir
	tardir=$(bindistdir) && $(am__tar) | GZIP=$(GZIP_ENV) gzip -c >$(bindistdir).tar.gz
	rm -rf $(bindistdir)
bindist-bzip2: bindistdir
	tardir=$(bindistdir) && $(am__tar) | BZIP2=$${BZIP2--9} bzip2 -c >$(bindistdir).tar.bz2
	rm -rf $(bindistdir)


.PHONY: bindistdir bindist-gzip bindist-bzip2


### Unit Tests ###

test_libs = \
    $(sina_libs) \
    $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
    $(BOOST_SYSTEM_LIB)

TESTS =  $(check_PROGRAMS)

export ARBHOME
LOG_COMPILER = $(top_srcdir)/tools/test_driver.sh

check_PROGRAMS = \
    src/unit_tests/cseq_test \
    src/unit_tests/cseq_comparator_test \
    src/unit_tests/idset_test \
    src/unit_tests/kmer_test \
    src/unit_tests/kmer_search_test

check_DATA = \
    test_data/ltp_reduced.arb


test_data/%: test_data/%.xz
	mkdir -p $$(dirname $@)
	unxz -c $< > $@

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

clean-local: clean-local-check
.PHONY: clean-local-check
clean-local-check:
	rm -f test.*.xml

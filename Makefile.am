# Put objects into subdirs
AUTOMAKE_OPTIONS = subdir-objects

# Some extra files to put into source tarfile
EXTRA_DIST=LICENSE CHANGELOG README.md arb_LICENSE.txt arb_PUBLIB_LICENSE.txt

# Default extension to use when no bla_SOURCES exists.
AM_DEFAULT_SOURCE_EXT = .cpp

# Macros (like boost) are in the m4 folder
ACLOCAL_AMFLAGS = -I m4

# Common flags for preprocessor and linker
AM_CPPFLAGS = \
    -DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN \
    $(BOOST_CPPFLAGS) \
    $(ARB_CPPFLAGS)
AM_LDFLAGS = \
    $(BOOST_LDFLAGS) \
    $(BOOST_THREAD_LDFLAGS) \
    $(BOOST_PROGRAM_OPTIONS_LDFLAGS) \
    $(ARB_LDFLAGS)

# calling pkglibdir our ARBHOME (.../lib/sina), we install 
# internal libs and binaries to .../bin and .../lib subdirs
arbbindir = $(pkglibdir)/bin
arblibdir = $(pkglibdir)/lib


### TARGETS ###

bin_SCRIPTS = src/sina
arbbin_PROGRAMS = src/sina.real arb_pt_server arb_message
arblib_DATA = $(ARB_FILES) arb_tcp.dat
noinst_PROGRAMS = src/sina_probematch src/kmer_test
noinst_LTLIBRARIES = src/libsina.la

arb_message_SOURCES=
arb_pt_server_SOURCES=

CLEANFILES    = $(bin_SCRIPTS)

### TARGET CONFIGS ###

src_libsina_la_SOURCES = \
    src/aligned_base.cpp \
    src/aligned_base.h \
    src/cseq.cpp \
    src/cseq.h \
    src/cseq_comparator.cpp \
    src/cseq_comparator.h \
    src/mseq.cpp \
    src/mseq.h \
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
    src/align.cpp \
    src/align.h \
    src/search_filter.cpp \
    src/search_filter.h \
    src/tray.cpp \
    src/tray.h \
    src/alignment_stats.cpp \
    src/alignment_stats.h \
    src/kmer_search.cpp \
    src/kmer_search.h

src_libsina_la_LDFLAGS = -static $(ARB_LDFLAGS)

sina_libs = \
    src/libsina.la \
    $(Z_LIB) \
    $(ARB_LIBS) \
    $(ARB_PROBE_LIBS) \
    $(ARB_HELIX_LIBS) \
    $(BOOST_PROGRAM_OPTIONS_LIBS) \
    $(BOOST_THREAD_LIBS) \
    $(BOOST_SERIALIZATION_LIBS)

src_sina_real_SOURCES = src/sina.cpp 
src_sina_real_LDADD = $(sina_libs)


src_sina_probematch_LDADD = $(sina_libs)
src_kmer_test_LDADD = $(sina_libs)

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
doc_DATA_DIST = doc/sina-manual.text


CLEANFILES += sina-manual.log 
LATEXTOMAN=$(top_srcdir)/tools/latex2man
export TEXINPUTS=$(top_srcdir):$(top_srcdir)/tools/:

CLEANFILES += $(dist_doc_DATA) $(man1_MANS) doc/sina-manual.log  doc/sina-manual.toc doc/sina-manual.aux doc/sina-manual.out


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

BINDISTDIR = sina-$(VERSION)-$(OSNAME)

bindist: $(BINDISTDIR)

$(BINDISTDIR): all
	rm -rf $(BINDISTDIR)
	$(MAKE) DESTDIR=$(abspath $(BINDISTDIR)) install
	mv -f $(BINDISTDIR)/$(prefix)/* $(BINDISTDIR)
	cd $(BINDISTDIR); for f in bin/sina share/doc/sina*; do \
	  test -f $$f && $(LN_S) $$f; done || true
	find $(BINDISTDIR) -depth -type d -empty -exec rmdir {} \;
	cp $(EXTRA_DIST) $(BINDISTDIR)

bindist-gzip: $(BINDISTDIR).tar.gz

$(BINDISTDIR).tar.gz: $(BINDISTDIR)
	tar -cv $^ | gzip --best > $@

.PHONY: $(BINDISTDIR) $(BINDISTDIR).tar.gz
### Unit Tests ###

test_libs = $(sina_libs) \
    $(BOOST_UNIT_TEST_FRAMEWORK_LIBS) \
    $(BOOST_SYSTEM_LIBS)

TESTS =  $(check_PROGRAMS)

export ARBHOME
LOG_COMPILER = $(top_srcdir)/tools/test_driver.sh

check_PROGRAMS = \
    src/unit_tests/cseq_test \
    src/unit_tests/cseq_comparator_test

src_unit_tests_cseq_test_LDADD = $(test_libs)
src_unit_tests_cseq_test_LDFLAGS = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS) $(ARB_LDFLAGS)
src_unit_tests_cseq_comparator_test_LDADD = $(test_libs)
src_unit_tests_cseq_comparator_test_LDFLAGS = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS) $(ARB_LDFLAGS)

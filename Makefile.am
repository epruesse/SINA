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
    $(ARB_CPPFLAGS)
AM_LDFLAGS = \
    $(BOOST_LDFLAGS) \
    $(BOOST_THREAD_LDFLAGS) \
    $(BOOST_PROGRAM_OPTIONS_LDFLAGS) 

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

CLEANFILES    = $(bin_SCRIPTS)

### TARGET CONFIGS ###

src_libsina_la_SOURCES = \
    src/aligned_base.cpp \
    src/cseq.cpp \
    src/cseq_comparator.cpp \
    src/mseq.cpp \
    src/pseq.cpp \
    src/query_arb.cpp \
    src/query_pt.cpp \
    src/rw_fasta.cpp \
    src/rw_arb.cpp \
    src/log.cpp \
    src/align.cpp \
    src/search_filter.cpp \
    src/tray.cpp \
    src/alignment_stats.cpp \
    src/kmer_search.cpp
src_libsina_la_LDFLAGS = -static

sina_libs = \
    src/libsina.la \
    $(Z_LIB) \
    $(ARB_SLIBS) \
    $(ARB_LIBS) \
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
libARBDB.so: $(ARBHOME)/lib/libARBDB.so 
	cp $< $@
libCORE.so: $(ARBHOME)/lib/libCORE.so 
	cp $< $@
arb_tcp.dat: $(ARBHOME)/lib/arb_tcp.dat
	cp $< $@


### DOCUMENTATION ###

dist_doc_DATA = doc/sina-manual.pdf doc/sina-manual.html
man1_MANS = doc/sina.man
CLEANFILES += sina-manual.log 
LATEXTOMAN=$(top_srcdir)/tools/latex2man
export TEXINPUTS=$(top_srcdir):$(top_srcdir)/tools/:

CLEANFILES += $(dist_doc_DATA) $(man1_MANS) doc/sina-manual.log  doc/sina-manual.toc doc/sina-manual.aux doc/sina-manual.out


%.txt: sina_%.txt
	-@mkdir `dirname $@`
	groff -T utf8 $< > $@

%.man: %.tex
	-@mkdir `dirname $@`
	sed 's/{}//g; s/``/"/g; s/'\'\''/"/g;' $< | ${LATEXTOMAN} - $@

%-manual.tex: %.tex
	-@mkdir `dirname $@`
	${LATEXTOMAN} -L $< $@

%.pdf: %.tex
	-@mkdir `dirname $@`
	pdflatex -output-directory `dirname $@` $< 
	pdflatex -output-directory `dirname $@` $< 

%-manual.html: %.tex
	-@mkdir `dirname $@`
	sed 's/{}//g; s/``/"/g; s/'\'\''/"/g;' $< | ${LATEXTOMAN} -H - $@


### PACKAGE ###

tgz: sina-$(VERSION).tgz

# create self-contained tgz
RELEASEDIR := sina-$(VERSION)
sina-$(VERSION).tgz: 
	TMPDIR=$(shell mktemp -dt sina-buildpkg-XXXX);\
	$(MAKE) DESTDIR="$$TMPDIR" install;\
	mv "$$TMPDIR$(prefix)" "$$TMPDIR/$(RELEASEDIR)";\
	(cd "$$TMPDIR/$(RELEASEDIR)"; ln -s -t . bin/sina share/doc/sina/*;); \
	find "$$TMPDIR" -depth -type d -empty -exec rmdir {} \; ;\
	tar -czvf sina-$(VERSION).tgz -C "$$TMPDIR" . ;\
	rm -rf "$$TMPDIR"

.PHONY: sina-$(VERSION).tgz

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
src_unit_tests_cseq_test_LDFLAGS = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)
src_unit_tests_cseq_comparator_test_LDADD = $(test_libs)
src_unit_tests_cseq_comparator_test_LDFLAGS = $(BOOST_UNIT_TEST_FRAMEWORK_LDFLAGS)

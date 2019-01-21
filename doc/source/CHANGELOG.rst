.. program:: sina

Changelog
=========

Version 1.4.0:
--------------

 - process sequences in parallel (:issue:`17`, :issue:`31`)
 - add support for gzipped read/write (:issue:`29`)
 - add support for "-" to read/write using pipes
 - remove internal pipeline in favor of TBB
 - add :option:`--add-relatives`; adding search result to output (:issue:`19`)
 - add logging with variable verbosity (:issue:`14`)
 - be smart about locating arb_pt_server binary (:issue:`30`)

Version 1.3.5:
--------------
 - report number of references discarded due to configured constraints
 - fix crash if no acceptable references found for a query
 - fix :option:`--search` causes a program option error (:issue:`28`)
 - fix race condition in terminating PT server

Version 1.3.4:
--------------
 - build binary releases for macOS and Linux (:issue:`26`)
 - fix "search.h" missing in source tar ball (:issue:`27`)

Version 1.3.3:
--------------
 - add option :option:`--fasta-write-dots`; writes dots on edges
 - add option :option:`--fasta-write-dna`; writes T/t instead of U/u (:issue:`24`)
 - fix PT server fails to build if ARBHOME not set (:issue:`15`)
 - fix psina not installed to $bindir
 - fix tab character in sequence causes sequence to be skipped (:issue:`21`)
 - fix last line of input FASTA ignored if missing newline (:issue:`16`)
 - fix :option:`--db` parameter demanded even if not required due to use of :option:`--prealigned`
 - fix SIGPIPE race on PT server shutdown (:issue:`11`)

Version 1.3.2:
--------------
 - split :option:`--help` into "common" and advanced options (:option:`--help-all`)
 - add psina wrapper script (runs parallel instances of SINA to align a
   single FASTA file)
 - fix memory access failure in cseq
 - fix memory access failure in mseq
 - fix crash on all references removed by filters
 - don't exit(1) on :option:`--help` (:issue:`9`)
 - added README.md (:issue:`5`)

Version 1.3.1:
--------------
 - add OSX support
 - change license to GPL
 - remove limitation on ARB integration mode
 - move revisioning to git
 - fix compilation with CLANG

Version 1.3.0:
--------------
 - dropped support for ARB 5.x

Version 1.2.13:
---------------
 - uppercase aligned bases if lowercase=unaligned
 - fix manual typos (thx to Mohamed El-hadidi)
 - search-db defaults to pt-db
 - search-port defaults to pt-port if search/align DBs are identical
   fixes unnecessary start of two PT servers (thx to Christian Wurzbacher)
 - change default for lca-quorum to 0.7
 - change default for search-min-sim to 0.7
 - be smarter about recoginizing FASTA format files and creating output FASTA name
   (".frn", ".fna", ".fas", "/dev/stdin" as input, ".fasta.aligned" and
   "/dev/stdout" as output)
 - write sequence ID in first column of CSV output
 - add fasta-block and fasta-idx options allowing to process only specific
   smaller blocks of larger fasta files (for parallelization)

Version 1.2.12:
---------------
 - use same ARB field type for align_ident_slv as SILVA uses
 - skip sequences with non-IUPAC characters when building reference
   and when loading sequences to be aligned from ARB file
   (complaint is issued on stderr)

Version 1.2.11:
---------------
 - fix :option:`--fs-req` was ignored
 - added option :option:`--calc-idty`
   Computes the minimum identity of the aligned query sequence with any
   of the reference sequences used for alignment. The value is exported
   in align_slv_idty.
 - added option :option:`--min-idty` IDTY
   Excludes sequences with align_slv_idty < IDTY from FASTA output.
   Implies :option:`--calc-idty`.

Version 1.2.10:
---------------
 - added option :option:`--fs-no-graph`
   Uses a column profile with PSP score as template (instead of the POA method)
   This feature is merely for completeness sake and evaluation. With SILVA
   SSU the POA based method is much more accurate.
 - changed default for :option:`--fs-cover-gene` to 0 (faster)
   The cover-gene feature only makes sense if `:option:`--gene-start` and :option:`--gene-end` are
   set such that the reference actually contains sequences touching these
   boundaries. If this is not the case, the reference selection algorithm
   wastes time with a futile search.
 - use unix socket as default for :option:`--ptport` and :option:`--search-port`
   Using "/tmp/sina_<PID>.socket" is a more suitable default than
   "localhost:4040", as it runs less risk of accessing a different PT server
   than intended.
 - fix inconsistencies in generated meta data fields and log output
 - updated ARB components to SVN revision 8225
 - added option :option:`--write-used-rels`
   The field used_rels is interpreted by ARB as the field containing the
   reference sequences that were used during alignment.
 - no longer write full_name content when exporting meta data encoded in
   the FASTA header
 - re-add clamped align_quality_slv
 - fix score normalization (scores > 1 were possible when fs-weight > 0)
 - fix calculation of bp score when orig-db no set (default ptdb)
 - added option :option:`--fs-req-gaps` n
   Ignores reference sequences having less than n gaps before the last base.
   I.e.: Ignores "unaligned" sequences. This is useful when running SINA
   out of ARB to prevent accidental alignment against unaligned sequences.
 - added options :option:`--search-iupac`, :option:`--search-correction` and :option:`--search-cover`
   These options configure how the "distance" (identity, similarity, ...)
   is calculated.
 - skip FASTA input sequences that contain invalid characters
   (i.e. not IUPAC encoded bases, '.', '-' or white space)

Version 1.2.9:
--------------
 - fixed sequence not filled with gap characters after copying full alignment

Version 1.2.8:
--------------
 - made --extra-fields actually load multiple fields from arb file
 - fixed sequence not filled with gap characters after copying subalignment
 - updated ARB components to SVN revision 7985
 - added changelog :)
Commandline Reference
=====================

SINA adds sequences to an existing multiple sequence alignment
(MSA). It can also execute a homology search based on the computed
alignment and generate a per sequence classifications from the search
results.


Synopsis
--------

.. program:: sina

**sina** [*options*] :option:`-i` <*unaligned*> :option:`-r` <*reference*> :option:`-o` <*aligned*>

**sina** [ :option:`-h` | :option:`--help` | :option:`--help-all` | :option:`--version` | :option:`--has-cli-vers` ]


			   
Description
-----------

:program:`sina` aligns the sequences in the file ``<unaligned>`` to match the alignment in ``<reference>`` and places the aligned sequences in the file ``<aligned>``.


General Options
---------------

.. program:: sina

.. option:: -h, --help

   Displays brief command line description.

.. option:: -H, --help-all

   Displays full command line description. If in doubt, refer to this
   as it will always be in sync with your installation of SINA.

.. option:: -V, --version

   Shows the version of SINA.

.. option:: -i filename, --in=filename (-)

   Specifies the file containing the input sequences. Allowable file
   formats are ARB and FASTA. Using "**-**" will read sequences from
   standard input. Using "**:**" *when running from within an ARB
   terminal** will read sequences from the current ARB database.

.. option:: -o filename, --out=filename (-)

   Specifies the file to which the aligned sequences will be
   written. Allowable file formats are ARB and FASTA. Using "**-**"
   will write sequences to standard output. Using "**:**" *when
   running from within an ARB terminal** will read sequences from the
   current ARB database.

.. option:: -r filename, --db=filename

   Specifies the file containing the reference alignment. This file
   must be in ARB format. To convert a reference alignment from FASTA
   to ARB format, run:

   ``sina -i reference.fasta --prealigned -o reference.arb``

.. option:: -t [all], --turn [=all]

   Enables turn check stage: Sequences not oriented in accordance with
   the reference database will be reverse complemented as needed.

   If *all* is specified, sequences will also be tested for only
   reversal or only complemented (this should only be necessary if
   your data was mishandled).

.. option:: -S, --search

   Enables the search stage. See `Search & Classify`_ below for more
   information.

.. option:: -P, --prealigned

   Disables the alignment stage. This is useful if you have already
   aligned sequences you wish to pass directly into the search stage,
   or if you want to use SINA to convert between any of its supported
   file formats.

.. option:: -v, --verbose

   Increase logging verbosity. Can be specified multiple times.

.. option:: -q, --quiet

   Decrease logging verbosity. Can be specified multiple times.

.. option:: --log-file=filename

   Specify log file. The output written to the log file will always be
   verbose and is not affected by using :option:`-v` or :option:`-q`.

.. option:: --meta-fmt=[none|header|comment|csv]

   Configures how meta data (such as alignment score or sequence classification results) are to be exported.

   **none**
     No output other than in the log is generated.

   **header**
     Appends meta data as ``[key=value]`` pairs to the FASTA header line

   **comment**
     Appends meta data as ``; key: value`` lines between the
     FASTA header and the sequence data.

   **csv**
     Writes meta data into a CSV side car file.

.. option:: -p, --threads (automatic)

   Override automatic detection of the number of threads used by
   SINA. This is usually only necessary if you need to constrain SINA
   to a lower number of threads. According to the Intel engineers
   whose *Threaded Building Blocks* library does the thread number
   detection for SINA, the only reason to use this parameter should be
   scalability testing.

.. option:: --num-pts (1)

   Set the maximum number of ARB PT server instances used by SINA. See
   also :option:`--fs-engine` below. If you are using the
   **pt-server** engine, this setting will be the limiting factor in
   your throughput. Be aware, however, that each PT server will occupy
   additional system memory. Choosing a too high value may cause SINA
   to fail with out-of-memory errors.

.. option:: --add-relatives=n (0)

   Add up to *n* reference sequences for each query sequence to the
   output file. If `Search & Classify`_ is enabled via :option:`--search`, the
   reference sequences are selected from the search result. Otherwise,
   they are selected from the query's alignment reference set.

   If the source set is smaller than *n*, no further sequences are
   added to the output. Sequences already included are skipped, but
   count towards the *n* of the query sequence.


   
Reference Selection Options
---------------------------

These options configure how the set of reference sequences used during alignment is selected from the configured reference database.

.. program:: sina

.. option:: --fs-engine=[internal|pt-server]

   Selects the search engine used to find closely related reference
   sequences for the alignment stage.

   **pt-server**
     Uses the ARB PT server to execute the k-mer search. The ARB
     PT server is a truncated suffix trie implementation implemented
     as part of the ARB package.

   **internal**
     Uses an internal k-mer search implementation.

.. option:: --fs-kmer-len=k (10)

   Set the size of *k* for the reference search. For SSU rRNA
   sequences, the default of 10 is a good value. For different
   sequence types, different values may perform better. For 5S, for
   example, 6 has shown to be more effective.

.. option:: --fs-min=n (15)

   Set the minimum number of reference sequences used for each query.

.. option:: --fs-max=n (40)

   Set the maximum number of reference sequences used for each query.

.. option:: --fs-msc=n (0.7)

   Set the minimum similarity reference sequences are required to have
   with the query sequence. This affects the range between
   :option:`--fs-min` and :option:`--fs-max`.

.. option:: --fs-req=n (1)

   Set the minimum number of reference sequences that must be found in
   order to attempt alignment. If fewer sequences than indicated here
   are found, the respective query sequence will be discarded.

.. option:: --fs-req-full=n (1)

   Set the minimum number of *full length* (see
   :option:`--fs-full-len`) reference sequences that must be included
   in the selected reference set. The search will proceed regardless
   of other settings until this setting has been satisfied. If it
   cannot be satisfied by any sequence in the reference database, the
   query sequence will be discarded.

   This setting exists to ensure that the entire length of the query
   sequence will be covered in the presence of partial sequences
   contained within your reference database.

   **Note:**
     If you are working with sequences other than 16S, you need to
     adjust this value or the value of :option:`--fs-full-len`
     accordingly. In particular when working with short reference
     sequences, this setting may prevent any acceptable reference
     sequences from being found, leading to no sequences being aligned.

.. option:: --fs-full-len=n (1400)

   Set the minimum length a sequence is required to have to be
   considered *full length*.

.. option:: --fs-req-gaps=n (10)

   Set the minimum number of gaps a reference sequence is required to
   contain to be considered. This setting ensures that unaligned
   sequences contained within the reference database are not used as
   reference (this may happen when SINA is used from within ARB).

.. option::  --fs-min-len=n (150)

   Set the minimum length reference sequences are required to
   have. Sequences shorter than this will not be included in the
   selection.

   **Note:**
     If you are working with particularly short reference sequences,
     you will need to lower this settings to allow any reference
     sequences to be found.
   

.. _`Search & Classify`:

Search & Classify Options
-------------------------

When enabled via :option:`--search`, SINA will execute a homology
search. Unlike most homology search tools, SINA uses the inferred
multiple sequence alignment to determine the similarity of each query
with the reference sequences, rather than computing pairwise optimal
alignments. **The similarity values will therefore be generally lower
than the results of a pairwise alignment based homology search**.

Based on the search results, SINA can be instructed to compute a
lowest common ancestor (LCA) based classification of the input
sequences. For this, your reference database must include a field
containing taxonomic classifications for each reference sequence. The
field contents must be in the format t *Domain;Phylum;...*. SINA will
compute query classifications as the deepest classification shared by
at least the fraction :option:`--lca-quorum` of the search result.

.. program:: sina

.. option:: --search-db=filename (=db)

   Specify an alternate reference database to use for search and
   classify. This can be useful if you have a specially curated
   alignment reference, but wish to search a larger set of sequences
   for classification purposes.

.. option:: --search-engine=[internal|pt-server]

   Override the value of :option:`--fs-engine` for use within the
   search module.

.. option:: --search-min-sim=id (0.7)

   The minimum fractional identity each result sequence must have with
   the query.
   
.. option:: --search-max-result=n (10)

   The maximum number of search results to return for each query sequence.

.. option:: --lca-fields=names

   Enables the classification stage. The parameter *name* must be a
   colon or comma separated list of field names in the search database
   containing the classification reference data. When using a SILVA
   ARB database as reference, the fields `tax_slv`, `tax_embl` and
   `tax_ltp` contain the reference classifications according to the
   SILVA, EMBL-EBI/ENA and LTP taxonomies, respectively. When using a
   SILVA SSU ARB database, the fields `tax_gg` and `tax_rdp` are
   available additionally, containing the reference classifications
   according to RDP II and Greengenes, respectively.

.. option:: --lca-quorum=fraction (0.7)

   Sets the fraction of the search result that must share the same
   classification. Using the default parameters
   :option:`--search-max-result`\=10 and :option:`--lca-quorum`\=0.7, this
   means that the deepest classification shared by 7 out of the top 10
   search results is chosen for the query sequence.


Advanced Options
----------------

.. option:: --show-conf

   Print the values of all configuration options (including defaults)
   at startup.

.. option::  --intype=[auto|arb|fasta] (auto)q

   Set the file format for :option:`--in`. If set to *auto* (default),
   the type is selected based on the file extension.
	     
.. option::  --outtype=[auto|arb|fasta] (auto)

   Set the file format for :option:`--out`. If set to *auto* (default),
   the type is selected based on the file extension.
	     
.. option::  --preserve-order

   Preserve the order of the input sequences in the output.

.. option::  --max-in-flight=n (2 * number of CPU cores)

   Set the maximum number of sequences "in flight", i.e. processed in
   parallel.
   
.. option::  --has-cli-vers=cliversion

   Verify that this version of SINA supports the CLI version
   **cliversion**. Exits immediately with exit code 0 if true and 1 if
   false.
   
.. option::  --no-align

   Backwards compatibility alias for :option:`--prealigned`.

.. option::  -f fields, --fields-fields

   Configures the set of fields written to the output file.

   .. todo:: reference description of SINA generated fields

Logging Options
---------------

.. option:: --show-diff

   Show differences between the inferred alignment and the original
   alignment. Requires either aligned sequences to be passed into sina
   via :option:`--in` or that a database with matching names is
   specified using :option:`--orig-db`.

   .. todo:: Fix text below

   This flag enables visualization of alignment differences. This
   feature allows you to quickly assess where your alignment differs
   from the one SINA computed. By also showing you the alignment of
   the reference sequences used for aligning the sequence, you can get
   an idea of why SINA came to its conclusions.  Many cases of
   "sub-optimal" alignment can be attributed to inconsistent
   alignment of the reference sequences.  To fix such problems, you
   could either correct the alignment of the reference sequences or
   add your corrected sequence to the reference alignment.

   Alignment difference visualization requires that the input
   sequences be already aligned in a way compatible with the used
   reference alignment. For positions at which the original alignment
   and the alignment computed by SINA differ, output as shown below
   will be printed to the log::

     Dumping pos 1121 through 1141:
     ---------  4 14 16-17 21 24
     G-C-AGUC-  40 <---(%% ORIG %%)
     GCA--GUC-  41 <---(## NEW ##)
     GCA-AGUC-  0-3 5-13 15 18-20 22-23 25-27 29-39
     GCAA-GUC-  28


   In this case, the bases '\texttt{C}' and '\texttt{A}' where placed
   in other columns than as per the original alignment. The original
   alignment is marked with \texttt{<-{}--(\%\% ORIG \%\%)}. The new
   alignment is marked with \texttt{<-{}--(\#\# NEW \#\#)}. The
   numbers to the right of the alignment excerpt indicate the indices
   of the sequences in the alignment reference (field

.. option:: --show-dist

   Show distance to original alignment

   .. todo:: describe values generated
	     
.. option::  --orig-db=filename

   Specify a database containing the original alignments for use with
   :option:`--show-dist` and :option:`--show-diff`. The sequence names
   in the input file and in the reference database must match exactly.

.. option::  --colors

   Use ANSI codes to show alignments dumped by :option:`--show-diff`
   in color.

ARB I/O Options
---------------

These options configure behavior supported only by the ARB backend for
input and output sequences.

.. option:: --markcopied

   Set *Mark* on sequences copied from the reference.

   .. todo: This feature is broken after reimplementing copy
	    
.. option:: --markaligned

   Set *Mark* on sequences updated or added by alignment stage.
	     
.. option:: --prot-level=n (4)

   Set the *protection level* to use when writing sequences to the
   output database.

.. option:: --select-file=filename

   Instead of iterating over the entire input database, process only
   the sequences listed in *filename*. The names must match the ARB
   *name* field and be separated by newlines. Use "**-**" to read from
   standard input.
   
.. option:: --select-step=n (1)

   Process only every *n*\th sequence. Can be combined with
   :option:`--select-file` and :option:`--select-skip`.

.. option:: --select-skip=n (0)

   Do not process the first *n* sequences. Can be combined with
   :option:`--select-file` and :option:`--select-step`.

FASTA I/O Options
-----------------

These options configure behavior supported only by the FASTA backend
for input and output sequences.

.. option:: --line-length=n (0)

   Output sequences using at most *n* characters per line. Set to 0 to
   place the entire output sequence on one line.
	    
.. option:: --min-idty=id

   Exclude sequences sharing less than *id* fractional identity with
   any of the alignment reference sequences from the output. Implies
   :option:`--calc-idty`.

.. option:: --fasta-write-dna

   Write output sequences as DNA, rather than the default
   RNA. (I.e. use T and t rather than U and U).

.. option:: --fasta-write-dots

   Use dots (".") rather than dashes ("-") for gaps that indicate
   missing data rather than an actual insertion/deletion. Most often,
   those are only the terminal gaps at the ends of the alignment.

   .. todo:: Check whether internal dots are handled correctly.

.. option:: --fasta-idx=n

   Only process sequences starting withing the *n*\th block of bytes
   within the input FASTA file.

   .. deprecated:: 1.4
      This feature was superseded by the built-in parallelization.

.. option:: --fasta-block=size

   Sets the size in bytes for the blocks used by :option:`--fasta-idx`.

   .. deprecated:: 1.4
      This feature was superseded by the built-in parallelization.

Alignment Options
-----------------

.. option:: --realign

   Forces computing the alignment of query sequences even if a
   reference sequence containing the exact sequence was found. Without
   this flag, SINA will copy the alignment from the reference
   sequence.

.. option:: --overhang=[attach|remove|edge] (attach)

   Configures how unaligned bases at the edge of the alignment (overhanging bases) should be handled.

   **attach**
     Overhang bases will be placed next to the last aligned base consecutively.

   **remove**
     Overhang bases will be deleted.
     
     .. todo:: This feature appears to be broken.

   **edge**
     Overhang bases will be placed next to the outer edge of the alignment.

.. option:: --lowercase=[none|original|unaligned]

   Configures which bases should be written using lower case characters.

   **none**
     All bases will use upper case characters

   **original**
     All bases will be written using the case they had in the input data.

   **unaligned**
     Aligned bases will be written in upper case; unaligned bases will
     be written in lower case. This serves to mark sections of the
     query sequences that could not be aligned because they were
     insertions (internal or edge) with respect to any of the
     reference sequences.

.. option:: --insertion=[shift|forbid|remove]

   Configures how the alignment width is preserved.

   **shift**
     The alignment is executed without constraining insertion
     sizes. Insertions for which insufficient columns exist between
     the adjoining aligned bases are force fitted into the alignment
     using NAST. That is, the minimum number of aligned bases to the
     left and right of the insertion are moved to accommodate the
     insertion.

     This mode will add warnings to the log for each sequence in which
     aligned bases had to be moved.

   **forbid**
     The alignment is executed using a scoring scheme disallowing
     insertions for which insufficient columns exist in the alignment.

     This mode causes less "misalignments" than the **shift** mode as
     it computes the best alignment under the constraint that no
     columns may be added to the alignment. However, it will not show
     if the computed alignment suffered from a lack of empty columns.

   **remove**
     The alignment is executed without constraining insertion
     sizes. Insertions larger than the number of columns between the
     adjoining aligned bases are truncated.

     While this mode yields the most accurate alignment for sequences
     with large insertions, it should be used with care as it modifies
     the original sequence.
   
.. option:: --fs-no-graph

   Instructs SINA to use a profile vector instead of a DAG to perform
   the alignment. That is, the base frequencies for all selected
   reference sequences are collected into a vector and the query is
   aligned to this vector weighting the alignment scores according to
   the respective frequencies.

   This feature was added in response to the requests of a reviewer of
   the original SINA publication and only intended to demonstrate that
   the DAG/POA approach is superior to the profile vector approach. Do
   not use this other than for testing.

.. option:: --fs-weight=weight (1)

   Adjust the weight factor for the frequency at which a node was
   observed in the reference alignment. Use 0 to disable weighting.

   This feature prefers the more common placement for bases with
   inconsistent alignment in the reference database.

.. option:: --match-score=n (2)

   Configures the score given for a match (should be positive).

.. option:: --mismatch-score=n (-1)

   Configures the score given for a mismatch (should be negative).
   
.. option:: --pen-gap=n                 gap open penalty (5)

   Configures the penalty subtracted from the score for opening a gap
   (should be positive).
   
.. option:: --pen-gapext=n

   Configures the penalty subtracted from the score for extending a
   gap (should be positive).
   
.. option:: --debug-graph

   Writes the DAG computed from the reference sequences for each query
   sequences to disk in dot format.
   
.. option:: --use-subst-matrix

   Weights the match and mismatch scores according to the overall base
   frequencies observed in the database.

   This feature is experimental and does not currently improve the
   results.
   
.. option:: --write-used-rels

   Writes the names of the alignment reference sequences into the
   field `used_rels`. This option allows using the ARB *mark used
   rels* feature to highlight the reference sequences used to align a
   given query sequence.

.. option:: --calc-idty

   Computes the highest similarity the aligned query sequence has with
   any of the sequences in the alignment reference set. The value is
   written to the field *align_ident_slv*.

Advanced Reference Selection Options
------------------------------------

.. option:: --ptdb=filename

   Alias of :option:`--db` for backwards compatibility.

.. option:: --ptport=port_or_socket (:/tmp/sina_pt_<pid>)

   Configures the port or socket on which the ARB PT server for the
   reference alignment is expected or started. To use a TCP port,
   specify *<hostname>*:*<port>*. If *<hostname>* is not `localhost`,
   the PT server must be launched externally. To use a Unix socket,
   specify `:<filename>`.

   When :option:`--num-pts` is greater than 1, the additional PT
   servers port names are generated by appending the respective
   number. Using port numbers greater of equal to 10000 will therefore
   not work.

   By default, the file `/tmp/sina_pt_<pid>` is used, where `<pid>` is
   replaced by the process ID of the SINA instance.
   
.. option:: --fs-kmer-no-fast

   Use all k-mers occurring in the query sequence in the search. By
   default, only k-mers starting with an A are used for extra
   performance.
   
.. option:: --fs-kmer-mm=n (0)

   Allow k-mer matches in the reference database to contain *n*
   mismatches. This feature is only supported by the **pt-server**
   search engine and requires substantial additional compute time (in
   particular for *n* > 1).

.. option:: --fs-kmer-norel

   Use absolute (number of shared k-mers) match scores in the kmer
   search rather than relative (number or shared k-mers divided by
   length of reference sequence) match scores.

.. option:: --fs-msc-max=id (2)

   Overriding all other options, reference sequences having a
   similarity with the query higher than this value are excluded from
   the alignment reference.

   This option artificially increases the difficulty of the alignment
   by increasing the distance of a query to any reference found in the
   database. It's purpose of this option is to generate a sufficiently
   large *N* of test cases for statistical analysis of SINA's accuracy
   for sequences distant to the reference alignment.
   
.. option:: --fs-leave-query-out

   Excludes sequences from the alignment reference sharing the same
   name as the respective query sequence. (For testing and evaluation).

.. option:: --gene-start=n (0)

   Sets the beginning of the gene within the reference alignment. See
   :option:`--fs-cover-gene`.

.. option:: --gene-end=n (0)

   Sets the end of the gene within the reference alignment. See
   :option:`--fs-cover-gene`.

.. option:: --fs-cover-gene=n (0)

   Similar to :option:`--fs-req-full`, this option requires a total of
   *n* sequences to cover each the beginning and the end of the gene
   within the alignment. This option is more precise than
   :option:`--fs-req-full`, but requires that the column numbers for
   the range in which the full gene is expected be specified via
   :option:`--gene-start` and :option:`--gene-end`.

.. option:: --filter=name

   Chooses an *ARB posvar filter* to use for weighting alignment
   positions by their variability.
	     
.. option:: --auto-filter-field=name

   Configures a database field using which the value of
   :option:`--filter` is determined by majority vote from the selected
   reference sequences. Since the filters are usually computed at
   domain level, this approach is usually sufficient to select an
   appropriate filter. For SILVA database, the field `tax_slv` contains
   appropriate data.
   
.. option:: --auto-filter-threshold arg
	    
   Sets the minimum quorum required for automatic filter
   selection. See :option:`--lca-quorum` for information on how the
   value is interpreted.

.. option:: --fs-oldmatch

   Use the pre-1.6.0 implementation for composing the alignment
   family. Requires :option:`--fs-engine` = ``pt-server``.

Search & Classify Options
-------------------------

.. option:: --search-port=port_or_socket (:/tmp/sina_pt2_<pid>)

   See :option:`--ptport`. This option sets the port for the search
   database. It is only used if :option:`--search-db` is specified and
   its value differs from the one given by :option:`--db`.
   
.. option:: --search-all

   Calculate the similarity of the query sequences with **all**
   reference sequences. Normally, SINA will only calculate the
   similarity for the sequences returned by a k-mer based similarity
   search. See also :option:`--search-kmer-candidates`.
   
.. option:: --search-no-fast              don't use fast family search

   See :option:`--fs-kmer-no-fast`. This option configures the same
   behavior for the search stage.
   
.. option:: --search-kmer-candidates=n (1000)

   Configures the number of candidate reference sequences retrieved
   from the k-mer based search. For each candidate, the MSA based
   similarity is calculated and the search result based on these
   numbers. A value for *n* one or two orders larger than
   :option:`--search-max-result` is usually quite sufficient.

.. option:: --search-kmer-len=n (10)

   See :option:`--fs-kmer-len`. Sets *k* for the kmer based candidate
   search.

.. option:: --search-kmer-mm arg

   See :option:`--fs-kmer-mm`. Sets the number of allowed mismatches
   within each kmer. Only available with the **pt-server** search
   engine.

.. option:: --search-kmer-norel

   See :option:`--fs-kmer-norel`. Configures the candidate search to
   use absolute rather than length-relative scores for ordering the
   results.

.. option:: --search-ignore-super

   Omit reference sequences of which the query is an exact sub-string
   from the result. Useful for testing and evaluation of the
   classification feature.

.. option:: --search-copy-fields=fields

   Specifies a (colon or comma separated) list of meta-data fields to
   be copied from each search result sequence into the output
   sequence. In the output sequence, the field names will each be
   prefixed with `copy_<acc>_` where `<acc>` is the value of the *acc*
   field in the reference.
   
.. option:: --search-iupac=[pessimistic|*optimistic|exact] (optimistic)

   Configures how ambiguous bases are matched when computing the
   scores for the search results.

   **pessimistic**
    Ambiguous bases do not match anything because they *could* always
    be a mismatch.

   **optimistic**
     Ambiguous bases are considered matches if a match with the other
     (potentially also ambiguous base) is possible. That is, `N` will
     match everything, including `Y`.

   **exact**
     Matches on character level. `N` matches exactly `N`.
     
.. option:: --search-correction=[none|jc] (none)

   Apply distance correction to search result scores.

   **none**
     Leave score unmodified.

   **jc**
     Apply Jukes-Cantor correction.

.. option:: --search-cover=[abs|query|target|min|max|avg|overlap|all|nogap] (query)

   Compute sequence similarity as the fraction of the number of matches and

   **abs**
     the number 1: yields the absolute number of matching bases

   **query**
     the length of the query sequence. Yields the fraction of the
     query covered by the reference sequence.

   **target**
     the length of the target sequence. Yields the fraction of the
     result sequence covered by the query sequence.

   **min**
     the length of the shorter of the sequences compared.

   **max**
     the length of the longer of the sequences compared.

   **avg**
     the average length of the two sequences compared.

   **nogap**
     the number of columns in which both sequences have bases. Yields
     the equivalent of *matches / (matches+mismatches)*.

   **all**
     the number of columns in which either sequence has a
     bases. Similar to **nogap**, but does not ignore indel events.

   **overlap**
     the length of the overlapping portion of the two sequences.
                                
.. option:: --search-filter-lowercase

   Ignore lowercase bases when scoring result sequences. This can be
   used in conjunction with :option:`--lowercase`\=unaligned to ignore
   unaligned bases during the search and classification stage.


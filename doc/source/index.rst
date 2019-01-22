Reference based multiple sequence alignment using SINA
======================================================

SINA allows incorporating additional sequences into an existing
multiple sequence alignment (MSA) without modifying the original
alignment. While adding sequences to an MSA with SINA is usually
faster than re-computing the entire MSA from an augmented set of
unaligned sequences, the primary benefit lies in protecting
investments made into the original MSA such as manual curation of the
alignment, compute intensive phylogenetic tree reconstruction and
taxonomic annotation of the resulting phylogeny.

Additionally, SINA includes a homology search which uses the
previously computed alignment to determine the most similar
sequences. Based on the search results, a LCA based classification of
the query sequence can be computed using taxonomic classifications
assigned to the sequences comprising the reference MSA.

SINA is used to compute the small and large subunit ribosomal RNA
alignments provided by the SILVA_ project and is able to use the ARB_
format reference databases released by the project `here
<https://www.arb-silva.de/download/arb-files/>`_.

An `online version of SINA <https://www.arb-silva.de/aligner/>`_ is provided
by the SILVA_ project.

Publication
~~~~~~~~~~~

If you use SINA in your work, please cite:

Pruesse E, Peplies J, Glöckner FO. SINA: accurate high-throughput
multiple sequence alignment of ribosomal RNA
genes. *Bioinformatics*. 2012;28(14):1823-9.
doi::doi:`10.1093/bioinformatics/bts252`

.. _SILVA: https://www.arb-silva.de
.. _ARB: https://www.arb-home.de



.. toctree::
   :maxdepth: 2
   :hidden:

   install
   commandline
   CHANGELOG
   GitHub Repo<https://www.github.com/epruesse/SINA>
   Report Bug <https://www.github.com/epruesse/SINA/issues/new>

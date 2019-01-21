SINA - reference based multiple sequence alignment
==================================================

|latest| |release| |Bioconda| |TravisCI| |CircleCI| |Read the Docs| |Codecov|

.. |latest| image:: https://img.shields.io/github/release/epruesse/SINA/all.svg?label=latest
.. |release| image:: https://img.shields.io/github/release/epruesse/SINA.svg
.. |Bioconda| image:: https://img.shields.io/conda/vn/Bioconda/sina.svg
   :target: https://bioconda.github.io/recipes/sina/README.html
.. |TravisCI| image:: https://img.shields.io/travis/epruesse/SINA.svg?label=build%20(TravisCI)
   :target: https://travis-ci.org/epruesse/SINA
.. |CircleCI| image:: https://img.shields.io/circleci/project/github/epruesse/SINA.svg?label=build%20(CircleCI)
   :target: https://circleci.com/gh/epruesse/SINA
.. |Codecov| image:: https://img.shields.io/codecov/c/github/epruesse/sina.svg
   :target: https://codecov.io/gh/epruesse/SINA
.. |Read the Docs| image:: https://img.shields.io/readthedocs/sina/latest.svg

SINA is a tool to add sequences to an existing multiple sequence
alignment. It needs about 1 second on a single core to add one 16S
full length sequence (about 100k/h on a 32-core workstation). It was
developed to create the multi-million sequence alignment that is the
core of the SILVA SSU and LSU rRNA databases.

Installation
------------

- Use the `online <https://www.arb-silva.de/aligner>`_ version hosted
  by the SILVA project to align small batches of LSU and SSU
  sequences to their databses.
- The preferred way to install SINA locally is via Bioconda
  (`full instructions <https://sina.readthedocs.io/en/latest/install.html>`_)::

    conda create -n sina sina
    conda activate sina

Documentation
-------------

- The `SINA manual <https://sina.readthedocs.io>`_ is hosted at readthedocs.

- Please refer to the `publication in
  Bioinformatics <https://doi.org/10.1093/bioinformatics/bts252>`_ for
  a description of the algorithm.

  If you use SINA in your research, please don't forget to cite us:

  Elmar Pruesse, Jörg Peplies, Frank Oliver Glöckner; *SINA: Accurate
  high-throughput multiple sequence alignment of ribosomal RNA
  genes.* Bioinformatics 2012; 28 (14): 1823-1829.
  doi:10.1093/bioinformatics/bts252

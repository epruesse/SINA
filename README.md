# SINA - reference based multiple sequence alignment
[![GitHub (pre-)release](https://img.shields.io/github/release/epruesse/SINA/all.svg?label=latest)]()
[![GitHub release](https://img.shields.io/github/release/epruesse/SINA.svg)]()
[![Bioconda](https://img.shields.io/conda/vn/Bioconda/sina.svg)](https://bioconda.github.io/recipes/sina/README.html)
[![Build Status Travis](https://img.shields.io/travis/epruesse/SINA.svg?label=build%20(TravisCI))](https://travis-ci.org/epruesse/SINA)
[![Build Status CircleCI](https://img.shields.io/circleci/project/github/epruesse/SINA.svg?label=build%20(CircleCI))](https://circleci.com/gh/epruesse/SINA)

SINA is a tool to add sequences to an existing multiple sequence alignment. It needs about 1 second on a single core to add one 16S full length sequence (about 100k/h on a 32-core workstation). It was developed to create the multi-million sequence alignment that is the core of the SILVA SSU and LSU rRNA databases.

# Installation

 - Use the [online](https://www.arb-silva.de/aligner) version hosted by the SILVA project to align small batches of LSU and SSU sequences to their databses.
 - The preferred way to install SINA locally is via Bioconda ([full instructions](https://github.com/epruesse/SINA/wiki/Installation#using-bioconda))
   ```
   conda create -n sina sina 
   source activate sina
   ```
# Documentation

 - Please refer to the [publication in Bioinformatics](https://doi.org/10.1093/bioinformatics/bts252)
   for a description of the algorithm. 

   If you use SINA in your research, please don't forget to cite us:

   Elmar Pruesse, Jörg Peplies, Frank Oliver Glöckner; *SINA: Accurate high-throughput multiple
   sequence alignment of ribosomal RNA genes.* Bioinformatics 2012; 28 (14): 1823-1829. doi: 10.1093/bioinformatics/bts252
   
 - You can find more information in the [wiki](https://github.com/epruesse/SINA/wiki) 
 
 - or the online [manual](https://github.com/epruesse/SINA/blob/master/doc/man.md).

<!---
[pubmed](https://www.ncbi.nlm.nih.gov/pubmed/22556368)
[bioinformatics](https://academic.oup.com/bioinformatics/article-lookup/doi/10.1093/bioinformatics/bts252)
[pmc](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3389763/)

[![Altimetric Donut](http://api.altmetric.com/donut/727541_100x100.png)](https://www.altmetric.com/details/727541)
---> 

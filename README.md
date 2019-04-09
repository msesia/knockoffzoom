# KnockoffZoom

A powerful and versatile framework for the statistical analysis of genome-wide association studies.

## Overview

[Accompanying paper.](???)

This repository contains a pipeline for the analysis large GWAS data sets, especially designed for the [UK Biobank](https://www.ukbiobank.ac.uk/).
The code for the data analysis is provided in the form of Bash and R scripts, while the core algorithms are implemented in the R package [SNPknock](https://bitbucket.org/msesia/snpknock). 

The *KnockoffZoom* methodology consists of different modules, each corresponding to a subset of scripts, as summarized in the following flowchart.

![KnockoffZoom flowchart](misc/flowchart.png "KnockoffZoom flowchart")

## Software dependencies

The following software is required:

   - [PLINK 1.9](https://www.cog-genomics.org/plink/1.9/)
   - [PLINK 2.0](https://www.cog-genomics.org/plink/2.0/) alpha
   - [fastPHASE](http://scheet.org/software.html) 1.4.8
   - [GNU datamash](https://www.gnu.org/software/datamash/) 1.3
   - [GNU awk](https://github.com/onetrueawk/awk) 4.0.2

The [R](https://www.r-project.org/) (version 3.5.1) code was tested on the following configuration:

   - [SNPknock](https://bitbucket.org/msesia/snpknock) 0.8.0
   - [adjclust](https://CRAN.R-project.org/package=adjclust ) 0.5.6
   - [bigsnpr](https://privefl.github.io/bigsnpr/) 0.9.1
   - [bigstatsr](https://privefl.github.io/bigstatsr/) 0.8.4
   - [snpStats](https://doi.org/doi:10.18129/B9.bioc.snpStats) 1.32
   - [Matrix](https://CRAN.R-project.org/package=Matrix ) 1.2.15
   - [data.table](https://CRAN.R-project.org/package=data.table) 1.12.0
   - [tidyverse](https://www.tidyverse.org/) 1.2.1
   - [devtools](https://CRAN.R-project.org/package=devtools) 1.13.6

## Authors

   - [Matteo Sesia](http://web.stanford.edu/~msesia/) [(Stanford University)](https://www.stanford.edu/)

See also the list of [contributors](misc/contributors.txt) who participated in this project.

## Acknowledgments

- [Christopher Chang](https://github.com/chrchang), for developing [PLINK](https://www.cog-genomics.org/plink) and offering valuable technical support.
- [Florian Privé](https://privefl.github.io/), for developing [bigsnpr](https://privefl.github.io/bigsnpr/) and promptly replying to related inquiries.

## References

This methodology was developed within the broader scope of [knockoffs](https://web.stanford.edu/group/candes/knockoffs/).

## License

This software is distributed under the [GPLv3 license](https://www.gnu.org/licenses/gpl-3.0.en.html) and it comes with ABSOLUTELY NO WARRANTY.

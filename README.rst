What is ClimDown?
=================

"ClimDown" is a Climate Downscaling package for the `R statistical
programming language`. It was written at the `Pacific Climate Impacts
Consortium`_ (PCIC) with support from `Environment and Climate Change
Canada`_.

The package provides routines for statistical downscaling of coarse
scale global climate model (GCM) output to a fine spatial resolution.

PCIC's suite of routines include several different (yet related)
downscaling techniques:

* Bias-Corrected Spatial Downscaling (BCSD)
* Bias-Correction Constructed Analogues (BCCA)
* Bias Correction/Climate Imprint (BCCI)
* Bias Correction/Constructed Analogues with Quantile mapping reordering (BCCAQ).

See the corresponding scientific literature for more details on the
pros and cons of each downscaling method.

  .. _R statistical programming language: http://www.r-project.org/
  .. _Pacific Climate Impacts Consortium: https://pacificclimate.org/
  .. _Environment and Climate Change Canada: http://ec.gc.ca/


Installation
------------

You can install the latest `ClimDown release from CRAN`_ using the R
interpreter: ``> install.packages('climdex.pcic')``

.. _ClimDown release from CRAN: http://cran.r-project.org/web/packages/ClimDown/index.html

If you are interested in a development version or a specific release
of ClimDown, you can use the `devtools` package as an installation
alternative.

     > install.packages('devtools')
     > devtools::install_github("pacificclimate/ClimDown", ref="release")
     # Or
     > devtools::install_github("pacificclimate/ClimDown", ref="1.0.1")

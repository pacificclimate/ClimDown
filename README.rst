What is ClimDown?
=================

"ClimDown" is a Climate Downscaling package for the `R statistical
programming language`. It was written at the `Pacific Climate Impacts
Consortium`_ (PCIC) with support from `Environment and Climate Change
Canada`_.

The package provides routines for statistical downscaling of coarse
scale global climate model (GCM) output to a fine spatial resolution.

PCIC's suite of routines include several different (yet related)
downscaling techniques. The entire process is named Bias
Correction/Constructed Analogues with Quantile mapping reordering
(BCCAQ) and is composed of the following steps.

* Constructed Analogues (CA)
* Climate Imprint (CI)
* Quantile Delta Mapping (QDM)
* Rerank

See refer to the package documentation for details on each step and
references to the corresponding scientific literature.

  .. _R statistical programming language: http://www.r-project.org/
  .. _Pacific Climate Impacts Consortium: https://pacificclimate.org/
  .. _Environment and Climate Change Canada: http://ec.gc.ca/

Climate Downscaling: What and Why?
==================================

Changes in global climate have widespread impacts on the environment,
economic activity, and human health, especially in high latitudes
where warming is proceeding more rapidly and where ecosystems and
traditional lifestyles are particularly sensitive to the impacts of
warming.

Planning for adapting to climate change requires scientifically sound
information about the future climate. Global climate models (GCMs)
simulate future climate under different emission scenarios. However,
GCMs simulate average conditions over large grid cells--typically on
the order of 10,000 square kilometers or more per cell--which is often
too coarse a resolution for regional and local applications. The use
of original GCM data is not always the best option to provide
adaptation-relevant information at the local scale.

Bias in model simulated local climate is of concern for many
applications. For example, compared with observations, the median
temperature simulated by GCMs from the Coupled Model Intercomparison
Project Phase 5 (CMIP5) shows biases relative to Climate Research Unit
high-resolution gridded dataset (`CRU TS3.10`_) ranging from -3° C to
1.5° C for seasonal and annual mean temperatures in 26 global land
areas (`Flato et al. 2013`_).  Precipitation simulated by CMIP5 models
is also biased relative to observations (`Flato et al. 2013`_). These
biases hinder the direct application of model simulated future climate
for impacts modelling and adaptation planning since climate impacts
are often related to certain physical or biophysical thresholds. As a
result, adaptation planning often uses model simulated future climate
information that has incorporated some sort of downscaling and bias
correction. Additionally, climate information is more useable and is
less prone to misinterpretation when presented in a manner specific to
the local community and/or impacts most relevant to a particular
sector. High-resolution future projections of impact-relevant climate
indices can be particularly useful in this regard.

ClimDown has been used to produce such bias corrected, downscaled GCMs
for current and future climates, and could be used to do so for
anywhere else in the world.

.. _Flato et al. 2013: http://www.ipcc.ch/pdf/assessment-report/ar5/wg1/WG1AR5_Chapter09_FINAL.pdf
.. _CRU TS3.10: http://dx.doi.org/10.1002/joc.3711

Installation
============

You can install the latest `ClimDown release from CRAN`_ using the R
interpreter: ``> install.packages('climdex.pcic')``

.. _ClimDown release from CRAN: http://cran.r-project.org/web/packages/ClimDown/index.html

If you are interested in a development version or a specific release
of ClimDown, you can use the `devtools` package as an installation
alternative.::

    > install.packages('devtools')
    > devtools::install_github("pacificclimate/ClimDown", ref="release")
    # Or
    > devtools::install_github("pacificclimate/ClimDown", ref="1.0.1")

System dependencies
-------------------

ClimDown reads all of its input and produces its output in `NetCDF
format`_ and manages numeric units using the `UDUNITS2 library`_. The
NetCDF and udunits2 libraries are system dependency of the
package. Ensure that the following packages are installed under
Debian/Ubuntu Linux systems: libnetcdf-dev, netcdf-bin, and
libudunits2-dev. For other systems, follow the installation `NetCDF
install instructions`_ and `UDUNITS2 install instructions`_ provided
by `Unidata`_.

.. _NetCDF format: https://www.unidata.ucar.edu/software/netcdf/docs/netcdf_introduction.html
.. _UDUNITS2 library: https://www.unidata.ucar.edu/software/udunits/udunits-current/doc/udunits/udunits2.html
.. _NetCDF install instructions: https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html
.. _UDUNITS2 install instructions: https://www.unidata.ucar.edu/software/udunits/udunits-current/doc/udunits/udunits2.html#Installation
.. _Unidata: https://www.unidata.ucar.edu/

Necessary Resources, Performance, and Platform
==============================================

The BCCAQ algorithm implemented by ClimDown is a complex, multi-stage
operation, and as such performance will vary widely depending on the
size of the input, the degree of parallelism selected by the user, and
the performance characteristics of user's system (CPU speed, available
RAM, I/O speed).

We have `previously written`_ about some of the `complexities involved
in downscaling performance`_, but it remains an area of active study.

Consider our experience as a matter of anecdote. We typically run
ClimDown for downscaling 150 year, daily GCM simulations to a
Canada-wide ANUSPLIN grid (approximately 10km resolution, 1068 by 510
cells). On our Linux systems, such runs can take up to 7 days to
complete. However each of the different downscaling steps has
different opportunities for parallelism and different performance
characteristics. Typical of our runs is something like this:

* CI: 1 core, 10 GB RAM, Run time ~ 7 hours
* CA: 8 cores, 10 GB RAM, Run time ~ 2 hours
* QDM: 1 core, 36 GB RAM, Run time ~ 1.5 days
* rerank: 4 cores, 8 GB RAM, Run time ~ 4 days

In general, this downscaling technique is very expensive for large
spatiotemporal domains. The more you can limit your domain, the faster
your runtime will be. For small domains, it may be possible to run
ClimDown on a typical workstation, but in general we do all of our
production runs on rack-mounted supercomputers.

Though Windows binaries for ClimDown are available `from CRAN`_, no
effort has been made to optimize this package for Windows and your
mileage may vary.

.. _previously written: http://james.hiebert.name/blog/work/2016/04/26/BCCA/
.. _complexities involved in downscaling performance: https://github.com/pacificclimate/ClimDown/blob/doc/doc/report.md#rewriting-numerous-algorithms
.. _from CRAN: https://cran.r-project.org/web/packages/ClimDown/index.html

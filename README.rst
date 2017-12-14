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

---
title: 'ClimDown: Climate Downscaling in R'
tags:
  - climate
  - downscaling
  - spatiotemporal
  - R
authors:
  - name: James Hiebert
    orcid: 0000-0002-4171-9586
    affiliation: 1
  - name: Alex Cannon
    affiliation: 2
  - name: Trevor Murdock
    affiliation: 1
  - name: Stephen Sobie
    affiliation: 1
  - name: Arelia Werner
    affiliation: 1
affiliations:
  - name: Pacific Climate Impacts Consortium
    index: 1
  - name: Environment and Climate Change Canada
    index: 2
date: 7 June 2017
---

# Summary

The ClimDown R package publishes the routines and techniques of the
[Pacific Climate Impacts Consortium](https://pacificclimate.org/)
(PCIC) for downscaling coarse scale Global Climate Models (GCMs) to
fine scale spatial resolution.

PCIC's overall downscaling algorithm is named Bias-corrected
constructed analogues with quantile mapping (BCCAQ). BCCAQ is a hybrid
downscaling method that combines outputs from Climate Analogues (CA)
and quantile mapping at the fine-scale resolution.  First, the CA and
Climate Imprint (CI) plus quantile delta mapping (QDM) algorithms are
run independently. BCCAQ then combines outputs from the two by taking
the daily QDM outputs at each fine-scale grid point and reordering
them within a given month according to the daily CA ranks, i.e., using
a form of Empirical Copula Coupling.

The package exports high-level wrapper functions that perform each of
three downscaling steps: CI, CA, and QDM, as well as one wrapper that
runs the entire BCCAQ pipeline.

# References

Cannon, A. J., Sobie, S. R., & Murdock, T. Q. (2015). Bias Correction of GCM Precipitation by Quantile Mapping: How Well Do Methods Preserve Changes in Quantiles and Extremes?. Journal of Climate, 28(17), 6938-6959. doi: 10.1175/JCLI-D-14-00754.1

Hunter, R. D., & Meentemeyer, R. K. (2005). Climatologically aided mapping of daily precipitation and temperature. Journal of Applied Meteorology, 44(10), 1501-1510.

Maurer, E. P., Hidalgo, H. G., Das, T., Dettinger, M. D., & Cayan, D. R. (2010). The utility of daily large-scale climate data in the assessment of climate change impacts on daily streamflow in California. Hydrology and Earth System Sciences, 14(6), 1125-1138.

Schefzik, R., Thorarinsdottir, T. L., & Gneiting, T. (2013). Uncertainty quantification in complex simulation models using ensemble copula coupling. Statistical Science, 28(4), 616-640.

Werner, A. T., & Cannon, A. J. (2016). Hydrologic extremes - an intercomparison of multiple gridded statistical downscaling methods. Hydrology and Earth System Sciences, 20(4), 1483-1508. doi: 10.5194/hess-20-1483-2016

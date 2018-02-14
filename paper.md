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
date: 1 February 2018
bibliography: paper.bib
---

# Summary

The ClimDown R package publishes the routines and techniques of the
[Pacific Climate Impacts Consortium](https://pacificclimate.org/)
(PCIC) for downscaling coarse scale Global Climate Models (GCMs) to
fine scale spatial resolution.

PCIC's overall downscaling algorithm is named Bias-corrected
constructed analogues with quantile mapping (BCCAQ)
[@cannon15; @werner16]. BCCAQ is a hybrid downscaling method that
combines outputs from Climate Analogues (CA) [@maurer10] and quantile
mapping at the fine-scale resolution.  First, the CA and Climate
Imprint (CI) [@hunter05] plus quantile delta mapping (QDM) [@cannon15]
algorithms are run independently. BCCAQ then combines outputs from the
two by taking the daily QDM outputs at each fine-scale grid point and
reordering them within a given month according to the daily CA ranks,
i.e., using a form of Empirical Copula Coupling [@schefzik13].

The package exports high-level wrapper functions that perform each of
three downscaling steps: CI, CA, and QDM, as well as one wrapper that
runs the entire BCCAQ pipeline.

# References

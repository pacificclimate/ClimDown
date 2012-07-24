#!/bin/bash

##******************************************************************************
# Bias Corrected Spatial Disaggregation (BCSD) downscaling algorithm
# Alex Cannon (acannon@uvic.ca)
##******************************************************************************

Rscript BCSD.1.aggregate_obs.R $1
Rscript BCSD.2.cdo_extract_clim.R $1
Rscript BCSD.3.bias_correct_dqm.R $1
Rscript BCSD.4.resample_rescale.R $1

##******************************************************************************

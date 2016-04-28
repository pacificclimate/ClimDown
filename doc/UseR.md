High performance climate downscaling in R

Global Climate Models (GCMs) can be used to assess the impacts of future climate change on particular regions of interest, municipalities or pieces of infrastructure. However, the coarse spatial scale of GCM grids (50km or more) can be problematic to engineers or others interested in more localized conditions. This is particularly true in areas with high topographic relief and/or substantial climate heterogeneity. A technique in climate statistics known as "downscaling" exists to map coarse scale climate quantities to finer scales.

Two main challenges are posed to potential climate downscalers: there exist few open source implementations of proper downscaling methods and many downscaling methods necessarily require information across both time and space. This requires high computational complexity to implement and substantial computational resources to execute.

The Pacific Climate Impacts Consortium in Victoria, BC, Canada has written a high-performance implementation of the detrended quantile mapping (QDM) method for downscaling climate variables. QDM is a proven climate downscaling technique that has been shown to preserve the relative changes in both the climate means and the extremes. We will release this software named "ClimDown" to CRAN under an open source license. Using ClimDown, researchers can downscale spatially coarse global climate models to an arbitrarily fine resolution for which gridded observations exist.

Our proof-of-concept for ClimDown has been to downscale all of Canada at 10km resolution for models and scenarios from the IPCC's Coupled Model Intercomparison Project. We will present our results and performance metrics from this exercize.

[1] http://journals.ametsoc.org/doi/abs/10.1175/JCLI-D-14-00754.1

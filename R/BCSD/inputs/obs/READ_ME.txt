forcings_symap_6_BC_29AUG2011_1950-2006
-this uses the same symap algorithm as below, but uses 6 nearest stations in the algorithm instead of 4.

forcings_new_symap_BC_22AUG2011_1950-2006
-this was run the same as the one below this, but with a re-write of symap in R, that implements symap as documented (REF).

forcings_upmet_BC_01JUN2011_1950-2006
-this dataset includes the VIC forcings for BC run including BCHydro stations and the original PRISM created just for BC, not run through ClimateWNA
-the south of BC portion of Washington is not included
-this should be pre-any changes to the code, it should be a functioning copy with no -99
-this should be distributed over the version below when people are interested in gridded T and P over the province, except it should be quality
controlled further first. This will be addressed in our validation plan.

forcings_1950-2006_10JUN2010
-this dataset includes BCH stations, EC stations
-this is the upmet version used to calibrate the Peace, Campbell and Columbia basins
-there may be some -99 in the southwest corner of the province
-this will include a portion that is south of the BC border in Washington to about 47 degrees N
-there were problems with the DEM in this southern portion. It was aggregated from 3 arc secs to 0.0625 degrees incorrectly, in a way that can't be reproduced.
-ClimateWNA was forced with this 'bad' DEM (SRTM v3) to get the PRISM Tmin, Tmax, Prec monthly climatologies and these are therefore not in good shape in Washington.
-those interested in replicating our hydrologic modelling work might be interested in this version over the one above.




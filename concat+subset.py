#!/usr/bin/python

import os

center = 'BNU'
gcm = 'BNU-ESM'
gcmdir = '../models/' + center + '/' + gcm + '/'
writedir = '../output/' + gcm + '/'
enss = ['r1i1p1']
scens = ['rcp26', 'rcp45', 'rcp85'] # , 'rcp26']
vars = ['pr', 'tasmin', 'tasmax']
h1 = '19500101'
h2 = '20051231'
f1 = '20060101'
f2 = '21001231'
s1 = '19500101'
sel1 = '1950-01-01T00:00'
sel2 = '2100-12-31T23:59'
for var in vars:
    for scen in scens:
        for ens in enss:
            vers1 = ''.join(os.listdir(gcmdir + 'historical/day/atmos/day/' + ens + '/'))
            vers2 = ''.join(os.listdir(gcmdir + scen + '/day/atmos/day/' + ens + '/'))
            #print vers1
            #print vers2
            concat = 'ionice -c3 cdo cat ' + \
                     gcmdir + 'historical/day/atmos/day/' + ens + '/' + vers1 + '/' + var + '/' + \
                     var + '_day_' + gcm + '_historical_' + ens + '_' + h1 + '-' + h2 + '.nc ' + \
                     gcmdir + scen + '/day/atmos/day/' + ens + '/' + vers2 +'/' + var + '/' + \
                     var + '_day_' + gcm + '_' + scen + '_' + ens + '_' + f1 + '-' + f2 + '.nc ' + \
                     writedir + var + '_day_' + gcm + '_historical+' + scen + '_' + ens + '_' + h1 + '-' + f2 + '.nc'
            comb1 = writedir +  var + '_day_' + gcm + '_historical+' + scen + '_' + ens + '_' + h1 + '-' + f2 + '.nc '
            comb2 = writedir +  var + '_day_' + gcm + '_historical+' + scen + '_' + ens + '_' + s1 + '-' + f2 + '.nc '
            short = 'ionice -c3 cdo seldate,' + sel1 + ',' + sel2 + ' ' + \
                    writedir +  var + '_day_' + gcm + '_historical+' + scen + '_' + ens + '_' + h1 + '-' + f2 + '.nc ' + \
                    writedir +  var + '_day_' + gcm + '_historical+' + scen + '_' + ens + '_' + s1 + '-' + f2 + '.nc'
            #print concat
            #print short
            os.system(concat)
            if comb1!=comb2:
                os.system(short)           
            
            

os.system('ionice -c3 ./global2local.py ' + writedir)



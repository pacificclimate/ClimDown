#!/usr/bin/python

import os

gcm = 'HadGEM2-CC'
gcmdir = 'models/' + gcm + '/'
enss = ['r1i1p1']
scens = ['rcp45', 'rcp85']
vars = ['pr', 'tasmin', 'tasmax']
h1 = '18600101'
h2 = '20051230'
f1 = '20060101'
f2 = '20991230'
s1 = '19500101'
sel1 = '1950-01-01T00:00'
sel2 = '2100-12-30T23:59'

for var in vars:
    for scen in scens:
        for ens in enss:
            concat = 'ionice -c3 cdo cat ' + gcmdir + 'historical/day/atmos/' + var + \
                     '/' + ens + '/' + var + '_day_' + gcm + \
                     '_historical_' + ens + '_' + h1 + '-' + h2 + '.nc ' + \
                     gcmdir + scen + '/day/atmos/' + var + '/' + ens + '/' + var + \
                     '_day_' + gcm + '_' + scen + '_' + ens + '_' + f1 + '-' + \
                     f2 + '.nc ' + var + '_day_' + gcm + '_historical+' + scen + \
                     '_' + ens + '_' + h1 + '-' + f2 + '.nc'
            short = 'ionice -c3 cdo seldate,' + sel1 + ',' + sel2 + ' ' + var + '_day_' + gcm + \
                    '_historical+' + scen + '_' + ens + '_' + h1 + '-' + f2 + \
                    '.nc ' + var + '_day_' + gcm + '_historical+' + scen + \
                    '_' + ens + '_' + s1 + '-' + f2 + '.nc'
            print concat
            print short
            os.system(concat)
            os.system(short)

os.system('ionice -c3 ./global2local.py')
os.system('mv subset.* subset/')
os.system('rm *.nc')
os.system('mv subset/subset.* .')


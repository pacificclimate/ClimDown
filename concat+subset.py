#!/usr/bin/python

import os

center = 'BCC'
gcm = 'bcc-csm1-1-m'
gcmdir = '../models/' + center + '/' + gcm + '/'
##writedir = '../output/' + gcm + '/'
writedir = '../output/' + gcm + '/subset/'
enss = ['r1i1p1']
scens = ['rcp85'] #, 'rcp45', 'rcp85'] # , 'rcp26']
vars = ['tasmin', 'tasmax']
h1 = '18500101'
h2 = '20121231'
f1 = '20060101'
f2 = '20991231'

##Fixed File descripters
hsub1 = '19500101'
hsub2 = '20051231'
fsub1 = '20060101'
fsub2 = '21001231'

##Fixed values to create a common subset
w1 = '19500101'
w2 = '21001231'
sel1 = '1950-01-01T00:00'
sel2 = '2100-12-31T23:59'

for var in vars:
    for scen in scens:
        for ens in enss:
            vers1 = ''.join(os.listdir(gcmdir + 'historical/day/atmos/day/' + ens + '/'))
            vers2 = 'v20120719'##''.join(os.listdir(gcmdir + scen + '/day/atmos/day/' + ens + '/'))
            print vers1
            print vers2
            ##Fix overlapping dates
            hist = 'ionice -c3 cdo seldate,1950-01-01T00:00,2005-12-31T23:59 ' + \
                   gcmdir + 'historical/day/atmos/day/' + ens + '/' + vers1 + '/' + var + '/' + \
                   var + '_day_' + gcm + '_historical_' + ens + '_' + h1 + '-' + h2 + '.nc ' + \
                   writedir + var + '_day_' + gcm + '_historical_' + ens + '_' + hsub1 + '-' + hsub2 + '.nc'
            os.system(hist)
            proj = 'ionice -c3 cdo seldate,2006-01-01T00:00,2100-12-31T23:59 ' + \
                gcmdir + scen + '/day/atmos/day/' + ens + '/' + vers2 + '/' + var + '/' + \
                var + '_day_' + gcm + '_' + scen + '_' + ens + '_' + f1 + '-' + f2 + '.nc ' + \
                writedir +  var + '_day_' + gcm + '_' + scen + '_' + ens + '_' + fsub1 + '-' + fsub2 + '.nc' 
            print proj
            os.system(proj)
            concat = 'ionice -c3 cdo cat ' + \
                  writedir +  var + '_day_' + gcm + '_historical_' + ens + '_' + hsub1 + '-' + hsub2 + '.nc ' + \
                  writedir +  var + '_day_' + gcm + '_' + scen + '_' + ens + '_' + fsub1 + '-' + fsub2 + '.nc ' + \
                  writedir +  var + '_day_' + gcm + '_historical+' + scen + '_' + ens + '_' + w1 + '-' + w2 + '.nc'
            
            #concat = 'ionice -c3 cdo cat ' + \
            #         gcmdir + 'historical/day/atmos/day/' + ens + '/' + vers1 + '/' + var + '/' + \
            #         var + '_day_' + gcm + '_historical_' + ens + '_' + h1 + '-' + h2 + '.nc ' + \
            #         gcmdir + scen + '/day/atmos/day/' + ens + '/' + vers2 +'/' + var + '/' + \
             #        var + '_day_' + gcm + '_' + scen + '_' + ens + '_' + f1 + '-' + f2 + '.nc ' + \
            #         writedir + var + '_day_' + gcm + '_historical+' + scen + '_' + ens + '_' + h1 + '-' + f2 + '.nc'
            print concat
            #comb1 = writedir +  var + '_day_' + gcm + '_historical+' + scen + '_' + ens + '_' + h1 + '-' + f2 + '.nc '
            #comb2 = writedir +  var + '_day_' + gcm + '_historical+' + scen + '_' + ens + '_' + s1 + '-' + s2 + '.nc '
            #short = 'ionice -c3 cdo seldate,' + sel1 + ',' + sel2 + ' ' + \
            #        writedir +  var + '_day_' + gcm + '_historical+' + scen + '_' + ens + '_' + h1 + '-' + f2 + '.nc ' + \
            #        writedir +  var + '_day_' + gcm + '_historical+' + scen + '_' + ens + '_' + w1 + '-' + w2 + '.nc'
            #print concat
            #print short
            os.system(concat)
            #if comb1!=comb2:
            #    os.system(short)           
            
            

os.system('ionice -c3 ./global2local.py ' + writedir)



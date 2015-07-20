#!/usr/bin/python

import os, glob, sys


writedir = str(sys.argv[1])
#files = glob.glob(writedir + '*.nc')
files = [os.path.basename(x) for x in glob.glob(writedir + '*.nc')]

files.sort()
print files
for file in files:
    file_out = 'subset.' + file
    print(file_out)
    os.system('ncks -d lon,180.,320. -d lat,5.,85. ' + writedir + file + ' ' + writedir + file_out)
    os.system('mv ' + writedir + 'subset.* ' + writedir + 'subset/')
    ##os.system('rm ' + writedir + '*.nc')
    ##os.system('mv ' + writedir + 'subset/subset.* ' + writedir + '.')

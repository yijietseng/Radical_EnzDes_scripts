from time import time
from sys import argv
start = time()

from sys import argv
from time import sleep, time
import os,json
from designClass_v4 import WT_validation as wt
import numpy as np
import pandas as pd



# load parameter file
with open(argv[1]) as jsonfile:
	xParams = json.load(jsonfile)

# input a file_increment, or it will grab the folder with the largest increment to start with
if len(argv) > 2:
	file_incre = "{:02d}".format(int(argv[2]))
else:
	for i in itertools.count(start=1):
		if not os.path.exists('../'+xParams['PDBID']+'/CLUSTERX_'+"{:02d}".format(i)+'/'):
			file_incre = "{:02d}".format(i-1)
			break

# making necessary folders
try:
	os.mkdir('../'+xParams['PDBID']+'/CLUSTERW_'+file_incre+'/outputstruct')
	os.mkdir('../'+xParams['PDBID']+'/CLUSTERW_'+file_incre+'/logs')

	print('\n\n\nFolders generated\n\n\n')
except:
	print('\n\n\nFolders exist\n\n\n')




# Reading data frame
pdbFN = '../'+xParams['PDBID']+'/CLUSTERX_W'+file_incre+'/WT/'+xParams['PDBID']+'_X.pdb'


proc = wt(**xParams)
proc.runProcesses()


end = time()

d = end - start

if d >= 3600:
	print('Total run time was', d/3600, 'hr')
if 60 <= d < 3600:
	print('Total run time was', d/60, 'min')
if d < 60:
	print('Total run time was', d, 'sec')




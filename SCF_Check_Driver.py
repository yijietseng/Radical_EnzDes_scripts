from designClass_v4 import enzymeRot
from sys import argv
import os, json, itertools, glob
import pandas as pd

def runDummy(rotNo,file_incre):
	#rot = int(PDBFILE.rsplit('/',1)[1][7:-6])
	rot1 = enzymeRot(rotamer=rotNo,designStage='ClashChecked/',designSuffix='_C.pdb',fileIncre=file_incre,dRes=True,**xParams)
	rot1.runDummyDesign()

# Input params
with open(argv[1]) as jsonfile:
	xParams = json.load(jsonfile)

# input a file_increment, or it will grab the folder with the largest increment to start with
for i in itertools.count(start=1):
	if not os.path.exists('../'+xParams['PDBID']+'/CLUSTERX_'+"{:02d}".format(i)+'/'):
		file_incre = "{:02d}".format(i-1)
		break

#PDBfile = argv[2]

rotNo = int(argv[2])

runDummy(rotNo, file_incre)

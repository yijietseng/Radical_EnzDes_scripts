from designClass_v4 import enzymeRot
from sys import argv
from time import sleep, time
import numpy as np
import pandas as pd
import json, glob, os, itertools


start = time()

# Input params
with open(argv[1]) as jsonfile:
	xParams = json.load(jsonfile)

# load in the passed log file
df = pd.read_csv(argv[2],delimiter='\t')

designStep = argv[3] # input hb for hbDesign or hp for hpDesign

rotNo = int(argv[4]) # indicate rotamer number for hbDesign or the row number for hpDesign

# input a file_increment, or it will grab the folder with the largest increment to start with
if len(argv) > 5:
	file_incre = "{:02d}".format(int(argv[5]))
else:
	for i in itertools.count(start=1):
		if not os.path.exists('../'+xParams['PDBID']+'/CLUSTERX_'+"{:02d}".format(i)+'/'):
			file_incre = "{:02d}".format(i-1)
			break

#PDB = glob.glob('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/'+xParams["hbSaveFile"]+'/*.pdb')

def runHBDesign(i):
	''' A function to run the hydrogen bonding design, will eventually turn into a function inside the design class, in here to make sure that it properly coordinates with the multithreading module for now.'''
	rot1 = enzymeRot(rotamer=i,designStage='ClashChecked/',designSuffix='_C.pdb',fileIncre=file_incre,dRes=True,**xParams)
	rot1.runDesign()

def runHPDesign(i):
	suffix = df.iloc[i]['pdbFN'][-9:]
	rot = int(df.iloc[i]['pdbFN'][-12:-9])
	rot1 = enzymeRot(rotamer=df.iloc[i]['Rotamer'],designStage='hbDesigned/',designSuffix=suffix,fileIncre=file_incre,dRes=True,**xParams)
	rot1.runDesign()

#------------------------------------------------------------------------------------------------------------------------------------------


if designStep == 'hb':
	# making necessary folders
	try:
		os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hbDesigned/')
		os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hbDesigned/logs')

		print('\n\n\nHB folders generated\n\n\n')
	except:
		print('\n\n\nHB folders exist\n\n\n')
	

	runHBDesign(rotNo)
	
	
	

if designStep == 'hp':
	# making necessary folders
	try:
		os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hpDesigned/')
		os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hpDesigned/logs')

		print('\n\n\nHP folders generated\n\n\n')
	except:
		print('\n\n\nHP folders exist\n\n\n')
	

	runHPDesign(rotNo)
	

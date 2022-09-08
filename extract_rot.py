from time import time, sleep
from sys import argv
import os, json, glob
import numpy as np
import pandas as pd

def add_rot(PDBID, file_incre):
	logLS = glob.glob('../'+PDBID+'/CLUSTERX_'+file_incre+'/Scaffold_Check/logs/*.txt')

	for i in logLS:
		df = pd.read_csv(i,delimiter='\t')
		rotNo = i[-7:-4]
		df['Rot'] = rotNo
		df.to_csv(i,sep='\t',index=False)

def extract_rot(PDBID, df, rep=10, SASA=30):
	rotLS = []
	repFilter = df['Gly_rep'] <= rep
	SASAFilter = df['Des_SASA'] <= SASA
	filtered = repFilter & SASAFilter
	passedDF = df[filtered]
	passedDF.to_csv('../'+PDBID+'/RotList2Des.txt',sep='\t',index=False)
	
def process(dir_list, file_incre):
	
	if not os.path.exists('../SCF_master_result/'):
		os.mkdir('../SCF_master_result/')
	
	for i in dir_list:
		# Add rotamer numbers to their corresponding data
		add_rot(i, file_incre)
		
		# Start combining into a master file
		logFolder = '../'+i+'/CLUSTERX_'+file_incre+'/Scaffold_Check/logs/*.txt'
		saveFN = '../SCF_master_result/'+i+'_master.txt'
		# Make a glob that has all the log files
		logFiles = glob.glob(logFolder)

		# Initialize a dataframe to hold all the information
		df = pd.DataFrame()

		for fn in logFiles:
			logDF = pd.read_csv(fn,delimiter='\t')
			df = df.append(logDF)
		df.reset_index(inplace=True)
		df.to_csv(saveFN,sep='\t',index=False)

		# filter through the rotamer that satisfied the filter
		extract_rot(i,df)

def remove_extra(dir_list):
	for i in dir_list:
		# Building passed rotmater list
		df_rot = pd.read_csv('../'+i+'/RotList2Des.txt')
		rotls = []
		for k in df_rot.Rot:
			rotls.append("{:03d}".format(k))
		# Deleting extra rotamers 
		for j in range(1,9):
			incre = "{:02d}".format(j)
			rot_new = []
			for y in rotls:
				rot_new.append('../'+i+'/CLUSTERX_'+incre+'/init_PS/'+i+'_PS'+y+'.pdb')
			rot_full = glob.glob('../'+i+'/CLUSTERX_'+incre+'/init_PS/*.pdb')
			for o in rot_full:
				if o not in rot_new:
					os.remove(o)

			
# ls the directory
alist = os.listdir('../') 
dir_list = []

# Pull the scaffold directories into a list
for i in alist:
	if len(i) == 4:
		dir_list.append(i) 
for f in dir_list:
	if os.path.isfile == True:
		dir_list.remove(f)
dir_list.sort()

# Identify file_incre
if len(argv) > 1:
		file_incre = "{:02d}".format(int(argv[1]))
else:
	file_incre = "01"

process(dir_list, file_incre)

remove_extra(dir_list)``

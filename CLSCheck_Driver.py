'''
*********************************************************************************************************************
This script will carry out the clash check of the initial structures. It will first mutate the residues into gly 
according to the resfile. It will then do a fixed-backbone relax to optimize the score; then fa_rep will be measured.
The gly poses will be saved and a log file will be outputted as well. The script will do filter the structrues based
on the fa_rep and fa_intra_rep, and it will output a master txt file of the passed structrues.
******You can instruct a specific file increment or let the program figure out the most recent file increment.******
*********************************************************************************************************************
'''

from time import time, sleep
from sys import argv
from designClass_v4 import enzymeRot
from multiprocessing import Pool
import os, json, glob, itertools
import numpy as np
import pandas as pd
start = time()

# loading parameter file
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
	# Getting the most recent design folder 
	os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/ClashChecked/')
	os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/ClashChecked/logs')
	os.mkdir('../'+xParams['PDBID']+'/passedCLS')
	print('\n\n\nFolders generated\n\n\n')
except:
	print('\n\n\nFolders exist\n\n\n')

pdbList = glob.glob('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/init_PS/*.pdb')
constraint = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/RC/'+xParams['PDBID']+'.cst'



#-------------------------------------------------------Below is a section for functions----------------------------------------------------
def runClashCheck(i):
	rot_No = int(pdbList[i][-7:-4])
	rot1 = enzymeRot(rotamer=rot_No, fileIncre=file_incre, glyRes = True,**xParams)
	rot1.gen_gly_pose()
	print('\n\n\nStarted on rotamer: ',pdbList[i-1][-14:],'\n\n\n')

def runClashFilter():
	# Start filtering  
	print('\n\nInitiating filtering process....')
	logFolder = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/ClashChecked/logs/*.txt'
	saveFN = '../'+xParams['PDBID']+'/passedCLS/'+xParams['PDBID']+'_CLUSTERX_'+file_incre+'_passedCLS.txt'
	# Make a glob that has all the log files
	logFiles = glob.glob(logFolder)

	# Initialize a dataframe to hold all the information
	df = pd.DataFrame()

	for fn in logFiles:
		logDF = pd.read_csv(fn,delimiter='\t')
		df = df.append(logDF)
	df.reset_index(inplace=True)
	#print(df.keys())
	# filter fa_rep and intra_fa_rep
	clashFilter = df['resFaRep'] <= xParams['clashMax']
	intraFilter = df['resFaIntraRep'] <= xParams['intraMax']

	combinedFilter = clashFilter & intraFilter
	combinedFilter = intraFilter
	goodCLS = df[combinedFilter]
	goodCLS.reset_index(inplace=True)
	goodCLS.to_csv(saveFN,sep='\t')

	print('\n',len(goodCLS),'structures passed clash check')

	print('\n Filtering process completed!!!! Clash check completed!!!!')

#--------------------------------------------Below section is for multitreading and runing functions--------------------------------------------

#myPool = Pool(processes=xParams['processes'])

myPool = Pool(processes=30)



# Initiate multithreading
for i in myPool.imap_unordered(runClashCheck,range(len(pdbList))):
	print('\nInitiating Clash Check....')


#runClashCheck(1) # This line is for single structure testing purpose


# Filter the structures that passed the fa_rep and fa_intra_rep filters
runClashFilter()


end = time()
d = end - start

if d >= 3600:
	print('\nTotal run time was', d/3600, 'hr')
if 60 <= d < 3600:
	print('\nTotal run time was', d/60, 'min')
if d < 60:
	print('\nTotal run time was', d, 'sec')
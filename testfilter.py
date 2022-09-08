from sys import argv
from time import sleep, time
import os,json,glob
from postDesignClass_v3 import postDesign
import numpy as np
import pandas as pd
from glob import glob
import os, json, shutil, itertools



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


# Set paths for log files and result file
logFolder = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/PFCheck/logs/*'
saveFN = pass_log_folder+xParams['PDBID']+'_CLUSTERX_'+file_incre+'_passedPFC.txt'

# Make a glob that has all the log files
logFiles = glob(logFolder)

# Initialize a dataframe to hold all the information
df = pd.DataFrame(columns=['pdbFN_CUT','pdbFN_MIN','initial_score','bound(minimized)','cut_score',
                            'rlx_test_score','translated(trans-preMin)','unbound(trans-Min)','numHB','total_hEn',
                            'SASA','SC','Bunsat','fa_rep','intra',
                            'fa_atr','elec','IE1_2','IE1_3','RMSD_BS','RMSD_ETE',
                            'Rotamer','FASTA'])


for fn in logFiles:
    logDF = pd.read_csv(fn,delimiter='\t')
    df = df.append(logDF,ignore_index=True)

# Set filter and output to a file
df.reset_index(inplace=True,drop=True)

boundhdxFilter = df['numHB'] >= df.numHB.quantile(0.9)
totalHBEnFilter = df['total_hEn'] <= df.total_hEn.quantile(0.1)
clsFilter = df['fa_rep'] <= df.fa_rep.quantile(0.25)
atrFilter = df['fa_atr'] <= df.fa_atr.quantile(0.25)
intraFilter = df['intra'] <= df.intra.quantile(0.25)
#elecFilter = df['elec'] <= xParams['elecMax']
SASAFilter = df['SASA'] <= df.SASA.quantile(0.25)
SCfilter = df['SC'] >= xParams['SCMin']


combinedFilter = boundhdxFilter & totalHBEnFilter & clsFilter & atrFilter & intraFilter & SASAFilter & SCfilter

good_PFC = df[combinedFilter]
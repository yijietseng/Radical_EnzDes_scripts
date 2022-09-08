'''
After all runs of the clash check and dummy design for the scaffold selection on Marylou have finished,
scores will be output as log files for each rotamer. Before running this analyzer script, download the 
directories from Marylou. This script will combine all log files of each scaffolds and will output a 
master file for plotting the result.
'''
from time import time, sleep
from sys import argv
import os, json, glob
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib as mlb
mlb.use('TkAgg')
import matplotlib.pyplot as plt


def save(arr,fileName):
	fileObject = open(fileName, 'wb')
	pkl.dump(arr, fileObject)
	fileObject.close()

def load(fileName):
	fileObject2 = open(fileName, 'rb')
	modelInput = pkl.load(fileObject2)
	fileObject2.close()
	return modelInput

def combinedata(dir_list, file_incre):
	if not os.path.exists('../SCF_master_result/'):
		os.mkdir('../SCF_master_result/')

	for i in dir_list:
		# Start combining  
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

def plot_result(data,file_incre,FigSize=(15,15)):
	fig,axes = plt.subplots(1,1,figsize=FigSize,sharey=True)
	for i in range(len(data)):
		mediansC = []
		#ave_per_C = []
		avesC = []
		mediansSA = []
		avesSA = []
		#ave_per_SA = []
		for key in data:
			df = data[key]
			clash = df.values[:,0]
			medianC = np.median(clash)
			mediansC.append(medianC)
			aveC = np.average(clash)
			avesC.append(aveC)
			#upperC = np.percentile(clash,75)-medianC
			#lowerC = medianC-np.percentile(clash,25)
			#ave_per_C.append(np.average([upperC,lowerC]))
				
			SASA = df.values[:,1]
			medianSA = np.median(SASA)
			mediansSA.append(medianSA)
			aveSA = np.average(SASA)
			avesSA.append(aveSA)
			#upperSA = np.percentile(SASA,75)-medianSA
			#lowerSA = medianSA-np.percentile(SASA,25)
			#ave_per_SA.append(np.average([upperSA,lowerSA]))
			
		axes.semilogx(mediansC, mediansSA, '.r')    
		#axes.semilogx(avesC, avesSA,'.r')   
			
		
		#axes.errorbar(mediansC,mediansSA, yerr=ave_per_SA, xerr=ave_per_C,capsize=2,fmt='.b')

		for j, txt in enumerate(data.keys()):
			axes.annotate(txt, (mediansC[j]+0.5, mediansSA[j]-1), fontsize=14)
			#axes.annotate(txt, (avesC[j]+0.01, avesSA[j]+2), fontsize=14)
			'''
			axes.annotate(txt, (mediansC[j]+0.01, mediansSA[j]+2), fontsize=14,
				textcoords='offset points', ha='center', va='bottom',
					bbox=dict(boxstyle='round,pad=0.2', fc='yellow', alpha=0.1),
						arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='red'))
			'''
		
	axes.set_title('Clash VS. SASA (solvent accessible suface area) medians of all scaffolds',fontsize= 20, pad=20)
		
	plt.xlabel('Gly Clash score (log)',labelpad=20, fontsize=18)
	plt.ylabel('Designed SASA score', labelpad=30, fontsize=18)
	axes.tick_params(axis='both', labelsize= 15)
	plt.xlim(1,1000)

	final_plt = '../SCF_master_result/clash_SASA'+file_incre+'.pdf'
	plt.savefig(final_plt,dpi = 300, bbox_inches = "tight")
	
	print('The result plot is generated:'+final_plt)


#------------------------------------------------------------------------------------------------------

alist = os.listdir('../') # ls the directory
dir_list = []

# Pull the scaffold directories into a list
for i in alist:
	if len(i) == 4:
		dir_list.append(i)
for f in dir_list:
	if os.path.isfile == True:
		dir_list.remove(f)
dir_list.sort()

'''
#Grab the folder with the largest increment to start with
for i in itertools.count(start=1):
	if not os.path.exists('../'++'/CLUSTERX_'+"{:02d}".format(i)+'/'):
		file_incre = "{:02d}".format(i-1)
		break
'''

# Identify file_incre
if len(argv) > 1:
		file_incre = "{:02d}".format(int(argv[1]))
else:
	file_incre = "01"


# Create the master result files
combinedata(dir_list, file_incre)

# Combine all the master result files into a super mater pickle file using a dictionary
d1 = {} 
for x in dir_list:
	d1[x] = pd.read_csv('../SCF_master_result/'+x+'_master.txt',delimiter='\t',header=None,skiprows=1,names=['Gly_rep','Des_SASA'],verbose=True,usecols=[1,2])
saveFN = '../SCF_master_result/super_masterLog.pkl'
save(d1, saveFN)

# Startloading
data = load(saveFN)
plot_result(data,file_incre)




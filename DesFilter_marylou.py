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

designStep = argv[2] # input hb for hbDesign or hp for hpDesign

# input a file_increment, or it will grab the folder with the largest increment to start with
if len(argv) > 3:
	file_incre = "{:02d}".format(int(argv[3]))
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
	rot1 = enzymeRot(rotamer=df.iloc[i]['Rotamer'],designStage='hbDesigned/',designSuffix=suffix,fileIncre=file_incre,dRes=True,**xParams)
	rot1.runDesign()

def testHpDesign(i):
	suffix = PDB[i][-9:]
	rot1 = enzymeRot(rotamer=363,designStage='hbDesigned/',designSuffix=suffix,dRes=True,**xParams)
	rot1.runDesign()

def testFunc(i):
	''' This is intended to quickly check that pool is running correctly, eventually will change a bit to make an actual test.'''
	dummy=np.linspace(0,i,10)
	print(dummy)

def runDesFilter():
	
	print('Initiating Design filters....')
	# Currently just manually put in some cutoffs, eventually read from file?
	hbCountMin = xParams['hbCountMin']
	clashMax = xParams['clashMax']
	hbEnMax = xParams['hbEnMax']
	intraMax = xParams['intraMax']
	SASAMax = xParams['SASAMax']
	atrMax = xParams['atrMax']

	if designStep == 'hb':
		try:
			os.mkdir('../'+xParams['PDBID']+'/passedHB')
		except:
			print('\nFolder exists:../'+xParams['PDBID']+'/passedHB')

		logFolder = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hbDesigned/logs/*.txt'
		saveFN = '../'+xParams['PDBID']+'/passedHB/'+xParams['PDBID']+'_CLUSTERX_'+file_incre+'_passedHB.txt'

		# Make a glob that has all the log files
		logFiles = glob.glob(logFolder)

		# Initialize a dataframe to hold all the information
		df = pd.DataFrame()

		for fn in logFiles:
			logDF = pd.read_csv(fn,delimiter='\t')
			df = df.append(logDF)
		df.reset_index(inplace=True)
		clashFilter = df['resFaRep'] <= clashMax
		hbCountFilter = df['resNumHB'] >= hbCountMin
		hbEnFilter = df['reshydroHEn'] <= hbEnMax
		intraFilter = df['resFaIntraRep'] <= intraMax
		#BUnsatFilter = df['resBUnSat'] <= BUnsatMax
		combinedFilter = clashFilter & hbCountFilter & hbEnFilter & intraFilter
		#combinedFilter = clashFilter & hbCountFilter & hbEnFilter & intraFilter & BUnsatFilter

	elif designStep == 'hp':
		try:
			os.mkdir('../'+xParams['PDBID']+'/passedHP')
		except:
			print('Folder exists!!!')

		logFolder = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hpDesigned/logs/*.txt'
		saveFN = '../'+xParams['PDBID']+'/passedHP/'+xParams['PDBID']+'_CLUSTERX_'+file_incre+'_passedHP.txt'

		# Make a glob that has all the log files
		logFiles = glob.glob(logFolder)

		# Initialize a dataframe to hold all the information
		df = pd.DataFrame()

		for fn in logFiles:
			logDF = pd.read_csv(fn,delimiter='\t')
			df = df.append(logDF)
		df.reset_index(inplace=True)
		clashFilter = df['resFaRep'] <= clashMax
		hbCountFilter = df['resNumHB'] >= hbCountMin
		hbEnFilter = df['reshydroHEn'] <= hbEnMax
		intraFilter = df['resFaIntraRep'] <= intraMax
		SASAfilter = df['resSASA'] <= SASAMax
		atrfilter = df['resFaAtr'] <= atrMax
		#BUnsatFilter = df['resBUnSat'] <= BUnsatMax
		combinedFilter = clashFilter & hbCountFilter & hbEnFilter & intraFilter & SASAfilter & atrfilter
		#combinedFilter = clashFilter & hbCountFilter & hbEnFilter & intraFilter & SASAfilter & atrfilter & BUnsatFilter



	good = df[combinedFilter]

	'''
	rotList = []
	for row in range(len(goodHB)):
		rotList.append(goodHB['pdbFN'].iloc[row][-12:-9])

	goodHB['rot'] = rotList
	'''

	uniFilter = []
	FASTA_uni = []

	for row in range(len(good)):
		j=good['Rotamer'].iloc[row],good['FASTA'].iloc[row]
		if (j not in FASTA_uni):
			FASTA_uni.append(j)
			uniFilter.append(True)
		else:
			uniFilter.append(False)
	good = good[uniFilter]
	good.reset_index(inplace=True)

	print(len(good),'structure passed filter')
	good.to_csv(saveFN,sep='\t')
	print('Design filter completed!!!!')

#------------------------------------------------------------------------------------------------------------------------------------------

runDesFilter()
	

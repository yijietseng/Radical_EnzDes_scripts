from sys import argv
import os,json,glob
from postDesignClass_v3 import postDesign
import numpy as np
import pandas as pd
from glob import glob
import os, json, shutil, itertools


# load parameter file
with open(argv[1]) as jsonfile:
	xParams = json.load(jsonfile)

stage = argv[2] # stage options are PFC or test
stage = stage.lower()

df = pd.read_csv(argv[3], delimiter='\t')

rowNo = int(argv[4])

# input a file_increment, or it will grab the folder with the largest increment to start with
if len(argv) > 5:
	file_incre = "{:02d}".format(int(argv[5]))
else:
	for i in itertools.count(start=1):
		if not os.path.exists('../'+xParams['PDBID']+'/CLUSTERX_'+"{:02d}".format(i)+'/'):
			file_incre = "{:02d}".format(i-1)
			break

# making necessary folders
try:
	os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/PFCheck')
	os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/PFCheck/logs')
	os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/PFCheck/Clv')
	os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/PFCheck/PFChecked')
	os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/PFCheck/Trans')
	os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/PFCheck/Min_des')
	os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/PFCheck/Rlx_test')
	print('\n\n\nFolders generated\n\n\n')
except:
	print('\n\n\nFolders exist\n\n\n')


#-------------------------------------------Below is the functions in this script---------------------------------
def runPostDesign(i):
	suffix = df.iloc[i]['pdbFN'][-12:]
	rot1 = postDesign(rotamer=df.iloc[i]['Rotamer'], designStage='hpDesigned', fileIncre=file_incre, designSuffix=suffix, **xParams)
	rot1.runPFC()
	print('\n\n\nStarted on rotamer: ',df['pdbFN'][i],'\n\n\n')

def runTestPFC(i,suffix):
	rot1 = postDesign(rotamer=i, designStage='hpDesigned', fileIncre=file_incre, designSuffix=suffix, **xParams)
	rot1.runPFC()

def runPreformFilter():

	# make directory for the result file
	try:
		os.mkdir('../'+xParams['PDBID']+'/passedPFC')
		os.mkdir('../'+xParams['PDBID']+'/passedPFC/result_plots')
		print('Folder created\n\n')
	except:
		print('Folder exists\n\n')

	pass_log_folder = '../'+xParams['PDBID']+'/passedPFC/'

	# Set paths for log files and result file
	logFolder = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/PFChecklogs/*'
	saveFN = pass_log_folder+xParams['PDBID']+'_CLUSTERX_'+file_incre+'_passedPFC.txt'

	# Make a glob that has all the log files
	logFiles = glob(logFolder)

	# Initialize a dataframe to hold all the information
	df = pd.DataFrame(columns=['pdbFN_CUT','pdbFN_MIN','total_ETE_score','initial_score',
							   'bound(minimized)','cut_score','rlx_test_score','translated(trans-preMin)',
							   'unbound(trans-Min)','numHB','total_hEn','SASA','SC','Bunsat','fa_rep','intra',
							   'fa_atr','IE1_2','IE1_3','RMSD_BS','RMSD_ETE','Rotamer','clusterCenter','FASTA'])
	
	for fn in logFiles:
		logDF = pd.read_csv(fn,delimiter='\t')
		df = df.append(logDF,ignore_index=True)

	'''
	# Currently just manually put in some cutoffs
	clashMax = 15
	atrMax = -15
	intra_Max = 20
	bound_hdx_Min = 3
	stg_hb_Min = 1
	hdx_hbNo_Min = 3
	per_hdx_hEnMax = -0.9
	SASA_Max = 8
	rmsd_Max = 1
	IE_Max = 10
	SC_Min = 0.6
	'''

	# Set filter and output to a file
	df.reset_index(inplace=True,drop=True)
	boundhdxFilter = df['numHB'] >= xPrarms['hbCountMin']
	totalHBEnFilter = df['total_hEn'] <= xPrarms['hbEnMax']
	clsFilter = df['fa_rep'] <= xPrarms['clashMax']
	atrFilter = df['fa_atr'] <= xPrarms['atrMax']
	intraFilter = df['intra'] <= xPrarms['intraMax']
	SASAFilter = df['SASA'] <= xPrarms['SASA_Max']
	SCfilter = df['SC'] >= xPrarms['SC_Min']
	#hdx_NoFilter = df['hdx_hbNo'] >= hdx_hbNo_Min
	#G_EnFilter = df['G_hEn'] <= per_hdx_hEnMax
	#D_EnFilter = df['D_hEn'] <= per_hdx_hEnMax
	#E_EnFilter = df['E_hEn'] <= per_hdx_hEnMax
	#IE1_2Filter = df['IE1_2'] <= IE_Max
	#IE1_3Filter = df['IE1_3'] <= IE_Max
	#rmsdBSFilter = df['RMSD_BS'] <= rmsd_Max
	#rmsdETEFilter = df['RMSD_ETE'] <= rmsd_Max
	


	combinedFilter = boundhdxFilter & totalHBEnFilter & clsFilter & atrFilter & intraFilter & SASAFilter & SCfilter
	#& G_EnFilter & D_EnFilter & E_EnFilter & IE1_2Filter & IE1_3Filter  & hdx_NoFilter & rmsdBSFilter & rmsdETEFilter

	good_PFC = df[combinedFilter]
	#good_PFC = good_PFC.sort_values(by=['IE1_3'],ascending=True)
	good_PFC.reset_index(inplace=True)
	print('\nPreformation check has Completed!!!!\n\n', len(good_PFC),'structures passed PFC!!!!')
	print('\nInitiating analysis and plotting processes....')
	designCode = file_incre.lstrip('0')

	good_PFC['short_name'] = [get_shorter_name(somepdb,designCode,False) for somepdb in good_PFC['pdbFN_CUT']]
	good_PFC['hb*SC'] = good_PFC['total_hEn']*good_PFC['SC']

	top10 = print_top(good_PFC, sort_by_colum='hb*SC', ascending=False)

	print('Top 10 structures are:\n',top10)

	good_PFC.to_csv(saveFN,sep='\t',index=False)

	plot_fig,plot_ax,plot_df = product_plot2D('SC','total_hEn','hb*SC',good_PFC)
	plot_fig.savefig('../'+xParams['PDBID']+'/passedPFC/result_plots/'+xParams['PDBID']+'_'+file_incre+'_HB&SC_Product.png',dpi=300)
	print('\nAnalysis and plotting are completed!!!!')

	# Initiate the process that pulls the candidates to a folder
	try:
		os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/top10/')
	except:
		print('Folder exists')

	top10folder = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/top10/'

	for i in range(0,10):
		shutil.copyfile(good_PFC['short_name'][i], top10folder+good_PFC['short_name'][i])
	print('\nThe top 10 structures are copied into:',top10folder)

#------------------------------------------------------------------------------------------------------------------


if stage == 'pfc':
	# Initiate multithreading for preformCheck
	runPostDesign(rowNo)
	print('\nPFC completed')
# For single test run 
elif stage == 'test':
	suffix = '_C_02_01.pdb'
	runTestPFC(2,suffix)
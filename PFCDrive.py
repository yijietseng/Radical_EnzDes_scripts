from sys import argv
from time import sleep, time
from multiprocessing import Pool
import os,json,glob
from postDesignClass_v2 import postDesign
import numpy as np
import pandas as pd
from glob import glob
import os, json, shutil, itertools


# load parameter file
with open(argv[1]) as jsonfile:
	xParams = json.load(jsonfile)

stage = argv[2] # stage options are PFC, RM, ordered or test
stage = stage.lower()
start = time()

if stage == 'ordered':
	PDB = glob.glob('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hpDesigned'+'/*.pdb')
else:
	# Reading data frame
	df = pd.read_csv(argv[3], delimiter='\t')


# input a file_increment, or it will grab the folder with the largest increment to start with
if len(argv) > 4:
	file_incre = "{:02d}".format(int(argv[4]))
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
def get_shorter_name(somepdb,des_code,cut):
	'''
	Nice function to change nasty pdb filenames into manageable codes that are still unique and identifying, based on current naming conventions (7-2-2021)
	Inputs: pdbname		str			Assumed to be formatted as usual in the scoring process
			des_code	str			The design number to be included in the short_name
			cut			bool		Whether this scorefile was taken from cleaved structures or not
	'''
	if cut:
		return 'd'+des_code+'_'+somepdb[-19:-15].lstrip('0')+somepdb[-13:-10].lstrip('0')+somepdb[-10:-8].lstrip('0')
	return 'd'+des_code+'_'+somepdb[-19:-15].lstrip('0')+somepdb[-13:-10].lstrip('0')+somepdb[-10:-8].lstrip('0')

def print_top(df,sort_by_colum,num_print=10,ascending=True):
	temp_df = df
	if ascending:
		temp_df = temp_df.sort_values(sort_by_colum,ignore_index=True, ascending=True)
	else:
		temp_df = temp_df.sort_values(sort_by_colum,ignore_index=True, ascending=False)

	return temp_df[sort_by_colum].iloc[:num_print]

def product_plot2D(term1,term2,prod_term,df,percent_list=[0.005,0.01,0.05,0.1,0.25],term1_lower_bound=False,term2_lower_bound=False,topnum=7,general_plot_marker='.',top_plot_marker='d',manual_plot_marker='x'):
	'''
	A function to make a 2D plot of 2 term product based metrics. Will serve as a blueprint for a 3D version. percent_list should list the percentages of data points (decimal form) that are as good or better than the line. Assumes that higher scores are preferred.
	INPUTS:		term1			str			key for first term in dataframe
				term2			str			key for second term in dataframe
				prod_term		str			key for product in dataframe
				df				DataFrame	dataframe of scores
				percent_list	list		list of float values determining dotted lines
	RETURNS:	(fig,ax)			pyplot figure and axis objects, dataframe of automatic selections
	'''
	# Simple scatter of all the points
	fig,ax = plt.subplots(1)
	ax.set_title(prod_term)
	df.plot.scatter(term1,term2,ax=ax,label='All sequences',marker=general_plot_marker)

	for some_cutoff in percent_list:
		some_index = int(len(df)*some_cutoff)
		prod_val = sorted_df.iloc[some_index][prod_term]
		x = np.linspace(df[term1].min(),df[term1].max(),100)
		ax.plot(x,prod_val/x,'--k')
	if term2_lower_bound:
		ax.set_ylim(term2_lower_bound,df[term2].max()*1.05)
	else:
		ax.set_ylim(df[term2].min(),df[term2].max()*1.05)
	if term1_lower_bound:
		ax.set_xlim(term1_lower_bound,x[-1])
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	return fig,ax

def runPostDesign(i):
	''' A function to run the hydrophobic design, will eventually turn into a function inside the design class, in here to make sure that it properly coordinates with the multithreading module for now.'''
	suffix = df.iloc[i]['pdbFN'][-12:]
	rot1 = postDesign(rotamer=df.iloc[i]['Rotamer'], designStage='hpDesigned', fileIncre=file_incre, designSuffix=suffix, **xParams)
	rot1.runPFC()
	print('\n\n\nStarted on rotamer: ',df['pdbFN'][i],'\n\n\n')

def runPFC_ordered(i):
	suffix = PDB[i][-12:]
	rot1 = postDesign(rotamer=int(PDB[i][-15:-12]), designStage='hpDesigned', fileIncre=file_incre, designSuffix=suffix, **xParams)
	rot1.runPFC()
	print('\n\n\nStarted on rotamer: ',PDB[i],'\n\n\n')

def runTestPFC(i,suffix):
	rot1 = postDesign(rotamer=i, designStage='hpDesigned', fileIncre=file_incre, designSuffix=suffix, **xParams)
	rot1.runPFC()

def runRevertMut(i):
	''' A function to run the hydrophobic design, will eventually turn into a function inside the design class, in here to make sure that it properly coordinates with the multithreading module for now.'''
	suffix = df.iloc[i]['pdbFN'][-12:]
	rot1 = postDesign(rotamer=df.iloc[i]['Rotamer'], designStage='hpDesigned', fileIncre=file_incre, designSuffix=suffix, **xParams)
	rot1.runPFC(revertMut=True)
	print('\n\n\nStarted on rotamer: ',df['pdbFN'][i],'\n\n\n')

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

myPool = Pool(processes=30)

if stage == 'pfc':
	# Initiate multithreading for preformCheck
	for i in myPool.imap_unordered(runPostDesign,range(len(df))):
		print('initiating preform check\n\n\n')
	print('\nPFC completed, initiating filter....')
	runPreformFilter()

elif stage == 'rm':
	# Initiate multithreading for RevertMut
	for i in myPool.imap_unordered(runRevertMut,range(len(df))):
		print('initiating preform check\n\n\n')
	print('\nPFC completed, initiating filter....')
	runPreformFilter()

elif stage == 'ordered':
	# Initiate multithreading for PFC for the ordered structrues
	for i in myPool.imap_unordered(runPFC_ordered,range(len(PDB))):
		print('initiating preform check\n\n\n')
	print('\nPFC completed, initiating filter....')
	runPreformFilter()

# For single test run 
elif stage == 'test':
	suffix = '_C_02_01.pdb'
	runTestPFC(362,suffix)


end = time()

d = end - start

if d >= 3600:
	print('\nTotal run time was', d/3600, 'hr')
if 60 <= d < 3600:
	print('\nTotal run time was', d/60, 'min')
if d < 60:
	print('\nTotal run time was', d, 'sec')
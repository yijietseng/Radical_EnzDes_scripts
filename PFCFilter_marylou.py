from sys import argv
from time import sleep, time
import os,json,glob
from postDesignClass_v3 import postDesign
import numpy as np
import pandas as pd
from glob import glob
import os, json, shutil, itertools
import matplotlib.pyplot as plt


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

def print_top(df,sort_by_colum,num_print=10,ascending=True, record='../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/top_struct/'+xParams['PDBID']+'_'+file_incre+'_topstruct.txt', collect=True):
	temp_df = df
	if ascending:
		temp_df = temp_df.sort_values(sort_by_colum, ascending=True)
		temp_df = temp_df.reset_index(drop=True)
	else:
		temp_df = temp_df.sort_values(sort_by_colum, ascending=False)
		temp_df = temp_df.reset_index(drop=True)

	temp2 = pd.DataFrame(columns=['short_name','hb*SC','numHB','total_hEn','SC','SASA','fa_rep','intra','fa_atr','elec','FASTA'])
	temp2['short_name'] = temp_df['short_name']
	temp2['hb*SC'] = temp_df['hb*SC']
	temp2['numHB'] = temp_df['numHB']
	temp2['total_hEn'] = temp_df['total_hEn']
	temp2['SC'] = temp_df['SC']
	temp2['SASA'] = temp_df['SASA']
	temp2['fa_rep'] = temp_df['fa_rep']
	temp2['intra'] = temp_df['intra']
	temp2['fa_atr'] = temp_df['fa_atr']
	temp2['elec'] = temp_df['elec']
	temp2['FASTA'] = temp_df['FASTA']

	topfolder = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/top_struct/'

	if record:
		temp2[:num_print].to_csv(record,sep='\t')

	if collect:
		if len(temp2) >= 10:
			for i in range(0,10):
				shutil.copyfile(temp_df['pdbFN_CUT'][i], topfolder+temp_df['pdbFN_CUT'][i][-26:])
				Rlxname = temp_df['pdbFN_CUT'][i][-26:-7]+'rlx.pdb'
				shutil.copyfile('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/PFCheck/Rlx_test/'+Rlxname,topfolder+Rlxname)
			print('\nThe top structures are copied into:',topfolder)
		if len(temp2) < 10:
			for i in range(len(temp2)):
				shutil.copyfile(temp_df['pdbFN_CUT'][i], topfolder+temp_df['pdbFN_CUT'][i][-26:])
				Rlxname = temp_df['pdbFN_CUT'][i][-26:-7]+'rlx.pdb'
				shutil.copyfile('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/PFCheck/Rlx_test/'+Rlxname,topfolder+Rlxname)
			print('\nThe top structures are copied into:',topfolder)

	return temp2[:num_print]

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
	#df.plot.scatter(term1,term2,ax=ax,label='All sequences',marker=general_plot_marker)
	plt.scatter(df[term1],df[term2],label='All sequences',marker=general_plot_marker)
	plt.legend(loc='lower right')

	for some_cutoff in percent_list:
		some_index = int(len(df)*some_cutoff)
		prod_val = df.iloc[some_index][prod_term]
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

def runPreformFilter(plot=False):

	# make directory for the result file
	try:
		os.mkdir('../'+xParams['PDBID']+'/passedPFC')
		os.mkdir('../'+xParams['PDBID']+'/passedPFC/result_plots')
		print('Folder created\n\n')
	except:
		print('Folder exists\n\n')

	pass_log_folder = '../'+xParams['PDBID']+'/passedPFC/'

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

	df.reset_index(inplace=True,drop=True)

	manual = input('\nWhich method you wish to use to set thresholds? (m for manul, p for parameter file or a for automatic): ')
	manual = manual.lower()

	if manual == 'm':
		# Set filter and output to a file
		MinNumHB = int(input('Set a value for number of hbonds:'))
		boundhdxFilter = df['numHB'] >= MinNumHB
		MaxHBEn = float(input('Set a value for HBEn:'))
		totalHBEnFilter = df['total_hEn'] <= MaxHBEn
		MaxRep = float(input('Set a value for fa_rep:'))
		clsFilter = df['fa_rep'] <= MaxRep
		MaxAtr = float(input('Set a value for fa_atr:'))
		atrFilter = df['fa_atr'] <= MaxAtr
		MaxIntra = float(input('Set a value for fa_intra:'))
		intraFilter = df['intra'] <= MaxIntra
		#elecFilter = df['elec'] <= xParams['elecMax']
		MaxSASA = float(input('Set a value for SASA:'))
		SASAFilter = df['SASA'] <= MaxSASA
		MinSC = float(input('Set a value for SC:'))
		SCfilter = df['SC'] >= MinSC	
		
		print('\nSetting numHB to:',MinNumHB)
		print('Setting total_hEn to:',MaxHBEn)
		print('Setting fa_rep to:',MaxRep)
		print('Setting fa_atr to:',MaxAtr)
		print('Setting intra to:',MaxIntra)
		print('Setting SASA to:',MaxSASA)
		print('Setting SC to:', MinSC)

	if manual == 'p':
		# Set filter and output to a file
		boundhdxFilter = df['numHB'] >= xParams["hbCountMin"]
		totalHBEnFilter = df['total_hEn'] <= xParams["hbEnMax"]
		clsFilter = df['fa_rep'] <= xParams["clashMax"]
		atrFilter = df['fa_atr'] <= xParams["atrMax"]
		intraFilter = df['intra'] <= xParams["intraMax"]
		#elecFilter = df['elec'] <= xParams['elecMax']
		SASAFilter = df['SASA'] <= xParams["SASAMax"]
		SCfilter = df['SC'] >= xParams['SCMin']	

		print('\nSetting thresholds based on the parameter file')

	if manual == 'a':
		print('\nSetting numHB to the 90% tile):',df.numHB.quantile(0.9))
		print('Setting total_hEn to the 10% tile:',round(df.total_hEn.quantile(0.1),2))
		print('Setting fa_rep to the median:',round(df.fa_rep.quantile(0.5),2))
		print('Setting fa_atr to the median:',round(df.fa_atr.quantile(0.5),2))
		print('Setting intra to the median:',round(df.intra.quantile(0.5),2))
		print('Setting SASA to the 25% tile):',round(df.SASA.quantile(0.25),2),'\n')
		boundhdxFilter = df['numHB'] >= df.numHB.quantile(0.9)
		totalHBEnFilter = df['total_hEn'] <= df.total_hEn.quantile(0.1)
		clsFilter = df['fa_rep'] <= df.fa_rep.quantile(0.5)
		atrFilter = df['fa_atr'] <= df.fa_atr.quantile(0.5)
		intraFilter = df['intra'] <= df.intra.quantile(0.5)
		#elecFilter = df['elec'] <= xParams['elecMax']
		SASAFilter = df['SASA'] <= df.SASA.quantile(0.25)
		SCfilter = df['SC'] >= xParams['SCMin']	
	
	combinedFilter = boundhdxFilter & totalHBEnFilter & clsFilter & atrFilter & intraFilter & SASAFilter & SCfilter
	#& G_EnFilter & D_EnFilter & E_EnFilter & IE1_2Filter & IE1_3Filter  & hdx_NoFilter & rmsdBSFilter & rmsdETEFilter

	good_PFC = df[combinedFilter]

	if len(good_PFC) != 0:
		if manual == 'm':
			#Recording the threshold change
			rec = open('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/threshold_record.txt','w+')
			rec.write('Final thresholds:\n')
			rec.write('numHB\t:\t'+str(MinNumHB)+'\n')
			rec.write('total_hEn\t:\t'+str(MaxHBEn)+'\n')
			rec.write('fa_rep\t:\t'+str(MaxRep)+'\n')
			rec.write('fa_atr\t:\t'+str(MaxAtr)+'\n')
			rec.write('intra\t:\t'+str(MaxIntra)+'\n')
			rec.write('SASA\t:\t'+str(MaxSASA)+'\n')
			rec.write('SC\t:\t'+str(MinSC)+'\n')
			rec.close()
		if manual == 'p':
			rec = open('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/threshold_record.txt','w+')
			rec.write('The Thresholds are based of the parameter file')
		if manual == 'a':
			#Recording the threshold change
			rec = open('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/threshold_record.txt','w+')
			rec.write('Final thresholds:\n')
			rec.write('numHB\t:\t'+str(df.numHB.quantile(0.9))+'\n')
			rec.write('total_hEn\t:\t'+str(round(df.total_hEn.quantile(0.1),2))+'\n')
			rec.write('fa_rep\t:\t'+str(round(df.fa_rep.quantile(0.5),2))+'\n')
			rec.write('fa_atr\t:\t'+str(round(df.fa_atr.quantile(0.5),2))+'\n')
			rec.write('intra\t:\t'+str(round(df.intra.quantile(0.5),2))+'\n')
			rec.write('SASA\t:\t'+str(round(df.SASA.quantile(0.25),2))+'\n')
			rec.close()

		#good_PFC = good_PFC.sort_values(by=['IE1_3'],ascending=True)
		good_PFC.reset_index(inplace=True)
		print('\nPreformation check has Completed!!!!\n\n', len(good_PFC),'structures passed PFC!!!!')
		print('\nInitiating analysis and plotting processes....')
		designCode = file_incre.lstrip('0')

		good_PFC['short_name'] = [get_shorter_name(somepdb,designCode,False) for somepdb in good_PFC['pdbFN_CUT']]
		good_PFC['hb*SC'] = abs(good_PFC['total_hEn'])*abs(good_PFC['SC'])

		# Initiate the process that pulls the candidates to a folder and record the top_struct stats
		try:
			os.mkdir('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/top_struct/')
		except:
			print('Folder exists')

		topfolder = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/top_struct/'

		if len(good_PFC) < 10:	
			top = print_top(good_PFC, sort_by_colum='hb*SC',num_print=len(good_PFC) ,ascending=False)
		else:
			top = print_top(good_PFC, sort_by_colum='hb*SC',ascending=False)
		#print('Top 10 structures are:\n',top_struct)

		good_PFC.to_csv(saveFN,sep='\t',index=False)
		if plot:
			plot_fig,plot_ax = product_plot2D(term1='SC',term2='total_hEn',prod_term='hb*SC',df=good_PFC)
			#plot_fig = product_plot2D(term1='SC',term2='total_hEn',prod_term='hb*SC',df=good_PFC)
			plot_fig.savefig('../'+xParams['PDBID']+'/passedPFC/result_plots/'+xParams['PDBID']+'_'+file_incre+'_HBSC_Product.png',dpi=300)
		print('\nAnalysis and plotting are completed!!!!')
	else:
		print('\nNo structrues passed PFC filter, please change the threshold of the filter!!!!\n')
	

#------------------------------------------------------------------------------------------------------------------


runPreformFilter()




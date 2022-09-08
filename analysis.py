'''
Script for analysis and plotting
Talmage Coates
'''
from typing_extensions import final
from matplotlib import pyplot as plt
from pandas import read_csv,concat
from sys import argv
from itertools import combinations
import numpy as np

SASA_cut=8

def print_top(df,sort_by,num_print=8,column='short_name',SASA_cut = False,ascending=True):
	temp_df = df
	if SASA_cut:
		temp_df = df[df.SASA <= SASA_cut]
	temp_df = temp_df.sort_values(sort_by,ignore_index=True,ascending=ascending)
	print(temp_df[column].iloc[:num_print])


def get_shorter_name(somepdb,des_code,cut):
	'''
	Nice function to change nasty pdb filenames into manageable codes that are still unique and identifying, based on current naming conventions (7-2-2021)
	Inputs: pdbname		str			Assumed to be formatted as usual in the scoring process
			des_num		str			The design number to be included in the short_name
			cut			bool		Whether this scorefile was taken from cleaved structures or not
	'''
	if cut:
		return des_code+'_'+somepdb[-19:-15]+somepdb[-13:-8]
	return des_code+'_'+somepdb[-15:-11]+somepdb[-9:-4]

# Next part turned nasty due to the unfortunate decision to name a run 5_6 instead of 6. Otherwise would have been a simple loop to load the dictionary with the dataframes and add needed columns to them
# Starting with the cleaved forms, will use that to filter down the hp dataframes
df_d2 = read_csv('des2_clv_combined.txt',delimiter='\t')
df_d2['short_name'] = [get_shorter_name(somepdb,'d2',True) for somepdb in df_d2['pdbFN']]

df_d5 = read_csv('des5_clv_combined.txt',delimiter='\t')
df_d5['short_name'] = [get_shorter_name(somepdb,'d5',True) for somepdb in df_d5['pdbFN']]

df_d5_6 = read_csv('des5_6_clv_combined.txt',delimiter='\t')
df_d5_6['short_name'] = [get_shorter_name(somepdb,'d5_6',True) for somepdb in df_d5_6['pdbFN']]

# Now using a temporary dataframe to read in the data, then copy the pertinent information to the cleaved dataframe. Repeated for each design run
keys_for_copying = ['total_score', 'hydroCount', 'fa_rep', 'intra', 'fa_atr', 'SASA', 
					'hydroHEn', 'totalEnzEnergy', 'active_siteSC', 'totalHB']

temp_df = read_csv('des2_hp_combined.txt',delimiter='\t')
temp_df['short_name'] = [get_shorter_name(somepdb,'d2',False) for somepdb in temp_df['pdbFN']]

# Printing out how many structures made it
print('Design run 2 had {:} sequences after hp and {:} went into PFC ({:.2f}%)'.format(len(temp_df),len(df_d2),100.0*len(df_d2)/len(temp_df)))
# Filtering out any structures that did not have a cleaved score
passed_short_names = df_d2.short_name.to_list()
temp_df = temp_df[temp_df.short_name.isin(passed_short_names)]

# Sorting and resetting indices to make the next part cleaner
temp_df = temp_df.sort_values('short_name',ignore_index=True)
df_d2 = df_d2.sort_values('short_name',ignore_index=True)

# Now copying over the values with the slightly modified column names
for somekey in keys_for_copying:
	newkey = 'fusion_'+somekey
	df_d2[newkey] = temp_df[somekey]


temp_df = read_csv('des5_hp_combined.txt',delimiter='\t')
temp_df['short_name'] = [get_shorter_name(somepdb,'d5',False) for somepdb in temp_df['pdbFN']]

# Printing out how many structures made it
print('Design run 5 had {:} sequences after hp and {:} went into PFC ({:.2f}%)'.format(len(temp_df),len(df_d5),100.0*len(df_d5)/len(temp_df)))

# Filtering out any structures that did not have a cleaved score
passed_short_names = df_d5.short_name.to_list()
temp_df = temp_df[temp_df.short_name.isin(passed_short_names)]

# Sorting and resetting indices to make the next part cleaner
temp_df = temp_df.sort_values('short_name',ignore_index=True)
df_d5 = df_d5.sort_values('short_name',ignore_index=True)

# Now copying over the values with the slightly modified column names
for somekey in keys_for_copying:
	newkey = 'fusion_'+somekey
	df_d5[newkey] = temp_df[somekey]


temp_df = read_csv('des5_6_hp_combined.txt',delimiter='\t')
temp_df['short_name'] = [get_shorter_name(somepdb,'d5_6',False) for somepdb in temp_df['pdbFN']]

# Printing out how many structures made it
print('Design run 5_6 had {:} sequences after hp and {:} went into PFC ({:.2f}%)'.format(len(temp_df),len(df_d5_6),100.0*len(df_d5_6)/len(temp_df)))

# Filtering out any structures that did not have a cleaved score
passed_short_names = df_d5_6.short_name.to_list()
temp_df = temp_df[temp_df.short_name.isin(passed_short_names)]

# Sorting and resetting indices to make the next part cleaner
temp_df = temp_df.sort_values('short_name',ignore_index=True)
df_d5_6 = df_d5_6.sort_values('short_name',ignore_index=True)

# Now copying over the values with the slightly modified column names
for somekey in keys_for_copying:
	newkey = 'fusion_'+somekey
	df_d5_6[newkey] = temp_df[somekey]

# Writing out these combined files in case they are desired
df_d2.to_csv('d2_combined.txt',index=False)
df_d5.to_csv('d5_combined.txt',index=False)
df_d5_6.to_csv('d5_6_combined.txt',index=False)

# Combining all 3 now
final_df = concat([df_d2,df_d5,df_d5_6],ignore_index=True)
# Making all values positive for later scoring purposes
for somekey in keys_for_copying:
	final_df[somekey] = np.abs(final_df[somekey])
	final_df['fusion_'+somekey] = np.abs(final_df['fusion_'+somekey])
'''
Putting info I had here from Josh about the previously selected sequences:
SASA		
4M7T/CLUSTERX_2/hpDesigned/4M7T_X_PS310_C_08_02.pdb       						11.46551776	large cavity	
4M7T/CLUSTERX_2/hpDesigned/4M7T_X_PS372_C_08_04.pdb         						7.013617235	Minimum cavity	
4M7T/CLUSTERX_2/hpDesigned/4M7T_X_PS155_C_02_08.pdb         						15.37265061	Large cavity 	
4M7T/CLUSTERX_2/hpDesigned/4M7T_X_PS302_C_01_09.pdb        						13.51254075	mid size cavity	
4M7T/CLUSTERX_2/hpDesigned/4M7T_X_PS245_C_07_08.pdb         						6.917810294	Minimum cavity	
4M7T/CLUSTERX_2/hpDesigned/4M7T_X_PS372_C_08_07.pdb        						6.917810294	Minimum cavity	
4M7T/CLUSTERX_2/hpDesigned/4M7T_X_PS314_C_08_06.pdb        						9.991870341	mid size cavity	
4M7T/CLUSTERX_2/hpDesigned/4M7T_X_PS169_C_02_08.pdb       						13.07063114	Large cavity
This leaves the following as 'bad'
d2_310_08_02
d2_155_02_08
d2_302_01_09
d2_314_08_06
d2_169_02_08

And the following as 'good'
d2_245_07_08
d2_372_08_07
d2_372_08_04


Now the previously selected sequence list, as none of them have holes that are too large:
d2_314_08_04
d2_371_04_04
d2_372_08_09
d2_373_06_04
d5_222_09_08
d5_363_00_04
d5_374_02_03
'''
# Turned that info into lists for better visualization and analysis
old_selections = ['d2_314_08_04','d2_371_04_04','d2_372_08_09','d2_373_06_04','d5_222_09_08','d5_363_00_04','d5_374_02_03']
good_hole_list = ['d2_245_07_08','d2_372_08_07','d2_372_08_04','d2_314_08_04','d2_371_04_04','d2_372_08_09','d2_373_06_04','d5_222_09_08','d5_363_00_04','d5_374_02_03']
bad_hole_list = ['d2_310_08_02','d2_155_02_08','d2_302_01_09','d2_314_08_06','d2_169_02_08']

good_hole_df = final_df[final_df['short_name'].isin(good_hole_list)]
bad_hole_df = final_df[final_df['short_name'].isin(bad_hole_list)]

# Pretty clear the fusion SASA is bad as a cutoff here
#print('Good hole fusion SASA:')
#print(good_hole_df.fusion_SASA)
#print('Bad hole fusion SASA')
#print(bad_hole_df.fusion_SASA)

print('Good hole cleaved SASA:')
print(good_hole_df.SASA)
print('Bad hole cleaved SASA')
print(bad_hole_df.SASA)

# Based on that output, looks like a cleaved SASA of about 8 or 9 should be effective

fig,ax = plt.subplots(3)
df_d2.hist('SASA',ax=ax[0],bins=30)
ax[0].set_xlim(0,25)
ax[0].axvline(x=8)
ax[0].set_title('D2')
df_d5.hist('SASA',ax=ax[1],bins=30)
ax[1].set_xlim(0,25)
ax[1].axvline(x=8)
ax[1].set_title('D5')
df_d5_6.hist('SASA',ax=ax[2],bins=30)
ax[2].set_xlim(0,25)
ax[2].axvline(x=8)
ax[2].set_title('D5-6')
fig.savefig('SASA_distributions.png',dpi=300)

# Making a new DF that takes all of those designs and applies the new cutoff
filtered_final_df = final_df[final_df.SASA <= 9]

# Now that that has cut it from 2135 to 1448, lets look at some rankings and scores
filtered_final_df['fusion_product'] = filtered_final_df.fusion_active_siteSC*filtered_final_df.fusion_totalHB # Metric proposed by Talmage, is focused on the fusion TSR, with total HB energy
filtered_final_df['ETE_product'] = filtered_final_df.active_siteSC*filtered_final_df.totalHB
#Modified version of the product, relies on cleaved scores instead
filtered_final_df['active_triple_product'] = filtered_final_df.active_siteSC/filtered_final_df.SASA*filtered_final_df.hydroHEn # Product proposed by Dr Moody, focuses on cleaved representation and hydroxyl hydrogen bond energy
filtered_final_df['macro_triple_product'] = filtered_final_df.wholeEnzSC/filtered_final_df.SASA*filtered_final_df.hydroHEn # Product proposed by Dr Moody, focuses on cleaved representation and hydroxyl hydrogen bond energy, but considering whole enzyme SC

# Sorting to show rankings
filtered_final_df = filtered_final_df.sort_values('fusion_product',ignore_index=True,ascending=False)
passed_new_df = filtered_final_df[filtered_final_df.short_name.isin(old_selections)]
print('Rankings for originally selected sequences by fusion product out of {:}'.format(len(filtered_final_df)))
print(passed_new_df[['short_name','fusion_product']])

# Sorting to show rankings
filtered_final_df = filtered_final_df.sort_values('ETE_product',ignore_index=True,ascending=False)
passed_new_df = filtered_final_df[filtered_final_df.short_name.isin(old_selections)]
print('Rankings for originally selected sequences by ETE product out of {:}'.format(len(filtered_final_df)))
print(passed_new_df[['short_name','ETE_product']])

# Now doing the same thing but with the active triple product
filtered_final_df = filtered_final_df.sort_values('active_triple_product',ignore_index=True,ascending=False)
passed_new_df = filtered_final_df[filtered_final_df.short_name.isin(old_selections)]
print('Rankings for originally selected sequences by active triple product out of {:}'.format(len(filtered_final_df)))
print(passed_new_df[['short_name','active_triple_product']])

# Now doing the same thing but with the macro triple product
filtered_final_df = filtered_final_df.sort_values('macro_triple_product',ignore_index=True,ascending=False)
passed_new_df = filtered_final_df[filtered_final_df.short_name.isin(old_selections)]
print('Rankings for originally selected sequences by whole enzyme triple product out of {:}'.format(len(filtered_final_df)))
print(passed_new_df[['short_name','macro_triple_product']])


# Making a list of score terms I want to see relationships between
cleaved_terms = ['SASA','fa_atr','active_siteSC','wholeEnzSC','totalHB','hydroHEn']
# And generating a list of the fusion counterparts
fusion_terms = []
for i in range(len(cleaved_terms)):
	fusion_terms.append('fusion_'+cleaved_terms[i])

'''
# Loop to make all the primitive 2D plots
for term1,term2 in combinations(cleaved_terms,2):
	fig,ax=plt.subplots(1)
	filtered_final_df.plot.scatter(term1,term2,ax=ax)
	fn = './primitive_plots/'+term1+'_'+term2+'.png'
	fig.savefig(fn,dpi=300)
'''
# now the product plots
def product_plot2D(term1,term2,prod_term,df,manual_names_list,percent_list=[0.005,0.01,0.05,0.1,0.25],term1_lower_bound=False,term2_lower_bound=False,topnum=7,general_plot_marker='.',top_plot_marker='d',manual_plot_marker='x'):
	'''
	A function to make a 2D plot of 2 term product based metrics. Will serve as a blueprint for a 3D version. percent_list should list the percentages of data points (decimal form) that are as good or better than the line. Assumes that higher scores are preferred.
	INPUTS:		term1			str			key for first term in dataframe
				term2			str			key for second term in dataframe
				prod_term		str			key for product in dataframe
				df				DataFrame	dataframe of scores
				percent_list	list		list of float values determining dotted lines
	RETURNS:	(fig,ax,df)			pyplot figure and axis objects, dataframe of automatic selections
	'''
	# Simple scatter of all the points
	fig,ax = plt.subplots(1)
	ax.set_title(prod_term)
	df.plot.scatter(term1,term2,ax=ax,label='All sequences',marker=general_plot_marker)
	# Scatter plot of manually selected points
	manual_df = df[df['short_name'].isin(manual_names_list)]
	print(manual_df.short_name.to_list())
	manual_df.plot.scatter(term1,term2,ax=ax,marker=manual_plot_marker,label='Manually Selected Sequences',color='r')
	# Scatter plot of top sequences
	sorted_df = df.sort_values(prod_term,ascending=False,ignore_index=True)
	top_df = sorted_df.iloc[:topnum]
	top_df.plot.scatter(term1,term2,ax=ax,marker=top_plot_marker,label='Automatic Selections',color='g')
	# Calculate and plot equipotentials
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
	return fig,ax,top_df

#nested_list = [['fusion_product','fusion_active_site_SC','fusion_totalHB'],['ETE_product','active_siteSC','totalHB']]
fusion_fig,fusion_ax,fusion_df = product_plot2D('fusion_active_siteSC','fusion_totalHB','fusion_product',filtered_final_df,old_selections)
cleaved_fig,cleaved_ax,cleaved_df = product_plot2D('active_siteSC','totalHB','ETE_product',filtered_final_df,old_selections)
fusion_fig.savefig('./Product_plots/fusion_prod_plot.png',dpi=300)
cleaved_fig.savefig('./Product_plots/cleaved_prod_plot.png',dpi=300)

# Making rank columns for each of the product values
# ranking_list = ['fusion_product','ETE_product',]


# Should make a function to fix the pdbFN column to match what Josh would have renamed them to

filtered_final_df.to_csv('filtered_all.csv',index=False)
'''
After generating the csvs, you can use something like the following to get rankings more easily:

import numpy as np
from pandas import read_csv
import matplotlib.pyplot as plt

df = read_csv('big_csv_Radical_SAM_designs.csv')
ranking_list = ['fusion_product','ETE_product','active_triple_product','macro_triple_product']

for someterm in ranking_list:
	new_rank = someterm+'_rank'
	df[new_rank] = df[someterm].rank()
'''
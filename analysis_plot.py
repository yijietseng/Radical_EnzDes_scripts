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
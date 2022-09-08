'''
Used to figure out Pyrosetta stuff
Talmage Coates
Often move end of long comments to just before functions when working with those
'''
import numpy as np
from pandas import read_csv, concat, DataFrame
from biopandas.pdb import PandasPdb
from scipy.spatial.distance import cdist
from scipy.optimize import curve_fit
from Bio.PDB.PDBParser import PDBParser
from glob import glob
from os import mkdir
from sys import argv
import subprocess as sp
from copy import deepcopy
import shlex

''' Commenting out Rosetta stuff for speed
# Initialize, no looking up PDB components just use database. Also, beta_nov16 sf
from pyrosetta import init as rosinit
rosinit('-load_PDB_components False','-beta_nov16 True')

# Other imports
from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
from pyrosetta import pose_from_pdb, dump_pdb, get_fa_scorefxn, Pose, standard_packer_task, MoveMap, create_score_function
from pyrosetta.rosetta.core.scoring.methods import EnergyMethodOptions
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, ResidueNameSelector
from pyrosetta.toolbox.mutants import mutate_residue
from pyrosetta.rosetta.core.select.residue_selector import NeighborhoodResidueSelector
from pyrosetta.rosetta.core.scoring import atom_pair_constraint, fa_atr, fa_rep,hbonds, hbond_sc, hbond_bb_sc, cart_bonded,fa_intra_rep,approximate_buried_unsat_penalty
# Watch stuff with PyMOL
from pyrosetta.rosetta.protocols.moves import PyMOLMover
#pyMove = PyMOLMover()

# Rosetta vectors so I can do things
from pyrosetta.rosetta.utility import vector1_std_string, vector1_unsigned_long

# Scorefunction stuff
sf = create_score_function("beta_nov16")

# Some FR stuff
from pyrosetta.rosetta.protocols.relax import FastRelax

# Setting up all the functions
def setupFR():
	FR = FastRelax()
	mMap = MoveMap()
	mMap.set_chi(True)
	mMap.set_bb(False)
	FR.set_movemap(mMap)
	FR.set_scorefxn(sf)
'''

def calc_scores(pdb,sf,pocket_res,resi=254):
	pose = pose_from_pdb(pdb)
	sf(pose)
	hbsc = pose.energies().residue_total_energies(resi)[hbond_sc]
	hbbb = pose.energies().residue_total_energies(resi)[hbond_bb_sc]
	total_HB = hbsc+hbbb
	sc_pocket_str = ''
	for i in range(len(pocket_res)):
		sc_pocket_str += str(pose.pdb_info().pdb2pose(pocket_res.iloc[i]['chain_id'],pocket_res.iloc[i]['residue_number']))+','
	sc_pocket_str = sc_pocket_str[:-4] # Currently assuming that last digits are the substrate, change after fixing selector
	sc = ShapeComplementarityFilter()
	sc.residues1(str(resi))
	sc.residues2(sc_pocket_str)
	sc_score = sc.score(pose)
	return [pdb,sc_score,total_HB]

def bio_dist_select(pdb,dist_cut=5,resi = 999):
	ppdb = PandasPdb()
	ppdb.read_pdb(pdb)
	atomdf = ppdb.df['ATOM']
	hetdf = ppdb.df["HETATM"]
	resi_df = hetdf[hetdf['residue_number'] == resi]
	otherlig_df = hetdf[hetdf['residue_number'] != resi]
	resi_coords = resi_df[['x_coord','y_coord','z_coord']].to_numpy()
	otherlig_coords = otherlig_df[['x_coord','y_coord','z_coord']].to_numpy()
	atom_coords = atomdf[['x_coord','y_coord','z_coord']].to_numpy()
	resi_atom_dists = cdist(atom_coords,resi_coords)
	resi_otherlig_dists = cdist(otherlig_coords,resi_coords)
	passed_atoms = atomdf[np.amin(resi_atom_dists,axis=1) <= dist_cut]
	passed_other = otherlig_df[np.amin(resi_otherlig_dists,axis=1) <= dist_cut]
	passed_resis = concat([passed_atoms[['residue_number','chain_id']],passed_other[['residue_number','chain_id']]]).drop_duplicates()
	passed_resis = passed_resis[passed_resis['residue_number']!=resi]
	return passed_resis

def find_sulfs(ppdb,cluster_resi,dist_max=3.5):
	# Start off by finding the distance to the sulfurs of all cysteine's
	cysdf = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name']=='CYS']
	sulf_df = cysdf[cysdf['atom_name']=='SG']
	sulf_coords = sulf_df[['x_coord','y_coord','z_coord']].to_numpy()
	clusterdf = ppdb.df['HETATM'][ppdb.df['HETATM']['residue_number']==cluster_resi]
	cluster_coords = clusterdf[['x_coord','y_coord','z_coord']].to_numpy()
	sulf_dists = cdist(sulf_coords,cluster_coords)
	passed_sulf = sulf_df[np.amin(sulf_dists,axis=1) <= dist_max]
	passed_resi = (passed_sulf['residue_number'].drop_duplicates()).to_list()
	# Account for floating methionine
	mlfdf = ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name']=='MLF']
	mlf_sulf_df = mlfdf[mlfdf['atom_name']=='SD'] 
	mlf_coords = mlf_sulf_df[['x_coord','y_coord','z_coord']].to_numpy()
	mlf_dists = cdist(mlf_coords,cluster_coords)
	# Need 0 because array is nested
	mlf_pass = (min(mlf_dists[0]) <= dist_max)
	if int(mlf_pass)+len(passed_resi) != 4:
		print('Number of coordinating cystein residues found:'+str(len(passed_resi)))
		print('Was MLF coordinating?')
		print(mlf_pass)
		bad_sulf_coords = passed_sulf[['x_coord','y_coord','z_coord']].to_numpy()
		bad_sulf_dists = cdist(bad_sulf_coords,cluster_coords)
		print(np.min(bad_sulf_dists,axis=1))
		print(passed_sulf)
		raise ArithmeticError('Wrong number of sulfurs coordinating!')
	return passed_resi,mlf_pass

def fake_find_sulfs(ppdb):
	pass

def get_sulf_dists(ppdb,sulf_list,include_mlf,mlf_resi = 251):
	'''Uses a PandasPDB, a list of sulfur residue indices, a boolean for if the MLF is supposed to be considered, and an optional residue index for MLF'''
	print('Getting dists for sulf list:')
	coords_list = []
	resi_list = sulf_list
	resn_list = []
	# Add cys coords
	for resi in sulf_list:
		resn_list.append('CYS')
		cys_slice = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_number']==resi]
		sulf_slice = cys_slice[cys_slice['element_symbol']=='S']
		coords_list.append(sulf_slice[['x_coord','y_coord','z_coord']].to_numpy()[0])
	# Add MLF if coordinating
	if include_mlf:
		resn_list.append('MLF')
		resi_list.append(mlf_resi)
		mlf_slice = ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name']=='MLF']
		sulf_slice = mlf_slice[mlf_slice['element_symbol']=='S']
		coords_list.append(sulf_slice[['x_coord','y_coord','z_coord']].to_numpy()[0])
	# Calculate distances
	atom_coords = np.array(coords_list)
	dist_array = cdist(atom_coords,atom_coords)
	# Append results to dataframe, using a beast of a lambda function
	dist_df_list = []
	for i in range(len(resi_list)-1):
		dist_df = concat(DataFrame([[resi_list[i],resn_list[i],resi_list[j],resn_list[j],dist_array[i,j]]],columns=['resi1','resn1','resi2','resn2','dist'])for j in range(i+1,len(resi_list)))
		dist_df_list.append(dist_df)
	return concat(dist_df_list,ignore_index=True)

def gen_sulf_dist_constraints(ppdb):
	# Start by getting list of unique resis for SF4
	cluster_df = ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name']=='SF4']
	cluster_resi_list = (cluster_df['residue_number'].drop_duplicates()).to_list()
	#print(cluster_resi_list)
	# Initialize a dataframe for the distances
	dist_df_list = []
	''' Commenting out for now, had to hard code in the residues but may want to return to this in the future
	# Loop through each cluster
	for cluster_resi in cluster_resi_list:
		print('Finding coordinating sulfur for SF4 '+str(cluster_resi))
		cys_resis,mlf = find_sulfs(ppdb,cluster_resi,dist_max=3.7)
		print('Found:')
		print(cys_resis)
		# Append to the distance dataframe
		dist_df_list.append(get_sulf_dists(ppdb,cys_resis,mlf))
	'''
	dist_df_list.append(get_sulf_dists(ppdb,[16, 20, 23],True))
	dist_df_list.append(get_sulf_dists(ppdb,[169, 187, 232, 235],False))
	const_df = concat(dist_df_list,ignore_index=True)
	return const_df

def gen_CA_constraints(ppdb):
	ca_df = ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name']=='CA']
	ca_coords = ca_df[['x_coord','y_coord','z_coord']].to_numpy()
	ca_dists = cdist(ca_coords,ca_coords)
	const_df = DataFrame()
	for i in range(len(ca_df)-2):
		const_df = const_df.append(concat(DataFrame([[ca_df.iloc[i]['residue_number'],ca_df.iloc[i]['residue_name'],[ca_df.iloc[j]['residue_number'],ca_df.iloc[j]['residue_name'],ca_dists[i,j]]]],columns=['resi1','resn1','resi2','resn2','dist']) for j in range(i+1,len(ca_df))),ignore_index=True)
	return const_df

def remove_ETE(pdb,newpdb):
	# Read in PDB dataframe
	ppdb = PandasPdb().read_pdb(pdb)
	# Make filters
	ETE_atom = ['H21','H22','H23','H11','H12','OG1','HG2','OD1','HD2','OE1','HE2','HZ1','HZ2','HE1','HD1','HG1','CB','CA','O1','O2','C1','C2','HA','CG','CD','CE','CZ','HB','C']
	UNJ_filter = ppdb.df['HETATM'].residue_name=='5AU'
	ETE_filter = (ppdb.df['HETATM']['atom_name'].isin(ETE_atom)[-59:]) & UNJ_filter
	# Remove ETE
	ppdb.df['HETATM'] = ppdb.df['HETATM'][~ETE_filter]
	# Reset UNJ filter
	UNJ_filter = ppdb.df['HETATM'].residue_name=='5AU'
	# Rename to 5AD
	ppdb.df['HETATM']['residue_name'][UNJ_filter] = '5AD'
	# Override exiting PDB
	ppdb.to_pdb(path=newpdb,records=['ATOM','HETATM', 'OTHERS'],append_newline=True)
	return

def get_atom_no(filename,resn,resi,atom_name):
	'''Do not let this function be run with unsanitized inputs, it uses shell=True in a run command for now'''
	print('Looking for the {:} of {:} {:} in {:}'.format(atom_name,resn,resi,filename))
	char_pattern = str(resi)+resn+"{:>7}".format(atom_name)
	print(char_pattern)
	atom_no = sp.run("awk '/"+char_pattern+"/{print $3}' "+filename,shell=True,text=True,stdout=sp.PIPE).stdout
	print('Printing the atom number')
	print(atom_no)
	return int(atom_no)

def fix_rad(infile,outfile,resn,atomn):
	from pymol import cmd
	cmd.load(infile)
	num_sel = cmd.select('rad','resn '+resn+' and name '+atomn)
	cmd.h_add('rad')
	cmd.save(outfile)
	return num_sel

def gen_conf_gro(new_dir,newpdb,mlf_gro,ad_gro,ad_df,fixed_ad,mlf_pdb):
	# Use a bash script to make basic GROMACS files
	sp.run(['/home/byu.local/tlc246/Documents/basic_gro.sh',new_dir,newpdb])
	# Get number of atoms in enzyme
	ori_num_atoms = int(sp.run('head -n 2 '+new_dir+'prep.gro | tail -n 1',shell=True,stdout=sp.PIPE,text=True).stdout)
	# Get and apply shifts for MLF atom numbers
	mlf_shift = int(sp.run("tail -n 2 "+new_dir+"prep.gro | head -n1 |awk '{print $3}'",shell=True,text=True,stdout=sp.PIPE).stdout)
	# Have to read in the mlf.gro file each time to use the += setup for atom number
	mlf_df = read_csv(mlf_gro,delim_whitespace=True,header=None,names=['resin','atom_name','atom_number','x','y','z'])
	mlf_df['atom_number'] = mlf_df['atom_number']+mlf_shift
	# Read in fixed 5AD ppdb df
	ad_pdb = (PandasPdb().read_pdb(fixed_ad)).df['HETATM']
	# Correct atom coords for mlf and 5AD, as they are indexed differently in pdb requires some logic and list comprehension
	mlf_df['x'] = [(mlf_pdb[mlf_pdb['atom_name']==mlf_df.iloc[i]['atom_name']]['x_coord'].to_numpy())[0] for i in range(len(mlf_df))]
	mlf_df['x'] = mlf_df['x']/10
	mlf_df['y'] = [(mlf_pdb[mlf_pdb['atom_name']==mlf_df.iloc[i]['atom_name']]['y_coord'].to_numpy())[0] for i in range(len(mlf_df))]
	mlf_df['y'] = mlf_df['y']/10
	mlf_df['z'] = [(mlf_pdb[mlf_pdb['atom_name']==mlf_df.iloc[i]['atom_name']]['z_coord'].to_numpy())[0] for i in range(len(mlf_df))]
	mlf_df['z'] = mlf_df['z']/10
	ad_df['x'] = [(ad_pdb[ad_pdb['atom_name']==ad_df.iloc[i]['atom_name']]['x_coord'].to_numpy())[0] for i in range(len(ad_df))]
	ad_df['x'] = ad_df['x']/10
	ad_df['y'] = [(ad_pdb[ad_pdb['atom_name']==ad_df.iloc[i]['atom_name']]['y_coord'].to_numpy())[0] for i in range(len(ad_df))]
	ad_df['y'] = ad_df['y']/10
	ad_df['z'] = [(ad_pdb[ad_pdb['atom_name']==ad_df.iloc[i]['atom_name']]['z_coord'].to_numpy())[0] for i in range(len(ad_df))]
	ad_df['z'] = ad_df['z']/10

	# Add columns for atom numbers for sulfur distance restraints
	sulf_num_list = []
	# Write filled gro file
	atomno = ori_num_atoms+len(mlf_df)+len(ad_df)
	# Python to write number of atoms
	f = open(new_dir+'conf.gro','w+')
	f.write('Giant Rising Ordinary Mutants for A Clerical Setup\n'+str(atomno)+'\n')
	f.close() # Close to copy over the original atoms
	sp.run('head -n -1 '+new_dir+'prep.gro | tail -n +3 >> '+new_dir+'conf.gro',shell=True)
	f = open(new_dir+'conf.gro','a')
	for i in range(len(mlf_df)):
		f.write('{:>8}{:>7}{:>5}{:>8.3f}{:>8.3f}{:>8.3f}\n'.format(mlf_df.iloc[i]['resin'],mlf_df.iloc[i]['atom_name'],mlf_df.iloc[i]['atom_number'],mlf_df.iloc[i]['x'],mlf_df.iloc[i]['y'],mlf_df.iloc[i]['z']))
	for i in range(len(ad_df)):
		f.write('{:>5}{:>5}{:>5}{:>5}{:>8.3f}{:>8.3f}{:>8.3f}\n'.format(ad_df.iloc[i]['resi'],ad_df.iloc[i]['resn'],ad_df.iloc[i]['atom_name'],ad_df.iloc[i]['atom_number'],ad_df.iloc[i]['x'],ad_df.iloc[i]['y'],ad_df.iloc[i]['z']))
	f.close()
	# Append last line using bash
	sp.run(['tail','-n1',new_dir+'prep.gro','>>',new_dir+'conf.gro'])
	sp.run(['echo','-e',"\n",'>>',new_dir+'conf.gro'])
	return mlf_shift

def prep_MD_folders(pdb_glob,mlf_gro,ad_gro,outerpath = 'SAM_analysis/'):
	# Copy 5AD.gro with corrected coordinates into the file
	ad_df = read_csv(ad_gro,delim_whitespace=True,header=None,names=['resi','resn','atom_name','atom_number','x','y','z'])
	
	for start_pdb in pdb_glob:
		# Make directory to move to supercomputer
		new_dir = outerpath+'for_md/'+start_pdb[2]+start_pdb[5:7]+start_pdb[16:19]+start_pdb[21:27]+'/'
		mkdir(new_dir)
		# Read in PDB
		ppdb = PandasPdb().read_pdb(outerpath+start_pdb)
		# MLF slice extracted
		mlf_pdb = ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name']== 'MLF']
		# Filter down the main ppdb to just 5AD and write to pdb to fix the radical
		ad_ppdb = deepcopy(ppdb)
		ad_ppdb.df['HETATM'] = ppdb.df['HETATM'][ppdb.df['HETATM']['residue_name']== '5AD']
		temp_ad=new_dir+'temp5ad.pdb'
		fixed_ad=new_dir+'fixed_5ad.pdb'
		ad_ppdb.to_pdb(path=temp_ad,records=['HETATM'],append_newline=True)
		sel_check = fix_rad(temp_ad,fixed_ad,'5AD',"C5'")
		# Stop if the selection did not work properly
		if not sel_check:
			print(fixed_ad)
			raise AssertionError('The radical was not fixed! Try again!')
		
		# Save pdb without ligands for GROMACS
		newpdb = new_dir+start_pdb[2]+start_pdb[5:7]+start_pdb[16:19]+start_pdb[21:27]+'.pdb'
		ppdb.to_pdb(path=newpdb,records=['ATOM', 'OTHERS'],append_newline=True)

		# Generate the gro files
		mlf_shift=gen_conf_gro(new_dir,newpdb,mlf_gro,ad_gro,ad_df,fixed_ad,mlf_pdb)

		# Get sulfur dist consts
		try:
			sulf_df = gen_sulf_dist_constraints(ppdb)
		except ArithmeticError:
			print(start_pdb)
			raise AssertionError
		sulf_dict = {'MLF':'SD','CYS':'SG'} #Using this to tell get_atom_no which atom it is looking for in the list comprehension below
		sulf_df['atom_number1'] = [get_atom_no(new_dir+'conf.gro',sulf_df.iloc[i]['resn1'],sulf_df.iloc[i]['resi1'],sulf_dict[sulf_df.iloc[i]['resn1']]) for i in range(len(sulf_df))]
		sulf_df['atom_number2'] = [get_atom_no(new_dir+'conf.gro',sulf_df.iloc[i]['resn2'],sulf_df.iloc[i]['resi2'],sulf_dict[sulf_df.iloc[i]['resn2']]) for i in range(len(sulf_df))]
		
		# Start realTopol.top with mlf_head.top, using bash
		realTop = new_dir+'realTopol.top'
		oldTop = new_dir+'topol.top'
		sp.run(['cp','topol_head.top',realTop])
		with open(realTop,'a') as outf:
			with open(new_dir+'topol.top','r') as inf:
				in_lines = inf.readlines()
			# Copy over the atoms from topol.top
			atoms_line = int(sp.run("awk '/\[ atoms \]/{print NR}' "+oldTop,shell=True,stdout=sp.PIPE,text=True).stdout)
			bonds_line = int(sp.run("awk '/\[ bonds \]/{print NR}' "+oldTop,shell=True,stdout=sp.PIPE,text=True).stdout)
			outf.writelines(in_lines[atoms_line-1:bonds_line-2])

			# Copy atoms from mlf_atoms.top, fixing atom and residue indexing
			mlf_atoms_df = read_csv('./mlf_atoms.top',delim_whitespace=True,header=None,names=['atomno','atom_name','resi','resn','atomtype','atomno2','charge','mass'])
			mlf_atoms_df['resi'] = mlf_atoms_df['resi'] + 250
			mlf_atoms_df['atomno'] = mlf_atoms_df['atomno']+mlf_shift
			mlf_atoms_df['atomno2'] = mlf_atoms_df['atomno']
			outf.writelines(['{:>6}{:>11}{:>7}{:>7}{:>7}{:>7}{:>11}{:>11}\n'.format(mlf_atoms_df.iat[i,0],mlf_atoms_df.iat[i,1],mlf_atoms_df.iat[i,2],mlf_atoms_df.iat[i,3],mlf_atoms_df.iat[i,4],mlf_atoms_df.iat[i,5],mlf_atoms_df.iat[i,6],mlf_atoms_df.iat[i,7]) for i in range(len(mlf_atoms_df))])
			outf.write('\n')

			# Copy bonds from topol.top
			pairs_line = int(sp.run("awk '/\[ pairs \]/{print NR}' "+oldTop,shell=True,stdout=sp.PIPE,text=True).stdout)
			outf.writelines(in_lines[bonds_line-1:pairs_line-2])
			# Copy bonds from MLF
			mlf_bonds_df = read_csv('./mlf_bonds.top',delim_whitespace=True,header=None,names=['a1','a2','func','r','k'])
			mlf_bonds_df['a1'] = mlf_bonds_df['a1']+mlf_shift
			mlf_bonds_df['a2'] = mlf_bonds_df['a2']+mlf_shift
			outf.writelines(['{:>5}{:>6}{:>6}{:>11.5f}{:>12}\n'.format(mlf_bonds_df.iat[i,0],mlf_bonds_df.iat[i,1],mlf_bonds_df.iat[i,2],mlf_bonds_df.iat[i,3],mlf_bonds_df.iat[i,4]) for i in range(len(mlf_bonds_df))])
			# Write bonds from sulfur constraints
			outf.writelines(['{:>5}{:>6}     6{:>9.3f}      259408.0\n'.format(sulf_df.iloc[i]['atom_number1'],sulf_df.iloc[i]['atom_number2'],sulf_df.iloc[i]['dist']/10) for i in range(len(sulf_df))])
			outf.write('\n')

			# Write pairs segment
			angle_line = int(sp.run("awk '/\[ angles \]/{print NR}' "+oldTop,shell=True,stdout=sp.PIPE,text=True).stdout)
			outf.writelines(in_lines[pairs_line-1:angle_line-2])
			mlf_pairs_df = read_csv('./mlf_pairs.top',delim_whitespace=True,header=None,names=['a1','a2','func'])
			mlf_pairs_df['a1'] = mlf_pairs_df['a1']+mlf_shift
			mlf_pairs_df['a2'] = mlf_pairs_df['a2']+mlf_shift
			outf.writelines(['{:>5}{:>6}{:>6}\n'.format(mlf_pairs_df.iat[i,0],mlf_pairs_df.iat[i,1],mlf_pairs_df.iat[i,2]) for i in range(len(mlf_pairs_df))])
			outf.write('\n')

			# Write angles segment
			prop_di_line = int(sp.run("awk '/\[ dihedrals \]/{print NR;exit;}' "+oldTop,shell=True,stdout=sp.PIPE,text=True).stdout)
			outf.writelines(in_lines[angle_line-1:prop_di_line-2])
			# mlf angles
			mlf_angles_df = read_csv('./mlf_angles.top',delim_whitespace=True,header=None,names=['a1','a2','a3','func','theta','cth'])
			mlf_angles_df['a1'] = mlf_angles_df['a1']+mlf_shift
			mlf_angles_df['a2'] = mlf_angles_df['a2']+mlf_shift
			mlf_angles_df['a3'] = mlf_angles_df['a3']+mlf_shift
			outf.writelines(['{:>5}{:>6}{:>6}{:>6}{:>14}{:>14}\n'.format(mlf_angles_df.iat[i,0],mlf_angles_df.iat[i,1],mlf_angles_df.iat[i,2],mlf_angles_df.iat[i,3],mlf_angles_df.iat[i,4],mlf_angles_df.iat[i,5]) for i in range(len(mlf_angles_df))])
			outf.write('\n')
			
			# write proper dihedrals
			improp_di_line = int(sp.run("tail -n +"+str(prop_di_line+10)+" "+oldTop+" | awk '/\[ dihedrals \]/{print NR}'",shell=True,stdout=sp.PIPE,text=True).stdout)+10+prop_di_line
			outf.writelines(in_lines[prop_di_line-1:improp_di_line-3])
			mlf_di_df = read_csv('./mlf_di.top',delim_whitespace=True,header=None,names=['a1','a2','a3','a4','func','C0','C1','C2','C3','C4','C5'])
			mlf_di_df['a1'] = mlf_di_df['a1']+mlf_shift
			mlf_di_df['a2'] = mlf_di_df['a2']+mlf_shift
			mlf_di_df['a3'] = mlf_di_df['a3']+mlf_shift
			mlf_di_df['a4'] = mlf_di_df['a4']+mlf_shift
			outf.writelines(['{:>5}{:>6}{:>6}{:>6}{:>6}{:>11.5f}{:>11.5f}{:>11.5f}{:>11.5f}{:>11.5f}{:>11.5f}\n'.format(mlf_di_df.iat[i,0],mlf_di_df.iat[i,1],mlf_di_df.iat[i,2],mlf_di_df.iat[i,3],mlf_di_df.iat[i,4],mlf_di_df.iat[i,5],mlf_di_df.iat[i,6],mlf_di_df.iat[i,7],mlf_di_df.iat[i,8],mlf_di_df.iat[i,9],mlf_di_df.iat[i,10]) for i in range(len(mlf_di_df))])
			outf.write('\n')
			
			# write improper dihedrals and tail
			outf.writelines(in_lines[improp_di_line-2:improp_di_line])
			outf.write('{:>5}{:>6}{:>6}{:>6}     1   180.00   4.60240   2\n'.format(2+mlf_shift,4+mlf_shift,3+mlf_shift,9+mlf_shift))
			outf.writelines(in_lines[improp_di_line:])
			outf.write('5AD                 1\n')
		# Launch script to do the rest of the Gromacs prep
		sp.call(['./gro_prep.sh',new_dir[:-1]])
	return

# New script to investigate the translation as it goes into prep.gro
def get_trans(grofile,pdbfile,resn='ALA',resi=132):
	ppdb = PandasPdb().read_pdb(pdbfile)
	gro_df = read_csv(grofile,skiprows=2,names=['resin','atom_na','atomno','x','y','z'],skipfooter=1,delim_whitespace=True,engine='python')
	pdb_ala = ppdb.df['ATOM'][ppdb.df['ATOM']['residue_name']=='ALA']
	some_pdb_ala = pdb_ala[pdb_ala['residue_number']==resi]
	some_gro = gro_df[gro_df['resin']==str(resi)+resn]
	diff_list = []
	for i in range(len(some_pdb_ala)):
		pdbcoords = some_pdb_ala.iloc[i].loc[['x_coord','y_coord','z_coord']].to_numpy()
		gro_slice = some_gro[some_gro['atom_na']==some_pdb_ala.iloc[i].loc['atom_name']]
		gro_coords = gro_slice[['x','y','z']].to_numpy()
		if gro_coords.any():
			diff_list.append(pdbcoords-gro_coords)
	return diff_list

def gen_score_table(out_file,ros_score_file='all_emptied.txt',xvg_dir='/home/byu.local/tlc246/Documents/SAM_analysis/analysis/'):
	# First read in Rosetta scores
	ros_df = read_csv(ros_score_file)
	# Add field for the xvg names
	ros_df['xvg'] = pdb2xvg(ros_df['PDB'].to_list())
	# Loop through the xvg files and calculate mean and std
	means=[]
	stds=[]
	print('Loading MD data')
	for i,xvg_file in enumerate(ros_df['xvg']):
		#print(i) # This will show current progress, as this part may be a bit slow
		data = np.loadtxt(xvg_dir+xvg_file,skiprows=24)
		means.append(np.mean(data[:,1]))
		stds.append(np.std(data[:,1]))
	ros_df['MD_mean'] = means
	ros_df['MD_std'] = stds
	ros_df.to_csv(out_file,index=False)
	return ros_df

def pdb2xvg(hp_pdb_list):
	conv_dict={"X":"5","2":"2"}
	return ['d'+conv_dict[pdb[-37]]+'_'+pdb[-15:-11]+pdb[-9:-4]+'_energy.xvg' for pdb in hp_pdb_list]

def get_fasta(pdb):
	ppdb = PandasPdb().read_pdb(pdb)
	seqlist = ppdb.amino3to1().residue_name.to_list()
	return ''.join(seqlist)
	

'''
# Short snippet of code used to determine active site for PROSPER analysis
active_resi_list = []
WT_glob = glob('/home/byu.local/tlc246/enzymedesign2/4M7T/CLUSTERX/initial375/*.pdb')
for pdb_file in WT_glob:
	pass_df = bio_dist_select(pdb_file)
	pass_list = pass_df.residue_number.to_list()
	for resi in pass_list:
		if resi not in active_resi_list and resi <251:
			active_resi_list.append(resi)

Used this to get the sequences for the top scoring pdbs
top_list=['d2_310_08_02_energy.xvg','d2_302_01_09_energy.xvg','d2_371_01_08_energy.xvg','d2_372_08_04_energy.xvg','d2_371_02_03_energy.xvg','d2_247_06_02_energy.xvg','d2_302_03_03_energy.xvg']
top_df=df.loc[df.xvg.isin(top_list)]
FASTA_list = []
for somepdb in top_df.PDB:
	ppdb = PandasPdb().read_pdb(somepdb)
	sequence=(ppdb.amino3to1())['residue_name'].tolist()
	FASTA_list.append(''.join(sequence))
top_df['FASTA']=FASTA_list
top_df.to_csv('top7_with_FASTA.csv')

# Legacy functions
def test_a_pdb(df_pdb_index):
	# Run test of shape complimentarity (Legacy code)
	dists = np.arange(5.0,6.1,step=0.5)
	strlist = []
	scorelist = []
	for somedist in dists:
		somestr,somescore = test_sc(df_pdb_index = 5,somedist=somedist)
		strlist.append(somestr)
		scorelist.append(somescore)
	for i in range(len(dists)):
		print()
		print(dists[i])
		print(scorelist[i])
		print(strlist[i])
	# Looping over difference scores for max difference
	percent_diffs = []
	for i in range(1,len(scorelist)):
		for j in range(i):
			percent_diffs.append(abs((scorelist[i]-scorelist[j])/scorelist[i]))
	return max(percent_diffs)
	
def test_sc(df_pdb_index = 5,somedist=5,set_neighbor_atoms='all'):
	# Legacy code, was used to determine distance cuttoff
	ppdb = PandasPdb()
	somepdb = josh_df.iloc[df_pdb_index]['Test_PDB']
	ppdb.read_pdb(somepdb)
	# Selecting the interface residues as presently calculated
	sub_set = bio_dist_select(ppdb,somedist)

	pymol_string = 'sele dist'+str(somedist)+', '

	for i in range(len(sub_set)):
		pymol_string += 'resi '+str(sub_set.iloc[i]['residue_number'])+' or '
	# Truncate extra or
	pymol_string = pymol_string[:-4]
	# Still need to fix the sc portion again, will be a bit complicated without changing all the chains to match or something
	talmage_sc = calc_scores(somepdb,sf,sub_set)
	#josh_sc = calc_sc(pose,get_josh_pocket_resis(pose))
	#print('\nMy Method')
	#print(talmage_sc)
	return pymol_string,talmage_sc

def trans_func(x,r11,r12,r13,r21,r22,r23,r31,r32,r33,b1,b2,b3):
	print('Shape of x given: {:}'.format(x.shape))
	R = np.asarray([(r11,r12,r13),(r21,r22,r23),(r31,r32,r33)])
	print('R has shape: {:}'.format(R.shape))
	print(R)
	b = (np.array([b1,b2,b3])).reshape((3,1))
	print('built b with shape: {:}'.format(b.shape))
	print(b)
	y = []
	for somerow in x:
		new_x = somerow.reshape(3,1)
		y.append(np.matmul(R,new_x)+b)
	y = np.array(y)
	if len(y.shape) > 2:
		y = y.reshape(len(x),3)
	print('y has type: {:}'.format(type(y)))
	print('Shape of y right before return: {:}'.format(y.shape))
	print(y)
	return y

def gen_trans_func_params(grofile,pdbfile):
	ppdb = PandasPdb().read_pdb(pdbfile)
	gro_df = read_csv(grofile,skiprows=2,names=['resin','atom_na','atomno','x','y','z'],skipfooter=1,delim_whitespace=True,engine='python')
	ca_df = ppdb.df['ATOM'][ppdb.df['ATOM']['atom_name']=='CA']
	gro_ca = gro_df[gro_df['atom_na']=='CA']
	xdata,ydata = [],[]
	for i in range(len(ca_df)):
		xdata.append(ca_df.iloc[i][['x_coord','y_coord','z_coord']].to_numpy())
		ydata.append(gro_ca.iloc[i][['x','y','z']].to_numpy())
	# Build an initial guess
	b0 = ydata[0]-xdata[0]
	my_p0 = [1,0,0,0,1,0,0,0,1,b0[0],b0[1],b0[2]]
	popt,pcov = curve_fit(trans_func,xdata[0:12],ydata[0:12],p0=my_p0)
	return popt
# Legacy scripting
# Used to remove ETE and get all in the same folder, as prep for migration to mary lou
mydf = read_csv('./All_Passed_HP_1.txt')
new_dir = '/home/byu.local/tlc246/Documents/SAM_analysis/'
new_pdb_list = []
for i in range(len(mydf)):
	somepdb = mydf.iloc[i]['PDB']
	if somepdb[-37]	=='X':
		newpdb = new_dir+'des5_'+somepdb[-24:-4]+somepdb[-37]+'_empty.pdb'
	else:
		newpdb = new_dir+'des2_'+somepdb[-24:-4]+somepdb[-37]+'_empty.pdb'
	print(newpdb)
	new_pdb_list.append(newpdb)
	remove_ETE(somepdb,newpdb)

# Used to score all pdbs that passed the hp filter
# Load dataframes
des5_df = read_csv('/home/byu.local/tlc246/enzymedesign2/pipeLineScripts/passedHP/4M7T_CLUSTERX_passedHP_5.txt',sep='\t')
des2_df = read_csv('/home/byu.local/tlc246/enzymedesign2/pipeLineScripts/passedHP/4M7T_CLUSTERX_passedHP_2.txt',sep='\t')
des_dir5 = '/home/byu.local/tlc246/enzymedesign2/4M7T/CLUSTERX/hpDesigned/'
des_dir2 = '/home/byu.local/tlc246/enzymedesign2/4M7T/CLUSTERX_2/hpDesigned/'
test_pdb_list = []
for i in range(len(des5_df)):
	test_pdb_list.append(des_dir5+des5_df.iloc[i]['pdbFN'][-24:])
for i in range(len(des2_df)):
	test_pdb_list.append(des_dir2+des2_df.iloc[i]['pdbFN'][-24:])
sf = create_score_function("beta_nov16")
emo = EnergyMethodOptions()
emo.hbond_options().decompose_bb_hb_into_pair_energies( True )
sf.set_energy_method_options(emo)

score_sets = []

for i in test_pdb_list:
	subset = bio_dist_select(i)
	scores = calc_scores(i,sf,subset)
	score_sets.append(scores)

score_df = DataFrame(score_sets,columns=['PDB','SC','HB'])
score_df['Product'] = score_df['SC']*score_df['HB']
'''


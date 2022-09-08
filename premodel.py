from pyrosetta import *
from pyrosetta.toolbox import *
from pyrosetta.rosetta.core.scoring import hbonds, fa_rep, atom_pair_constraint, approximate_buried_unsat_penalty, fa_intra_rep
from pyrosetta.rosetta.core.pack.task.operation import ReadResfile, RestrictToRepacking
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.conformation import Conformation
from pyrosetta.rosetta.core.id import AtomID
import pymol
from pymol import cmd
from sys import argv
import sys, os, os.path, string, time, subprocess, shutil, itertools, math, random, glob, json
from itertools import combinations, filterfalse
import numpy as np
from numpy import pi, array
import pandas as pd
init('-load_PDB_components False','-beta_nov16 True')

print('\nPyRosetta initiated!!!!')


#--------------------------------------- Section for functions ---------------------------------------
def phase1N2(PDBID,raw_pdb_file):
	'''
	**************** Phase 1: Download PDB and clean out the unwanted organics ****************
	'''
	# Checking the the existance of the PREDES folder and its subfolder 
	if not os.path.exists('../'+PDBID+'/'):
		os.mkdir('../'+PDBID+'/')
	if not os.path.isdir('../'+PDBID+'/PREDES'):
		os.mkdir('../'+PDBID+'/PREDES')
	if not os.path.isdir('../'+PDBID+'/PREDES/INPUTS'):
		os.mkdir('../'+PDBID+'/PREDES/INPUTS')

	print('\nInitiating Phase 1: Generating PDBs....')
	# Please add more if ther is are more
	unwanted_organic = ['GOL','PEG','DTU','25W',
						'IGP','TRS','DTB','PLP',
						'LYS','MT2','PGE','2CN',
						'POP','OTT','2KA','PO4',
						'TRP','ALA','IPA','IMD',
						'X0K','CPS','5X8','41K',
					   '7C5'] 
	unwanted_inorganic = ['SO4','NA','CL']
	cmd.fetch(PDBID)
	# remove not necessary molecules and solvent
	cmd.remove('solvent')
	for i in unwanted_organic:
		cmd.remove('organic and resn '+i)
	for i in unwanted_inorganic:
		cmd.remove('inorganic and resn '+i)
	# Extract the first chain and save it as the raw file
	chain_ls = cmd.get_chains(PDBID)
	cmd.select('firstChain',PDBID+' and chain '+chain_ls[0])
	cmd.alter('firstChain','chain="A"')
	cmd.save(raw_pdb_file, 'firstChain')
	# Further clean the structure and save it as the cleaned file
	cmd.remove('organic')
	cmd.remove('inorganic')
	cmd.select('clean',PDBID+' and chain '+chain_ls[0])
	cmd.alter('clean','segi=""')
	cmd.save(cleaned_pdb,'clean')

	try:
		os.remove(PDBID+'.cif')
	except:
		PDBID_low = PDBID.lower()
		os.remove(PDBID_low+'.cif')

	print('\nPhase 1 completed!!!!')
	print('\nInitiating Phase 2: fixing atom names....')
	cmd.delete('*')

	'''
	**************** Phase 2: Generate raw pdb file ****************
	In this phase, the artificial names from the crystal solving process
	will be fixed. In addition, it will also check for the presence of 
	SF4, SAM, MET, and 5AD. If any one of them is not presence, it will
	stop the program and maunal modeling will be required. The script can
	be continued when manual modeling is done. 
	'''
	cmd.load(raw_pdb)
	# Calling phase 2 script
	subprocess.Popen(['./2_fix_cluster_names.sh', raw_pdb_file])

	# Checking for cofactor existance
	precut_co = ['SAM', 'EEM']
	cut_co = ['MET','5AD']
	cofactor_inor = ['SF4','F3S']
	# Checking if SF4 is there  
	if not bool(cmd.select(PDBID+' and inorganic and resn SF4')):
		print('\n'+cofactor_inor[0]+' is not found in the crystal structure!!!! ')
		print('Manual modeling is required for '+PDBID+'. Program paused!!!!')
		go_or_not = input('Are you done with manual modeling and wish to initiate phase 3? (y/n):')
		if go_or_not in N_ls:
			sys.exit('Program terminated!!!!')
	# Check if SAM, MET or 5AD is there
	else:
		not_there_ls = []
		for i in precut_co:
			# if SAM is in the structrue then break out
			if bool(cmd.select(PDBID+' and organic and resn '+i)):
				print('\nDetacted:'+i)
				break
			# if MET and 5AD are BOTH in the structure then break out
			elif bool(cmd.select(PDBID+' and organic and resn '+'+'.join(cut_co))):
				not_there_ls.append('SAM')
				for j in range(len(cut_co)):
					print('Detacted:'+cut_co[j])
				break
			# Otherwise manual modeling is required
			else:
				for i in range(len(cut_co)):
				   if not bool(cmd.select(PDBID+' and organic and resn '+'+'.join(cut_co))):
					   not_there_ls.append(cut_co[i])
				print('\n This crystal structure is missing:'+','.join(not_there_ls))
				print('Manual modeling is required for '+PDBID+'. Program paused!!!!')
				go_or_not = input('Are you done with manual modeling and wish to initiate phase 3? (y/n):')
				if go_or_not in N_ls:
					sys.exit('Program terminated!!!!')

	cmd.delete('*')
	
	cmd.load(raw_pdb)

	# Count and rename SF4's residue numeber and chain
	cmd.select('SF4',PDBID+' and inorganic and resn SF4')
	No_of_SF4 = len(cmd.get_model('SF4').get_residues())
	print('\n'+str(No_of_SF4)+' SF4 have been detected!!!!')
	pymol.stored.dict = {}
	cmd.iterate('(SF4)', 'stored.dict[(resi)]=1')
	SF4_resi_ls = list(pymol.stored.dict.keys())
	resi_Starts = 500
	for i in range(len(SF4_resi_ls)):
		cmd.select('SF4_'+str(i),raw_pdb[15:-4]+' and resi '+SF4_resi_ls[i])
		cmd.alter('SF4_'+str(i),"chain='B'")
		cmd.alter('SF4_'+str(i),'segi=""')
		cmd.alter('SF4_'+str(i),'resi='+str(resi_Starts))
		resi_Starts += 1

	# Crop and save SAM or 5AD and MET, SF4, and other cofactors
	for i in precut_co:
		if bool(cmd.select(PDBID+' and resn '+i)):
			cmd.select('SAM', PDBID+' and resn '+'+'.join(cofactor_inor)+'+'+i)
			cmd.create('SAM_PDB', i)
			cmd.save(raw_5ad,'SAM_PDB')
			break
	if bool(cmd.select(PDBID+' and organic and resn '+'+'.join(cut_co))):
		cmd.select('5AD', PDBID+' and organic and resn '+'+'.join(cut_co))
		cmd.create('5AD_PDB', '5AD')
		cmd.select('SF4',PDBID+' and inorganic and resn '+'+'.join(cofactor_inor))
		cmd.create('SF4_PDB', 'SF4')
		cmd.save(raw_5ad,'SF4_PDB, 5AD_PDB')

	print('\nraw PDB file for SAM or MET+5AD is saved!!!!')
	print('\nPhase 2 completed!!!!')
	cmd.delete("*")

def phase3(PDBID,AD_path):   
	'''
	**************** Phase 3: Modify SAM ****************
	'''
	print('\nInitiating Phase 3: Modifying SAM')
	precut_co = ['SAM', 'EEM']
	cut_co = ['MET','5AD']
	cofactor_inor = ['SF4','F3S']
	cmd.load(raw_pdb)

	# Check if SAM or MET+5AD are in the crystal structure and call phase 3 script
	for i in precut_co:
		if bool(cmd.select(PDBID+' and organic and resn '+i)):
			SAM_name = i
			subprocess.Popen(['./3_modify_SAM.sh', PDBID, SAM_name])
			break
		elif bool(cmd.select(PDBID+' and organic and resn '+'+'.join(cut_co))):
			MET_ori_name = cut_co[0]
			AD_ori_name = cut_co[1]
			subprocess.Popen(['./3_modify_SAM.sh', PDBID, MET_ori_name, AD_ori_name])
			break
	cmd.delete('*')
	print('\nName modified for: '+ad_path)
	# Pause for the file to be ready
	while not os.path.exists(ad_path):
		time.sleep(0.5)
	# Load back into pymol to cleavage the S-C bond
	cmd.load(AD_path)
	cmd.select('SD','name SD')
	cmd.select('C5',"name C5'")
	cmd.unbond('SD','C5')
	print('\nS-C cleaved!!!!')
	cmd.save(AD_path)
	print('\n'+AD_path+' is saved')
	cmd.load(cleaned_pdb)
	cmd.order(AD_path[15:-4], 'off', 'bottom')
	cmd.save(p_N_s_fixed)
	cmd.delete('*')
	print('\nPhase 3 completed!!!!')

def phase4(init_pdb, cleaned_pdb, init_cons_file):
	'''
	**************** Phase 4: generate constraints for 5AD+MET hbonding patterns and Ca of protein backbone ****************
	'''
	print('\nInitiating Phase 4: Generating constraint files...')
	#*************************** Get the initial Hbonding patterns for 5AD and MET *************************** 
	pose = pose_from_pdb(init_pdb)
	AD_resi = pose.pdb_info().pdb2pose(xParams['tsrCode'], xParams['tsrPDBNum'])
	MET_resi = pose.pdb_info().pdb2pose('C', 503)
	sf = create_score_function('beta_nov16')
	sf(pose)
	hs = hbonds.HBondSet()
	pose.update_residue_neighbors()
	hbonds.fill_hbond_set(pose,False,hs)

	#hbond_set = pose.get_hbonds()
	AD_hbonds = hs.residue_hbonds(AD_resi)
	MET_hbonds = hs.residue_hbonds(MET_resi)
	print('\n5AD has', len(AD_hbonds), 'hbonds')
	print('MET has', len(MET_hbonds), 'hbonds\n')


	f = open (init_cons_file,"w+")
	f.write('# Hbonding patterns for 5AD and MET\n')
	for i in range(1, len(AD_hbonds)+1):
		hbond_info = AD_hbonds[i]
		don_res = hbond_info.don_res()
		acc_res = hbond_info.acc_res()
		don_hatm = hbond_info.don_hatm()
		acc_atm = hbond_info.acc_atm()
		don_name = pose.residue(don_res).atom_name(don_hatm).strip()
		acc_name = pose.residue(acc_res).atom_name(acc_atm).strip()
		don_chain = pose.pdb_info().chain(don_res)
		acc_chain = pose.pdb_info().chain(acc_res)

		don_xyz = pose.residue(don_res).xyz(don_name)
		acc_xyz = pose.residue(acc_res).xyz(acc_name)
		dist = round((don_xyz - acc_xyz).length(),4)

		#print('AtomPair '+don_name+' '+str(don_res)+don_chain+' '+acc_name+' '+str(acc_res)+acc_chain+' FLAT_HARMONIC '+str(dist)+' 0.5 0.5\n')

		f.write('AtomPair '+don_name+' '+str(don_res)+' '+acc_name+' '+str(acc_res)+' FLAT_HARMONIC '+str(dist)+' 0.5 0.5\n')
	
	for i in range(1, len(MET_hbonds)+1):
		hbond_info = MET_hbonds[i]
		don_res = hbond_info.don_res()
		acc_res = hbond_info.acc_res()
		don_hatm = hbond_info.don_hatm()
		acc_atm = hbond_info.acc_atm()
		don_name = pose.residue(don_res).atom_name(don_hatm).strip()
		acc_name = pose.residue(acc_res).atom_name(acc_atm).strip()
		don_chain = pose.pdb_info().chain(don_res)
		acc_chain = pose.pdb_info().chain(acc_res)

		don_xyz = pose.residue(don_res).xyz(don_name)
		acc_xyz = pose.residue(acc_res).xyz(acc_name)
		dist = round((don_xyz - acc_xyz).length(),4)

		#print('AtomPair '+don_name+' '+str(don_res)+don_chain+' '+acc_name+' '+str(acc_res)+acc_chain+' FLAT_HARMONIC '+str(dist)+' 0.5 0.5\n')

		f.write('AtomPair '+don_name+' '+str(don_res)+' '+acc_name+' '+str(acc_res)+' FLAT_HARMONIC '+str(dist)+' 0.5 0.5\n')
	
	pose.dump_pdb(init_pdb)
	print('New INIT.pdb is saved as: '+init_pdb+'\n')
	'''
	*********************************************************************************
	The below section is to calculate the constraints of the Ca of the protein backbone.
	The rules are as follow:
	1. CA constraint pairs have to be at least 7 residues apart
	2. The distances between the Ca pair should be within 15 angstroms
	*********************************************************************************
	'''

	#------------------Get rid of combinations within 7 residues--------------------
	f.write('\n# Constraints between CA atoms of the protein backbone\n')
	pose2 = pose_from_pdb(cleaned_pdb)
	total = pose2.total_residue()
	comb = itertools.combinations(range(1, total+1), 2)
	filtered = filterfalse(lambda x: -7< (x[1]-x[0]) < 7, comb) 

	#------------------Generate CA combinations-----------------------------
	for combinations in filtered:
		res1 = combinations[0]
		res2 = combinations[1]
		
		CA_xyz_1 = pose2.residue(res1).xyz('CA')
		CA_xyz_2 = pose2.residue(res2).xyz('CA')
		dist = (CA_xyz_1 - CA_xyz_2).length()

	#-----------------Keeping the pair within 15 angstroms------------------    
		if dist <= 15:
			f.write('AtomPair CA '+str(res1)+' CA '+str(res2)+' FLAT_HARMONIC '+str(dist)+' 0.25 0.25 \n')
	f.close()

	print('\nConstraint file:'+init_cons_file+' has been generated!!!!')
	print('\nPhase 4 completed!!!!')
	return pose
	
def phase5(pose, init_cons_file, result_file):
	'''
	**************** Phase 5: Relax the initial structure ****************
	'''
	print('\nInitiating Phase 5: Relaxing structure....')
	start = time.time()
	# set up score function
	sf = create_score_function('beta_nov16')
	sf.set_weight(atom_pair_constraint, 1)
	score_before = sf(pose)

	# set up constraint
	cons = pyrosetta.rosetta.protocols.constraint_movers.ConstraintSetMover()
	cons.constraint_file(init_cons_file)
	cons.add_constraints(True)
	cons.apply(pose) 
	
	# set up movemap
	movemap = MoveMap()
	movemap.set_bb(True) 
	movemap.set_chi(True)
	
	# set up relax and initiate it
	relax = FastRelax()
	relax.set_scorefxn(sf)
	relax.set_movemap(movemap)
	relax.apply(pose)

	score_after = sf(pose)

	print ('\nThe initial score of this structure is:', score_before)
	print ('The score after relaxing is:', score_after)
	print ('\n','The distribution of scores is as follow:', '\n')
	sf.show(pose)

	pose.dump_pdb(result_file)
	print('\nRelaxed PDB is saved as: '+relaxed_pdb)

	end = time.time()
	d = end-start
	if d >= 3600:
		print ('\n'*2, 'The time spent on relaxing with constraints was:', d/3600, 'hours')
	if 60 <= d < 3600:
		print ('\n'*2, 'The time spent on relaxing with constraints was:', d/60, 'minutes')
	if d < 60:
		print ('\n'*2, 'The time spent on relaxing with constraints was:', d, 'seconds')
	
	print('\nPhase 5 completed!!!!')  

def phase6(AU_pdb, relaxed_pdb, TS_pdb):
	'''
	 **********************************************************************************************************
	 Phase 6: Using the relaxed 5AD as the reference point to build TS model into the scaffold. If there are
	 more than one ligands, make sure to name the second reactive carbon as "CR", residue# as 3, chain as D
	 and the residue name as "UL2". THis phase will build the sudo bonds into the model.
	 **********************************************************************************************************
	'''
	print('\nInitiating Phase 6: Generating TS fusion and model into scaffold...')
	Enz = relaxed_pdb
	enzyme = Enz[15:-4]
	TS = TS_pdb
	substrate = TS[13:-4]
	result_file = Enz[:-7]+'PS.pdb'
	AU_pdb = struc_folder+xParams['fusligName']+'.pdb'

	# Loading relaxed enzyme (including 5AD) and the transition state model
	cmd.load(Enz)
	cmd.load(TS)
	# Check how many residues there are in the TS pdb.
	TS_index = TS_pdb.rfind('/')+1
	TS_obj_name = TS_pdb[TS_index:-4]
	TS_resi_No = len(cmd.get_model(TS_obj_name).get_residues())
	print('\n'+str(TS_resi_No)+' residues detected in '+TS_pdb)
	

	# Define atom paths for pair fit to the 5AD
	C5_s = substrate+"//D/5AD`2/C5'"
	C5_e = enzyme+"//"+xParams['tsrCode']+"/5AD`"+str(xParams['tsrPDBNum'])+"/C5'"
	C4_s = substrate+"//D/5AD`2/C4'"
	C4_e = enzyme+"//"+xParams['tsrCode']+"/5AD`"+str(xParams['tsrPDBNum'])+"/C4'"
	C3_s = substrate+"//D/5AD`2/C3'"
	C3_e = enzyme+"//"+xParams['tsrCode']+"/5AD`"+str(xParams['tsrPDBNum'])+"/C3'"
	O4_s = substrate+"//D/5AD`2/O4'"
	O4_e = enzyme+"//"+xParams['tsrCode']+"/5AD`"+str(xParams['tsrPDBNum'])+"/O4'"

	# pair_fit the TS and the scaffold
	cmd.pair_fit(C5_s, C5_e, C4_s, C4_e, C3_s, C3_e, O4_s, O4_e)

	# Select H atoms of the C5'
	cmd.select('TS_H1', substrate+" and name H5'2")
	cmd.select('TS_H2', substrate+" and name H5'1")
	cmd.select('TS_HR', substrate+' and name HR')
	cmd.select('EZ_H1', enzyme+" and name H5'1")
	cmd.select('EZ_H2', enzyme+" and name H5'2")
	#cmd.select('EZ_H3', enzyme+" and name H5'1")

	# Extract coordinates for the hydrogens on the C5' of the transition state
	H1_coords = cmd.get_coords('TS_H1', 1).tolist()[0]
	H2_coords = cmd.get_coords('TS_H2', 1).tolist()[0]

	#Setting new coords to the H5's
	cmd.alter_state('1', 'EZ_H1', '(x,y,z) = ('+str(H1_coords[0])+','+str(H1_coords[1])+','+str(H1_coords[2])+')')
	cmd.alter_state('1', 'EZ_H2', '(x,y,z) = ('+str(H2_coords[0])+','+str(H2_coords[1])+','+str(H2_coords[2])+')')
	
	#Creating bond between C5' and Cz for 5AU file
	cmd.select('TS_5AD', substrate+' and resn 5AD')
	cmd.remove("TS_5AD") # Remove the original TS 5AD atoms
	cmd.select('HR', 'name HR*')
	cmd.remove('HR') # remove all HR
	
	cmd.select('5AF', 'resn 5AD+U*')
	cmd.create('5AU', '5AF')
	cmd.select('5AU_C5', "5AU and name C5'")
	cmd.select('5AU_CZ', '5AU and name CZ')
	cmd.bond('5AU_C5','5AU_CZ')
	if TS_resi_No > 2:
		CR_s2 = substrate+"//D/UL2`3/CR"
		cmd.select('5AU_CR','5AU and name CR')
		cmd.bond('5AU_CZ', '5AU_CR')
	cmd.alter('5AU', 'q="1.00"')
	cmd.alter('5AU', 'resn="'+xParams['fusligName']+'"')
	cmd.alter('5AU', 'resi="'+str(xParams['tsrPDBNum'])+'"')
	cmd.alter('5AU', 'chain="'+xParams['tsrCode']+'"')
	cmd.save(AU_pdb,"5AU")

	cmd.select('UNL','resn U*')
	cmd.select('5AD','resn 5AD')
	cmd.remove('UNL 5AD')
	cmd.save(result_file)
	cmd.delete("*")
	print('\nPhase 6 completed!!!!')

def phase7(struc_folder, AU_pdb, TS_pdb, intraMax=10):
	'''
	 **********************************************************************************************************
	 Phase 7: Generate Rotamers based on the psudo bond created in the previous phase
	 This phase has to be modified for two-ligand system
	 **********************************************************************************************************
	'''
	print('\nInitiating Phase 7: Generating rotameric states....')
	cmd.load(TS_pdb)
	# Check how many residues there are in the TS pdb.
	TS_index = TS_pdb.rfind('/')+1
	TS_obj_name = TS_pdb[TS_index:-4]
	TS_resi_No = len(cmd.get_model(TS_obj_name).get_residues())
	print('\n'+str(TS_resi_No-1)+' residues detected in '+TS_pdb+'\n')

	# Generating rotamer folders
	if not os.path.exists(struc_folder+'Fu_Ro/'):
		os.mkdir(struc_folder+'Fu_Ro/')
		os.mkdir(struc_folder+'Fu_Ro/good_rot/')
		print('\nRotamer folders generated!!!!')

	# Set up file path and fa_intra_rep max (by default is 20)
	out_rot = struc_folder+'Fu_Ro/'
	good_rot = struc_folder+'Fu_Ro/good_rot/'
	intra_max = intraMax

	# Load structure into PyRosetta and set score function
	pose = pose_from_pdb(AU_pdb)
	sf = create_score_function('beta_nov16')
	resi = pose.pdb_info().pdb2pose(xParams['tsrCode'],xParams['tsrPDBNum'])
	
	# define atoms for defined torsion angles
	# identify the atoms on 5AD
	C3_atom = AtomID(pose.residue(resi).atom_index("C3'"),resi)
	C4_atom = AtomID(pose.residue(resi).atom_index("C4'"),resi)
	C5_atom = AtomID(pose.residue(resi).atom_index("C5'"),resi)
	O4_atom = AtomID(pose.residue(resi).atom_index("O4'"),resi)
	# identify the atoms on TS
	Cz_atom = AtomID(pose.residue(resi).atom_index(xParams['torsion_atom_list'][0]),resi)
	Ce_atom = AtomID(pose.residue(resi).atom_index(xParams['torsion_atom_list'][1]),resi)
	Cd_atom = AtomID(pose.residue(resi).atom_index(xParams['torsion_atom_list'][2]),resi)
	
	
	# Get the default torsion angle defined by the atoms in degree
	rotate = pose.conformation()
	O4_C4_C5_Cz_old = (rotate.torsion_angle(O4_atom, C4_atom, C5_atom, Cz_atom))*180/pi
	#C4_C5_Cz_Ce_old = (rotate.torsion_angle(Ce_atom, Cz_atom, C5_atom, C4_atom))*180/pi
	#C5_Cz_Ce_Cd_old = (rotate.torsion_angle(C5_atom, Cz_atom, Ce_atom, Cd_atom))*180/pi
	
	# Generate angles
	p_list = [-10,-20,0,10,20]
	ang_Matrix = []
	ang_list = array([O4_C4_C5_Cz_old, 180, 180])

	for ang in ang_list:
		ang_Matrix.append((p_list+ang)*pi/180) 
	for ang in ang_Matrix[1]:
		ang_Matrix[1]=np.append(ang_Matrix[1], ang+2*pi/3)
		ang_Matrix[1]=np.append(ang_Matrix[1], ang+4*pi/3)
		
	# Set torsion angles and filter by intra_rep max
	i = 1
	for ang in ang_Matrix[0]:
		for ang2 in ang_Matrix[1]:
			for ang3 in ang_Matrix[2]:
				rotate.set_torsion_angle(O4_atom, C4_atom, C5_atom, Cz_atom, ang)
				rotate.set_torsion_angle(Ce_atom, Cz_atom, C5_atom, C4_atom, ang2)
				#rotate.set_torsion_angle(C5_atom, Cz_atom, Ce_atom, Cd_atom, ang3)			
				pose.dump_pdb(out_rot+'ROT'+"{:03d}".format(i)+'.pdb')
				sf(pose)
				intra = pose.energies().residue_total_energies(resi)[fa_intra_rep]
				if intra <= intra_max:
					pose.dump_pdb(good_rot+'ROT'+"{:03d}".format(i)+'.pdb')
				i = i+1
	cmd.delete('*')
	print('\nPhase 7 completed!!!!')
	
def phase8(PDBID, PS_pdb, fileIncre=False):
	'''
	 **********************************************************************************************************
	 Phase 8: Model the generated rotamers into the scaffold. 
	 **********************************************************************************************************
	'''
	print('\nInitiating Phase 8: Generating Protein-Substrate models....')
	# Creating design folders 
	if fileIncre:
		file_incre = "{:02d}".format(int(fileIncre))
	for i in itertools.count(start=1):
		if not os.path.exists('../'+PDBID+'/CLUSTERX_'+"{:02d}".format(i)+'/'):
			os.mkdir('../'+PDBID+'/CLUSTERX_'+"{:02d}".format(i))
			os.mkdir('../'+PDBID+'/CLUSTERX_'+"{:02d}".format(i)+'/init_PS')
			file_incre = "{:02d}".format(i)
			break
	
	log_of_all = '../'+PDBID+'/Description_of_all_design.txt'
	des = open(log_of_all, 'a+')
	des.write('CLUSTERX_'+file_incre+':\t'+xParams['DesignDescription']+'\n')
	des.close()


	rotList = glob.glob('../'+PDBID+'/PREDES/Fu_Ro/good_rot/*.pdb')
	struc_r = '../'+PDBID+'/PREDES/Fu_Ro/good_rot/'

	enzyme = PDBID+"_PS"
	PS_folder = '../'+PDBID+'/CLUSTERX_'+file_incre+'/init_PS/'
	TS_file = "ROT"
	
	for i in range (len(rotList)):
		rotamer = rotList[i][-10:]
		cmd.load(PS_pdb)
		cmd.load(struc_r+rotamer)

		C5_r = rotamer[:-4]+"//"+xParams['tsrCode']+"/"+xParams['fusligName']+"`"+str(xParams['tsrPDBNum'])+"/C5'"
		C5_e = enzyme+"//"+xParams['tsrCode']+"/"+xParams['fusligName']+"`"+str(xParams['tsrPDBNum'])+"/C5'"
		C4_r = rotamer[:-4]+"//"+xParams['tsrCode']+"/"+xParams['fusligName']+"`"+str(xParams['tsrPDBNum'])+"/C4'"
		C4_e = enzyme+"//"+xParams['tsrCode']+"/"+xParams['fusligName']+"`"+str(xParams['tsrPDBNum'])+"/C4'"
		C3_r = rotamer[:-4]+"//"+xParams['tsrCode']+"/"+xParams['fusligName']+"`"+str(xParams['tsrPDBNum'])+"/C3'"
		C3_e = enzyme+"//"+xParams['tsrCode']+"/"+xParams['fusligName']+"`"+str(xParams['tsrPDBNum'])+"/C3'"
		O4_r = rotamer[:-4]+"//"+xParams['tsrCode']+"/"+xParams['fusligName']+"`"+str(xParams['tsrPDBNum'])+"/O4'"
		O4_e = enzyme+"//"+xParams['tsrCode']+"/"+xParams['fusligName']+"`"+str(xParams['tsrPDBNum'])+"/O4'"

		cmd.pair_fit(C5_r, C5_e, C4_r, C4_e, C3_r, C3_e, O4_r, O4_e)	
		
		cmd.select("sele,"+ enzyme+" and chain X")
		cmd.remove("sele")
		cmd.save(PS_folder+enzyme+rotamer[-7:])
		cmd.delete(enzyme+ ' and '+ rotamer[:-4])
	print('\nPhase 8 completed!!!!')

def phase9(struc_folder, center_atom='ca', radius=14, min_angle=100, max_angle=180):
	'''
	-------------------------Read me before running the phase----------------------------------
	This section is encoded to find the residues that are pointing at the binding site. The code
	is able to get most of the residues in the binding site. However, manual addtion or deletion
	residues might be required. The program will be paused when the residues are generated, and
	will continue when upon instruction. The default center point is the Ca of the substrate,
	which the radius by default is 14 angstroms. 
	-------------------------------------------------------------------------------------------
	'''
	print('\nInitiating Phase 9: Evaluating binding site residues....')
	
	enzyme = PS_pdb
	result_file = struc_folder+'RESFILE_FORMAT.txt'
	BS_pdb = struc_folder+'BS.pdb'
	#set up threshold
	radius = 14 # the distance to Ca atom of the substrate
	cent_atom = center_atom.lower()
	# below are angles in degrees
	min_ang = min_angle 
	max_ang = max_angle

	# Generating head line for the resi_file formate. 
	r = open(result_file,"w+")
	r.write('resi chain\n')
	r.close()

	# Load enzyme to evaluate the surrounding residues
	cmd.load(enzyme)
	cmd.show('sticks', 'name ca+cb')
	cmd.set('cartoon_flat_sheet', '0')
	cmd.select('cent_atom', 'resi '+str(xParams['tsrPDBNum'])+' and name '+cent_atom)

	# Counting total residues and get their residue ID
	cmd.select('pp', 'polymer') # seleting protein backbone
	No_of_resi = len(cmd.get_model('pp').get_residues())
	pymol.stored.dict = {}
	cmd.iterate('(pp)', 'stored.dict[(resi)]=1')
	pp_resi_ls = list(pymol.stored.dict.keys())

	# Evaluate the angle and distances of the residues to the center atom
	BS_resi_ls = []
	for i in range (len(pp_resi_ls)):
		cmd.select('ca', 'resi '+pp_resi_ls[i]+' and name ca')
		cmd.select('cb', 'resi '+pp_resi_ls[i]+' and name cb')
		dis = round(cmd.distance('dis', 'cent_atom', 'cb'), 2)
		ang = round(cmd.angle('ang', 'cent_atom', 'cb', 'ca'),2)
		if dis <= radius:
			if min_ang < ang <= max_ang:
				cmd.color('red', 'ca + cb')
				BS_resi_ls.append(pp_resi_ls[i])
				r = open (result_file,"a+")	
				r.write(pp_resi_ls[i]+' A\n')
				r.close()
		cmd.delete('dis and ang')
	cmd.delete('pp and cent_atom and ca and cb')
	cmd.save(BS_pdb)
	resi_series = '+'.join(BS_resi_ls)
	os.system('pymol '+BS_pdb+' -d "show sticks,name ca+cb;color white,organic and ele c;color red, resi '+resi_series+'"'+'& code '+result_file+' &')
	print('\nRes_file format generated as:'+result_file)
	print('Residues generated and the PDB file is saved as:'+BS_pdb)
	print('\nPlease evaluate the selected binding site residues manually before continue!!!!')
	
	#print('When you load the BS.pdb into pymol, please copy the commands below and execute them in PyMOL:')
	#print('\nshow sticks, name ca+cb')
	#print('color red, resi '+'+'.join(BS_resi_ls))
	print('\nPhase 9 completed!!!!')

def phase10(PDBID, PS_pdb, struc_folder, SF4chain='B', SF4pdbi=500,MLFchain='C',MLFpdbi=503,TSchain='X',TSpdbi=999):
	'''
	**********************************************************************************************************
	This phase will take in a text file and a pdb, and create a both resfiles and a constraint file. The text
	file should just have an entry for each residue number (based on the pdb) that may hydrogen bond. This file
	should separate the entries with spaces or newline characters.
	**********************************************************************************************************
	'''
	print('\nInitiating phase 10: Generating resfiles and constraint file....')
	cmd.delete('*')

	for i in itertools.count(start=1):
		if not os.path.exists('../'+PDBID+'/CLUSTERX_'+"{:02d}".format(i)+'/'):
			if i > 1:
				file_incre = i-1
			else:
				file_incre = i
			break


	des_folder = '../'+PDBID+'/CLUSTERX_'+"{:02d}".format(file_incre)+'/'
	if not os.path.exists(des_folder):
		os.mkdir(des_folder)
	try:
		os.mkdir(des_folder+'RC/') 
	except:
		print(des_folder+'RC/ already exists!!!!')
	pdbfn = PS_pdb

	# Generating filenames
	resfn = struc_folder+'RESFILE_FORMAT.txt'		# File name for the interesting residues to read in
	resfile = des_folder+'RC/'+PDBID+'.resfile'			# File name for design resfile
	gresfile = des_folder+'RC/'+PDBID+'_GLY.resfile'		# File name for glycine resfile
	confile = des_folder+'RC/'+PDBID+'.cst'					# File Name for the constraint file


	# Gathering the information from the pose, using the pdb info
	# After getting each residue number, getting vectors representing the position of some of the atoms
	pose = pose_from_pdb(pdbfn)

	# Count and rename SF4's residue numeber and chain
	cmd.load(pdbfn)
	cmd.select('SF4',PDBID+' and inorganic and resn SF4')
	No_of_SF4 = len(cmd.get_model('SF4').get_residues())
	print('\n'+str(No_of_SF4)+' SF4 have been detected!!!!')
	pymol.stored.dict = {}
	cmd.iterate('(SF4)', 'stored.dict[(chain, resi)]=1')
	SF4_resi_ls = list(pymol.stored.dict.keys())
	cmd.delete('*')

	SF4resi = pose.pdb_info().pdb2pose(SF4chain,SF4pdbi)
	SF4_Fe4 = pose.residue(SF4resi).xyz("FE4")
	SF4_S = pose.residue(SF4resi).xyz("S3")
	MLFresi = pose.pdb_info().pdb2pose(MLFchain,MLFpdbi)
	MLF_N = pose.residue(MLFresi).xyz("N")
	MLF_O = pose.residue(MLFresi).xyz("OXT")
	MLF_S = pose.residue(MLFresi).xyz("SD")
	MLF_met = pose.residue(MLFresi).xyz("CE")
	TSresi = pose.pdb_info().pdb2pose(TSchain,TSpdbi)
	TS_N1 = pose.residue(TSresi).xyz("N1")
	TS_N6 = pose.residue(TSresi).xyz("N6")

	# Load in the pdb indices for the residues of interest
	resArray = pd.read_csv(resfn,delimiter=' ',dtype={'resi':int,'chain':str})

	# Start writing the resfiles and constraint file

	f = open(resfile,"w+")
	f.write('NATRO\nUSE_INPUT_SC\nstart\n')
	#f.write(str(TSpdbi)+' '+TSchain+' NATAA\n')
	g = open(gresfile,"w+")
	g.write('NATRO\nUSE_INPUT_SC\nstart\n')
	#g.write(str(TSpdbi)+' '+TSchain+' NATAA\n')
	c1 = open(confile,"w+")
	c1.write('# Starting with binding the TS residue to the CA of the residues of interest\n\n')

	# Bind potentially hydrogen bonding residues to the ligands, while finishing resfiles

	for i in range(len(resArray['resi'])):
		f.write(str(resArray['resi'][i])+' '+resArray['chain'][i]+' ALLAA\n')
		g.write(str(resArray['resi'][i])+' '+resArray['chain'][i]+' PIKAA'+' G\n')
		# Find distances to TS to constrain
		eXYZ = pose.residue(pose.pdb_info().pdb2pose(resArray['chain'][i],resArray['resi'][i])).xyz("CA")
		# Calculating displacement vectors for distance
		dVec1 = eXYZ - TS_N1
		dVec2 = eXYZ - TS_N6
		c1.write('AtomPair CA '+str(resArray['resi'][i])+resArray['chain'][i]+' N1 '+str(TSpdbi)+TSchain+' FLAT_HARMONIC '+str(dVec1.norm())+' 0.25 0.25\n')
		c1.write('AtomPair CA '+str(resArray['resi'][i])+resArray['chain'][i]+' N6 '+str(TSpdbi)+TSchain+' FLAT_HARMONIC '+str(dVec2.norm())+' 0.25 0.25\n')

	f.close()
	print('Resfile '+resfile+' written')
	g.close()
	print('Resfile '+gresfile+' written')


	# Bind the MLF to the SF4, now uses pdb indices
	c1.write('\n# Now binding the MLF and SF4 residues\n\n')
	dVec = MLF_N - SF4_Fe4
	c1.write('AtomPair N '+str(MLFpdbi)+MLFchain+' FE4 '+str(SF4pdbi)+SF4chain+' FLAT_HARMONIC '+str(dVec.norm())+' 0.25 0.25\n')
	dVec = MLF_O - SF4_Fe4
	c1.write('AtomPair OXT '+str(MLFpdbi)+MLFchain+' FE4 '+str(SF4pdbi)+SF4chain+' FLAT_HARMONIC '+str(dVec.norm())+' 0.25 0.25\n')
	dVec = MLF_S - SF4_Fe4
	c1.write('AtomPair SD '+str(MLFpdbi)+MLFchain+' FE4 '+str(SF4pdbi)+SF4chain+' FLAT_HARMONIC '+str(dVec.norm())+' 0.25 0.25\n')
	dVec = MLF_met - SF4_S
	c1.write('AtomPair CE '+str(MLFpdbi)+MLFchain+' S3 '+str(SF4pdbi)+SF4chain+' FLAT_HARMONIC '+str(dVec.norm())+' 0.25 0.25\n')

	'''
	**********************************************************************************************************
	Below section is only valid for ETE. The purpose of the below constraint is to fix the torsion angles of 
	ETE on the ribose region during design runs.
	**********************************************************************************************************
	
	c1.write('\n# Constrain the dihedrals on before the carbonyl\n')
	c1.write('Dihedral HB '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CB '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CA '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' HA '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CIRCULARHARMONIC 3.14 0.01\n')
	c1.write('Dihedral OE1 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CE '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CD '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' OD1 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CIRCULARHARMONIC 1.082 0.01\n')
	c1.write('Dihedral OD1 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CD '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CG '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' OG1 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CIRCULARHARMONIC -1.061 0.01\n')
	c1.write('Dihedral OG1 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CG '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CB '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' HB '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CIRCULARHARMONIC 1.286 0.01\n')
	c1.write('Dihedral HA '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CA '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' C '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' O1 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CIRCULARHARMONIC 0.026 0.1\n')
	c1.write('Dihedral O1 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' C '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' O2 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' C1 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CIRCULARHARMONIC 0.00349 0.2\n')
	c1.write('Dihedral C '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' O2 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' C1 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' C2 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CIRCULARHARMONIC 3.135 0.3\n')
	c1.write('Dihedral O2 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' C1 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' C2 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' H21 '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' CIRCULARHARMONIC 0.108 1.05\n')
	
	c1.write("\n# Distance constraints on "+xParams['cst_atoms'][0]+"-"+xParams['cst_atoms'][1]+', and '+xParams['cst_atoms'][0]+"-"+xParams['cst_atoms'][2]+'\n')
	atom_coord1 = pose.residue(pose.pdb_info().pdb2pose(xParams['tsrCode'],xParams['tsrPDBNum'])).xyz(xParams['cst_atoms'][0])
	atom_coord2 = pose.residue(pose.pdb_info().pdb2pose(xParams['tsrCode'],xParams['tsrPDBNum'])).xyz(xParams['cst_atoms'][1])
	atom_coord3 = pose.residue(pose.pdb_info().pdb2pose(xParams['tsrCode'],xParams['tsrPDBNum'])).xyz(xParams['cst_atoms'][2])

	vec1 = atom_coord1 - atom_coord2
	vec2 = atom_coord1 - atom_coord3

	dist1 = str(vec1.norm())
	dist2 = str(vec2.norm())

	c1.write('AtomPair '+xParams['cst_atoms'][0]+' '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' '+xParams['cst_atoms'][1]+' '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' FLAT_HARMONIC '+dist1+' 0.01 0.01\n')
	c1.write('AtomPair '+xParams['cst_atoms'][0]+' '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' '+xParams['cst_atoms'][2]+' '+str(xParams['tsrPDBNum'])+xParams['tsrCode']+' FLAT_HARMONIC '+dist2+' 0.01 0.01\n')
	
	******************************************* Section ends ****************************************************
	'''


	# Bind the CA of the backbone
	# List out the residue indices
	bkbone = np.arange(1,pose.total_residue()+1,dtype=int)
	# Filter out the ligand indices

	if No_of_SF4 > 1:
		for i in range(len(SF4_resi_ls)):
			pose_ID_SF4 = pose.pdb_info().pdb2pose(SF4_resi_ls[i][0],int(SF4_resi_ls[i][1]))
			bkbone[pose_ID_SF4-1] = 0
	bkbone[MLFresi-1] = 0
	bkbone[TSresi-1] = 0
	bkbone = bkbone[np.nonzero(bkbone)]


	# Loop through and only bind together residues that are at least 7 away in index, which typically follows the sequence
	c1.write('\n# Now to bind the CA of the backbone together\n\n')
	try:
		for index in bkbone:
			for ind in bkbone:
					if (index-ind > 7):
						dVec= pose.residue(index).xyz("CA")-pose.residue(ind).xyz("CA")
						# We will see if it allows me to remove the spaces from the output of pose2pdb this easily
						pdbind = pose.pdb_info().pose2pdb(ind).replace(' ','')
						pdbindex = pose.pdb_info().pose2pdb(index).replace(' ','')
						c1.write('AtomPair CA '+pdbind+' CA '+pdbindex+' FLAT_HARMONIC '+str(dVec.norm())+' 0.25 0.25\n')
	except:
		print('some residue has no CA')		
	print('Constraint file '+confile+' written')
	c1.close()
	print('\nPhase 10 completed!!!!')


#--------------------------------------- Section of parameters ---------------------------------------
cmd.delete("*") # clear anything that might be in pymol
N_ls = ['n','no']

ID = argv[1]
ID = ID.upper()

params = 'params4marylou/'+ID+'_params.json'

# load parameter file
with open(params) as jsonfile:
	xParams = json.load(jsonfile)

# Getting the scaffold ID and transition state PDB
PDBID = xParams['PDBID']
PDBID = PDBID.upper()
TS_pdb = '../structure/'+xParams['TSFileName']

# Checking the the existance of the PREDES folder and its subfolder 
if not os.path.isdir('../'+PDBID):
	os.mkdir('../'+PDBID)
if not os.path.isdir('../'+PDBID+'/PREDES'):
	os.mkdir('../'+PDBID+'/PREDES')
if not os.path.isdir('../'+PDBID+'/PREDES/INPUTS'):
	os.mkdir('../'+PDBID+'/PREDES/INPUTS')

# define all paths
struc_folder = '../'+PDBID+'/PREDES/'
struc_folder = struc_folder.upper()
input_folder = struc_folder+'INPUTS/'
init_cons_file = input_folder+'CONS_INIT.cst'
raw_pdb = struc_folder+PDBID+'_raw.pdb'
cleaned_pdb = struc_folder+PDBID+'_clean.pdb'
p_N_s_fixed = struc_folder+PDBID+'_INIT.pdb'
raw_5ad = struc_folder+'5AD_raw.pdb'
ad_path = struc_folder+'5AD.pdb'
relaxed_pdb = struc_folder+PDBID+'_RLX.pdb'
AU_pdb = struc_folder+xParams['fusligName']+'.pdb'
PS_pdb = struc_folder+PDBID+'_PS.pdb'


#--------------------------------------- Execution of sections ---------------------------------------
# print the phase descriptions
print('\n**************************************************************************************************************',
	'\nPhase 1 and 2: Download, clean and fix the names of residues, atoms, and chains',
	"\nPhase 3: Modify SAM including cleaving the bond between C5' and SD",
	'\nPhase 4: Generate initial constraint file for the protein backbone and the hbonding patterns for 5AD and MET',
	'\nPhase 5: Relax the initial structure',
	'\nPhase 6: Build the transition state model into the relaxed scaffold',
	'\nPhase 7: Generate rotameric states',
	'\nPhase 8: Model the rotameric states into the scaffold',
	'\nPhase 9: Find the possible binding site residues and generate the resfile format file',
	'\nPhase 10: Generate resfile and constraint file for design runs',
	'\n**************************************************************************************************************')
phase = input('Which phase do you want to run? (Please refer to the list above and enter a number from 1 to 10):')

if (phase == str(1) or phase == str(2)):
	phase1N2(PDBID, raw_pdb)
	phase3(PDBID,ad_path)
	phase4(p_N_s_fixed, cleaned_pdb, init_cons_file)
	pose = pose_from_pdb(p_N_s_fixed)
	phase5(pose, init_cons_file, relaxed_pdb)
	phase6(AU_pdb, relaxed_pdb, TS_pdb)
	phase7(struc_folder, AU_pdb, TS_pdb)
	if len(argv) > 2:
		file_incre = "{:02d}".format(int(argv[2]))
		phase8(PDBID, PS_pdb, file_incre)
	else:
		phase8(PDBID, PS_pdb)
	if os.path.exists('../'+xParams['PDBID']+'/PREDES/RESFILE_FORMAT.txt'):
		resi_format_check = input('Do you want to use existing RESFILE_FORMAT.txt: (y/n')
		while resi_format_check in N_ls:
			phase9(struc_folder)
			pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
			while pauseNcheck in N_ls:
				term = input('Do you want to terminate the program? (y/n):')
				if not term in N_ls:
					sys.exit('Program terminated!!!!')
				else:
					pauseNcheck = input('Are you done with manual evaluation? (y/n):')
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
		else:
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
	else:
		phase9(struc_folder)
		pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
		while pauseNcheck in N_ls:
			term = input('Do you want to terminate the program? (y/n):')
			if not term in N_ls:
				sys.exit('Program terminated!!!!')
			else:
				pauseNcheck = input('Are you done with manual evaluation? (y/n):')
		phase10(PDBID, PS_pdb, struc_folder)
	print('\n'+PDBID+' pre-design modeling completed!!!!')

if phase == str(3):
	phase3(PDBID,ad_path)
	phase4(p_N_s_fixed, cleaned_pdb, init_cons_file)
	pose = pose_from_pdb(p_N_s_fixed)
	phase5(pose, init_cons_file, relaxed_pdb)
	phase6(AU_pdb, relaxed_pdb, TS_pdb)
	phase7(struc_folder, AU_pdb, TS_pdb)
	if len(argv) > 2:
		file_incre = "{:02d}".format(int(argv[2]))
		phase8(PDBID, PS_pdb, file_incre)
	else:
		phase8(PDBID, PS_pdb)
	if os.path.exists('../'+xParams['PDBID']+'/PREDES/RESFILE_FORMAT.txt'):
		resi_format_check = input('Do you want to use existing RESFILE_FORMAT.txt: (y/n')
		while resi_format_check in N_ls:
			phase9(struc_folder)
			pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
			while pauseNcheck in N_ls:
				term = input('Do you want to terminate the program? (y/n):')
				if not term in N_ls:
					sys.exit('Program terminated!!!!')
				else:
					pauseNcheck = input('Are you done with manual evaluation? (y/n):')
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
		else:
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
	else:
		phase9(struc_folder)
		pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
		while pauseNcheck in N_ls:
			term = input('Do you want to terminate the program? (y/n):')
			if not term in N_ls:
				sys.exit('Program terminated!!!!')
			else:
				pauseNcheck = input('Are you done with manual evaluation? (y/n):')
		phase10(PDBID, PS_pdb, struc_folder)
	print('\n'+PDBID+' pre-design modeling completed!!!!')

if phase == str(4):
	phase4(p_N_s_fixed, cleaned_pdb, init_cons_file)
	pose = pose_from_pdb(p_N_s_fixed)
	phase5(pose, init_cons_file, relaxed_pdb)
	phase6(AU_pdb, relaxed_pdb, TS_pdb)
	phase7(struc_folder, AU_pdb, TS_pdb)
	if len(argv) > 2:
		file_incre = "{:02d}".format(int(argv[2]))
		phase8(PDBID, PS_pdb, file_incre)
	else:
		phase8(PDBID, PS_pdb)
	if os.path.exists('../'+xParams['PDBID']+'/PREDES/RESFILE_FORMAT.txt'):
		resi_format_check = input('Do you want to use existing RESFILE_FORMAT.txt: (y/n):')
		while resi_format_check in N_ls:
			phase9(struc_folder)
			pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
			while pauseNcheck in N_ls:
				term = input('Do you want to terminate the program? (y/n):')
				if not term in N_ls:
					sys.exit('Program terminated!!!!')
				else:
					pauseNcheck = input('Are you done with manual evaluation? (y/n):')
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
		else:
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
	else:
		phase9(struc_folder)
		pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
		while pauseNcheck in N_ls:
			term = input('Do you want to terminate the program? (y/n):')
			if not term in N_ls:
				sys.exit('Program terminated!!!!')
			else:
				pauseNcheck = input('Are you done with manual evaluation? (y/n):')
		phase10(PDBID, PS_pdb, struc_folder)
		print('\n'+PDBID+' pre-design modeling completed!!!!')

if phase == str(5):
	pose = pose_from_pdb(p_N_s_fixed)
	phase5(pose, init_cons_file, relaxed_pdb)
	phase6(AU_pdb, relaxed_pdb, TS_pdb)
	phase7(struc_folder, AU_pdb, TS_pdb)
	if len(argv) > 2:
		file_incre = "{:02d}".format(int(argv[2]))
		phase8(PDBID, PS_pdb, file_incre)
	else:
		phase8(PDBID, PS_pdb)
	if os.path.exists('../'+xParams['PDBID']+'/PREDES/RESFILE_FORMAT.txt'):
		resi_format_check = input('Do you want to use existing RESFILE_FORMAT.txt: (y/n):')
		while resi_format_check in N_ls:
			phase9(struc_folder)
			pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
			while pauseNcheck in N_ls:
				term = input('Do you want to terminate the program? (y/n):')
				if not term in N_ls:
					sys.exit('Program terminated!!!!')
				else:
					pauseNcheck = input('Are you done with manual evaluation? (y/n):')
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
		else:
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
	else:
		phase9(struc_folder)
		pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
		while pauseNcheck in N_ls:
			term = input('Do you want to terminate the program? (y/n):')
			if not term in N_ls:
				sys.exit('Program terminated!!!!')
			else:
				pauseNcheck = input('Are you done with manual evaluation? (y/n):')
		phase10(PDBID, PS_pdb, struc_folder)
		print('\n'+PDBID+' pre-design modeling completed!!!!')

if phase == str(6):
	phase6(AU_pdb, relaxed_pdb, TS_pdb)
	phase7(struc_folder, AU_pdb, TS_pdb)
	if len(argv) > 2:
		file_incre = "{:02d}".format(int(argv[2]))
		phase8(PDBID, PS_pdb, file_incre)
	else:
		phase8(PDBID, PS_pdb)
	if os.path.exists('../'+xParams['PDBID']+'/PREDES/RESFILE_FORMAT.txt'):
		resi_format_check = input('Do you want to use existing RESFILE_FORMAT.txt: (y/n):')
		while resi_format_check in N_ls:
			phase9(struc_folder)
			pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
			while pauseNcheck in N_ls:
				term = input('Do you want to terminate the program? (y/n):')
				if not term in N_ls:
					sys.exit('Program terminated!!!!')
				else:
					pauseNcheck = input('Are you done with manual evaluation? (y/n):')
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
		else:
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
	else:
		phase9(struc_folder)
		pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
		while pauseNcheck in N_ls:
			term = input('Do you want to terminate the program? (y/n):')
			if not term in N_ls:
				sys.exit('Program terminated!!!!')
			else:
				pauseNcheck = input('Are you done with manual evaluation? (y/n):')
		phase10(PDBID, PS_pdb, struc_folder)
		print('\n'+PDBID+' pre-design modeling completed!!!!')

if phase == str(7):
	phase7(struc_folder, AU_pdb, TS_pdb)
	if len(argv) > 2:
		file_incre = "{:02d}".format(int(argv[2]))
		phase8(PDBID, PS_pdb, file_incre)
	else:
		phase8(PDBID, PS_pdb)
	if os.path.exists('../'+xParams['PDBID']+'/PREDES/RESFILE_FORMAT.txt'):
		resi_format_check = input('Do you want to use existing RESFILE_FORMAT.txt: (y/n):')
		while resi_format_check in N_ls:
			phase9(struc_folder)
			pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
			while pauseNcheck in N_ls:
				term = input('Do you want to terminate the program? (y/n):')
				if not term in N_ls:
					sys.exit('Program terminated!!!!')
				else:
					pauseNcheck = input('Are you done with manual evaluation? (y/n):')
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
		else:
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
	else:
		phase9(struc_folder)
		pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
		while pauseNcheck in N_ls:
			term = input('Do you want to terminate the program? (y/n):')
			if not term in N_ls:
				sys.exit('Program terminated!!!!')
			else:
				pauseNcheck = input('Are you done with manual evaluation? (y/n):')
		phase10(PDBID, PS_pdb, struc_folder)
		print('\n'+PDBID+' pre-design modeling completed!!!!')

if phase == str(8):
	if len(argv) > 2:
		file_incre = "{:02d}".format(int(argv[2]))
		phase8(PDBID, PS_pdb, file_incre)
	else:
		phase8(PDBID, PS_pdb)
	if os.path.exists('../'+xParams['PDBID']+'/PREDES/RESFILE_FORMAT.txt'):
		resi_format_check = input('Do you want to use existing RESFILE_FORMAT.txt: (y/n):')
		while resi_format_check in N_ls:
			phase9(struc_folder)
			pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
			while pauseNcheck in N_ls:
				term = input('Do you want to terminate the program? (y/n):')
				if not term in N_ls:
					sys.exit('Program terminated!!!!')
				else:
					pauseNcheck = input('Are you done with manual evaluation? (y/n):')
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
		else:
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
	else:
		phase9(struc_folder)
		pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
		while pauseNcheck in N_ls:
			term = input('Do you want to terminate the program? (y/n):')
			if not term in N_ls:
				sys.exit('Program terminated!!!!')
			else:
				pauseNcheck = input('Are you done with manual evaluation? (y/n):')
		phase10(PDBID, PS_pdb, struc_folder)
		print('\n'+PDBID+' pre-design modeling completed!!!!')

if phase == str(9):
	if os.path.exists('../'+xParams['PDBID']+'/PREDES/RESFILE_FORMAT.txt'):
		resi_format_check = input('Do you want to use existing RESFILE_FORMAT.txt: (y/n):')
		while resi_format_check in N_ls:
			phase9(struc_folder)
			pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
			while pauseNcheck in N_ls:
				term = input('Do you want to terminate the program? (y/n):')
				if not term in N_ls:
					sys.exit('Program terminated!!!!')
				else:
					pauseNcheck = input('Are you done with manual evaluation? (y/n):')
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
		else:
			phase10(PDBID, PS_pdb, struc_folder)
			print('\n'+PDBID+' pre-design modeling completed!!!!')
	else:
		phase9(struc_folder)
		pauseNcheck = input('\nAre you done with manual evaluation? (y/n):')
		while pauseNcheck in N_ls:
			term = input('Do you want to terminate the program? (y/n):')
			if not term in N_ls:
				sys.exit('Program terminated!!!!')
			else:
				pauseNcheck = input('Are you done with manual evaluation? (y/n):')
		phase10(PDBID, PS_pdb, struc_folder)
		print('\n'+PDBID+' pre-design modeling completed!!!!')

if phase == str(10):
	phase10(PDBID, PS_pdb, struc_folder)
	print('\n'+PDBID+' pre-design modeling completed!!!!')
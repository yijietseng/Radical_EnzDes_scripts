from pyrosetta import init, dump_pdb, get_fa_scorefxn, Pose, standard_packer_task, MoveMap, create_score_function
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory, mm_enable
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.scoring import atom_pair_constraint, fa_atr, fa_rep,hbonds, hbond_sc, hbond_bb_sc, cart_bonded,fa_intra_rep,approximate_buried_unsat_penalty, fa_elec
from pyrosetta.rosetta.protocols.rigid import RigidBodyTransMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring.methods import EnergyMethodOptions
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, ResidueNameSelector
import pyrosetta.rosetta.core.simple_metrics.per_residue_metrics as sm
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.optimization import MinimizerOptions
from biopandas.pdb.pandas_pdb import PandasPdb
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric 
import pandas as pd
import sys, fileinput, re, operator
from pyrosetta.rosetta.core.pack.task import TaskFactory, parse_resfile
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.core.pack.task.operation import OperateOnResidueSubset, RestrictToRepacking,PreventRepackingRLT,RestrictToRepackingRLT
from pyrosetta.rosetta.protocols.moves import TrialMover,MonteCarlo,SequenceMover
from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
from pyro_scratch import bio_dist_select # This is a script that Talmage wrote
from typing_extensions import final
from matplotlib import pyplot as plt
from itertools import combinations
import numpy as np


#init(extra_options='-dunbrack_prob_buried 1.0 -dunbrack_prob_nonburied 1.0 -dunbrack_prob_buried_semi 1.0')
init('-load_PDB_components False','-beta_nov16 True')


class postDesign():
	def __init__(self, rotamer=1, designStage='hpDesigned/', fileIncre="01", designSuffix='.pdb',analysisFile=False, **kwargs):
		self.__dict__.update(kwargs)
		self.rotamer = rotamer
		self.designStage = designStage
		self.file_incre = fileIncre
		self.designSuffix = designSuffix
		self.revert_resfile = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/RC/'+self.revert_resfile
		self.SASA_cut= self.SASAMax
		self.designNo = fileIncre.lstrip('0') # remove the leading zeros
		if analysisFile:
			self.analysisFN = analysisFile

		# Try to automatically set self.pdbFile, warn if it fails
		try:
			self.setPDB()
		except:
			print('WARNING: PDBFILE NOT AUTOMATICALLY SET')

		try:
			self.setOutputs()
		except:
			print('WARNING: OUTPUT PATH NOT AUTOMATICALLY SET')

	def setPDB(self):
		'''This function sets self.pdbFile according to the current member variables, following the current file tree'''
		self.pdbFile = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.designStage+'/'+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix
		print('PDB file set to: '+self.pdbFile)
		return

	def setOutputs(self, outputStage='PFC',setClv='Clv/', setMin_Des='Min_des/', setTrans_PDB='Trans/', setTRlx_PDB='PFChecked/', setTRlx_log='logs/', setRlx_test='Rlx_test/'):
		'''outputStage options: PFC, RM'''
		if outputStage == 'PFC':
			if setClv:
				self.ClvPath = 'PFCheck/'+setClv
				self.ClvPDB = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.ClvPath+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix[:-4]+'_CUT.pdb'
				print('Ligand-Cleaved PDB path is set to:',self.ClvPDB)
			if setMin_Des:
				self.Min_Des_Path = 'PFCheck/'+setMin_Des
				self.MinPDB = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.Min_Des_Path+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix[:-4]+'_MIN.pdb'
				print('Minimized design PDB path is set to:',self.MinPDB)
			if setTrans_PDB:
				self.Trans_Path = 'PFCheck/'+setTrans_PDB
				self.TransPDB = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.Trans_Path+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix[:-4]+'_TRANS.pdb'
				print('Translated PDB path is set to:',self.TransPDB)
			if setTRlx_PDB:
				self.TRlx_Path = 'PFCheck/'+setTRlx_PDB
				self.TRlxPDB = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.TRlx_Path+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix[:-4]+'_TRlx.pdb'
				print('PreformChecked PDB path is set to:',self.TRlxPDB)
			if setTRlx_log:
				self.log_Path = 'PFCheck/'+setTRlx_log
				self.TRlx_log = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.log_Path+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix[:-4]+'_log.txt'
				print('PreformChecked log path is set to:',self.TRlx_log)
			if setRlx_test:
				self.Rlx_Path = 'PFCheck/'+setRlx_test
				self.Rlx_test = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.Rlx_Path+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix[:-4]+'_rlx.pdb'
				print('Relax test path is set to:',self.Rlx_test)
		if outputStage == 'RM':
			if setClv:
				self.ClvPath = 'RevertMut/'+setClv
				self.ClvPDB = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.ClvPath+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix[:-4]+'_CUT.pdb'
				print('Ligand-Cleaved PDB path is set to:',self.ClvPDB)
			if setMin_Des:
				self.Min_Des_Path = 'RevertMut/'+setMin_Des
				self.MinPDB = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.Min_Des_Path+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix[:-4]+'_MIN.pdb'
				print('Minimized design PDB path is set to:',self.MinPDB)
			if setTrans_PDB:
				self.Trans_Path = 'RevertMut/'+setTrans_PDB
				self.TransPDB = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.Trans_Path+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix[:-4]+'_TRANS.pdb'
				print('Translated PDB path is set to:',self.TransPDB)
			if setTRlx_PDB:
				self.TRlx_Path = 'RevertMut/'+setTRlx_PDB
				self.TRlxPDB = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.TRlx_Path+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix[:-4]+'_TRlx.pdb'
				print('PreformChecked PDB path is set to:',self.TRlxPDB)
			if setTRlx_log:
				self.log_Path = 'RevertMut/'+setTRlx_log
				self.TRlx_log = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.log_Path+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix[:-4]+'_log.txt'
				print('PreformChecked log path is set to:',self.TRlx_log)
			if setRlx_test:
				self.Rlx_Path = 'RevertMut/'+setRlx_test
				self.Rlx_test = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.Rlx_Path+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix[:-4]+'_rlx.pdb'
				print('Relax test path is set to:',self.Rlx_test)
		return
	
	def cut_lig(self, filename, load_pdb=True):

		#Generating all Ligand atom name
		pose2cut = pose_from_pdb(filename)
		AD_atom_names=[	" C5'"," C4'"," O4'"," C1'"," C2'",
 						" C3'"," O3'"," O2'",' N9 ',' C4 ',
 						' N3 ',' C2A',' N1 ',' C6 ',' C5 ',
						' N7 ',' C8 ',' N6 ',"H5'1","H5'2",
 						" H4'",' H1 '," H2'",' H31',' H32',
 						' H22',' H2A',' H8 ','HN61','HN62']

		Lig_resi = pose2cut.pdb_info().pdb2pose(self.tsrCode,self.tsrPDBNum)
		Lig_atom_No = pose2cut.residue(Lig_resi).natoms()
		Lig_atom_ls = []
		
		for i in range(1, Lig_atom_No+1):
			Lig_atom_ls.append(pose2cut.residue(Lig_resi).atom_name(i))
		
		# Take out the 5AD atom from the fusion ligand
		Lig_atom_ls = [x for x in Lig_atom_ls if x not in AD_atom_names]
		Lig_atom_ls = [s.strip() for s in Lig_atom_ls]
		
		# Read in PDB dataframe
		ppdb = PandasPdb().read_pdb(filename)

		# Change the residue names to separate ligand and 5AD
		Lig_filter = ppdb.df['HETATM']['atom_name'].isin(Lig_atom_ls)[-Lig_atom_No:]
		ppdb.df['HETATM']['residue_name'][-Lig_atom_No:].loc[Lig_filter] =  self.cutligName
		ppdb.df['HETATM']['chain_id'][-Lig_atom_No:].loc[Lig_filter] =  self.tsrCode2
		ppdb.df['HETATM']['residue_number'][-Lig_atom_No:].loc[Lig_filter] =  self.tsrPDBNum2

		# Override exiting PDB
		ppdb.to_pdb(path= self.ClvPDB,
					records=['ATOM','HETATM', 'OTHERS'],
					gz = False,
					append_newline=True)

		# Change all remain fused ligand to 5AD
		for line in fileinput.input(self.ClvPDB, inplace=True): 
			sys.stdout.write(line.replace(self.fusligName,'5AD'))
	
		# Section below is to change the CONACT info in the pdb
		a1 = re.compile(r".{4}  C5' 5AD "+self.tsrCode+" "+str(self.tsrPDBNum))
		a2 = re.compile(r".{4}  CZ  "+self.cutligName+" "+self.tsrCode2+" "+str(self.tsrPDBNum2))

		with open(self.ClvPDB) as f:
			for aline in f:
				if a1.search(aline):
					a1_No = aline[7:11]
				if a2.search(aline):
					a2_No = aline[7:11]
					
		c1 = re.compile(r"CONECT\s*"+a1_No)
		c2 = re.compile(r"CONECT\s*"+a2_No)

		with open(self.ClvPDB) as f:
			for cline in f:
				if c1.search(cline):
					c1_str = cline
				if c2.search(cline):
					c2_str = cline

		c1_cut = c1_str.replace(' '+a2_No,'')
		c2_cut = c2_str.replace(' '+a1_No,'')

		for cline in fileinput.input([self.ClvPDB], inplace=True):
			sys.stdout.write(cline.replace(c1_str,c1_cut))
		for cline in fileinput.input([self.ClvPDB], inplace=True):
			sys.stdout.write(cline.replace(c2_str,c2_cut))

		if load_pdb:
			pose = pose_from_pdb(self.ClvPDB)
			print('Cut PDB:', self.ClvPDB, 'loaded')
		
		return pose

	def getScoreString(self, dumpPDB=False):
		'''
		This method standardizes the way that scores are saved throughout the design and evaluation process.

		Input:
			dumpPDB 	<str>	If provided, will attempt to save a pdb with this filename
		
		Output:
			<str>				Tab separated string with the score-terms, in the same order as in self.headStr
		
		Sets:
			self.headStr		Header used in writing log files. Current format:
		'pdbFN\ttotal_score\tnumHB\tfa_rep\tintra\tfa_atr\tSASA\tRotamer\tclusterID\tFASTA\thydroHEn\n'
		'''
		# Set up the pose and the default scorefunction
		sf = create_score_function('beta_nov16')
		sf.set_weight(approximate_buried_unsat_penalty,1)
		sf.set_weight(fa_intra_rep, 0.2)
		# Load the pose from file only if a pose has not already been loaded and modified
		try:
			self.pose.sequence()
		except:
			self.pose = pose_from_pdb(self.pdbFile)     
		# Next three lines are there to make sure that hydrogen bonding energy is pairwise decomposed
		emo = EnergyMethodOptions()
		emo.hbond_options().decompose_bb_hb_into_pair_energies( True )
		sf.set_energy_method_options(emo)
		totalEnzEn = sf(self.pose)

		# Record the sequence
		self.FASTA = self.pose.sequence()
		hSet = hbonds.HBondSet()
		self.pose.update_residue_neighbors()
		hbonds.fill_hbond_set(self.pose,False,hSet)

		# Identify the number of the TSR
		resi = self.pose.pdb_info().pdb2pose(self.tsrCode,self.tsrPDBNum)
		resHSet = hSet.residue_hbonds(resi)
		self.hydroHEn = 0
		self.hydroCount = len(resHSet)
		for i in range(1, len(resHSet)+1):
			if resHSet[i].don_res() != resHSet[i].acc_res():
				self.hydroHEn += resHSet[i].energy()*resHSet[i].weight()     

		# Repulsive energy for TSR
		self.rep = self.pose.energies().residue_total_energies(resi)[fa_rep]
		# Intramolecular rep for TSR
		self.intra = self.pose.energies().residue_total_energies(resi)[fa_intra_rep]
		# Leonard-Jones attractive energy
		self.atr = self.pose.energies().residue_total_energies(resi)[fa_atr]
		# Total TSR energy
		self.resiEn = self.pose.energies().residue_total_energy(resi)
		
		# Some setup for SASA
		TSR_selector = ResidueNameSelector(self.fusligName)
		sasa_metric = sm.PerResidueSasaMetric()
		sasa_metric.set_residue_selector(TSR_selector)
		# Getting the SASA
		resi_sasa = sasa_metric.calculate(self.pose)
		# Formatting the SASA
		resi_sasa = sorted(resi_sasa.items(), key=operator.itemgetter(1), reverse=False)
		relist, salist = zip(*resi_sasa)
		self.SASA = salist[0]
		# Making and storing header string, in case that is also wanted
		self.headStr = 'pdbFN\ttotal_score\thydroCount\thydroHEn\ttotalEnzEnergy\tfa_rep\tintra\tfa_atr\tSASA\tRotamer\tclusterID\tFASTA\n'
		# Building the score string. This would be on one line of a text file in a log
		scoreStr = self.pdbFile+'\t'+str(self.resiEn)+'\t'+str(self.hydroCount)+'\t'+str(self.hydroHEn)+'\t'+str(totalEnzEn)+'\t'+str(self.rep)+'\t'+str(self.intra)+'\t'+str(self.atr)+'\t'+str(self.SASA)+'\t'+str(self.rotamer)+self.FASTA+'\n'
		if bool(dumpPDB):
			self.pose.dump_pdb(dumpPDB)
		return scoreStr	

	def runRelax(self,sf,pose,relBackBone = False, numRelax = 1,mMap = False,dumpTraj=False,dumpName='traj'):
		'''
		Used to run a FastRelax protocol.
		Inputs:
			sf 			<scfxn>		Scorefunction to use during the FastRelax
		Optional:
			relBackBone	<bool>		Allow backbone to move, default False
			numRelax	<int>		Number of times to repeat FastRelax protocol
			mMap 		<MoveMap>	MoveMap for FastRelax to use, will respect disabled chi
			dumpTraj 	<int>		If nonzero, dump the trajectory. Will set stride to match int passed in
		'''
		# Start setting up FastRelax
		relax = FastRelax()
		# Dump the trajectory if requested
		if dumpTraj:
			emo = EnergyMethodOptions()
			emo.dump_trajectory_stride(dumpTraj)
			emo.dump_trajectory_prefix(dumpName)
			sf.set_energy_method_options(emo)
			sf.set_weight(dump_trajectory,1)
		relax.set_scorefxn(sf)
		# If a moveMap was supplied, use that (including chi angles, unlike default)
		if mMap:
			relax.set_movemap(mMap)
			#relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
		# Otherwise, set up a default MoveMap
		else:
			movemap = MoveMap()
			movemap.set_bb(relBackBone)
			movemap.set_chi(True)
			relax.set_movemap(movemap)
		# Run the FastRelax the specified number of times
		for i in range(numRelax):
			relax.apply(pose)
			#print('passed relax')
		# undo the dump_traj setting in the sfxn
		if dumpTraj:
			sf.set_weight(dump_trajectory,0)
		return pose
	
	def runPac_Min(self, mMap, sf, pose,set_PackerFocus=False, Min_cycles=5, MC_cycles=20):
		# Restrict to the residues indicated in residueSelector
		if bool(set_PackerFocus):
			# Prevent other residues being repacked
			prevent = PreventRepackingRLT()
			restrict_to_focus = OperateOnResidueSubset(prevent, set_PackerFocus, True)
		
		# Creating TaskFactory and read in the above parameters
		tf = TaskFactory()

		tf.push_back(RestrictToRepacking()) # Restrict the focused residues to repack and turn off Design
		if bool(set_PackerFocus):
			tf.push_back(restrict_to_focus)
		#tf.create_task_and_apply_taskoperations(pose)
		packer = PackRotamersMover(sf)
		packer.task_factory(tf)
		mc = MonteCarlo(pose,sf,1.0)
		testMove = TrialMover(packer,mc)
		
		for i in range(MC_cycles):
			testMove.apply(pose)
			print('test')
		mc.recover_low(pose)

		for i in range(Min_cycles-1):
			min_mover = MinMover(mMap, sf, 'dfpmin', 0.01, True)
			min_mover.apply(pose)
		min_mover = MinMover(mMap, sf, 'dfpmin', 0.00001, True)
		min_mover.apply(pose)

		return pose
	
	def revertMut(self,mMap, sf, pose, resfile, Min_cycles=5, MC_cycles=20):
		tf = TaskFactory.create_packer_task(pose)
		parse_resfile(pose,tf,resfile)
		rmover = PackRotamersMover(sf,tf)
		rmover.apply(pose)

		return pose

	def RBTranslate(self, pose, step_size=500):
		# Count ligand No.
		ppdb2 = PandasPdb().read_pdb(self.ClvPDB)
		res_list = []
		self.lig_count = 0
		for i in range(len(ppdb2.df['HETATM'])): 
			if ppdb2.df['HETATM']['residue_number'][i] not in res_list: 
				res_list.append(ppdb2.df['HETATM']['residue_number'][i])
				self.lig_count += 1
		
		# lig_count is tell the mover the ligand to be translated
		
		trans_mover = RigidBodyTransMover(pose, self.lig_count) 
		trans_mover.step_size(step_size) # 500 Ang is the default
		trans_mover.apply(pose)
		return pose

	def get_whole_protein_resi_string(self, pose):
		resi_string = ''
		for i in range(1,pose.total_residue()-4):
			resi_string += str(i)+','
		return resi_string

	def get_sc_score(self, pose, ligand_resi_string, enzyme_resi_string):
		sc = ShapeComplementarityFilter()
		sc.residues1(str(ligand_resi_string))
		sc.residues2(enzyme_resi_string)
		sc_score = sc.score(pose)

		return sc_score

	def get_whole_protein_resi_string(self, pose):
		resi_string = ''
		for i in range(1,pose.total_residue()-4):
			resi_string += str(i)+','
		return resi_string

	def get_pocket_res(self, pose, pdbFile):
		sc_pocket_str = ''
		pocket_res = bio_dist_select(pdbFile)
		for i in range(len(pocket_res)):
			sc_pocket_str += str(pose.pdb_info().pdb2pose(pocket_res.iloc[i]['chain_id'],pocket_res.iloc[i]['residue_number']))+','
		sc_pocket_str = sc_pocket_str[:-4]
		return sc_pocket_str

	def runPFC(self, revertMut=False):

		# Load hpDesigned structure to be scored
		pose1_hpdesigned = pose_from_pdb(self.pdbFile)
		sf = create_score_function('beta_nov16')
		sf.set_weight(approximate_buried_unsat_penalty,1)
		sf.set_weight(fa_intra_rep,0.2)
		# Get FASTA sequence
		self.FASTA = pose1_hpdesigned.sequence()
		
		# Next three lines are there to make sure that hydrogen bonding energy is pairwise decomposed
		emo = EnergyMethodOptions()
		emo.hbond_options().decompose_bb_hb_into_pair_energies( True )
		sf.set_energy_method_options(emo)

		# Load in residue numbers indicated in resfile
		df = pd.read_csv('../'+self.PDBID+'/'+self.ResiList,delimiter=' ')
		resString = ''
		for i in range(len(df)): 
			resString +=(str(df.iloc[i]['resi'])+df.iloc[i]['chain'])+','
		resString = resString+str(self.tsrPDBNum)+self.tsrCode
		

		# Set up residue selectors and movemap
		resiSelector = ResidueIndexSelector(resString) # Select only designed residues
		mmf = MoveMapFactory()
		mmf.all_bb(False) # Turn off all backbone movement
		mmf.all_chi(False) # Turn off all sidechain movement
		mmf.add_chi_action(mm_enable, resiSelector)

		mm  = mmf.create_movemap_from_pose(pose1_hpdesigned)
		#print (mm)

		inital_score = sf(pose1_hpdesigned) 
		
		if revertMut:
			inital_score = sf(pose1_hpdesigned)
			# Revert certain mutations
			self.revertMut(mm,sf,pose1_hpdesigned,self.revert_resfile)
		# Repack-MinMover

		self.runPac_Min(mm,sf,pose1_hpdesigned,set_PackerFocus=resiSelector) 
		
		fusion = pose1_hpdesigned.pdb_info().pdb2pose(self.tsrCode,self.tsrPDBNum)

		sf(pose1_hpdesigned)
		hSet = hbonds.HBondSet()
		pose1_hpdesigned.update_residue_neighbors()
		hbonds.fill_hbond_set(pose1_hpdesigned,False,hSet)

		resHSet = hSet.residue_hbonds(resi)
		self.hydroHEn = 0
		self.hydroCount = len(resHSet)
		for i in range(1, len(resHSet)+1):
			if resHSet[i].don_res() != resHSet[i].acc_res():
				self.hydroHEn += resHSet[i].energy()*resHSet[i].weight()     
		
		# Repulsive energy for TSR
		self.rep = pose1_hpdesigned.energies().residue_total_energies(fusion)[fa_rep]
		# Intramolecular rep energy for TSR
		self.intra = pose1_hpdesigned.energies().residue_total_energies(fusion)[fa_intra_rep]
		# Modified Leonard-Jones attractive energy
		self.atr = pose1_hpdesigned.energies().residue_total_energies(fusion)[fa_atr]
		# Coulombic electrostatic potental with a distance-dependent dielectric
		self.elec = pose1_hpdesigned.energies().residue_total_energies(fusion)[fa_elec]
		# Total TSR energy
		self.resiEn = pose1_hpdesigned.energies().residue_total_energy(fusion)    
		
		pose1_hpdesigned.dump_pdb(self.MinPDB)

		print('\n\n\npassed first minMover\n\n\n')
		# Score the bound state
		bound_score = sf(pose1_hpdesigned)

		print('\n\n\nBound state distribution for',self.pdbFile)
		sf.show(pose1_hpdesigned)


		print('\n\n\nEvaluation on original strucutre completed\n\n\n')
		#Cut fused ligand from the minimized crystal structure
		pose2 = self.cut_lig(self.MinPDB) # Cut and load the cleaved PDB
		resi = pose2.pdb_info().pdb2pose(self.tsrCode2,self.tsrPDBNum2) # Define residue No for cleaved ETE.

		# Get SASA for the cut substrate
		TSR_selector = ResidueNameSelector(self.cutligName) 
		sasa_metric = sm.PerResidueSasaMetric()
		sasa_metric.set_residue_selector(TSR_selector)
		resi_sasa = sasa_metric.calculate(pose2)
		resi_sasa = sorted(resi_sasa.items(), key=operator.itemgetter(1), reverse=False)
		relist, salist = zip(*resi_sasa)
		self.SASA = salist[0]

		# Get Shape Complementarity
		p_string = self.get_whole_protein_resi_string(pose2)
		p_string = p_string+'503C,'+str(self.tsrPDBNum)+self.tsrCode
		lig_string = str(self.tsrPDBNum2)+self.tsrCode2
		#self.SC = self.get_sc_score(pose2,lig_string,p_string)
		# Active Site SC
		self.active_sc = self.get_sc_score(pose2, ligand_resi_string=resi, enzyme_resi_string = self.get_pocket_res(pose2, self.pdbFile))

		self.cut_score = sf(pose2)
		# Generate pose for RMSD and relax
		#pose2 = pose_from_pdb(self.ClvPDB)  # Load in the cleaved pdb
		
		self.Bunsat = pose2.energies().residue_total_energies(resi)[approximate_buried_unsat_penalty]

		
		# copy pose for RMSD
		rmsd_pose1_hp_min = Pose()
		rmsd_pose1_hp_min.assign(pose2)
		# Delete substrate for the sake of rmsd
		rmsd_pose1_hp_min.delete_polymer_residue(resi)
		
		#------------------------------------------- Fidelity test below-------------------------------------
		# Relax test (Cut and relax test)
		pose_4_rlx = Pose()
		pose_4_rlx.assign(pose2)
		self.runRelax(sf,pose_4_rlx)
		rlxString = resString+','+str(self.tsrPDBNum)+self.tsrCode+','+str(self.tsrPDBNum2)+self.tsrCode2

		# Set up residue selectors and movemap
		rlxSelector = ResidueIndexSelector(rlxString)

		mmf_rlx = MoveMapFactory()
		mmf_rlx.all_bb(False) # Turn off all backbone movement
		mmf_rlx.all_chi(False) # Turn off all sidechain movement
		mmf_rlx.add_chi_action(mm_enable, rlxSelector)

		mm_rlx  = mmf_rlx.create_movemap_from_pose(pose_4_rlx)
		self.runPac_Min(mm_rlx,sf,pose_4_rlx,set_PackerFocus=rlxSelector)

		self.rlx_test_score =sf(pose_4_rlx)

		ETE_selector = ResidueNameSelector(self.cutligName)
		rmsd_ETE = RMSDMetric()
		rmsd_ETE.set_residue_selector(ETE_selector) 
		rmsd_ETE.set_residue_selector_reference(ETE_selector)
		rmsd_ETE.set_comparison_pose(pose_4_rlx)
		self.RMSD_ETE = rmsd_ETE.calculate(pose2)
		# Save structure for relax test	
		pose_4_rlx.dump_pdb(self.Rlx_test)
		#------------------------------------------- End Fidelity test-------------------------------------


		#Setting up for translation
		trans_pose = Pose()
		trans_pose.assign(pose2)
		self.RBTranslate(trans_pose)

				# Score the translated-prerelax state
		self.translated = sf(trans_pose)

		print('\n\n\nTranslated-prerelax distribution for',self.pdbFile[-24:])
		sf.show(trans_pose)

		# Save the translated pose
		trans_pose.dump_pdb(self.TransPDB)
		
		# Set up movemap for trans_pose
		mm  = mmf.create_movemap_from_pose(trans_pose)
		#print(mm)

		# Initiate fix backbone relax
		#print('\n\n\nInitiating fix backbone relax!!!\n\n\n')
		self.runRelax(sf,trans_pose,mMap=mm)
		
		
		# Packer_MinMover
		self.runPac_Min(mm, sf, trans_pose, set_PackerFocus=resiSelector)
		


		print('\n\n\npassed the second minMover\n\n\n')
		# Score the translated-relaxed state
		self.unbound_score = sf(trans_pose)

		print('\n\n\nTranslated-postrelax distribution for',self.pdbFile[-24:])
		#sf.show(trans_pose)
		
		# Save relaxed structrue
		trans_pose.dump_pdb(self.TRlxPDB)

		# Remove substrate in pose2 
		rmsd_pose2 = trans_pose
		rmsd_pose2.delete_polymer_residue(resi)
		
		print('\n\n\nstart calculating IE and rmsd\n\n\n')
		# Calculate IE and sc_rmsd
		self.IE1_2 = bound_score - self.translated
		self.IE1_3 = bound_score - self.unbound_score
		#self.sc_rmsd = all_scatom_rmsd_nosuper(rmsd_pose1_hp_min, rmsd_pose2)
		rmsd = RMSDMetric()
		rmsd.set_residue_selector(resiSelector) 
		rmsd.set_residue_selector_reference(resiSelector)
		rmsd.set_comparison_pose(rmsd_pose2)
		self.RMSD_BS = rmsd.calculate(rmsd_pose1_hp_min)
		print('\n\n\n\nPreform check completed\n\n\n')

		# Saving output logs
		self.headers = 'pdbFN_CUT\tpdbFN_MIN\tinitial_score\tbound(minimized)\tcut_score\trlx_test_score\ttranslated(trans-preMin)\tunbound(trans-Min)\tnumHB\ttotal_hEn\tSASA\tSC\tBunsat\tfa_rep\tintra\tfa_atr\telec\tIE1_2\tIE1_3\tRMSD_BS\tRMSD_ETE\tRotamer\tFASTA\n'
		self.scores = self.ClvPDB+'\t'+self.MinPDB+'\t'+str(inital_score)+'\t'+str(bound_score)+'\t'+str(self.cut_score)+'\t'+str(self.rlx_test_score)+'\t'+str(self.translated)+'\t'+str(self.unbound_score)+'\t'+str(self.hydroCount)+'\t'+str(self.hydroHEn)+'\t'+str(self.SASA)+'\t'+str(self.active_sc)+'\t'+str(self.Bunsat)+'\t'+str(self.rep)+'\t'+str(self.intra)+'\t'+str(self.atr)+'\t'+str(self.elec)+'\t'+str(self.IE1_2)+'\t'+str(self.IE1_3)+'\t'+str(self.RMSD_BS)+'\t'+str(self.RMSD_ETE)+'\t'+str(self.rotamer)+'\t'+self.FASTA+'\n'

		
		f = open(self.TRlx_log,"w+")
		f.write(self.headers)
		f.write(self.scores)
		f.close()
		print('\n\n\n\nLog file for Preform check',self.pdbFile[-24:],'is saved\n\n\n')
		
		return


	
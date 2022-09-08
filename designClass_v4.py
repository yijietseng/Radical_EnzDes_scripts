from pyrosetta import init, Pose, MoveMap,create_score_function
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta.core.scoring import fa_rep,hbond_sc,hbond_bb_sc,atom_pair_constraint,dihedral_constraint,hbonds,fa_atr, fa_intra_rep,dump_trajectory, calpha_superimpose_pose, fa_elec
from pyrosetta.rosetta.core.scoring.methods import EnergyMethodOptions
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.select.residue_selector import ResidueNameSelector,ResidueIndexSelector
import pyrosetta.rosetta.core.simple_metrics.per_residue_metrics as sm
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory, mm_enable
from pyrosetta.rosetta.core.pack.task.operation import OperateOnResidueSubset, RestrictToRepacking,PreventRepackingRLT,RestrictToRepackingRLT,ReadResfile
from pyrosetta.rosetta.core.pack.task import TaskFactory, parse_resfile,operation
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.denovo_design.movers import FastDesign
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover,MinMover
from pyrosetta.rosetta.protocols.constraint_movers import ConstraintSetMover
from pyrosetta.rosetta.protocols.rigid import RigidBodyTransMover
from pyrosetta.rosetta.protocols.moves import SequenceMover,TrialMover,MonteCarlo
from pyrosetta.toolbox.mutants import mutate_residue
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, ResidueNameSelector
from pandas import read_csv
import operator, os

from pyrosetta.rosetta.std import vector_std_string
import pickle as pkl
from copy import deepcopy


#init(extra_options='-dunbrack_prob_buried 1.0 -dunbrack_prob_nonburied 1.0 -dunbrack_prob_buried_semi 1.0')
init('-load_PDB_components False','-beta_nov16 True')

class enzymeRot():
	'''
	This class holds information about enzymes and holds methods to score, mutate, or design them.

	Expected input:
		rotamer= int, refers to the rotamer of the TSR
		**kwargs: should contain dereferenced dictionary of variables, as in params.json
	Additional optional arguments:
		designStage:	refers to stage of design you are in, is used in setting the pdb file, defaults to 'init_PS'
		designSuffix:	defaults to '.pdb'
		glyRes:			boolean for whether or not to use the glycine resfile, defaults to False
		dRes:			boolean for whether or not to use the design resfile, defaults to False
	Methods:
		setPDB
		getScoreString
		getTSRMM
		addConst
		runRelax
		gen_gly_pose
	
	Last updated by Josh on 11-8-2021 at 11:43
	'''
	def __init__(self,rotamer=1,designStage='init_PS/',designSuffix='.pdb',fileIncre="01",glyRes=False,dRes=False,**kwargs):
		self.__dict__.update(kwargs)
		self.rotamer = rotamer
		self.designStage = designStage
		self.designSuffix = designSuffix
		# Getting the design number
		self.file_incre = fileIncre
		self.cnst = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/RC/'+self.PDBID+'.cst'
		if dRes:
			self.dRes = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/RC/'+self.PDBID+'.resfile'
		if glyRes:
			self.gres = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/RC/'+self.PDBID+'_GLY.resfile'
		# Try to automatically set self.pdbFile, warn if it fails
		try:
			self.setPDB()
		except:
			print('WARNING: PDBFILE NOT AUTOMATICALLY SET')

	def setPDB(self):
		'''This function sets self.pdbFile according to the current member variables, following the current file tree'''
		self.pdbFile = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.designStage+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix
		print('PDB file set to: '+self.pdbFile)
		return

	def getScoreString(self, dumpPDB=False,isRepeat=False,givenSF=False):
		'''
		This method standardizes the way that scores are saved throughout the design and evaluation process.

		Input:
			dumpPDB 	<str>	If provided, will attempt to save a pdb with this filename
		
		Output:
			<str>				Tab separated string with the score-terms, in the same order as in self.headStr
		
		Sets:
			self.headStr		Header used in writing log files. Current format:
		'pdbFN\ttotal_score\tnumHB\tfa_rep\tintra\tfa_atr\tSASA\t`Rotamer`\tclusterID\tFASTA\thydroHEn\n'
		'''
		# Setup default score function, then if one is given, match that
		sf = create_score_function('beta_nov16')
		if givenSF != False:
			sf.assign(givenSF)
		# Load the pose from file only if a pose has not already been loaded and modified
		try:
			self.pose.sequence()
		except:
			self.pose = pose_from_pdb(self.pdbFile)
		# Next three lines are there to make sure that hydrogen bonding energy is pairwise decomposed
		emo = EnergyMethodOptions()
		emo.hbond_options().decompose_bb_hb_into_pair_energies(True)
		sf.set_energy_method_options(emo)
		
		# Record the sequence
		self.FASTA = self.pose.sequence()
		# Identify the number of the TSR
		resi = self.pose.pdb_info().pdb2pose(self.tsrCode,self.tsrPDBNum)
		
		# Starting the strings for the score and header if the first run
		if not isRepeat:
			self.headStr = 'pdbFN\tFASTA\tRotamer'
			self.scoreStr = self.pdbFile+'\t'+self.FASTA+'\t'+"{:03d}".format(self.rotamer)

		# Adding total Enzyme Energy
		self.scoreStr+='\t'+str(sf(self.pose))
		if bool(isRepeat):
			self.headStr += '\t'+isRepeat+'totalEnzEn'
		else:
			self.headStr+='\t'+'totalEnzEn'

		'''		
		# Storing all the terms for the whole enzyme
		for scTerm in (sf.get_nonzero_weighted_scoretypes()):
			termStr = str(scTerm)[10:]
			self.scoreStr += '\t'+str(self.pose.energies().total_energies_array()[termStr][0]*sf.get_weight(scTerm))
			if bool(isRepeat):
				termStr = isRepeat+termStr
			self.headStr += '\t'+termStr
		'''
		# Setting up hbond information
		hSet = hbonds.HBondSet()
		self.pose.update_residue_neighbors()
		hbonds.fill_hbond_set(self.pose,False,hSet)

		resHSet = hSet.residue_hbonds(resi)
		# Loops through the hbonds and adds energy contribution if it involves one of the hydroxyls
		hydroHEn = 0
		hydroCount = len(resHSet)

		# Calculating the intra-hbonding score
		for i in range(1, len(resHSet)+1):
			if resHSet[i].don_res() != resHSet[i].acc_res():
				hydroHEn += resHSet[i].energy()*resHSet[i].weight()
		self.scoreStr+='\t'+str(hydroHEn)
		if bool(isRepeat):
			self.headStr += '\t'+isRepeat+'reshydroHEn'
		else:
			self.headStr+='\t'+'reshydroHEn'
		
		# Calculate SASA
		if bool(isRepeat):
			self.headStr += '\t'+isRepeat+'resSASA'
		else:
			self.headStr+='\t'+'resSASA'
	
		TSR_selector = ResidueNameSelector(self.fusligName) 
		sasa_metric = sm.PerResidueSasaMetric()
		sasa_metric.set_residue_selector(TSR_selector)
		# Getting the SASA
		resi_sasa = sasa_metric.calculate(self.pose)
		# Formatting the SASA
		resi_sasa = sorted(resi_sasa.items(), key=operator.itemgetter(1), reverse=False)
		relist, salist = zip(*resi_sasa)
		self.SASA = salist[0]
		self.scoreStr+='\t'+str(self.SASA)

		# Count the number of hydroxyls with hbonds
		self.scoreStr+='\t'+str(hydroCount)
		if bool(isRepeat):
			self.headStr += '\t'+isRepeat+'resNumHB'
		else:
			self.headStr+='\t'+'resNumHB'
		
		# Repulsive energy for TSR
		self.fa_rep = self.pose.energies().residue_total_energies(resi)[fa_rep]
		self.scoreStr+='\t'+str(self.fa_rep)
		if bool(isRepeat):
			self.headStr += '\t'+isRepeat+'resFaRep'
		else:
			self.headStr+='\t'+'resFaRep'
		# Intramolecular rep for TSR

		self.scoreStr+='\t'+str(self.pose.energies().residue_total_energies(resi)[fa_intra_rep])
		if bool(isRepeat):
			self.headStr += '\t'+isRepeat+'resFaIntraRep'
		else:
			self.headStr+='\t'+'resFaIntraRep'

		# Leonard-Jones attractive energy for the residue
		self.scoreStr+='\t'+str(self.pose.energies().residue_total_energies(resi)[fa_atr])
		if bool(isRepeat):
			self.headStr += '\t'+isRepeat+'resFaAtr'
		else:
			self.headStr+='\t'+'resFaAtr'

		# Coulombic electrostatic potental with a distance-dependent dielectric
		self.scoreStr+='\t'+str(self.pose.energies().residue_total_energies(resi)[fa_elec])
		if bool(isRepeat):
			self.headStr += '\t'+isRepeat+'resElec'
		else:
			self.headStr+='\t'+'resElec'

		# Total TSR energy
		self.scoreStr+='\t'+str(self.pose.energies().residue_total_energy(resi))+'\n'
		if bool(isRepeat):
			self.headStr += '\t'+isRepeat+'resTotal\n'
		else:
			self.headStr+='\t'+'resTotal\n'

		if bool(dumpPDB):
			self.pose.dump_pdb(dumpPDB)


			
		return self.scoreStr
	
	def getTSRMM(self):
		'''Will return a MoveMap that only allows the TSR to be moved based on default member variables'''
		resi = self.pose.pdb_info().pdb2pose(self.tsrCode,self.tsrPDBNum)
		mMap = MoveMap()
		mMap.set_bb(False)
		mMap.set_chi(False)
		mMap.set_chi(resi,True)
		mMap.set_bb(resi,True)
		return mMap
	
	def addConst(self,sf):
		'''
		Sets up constraints using the file previously attached to the object, and sets weights for the score function that is passed in.
		Input:
			sf 		<scfxn>		Scorefunction to set nonzero weights on constraints
		'''
		cnst=ConstraintSetMover()
		cnst.constraint_file(self.cnst)
		cnst.add_constraints(True)
		cnst.apply(self.pose)		
		sf.set_weight(atom_pair_constraint,self.distConst)
		sf.set_weight(dihedral_constraint,self.angleConst)
	
	def runRelax(self,sf,relBackBone = False, numRelax = 1,mMap = False,dumpTraj=False,dumpName='traj'):
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
			relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
		# Otherwise, set up a default MoveMap
		else:
			movemap = MoveMap()
			movemap.set_bb(relBackBone)
			movemap.set_chi(True)
			relax.set_movemap(movemap)
		# Run the FastRelax the specified number of times
		for i in range(numRelax):
			relax.apply(self.pose)
			print('run relax')
		# undo the dump_traj setting in the sfxn
		if dumpTraj:
			sf.set_weight(dump_trajectory,0)
		return
	
	def runPac_Min(self, mMap, sf,set_PackerFocus=False, MC_cycles=10):
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
		packer = PackRotamersMover(sf)
		packer.task_factory(tf)
		mc = MonteCarlo(self.pose,sf,1.0)
		min_mover = MinMover(mMap, sf, 'dfpmin', 0.01, True)
		totalMove = SequenceMover(packer,min_mover)
		testMove = TrialMover(totalMove,mc)
		for i in range(MC_cycles):
			testMove.apply(self.pose)
		mc.recover_low(self.pose)
		return

	def gen_gly_pose(self,relBackBone = False):	
		'''
		*******************************************************************************************************************
		This function is used to check the clashing. First it will mutate the identified residues into gly. It will then do
		a fixed-backbone relax and measure the fa_rep. It outputs a pdbfile of the glypose and a txt file with the scores.
		*******************************************************************************************************************
		'''	
		# Set up the pose
		self.pose = pose_from_pdb(self.pdbFile)
		# Set up the a default scorefunction
		sf = create_score_function('beta_nov16')
		# set up and apply the constraints
		self.addConst(sf)
		# Make it into a glypose
		gtask=TaskFactory.create_packer_task(self.pose)
		parse_resfile(self.pose,gtask,self.gres)
		gmover = PackRotamersMover(sf,gtask)
		gmover.apply(self.pose)
		self.addConst(sf)
		sf.set_weight(fa_intra_rep,0.55)
		# Set up and run relax, currently fixed backbone, consider relBackBone=True
		self.runRelax(sf,relBackBone=True)
		# Relax just tsr
		resi = self.pose.pdb_info().pdb2pose(self.tsrCode,self.tsrPDBNum)
		# Relax just tsr
		#mMap = self.getTSRMM()
		#self.runRelax(sf,mMap = mMap)
		# Generate the new pdb file name
		self.designStage='ClashChecked/'
		self.designSuffix='_C.pdb'
		self.setPDB()
		# Score and store the pose
		scoreStr = self.getScoreString(dumpPDB=self.pdbFile)
		# Write the log file
		f = open('../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/ClashChecked/logs/'+self.PDBID+'_PS'+str(self.rotamer)+'_log'+self.file_incre+'.txt','w+')
		f.write(self.headStr)
		f.write(scoreStr)
		f.close()
		return

	def runDesign(self, designing=True, setPose=True):
		'''Should run the appropriate design based on the design stage member variable'''

		# Define function to load pickle files
		def load(fileName):
			fileObject2 = open(fileName, 'rb')
			modelInput = pkl.load(fileObject2)
			fileObject2.close()
			return modelInput

		# Load pose if not already loaded
		sf = create_score_function('beta_nov16')
		# set weights appropriately (set in the relax script)
		if self.designStage == 'ClashChecked/':
			designRun = self.hbDesignRuns
			designStr = load(self.hbRelaxScript)
			currentstage = 'hbDesign'
			self.outputPDBBase = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/hbDesigned/'+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix
			svLog = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/hbDesigned/'+'logs/'+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix[:-4]+'.txt'
		elif self.designStage == 'hbDesigned/':
			designRun = self.hpDesignRuns
			designStr = load(self.hpRelaxScript)
			currentstage = 'hpDesign'
			self.outputPDBBase = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/hpDesigned/'+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix
			svLog = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/hpDesigned/logs/'+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix[:-4]+'.txt'
		
		if setPose:
			self.pose = pose_from_pdb(self.pdbFile)
		# Set up the a default scorefunction
		sf = create_score_function('beta_nov16')
		self.addConst(sf)
		self.getScoreString()

		f = open(svLog,"w+")
		f.write(self.headStr)
		f.close()
		basePose = self.pose.clone()
		if not designing:
			lowScore = sf(self.pose)
		
		self.designSuffix_ori = self.designSuffix # To set original suffix
		for i in range(1,designRun+1):
			# Design a copy of the pose
			self.pose.assign(basePose)

			mmp = MoveMap()
			mmp.set_bb(True)
			mmp.set_chi(True)

			design = FastRelax()
			if designing:
				design.set_enable_design(True)
				task_design = TaskFactory()
				read = ReadResfile(self.dRes)
				task_design.push_back(read)
				design.set_task_factory(task_design)
			# If a custom ramping scheme is given, apply it here
			relax_script = vector_std_string(designStr)
			design.set_script_from_lines(relax_script)
			print('\n\n\n\nCustom relax script loaded:\n',designStr,'\n\n\n')
			design.set_scorefxn(sf)
			design.set_movemap(mmp)
			design.apply(self.pose)
			print('passed designmover')
			# Saving PDB
			print('\n\n\n\n', currentstage,' completed') 
			self.pose.dump_pdb(self.outputPDBBase[:-4]+"_{:02d}".format(i)+'.pdb')

			# Saving output logs
			print('\n\n\n\nDesign run',self.outputPDBBase[:-4]+"_{:02d}".format(i)+'.pdb','is saved')
				
			
			# For hbDesign
			if self.designStage == 'ClashChecked/':
				# Set new path for output scorestring
				self.designSuffix = self.designSuffix_ori[:-4]+"_{:02d}".format(i)+'.pdb'
				self.designStage = 'hbDesigned/'
				self.setPDB()

				score_string = self.getScoreString()
				f = open(svLog,"a+")
				f.write(score_string)
				f.close()
				print('\n\n\n\nOutput file recorded')
				
				# reset parameters for next run
				
				self.designStage = 'ClashChecked/'
			
			elif self.designStage == 'hbDesigned/':
				# Set new path for output scorestring
				self.designSuffix = self.designSuffix_ori[:-4]+"_{:02d}".format(i)+'.pdb'
				self.designStage = 'hpDesigned/'
				self.setPDB()

				score_string = self.getScoreString()
				print(score_string)
				f = open(svLog,"a+")
				f.write(score_string)
				f.close()
				print('\n\n\n\nOutput file recorded')
				
				# reset parameters for next run
				self.designStage = 'hbDesigned/'
			if not designing:
				if (sf(self.pose)<=lowScore):
					lowScore=sf(self.pose)
					lowPDB=deepcopy(self.pdbFile)
					return lowPDB

	def runDummyDesign(self):
		'''
		This function will carry out the dummy designs for the scaffold check. This function is a standard FastDesign with fixed backbones.
		'''

		if not os.path.exists('../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.ScaffoldCheckFolder):
			os.mkdir('../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.ScaffoldCheckFolder)
		if not os.path.exists('../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.ScaffoldCheckFolder+'/logs'):
			os.mkdir('../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.ScaffoldCheckFolder+'/logs')

		# Load in PDB
		self.pose = pose_from_pdb(self.pdbFile)

		# Set up the scorefunction
		sf = create_score_function('beta_nov16')
		emo = EnergyMethodOptions()
		emo.hbond_options().decompose_bb_hb_into_pair_energies(True)
		sf.set_energy_method_options(emo)
		sf(self.pose)

		# Find the TSR
		resi = self.pose.pdb_info().pdb2pose(self.tsrCode,self.tsrPDBNum)

		# Apply all the constraints
		self.addConst(sf)
		self.getScoreString()
		self.dummyRep = str(self.fa_rep)
		sf.set_weight(atom_pair_constraint,1)
		#sf.set_weight(dihedral_constraint,1)

		# Set up MoveMap and FastDesign mover
		mmp = MoveMap()
		mmp.set_bb(False)
		mmp.set_chi(True)

		# Reading resfile
		task_design = TaskFactory()
		read = ReadResfile(self.dRes)
		task_design.push_back(read)

		# Setting up and executing FastDesign Mover
		dummy = FastDesign(sf)
		dummy.set_movemap(mmp)
		dummy.set_task_factory(task_design)
		dummy.apply(self.pose)

		# Set up outputs
		suffix = self.pdbFile.rsplit('/',1)[1][5:-6]
		newPDB = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.ScaffoldCheckFolder+'/'+suffix+'.pdb'
		logFn ='../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.ScaffoldCheckFolder+'/logs/'+suffix+'.txt'
		self.pose.dump_pdb(newPDB)

		# Set up output file
		self.pdbFile = newPDB
		self.getScoreString()
		self.dummySASA = str(self.SASA)
		log_head = 'Gly_rep\tDes_SASA\n'
		
		f = open(logFn,'w+')
		f.write(log_head)
		f.write(self.dummyRep+'\t'+self.dummySASA+'\n')
		f.close()


class WT_validation():
	'''
	This class holds information about enzymes and holds methods to simulate the same design processes for WT structures.

	Expected input:
		rotamer= int, refers to the rotamer of the TSR
		**kwargs: should contain dereferenced dictionary of variables, as in params.json
	Additional optional arguments:
		designStage:	refers to stage of design you are in, is used in setting the pdb file, defaults to 'init_PS'
		designSuffix:	defaults to '.pdb'
		glyRes:			boolean for whether or not to use the glycine resfile, defaults to False
		dRes:			boolean for whether or not to use the design resfile, defaults to False
	Methods:
		setPDB
		getTSRMM
		addConst
		runRelax
	
	Last updated by Josh on 1-26-2022 at 11:59
	'''
	def __init__(self,rotamer=1,designStage='WT/',designSuffix='.pdb',fileIncre="01",glyRes=False,dRes=False,**kwargs):
		self.__dict__.update(kwargs)
		self.rotamer = rotamer
		self.designStage = designStage
		self.designSuffix = designSuffix
		# Getting the design number
		self.file_incre = fileIncre
		self.cnst = '../'+self.PDBID+'/CLUSTERW_'+self.file_incre+'/RC/'+self.PDBID+'.cst'
		if dRes:
			self.dRes = '../'+self.PDBID+'/CLUSTERW_'+self.file_incre+'/RC/'+self.PDBID+'.resfile'
		if glyRes:
			self.gres = '../'+self.PDBID+'/CLUSTERW_'+self.file_incre+'/RC/'+self.PDBID+'_GLY.resfile'
		self.output = '../'+self.PDBID+'/CLUSTERW_'+self.file_incre+'/outputstruct/'+self.PDBID
		self.saveFN = '../'+self.PDBID+'/CLUSTERW_'+self.file_incre+'/logs/'+self.PDBID+'_log.txt'
		
		# Try to automatically set self.pdbFile, warn if it fails

		try:
			self.setPDB()
		except:
			print('WARNING: PDBFILE NOT AUTOMATICALLY SET')

	def setPDB(self):
		'''This function sets self.pdbFile according to the current member variables, following the current file tree'''
		self.pdbFile = '../'+self.PDBID+'/CLUSTERX_'+self.file_incre+'/'+self.designStage+self.PDBID+'_PS'+"{:03d}".format(self.rotamer)+self.designSuffix
		print('PDB file set to: '+self.pdbFile)
		return

	def getTSRMM(self):
		'''Will return a MoveMap that only allows the TSR to be moved based on default member variables'''
		resi = self.pose.pdb_info().pdb2pose(self.tsrCode,self.tsrPDBNum)
		mMap = MoveMap()
		mMap.set_bb(False)
		mMap.set_chi(False)
		mMap.set_chi(resi,True)
		mMap.set_bb(resi,True)
		return mMap
	
	def addConst(self,sf):
		'''
		Sets up constraints using the file previously attached to the object, and sets weights for the score function that is passed in.
		Input:
			sf 		<scfxn>		Scorefunction to set nonzero weights on constraints
		'''
		cnst=ConstraintSetMover()
		cnst.constraint_file(self.cnst)
		cnst.add_constraints(True)
		cnst.apply(self.pose)		
		sf.set_weight(atom_pair_constraint,self.distConst)
		sf.set_weight(dihedral_constraint,self.angleConst)
	
	def runRelax(self,sf,relBackBone = False, numRelax = 1,mMap = False,dumpTraj=False,dumpName='traj'):
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
			relax.set_movemap_disables_packing_of_fixed_chi_positions(True)
		# Otherwise, set up a default MoveMap
		else:
			movemap = MoveMap()
			movemap.set_bb(relBackBone)
			movemap.set_chi(True)
			relax.set_movemap(movemap)
		# Run the FastRelax the specified number of times
		for i in range(numRelax):
			relax.apply(self.pose)
			print('run relax')
		# undo the dump_traj setting in the sfxn
		if dumpTraj:
			sf.set_weight(dump_trajectory,0)
		return
	
	def runPac_Min(self, mMap, sf,set_PackerFocus=False, MC_cycles=10):
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
		packer = PackRotamersMover(sf)
		packer.task_factory(tf)
		mc = MonteCarlo(self.pose,sf,1.0)
		min_mover = MinMover(mMap, sf, 'dfpmin', 0.01, True)
		totalMove = SequenceMover(packer,min_mover)
		testMove = TrialMover(totalMove,mc)
		for i in range(MC_cycles):
			testMove.apply(self.pose)
		mc.recover_low(self.pose)
		return
	
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

	def runProcesses(self):
		# Define function to load pickle files
		def load(fileName):
			fileObject2 = open(fileName, 'rb')
			modelInput = pkl.load(fileObject2)
			fileObject2.close()
			return modelInput

		# Load WT PDB
		self.pose = pose_from_pdb(self.pdbFile)

		# Set up the a default scorefunction
		#sf = get_fa_scorefxn()
		sf = create_score_function('beta_nov16')
		sf.set_weight(approximate_buried_unsat_penalty,1)
		sf.set_weight(fa_intra_rep,0.2)
		emo = EnergyMethodOptions()
		emo.hbond_options().decompose_bb_hb_into_pair_energies( True )
		sf.set_energy_method_options(emo)
		self.initial_score = sf(self.pose)
		print('\n\n\nInitiail score = ',self.initial_score,'\n\n\n')

		# Simulating CLS, HB, and HP process
		self.addConst(self.pose, sf, self.cnst)
		print('\n\n\nInitiating first flexible bb relax process\n\n\n')
		for i in range(3):
			self.runRelax(sf,relBackBone=True)
			self.pose.dump_pdb(self.output+'_RLX'+"{:02d}".format(i)+self.designSuffix)
	
		# Load in residue numbers indicated in resfile
		df = pd.read_csv('../'+self.PDBID+'/'+self.ResiList,delimiter=' ')
		resString = ''
		for i in range(len(df)): 
			resString +=(str(df.iloc[i]['resi'])+df.iloc[i]['chain'])+','
		resString = resString[:-1]
		sf(self.pose) # score before the selector for maximum performance

		# Set up residue selectors and movemap
		resiSelector = ResidueIndexSelector(resString) # Select only designed residues
		mmf = MoveMapFactory()
		mmf.all_bb(False) # Turn off all backbone movement
		mmf.all_chi(False) # Turn off all sidechain movement
		mmf.add_chi_action(mm_enable, resiSelector)

		mm  = mmf.create_movemap_from_pose(self.pose)

		# Initiate first Repack-Min to minimize BS residues
		pose_4_rlx = Pose()
		pose_4_rlx.assign(self.pose)
		Rlx_resiString = resString+',997Z'
		Rlx_resiSelector = ResidueIndexSelector(resString)
		self.runPac_Min(mm,sf,self.pose,set_PackerFocus=Rlx_resiSelector)
		print('\n\n\nFirst repack min=',sf(self.pose),'\n\n\n')
		Sub_selector = ResidueNameSelector(self.cutligName)
		
		# Initialte Relax test and calculate RMSD
		print('\n\n\ninitiating fix bb relax\n\n\n')

		self.runRelax(pose_4_rlx,sf, relBackBone=False)
		print('\n\n\nscore after relax=',sf(pose_4_rlx),'\n\n\n')
		self.RMSD_Sub = self.calRMSD(self.pose, pose_4_rlx, resiSelector=Sub_selector)
		
		self.bound = sf(self.pose)
		self.pose.dump_pdb(self.output+'_min'+self.designSuffix)

		resi = self.pose.pdb_info().pdb2pose(self.tsrCode2,self.tsrPDBNum2)

		# Calulate hbEn
		hSet = hbonds.HBondSet()
		self.pose.update_residue_neighbors()
		hbonds.fill_hbond_set(self.pose,False,hSet)

		resHSet = hSet.residue_hbonds(resi)
		self.total_hb = len(resHSet)

		self.hydroHEn = 0
		hydroCount = len(resHSet)
		# Calculating the intra-hbonding score
		for i in range(1, len(resHSet)+1):
			if resHSet[i].don_res() != resHSet[i].acc_res():
				self.hydroHEn += resHSet[i].energy()*resHSet[i].weight()

		# Repulsive energy for TSR
		self.rep = self.pose.energies().residue_total_energies(resi)[fa_rep]
		# Intramolecular rep energy for TSR
		self.intra = self.pose.energies().residue_total_energies(resi)[fa_intra_rep]
		# Modified Leonard-Jones attractive energy
		self.atr = self.pose.energies().residue_total_energies(resi)[fa_atr]
		# Total TSR energy
		self.resiEn = self.pose.energies().residue_total_energy(resi)
		# Calculating SASA
		TSR_selector = ResidueNameSelector(self.lig_name)
		self.SASA = self.get_SASA(self.pose,resiSelector=TSR_selector)
		# approximate_buried_unsat_penalty
		self.Bunsat = self.pose.energies().residue_total_energies(resi)[approximate_buried_unsat_penalty]


		# Initiating translation
		self.pose2 = Pose()
		self.pose2.assign(self.pose)
		self.RBTranslate(self.pose2)
		print('\n\n\nLigand translated\n\n\n')
		self.pose2.dump_pdb(self.output+'_trans'+self.designSuffix)
		self.unbound2 =sf(self.pose2)
		print('\n\n\nMinimizing\n\n\n')
		self.runPac_Min(mm,sf,self.pose2,set_PackerFocus=resiSelector)

		self.unbound3 = sf(self.pose2)

		self.pose2.dump_pdb(self.output+'_pfc'+self.designSuffix)

		#Calculating IE and RMSD
		self.IE1_2 = self.bound - self.unbound2
		self.IE1_3 = self.bound - self.unbound3
		self.RMDS_BS = self.calRMSD(self.pose,self.pose2,resiSelector=resiSelector)

		# Get Shape Complementarity
		p_string = self.get_whole_protein_resi_string(self.pose)
		p_string = p_string+'503C,'+self.tsrPDBNum+self.tsrCode
		lig_string = self.tsrPDBNum2+self.tsrCode2

		# Active Site SC
		self.active_sc = self.get_sc_score(self.pose, ligand_resi_string=resi, enzyme_resi_string = self.get_pocket_res(self.pose, self.pdbFile))



		# Set up output
		self.headStr = 'pdbFN\ttotalEn\tnumHB\thbEn\tSC\tSASA\tBunsat\tfa_rep\tintra\tfa_atr\tIE1_2\tIE1_3\tRMSD_BS\tRMSD_Sub\n'
		self.scoreStr = self.pdbFile+'\t'+str(self.resiEn)+'\t'+str(self.total_hb)+'\t'+str(self.hydroHEn)+'\t'+str(self.active_sc)+'\t'+str(self.SASA)+'\t'+str(self.Bunsat)+'\t'+str(self.rep)+'\t'+str(self.intra)+'\t'+str(self.atr)+'\t'+str(self.IE1_2)+'\t'+str(self.IE1_3)+'\t'+str(self.RMDS_BS)+'\t'+str(self.RMSD_Sub)+'\n'

		f = open(self.saveFN, "w+")
		f.write(self.headStr)
		f.write(self.scoreStr)
		f.close()
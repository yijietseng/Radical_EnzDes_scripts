'''
This Class is to run the same design process but on the WT crystal structure
'''
from pyrosetta import init, dump_pdb, get_fa_scorefxn, Pose, standard_packer_task, MoveMap, create_score_function
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory, mm_enable
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.scoring import atom_pair_constraint, fa_atr, fa_rep,hbonds, hbond_sc, hbond_bb_sc, cart_bonded,fa_intra_rep,approximate_buried_unsat_penalty
from pyrosetta.rosetta.protocols.rigid import RigidBodyTransMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.denovo_design.movers import FastDesign
from pyrosetta.rosetta.core.scoring.methods import EnergyMethodOptions
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, ResidueNameSelector
import pyrosetta.rosetta.core.simple_metrics.per_residue_metrics as sm
from pyrosetta.rosetta.protocols.minimization_packing import MinMover
from pyrosetta.rosetta.core.optimization import MinimizerOptions
import operator
from pyrosetta.rosetta.protocols.constraint_movers import ConstraintSetMover
from biopandas.pdb.pandas_pdb import PandasPdb
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric 
import pandas as pd
import sys, fileinput, re
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.core.pack.task import operation
from pyrosetta.rosetta.core.pack.task.operation import OperateOnResidueSubset, RestrictToRepacking,PreventRepackingRLT,RestrictToRepackingRLT, ReadResfile
from pyrosetta.rosetta.protocols.moves import TrialMover,MonteCarlo 
from pyrosetta.rosetta.std import vector_std_string
import pickle as pkl
from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
init('-load_PDB_components False','-beta_nov16 True')

class WT_validation():
	def __init__(self, designStage='WT/', designSuffix='.pdb', lig_name='25W', **kwargs):
		self.__dict__.update(kwargs)
		self.designStage = designStage
		self.designSuffix = designSuffix
		self.dRes = self.topPath+self.enzymeName+'/CLUSTERW'+'/RC/'+self.d_resfile
		self.lig_name = lig_name
		self.output = self.topPath+self.enzymeName+'/CLUSTERW/'+self.output_folder+self.enzymeName+'_'+self.clusterID
		self.saveFN = self.topPath+self.enzymeName+'/CLUSTERW'+'/logs/'+self.enzymeName+'_'+self.clusterID+'_log.txt'

		# Try to automatically set self.pdbFile, warn if it fails
		try:
			self.setPDB()
		except:
			print('WARNING: PDBFILE NOT AUTOMATICALLY SET')

	def setPDB(self):
		'''This function sets self.pdbFile according to the current member variables, following the current file tree'''
		self.pdbFile = self.topPath+self.enzymeName+'/CLUSTERW/'+self.designStage+self.enzymeName+'_'+self.clusterID+self.designSuffix
		print('PDB file set to: '+self.pdbFile)
		return

	def addConst(self, pose, sf, cnstfile):
		'''
		Sets up constraints using the file previously attached to the object, and sets weights for the score function that is passed in.
		Input:
			sf 		<scfxn>		Scorefunction to set nonzero weights on constraints
		'''
		cnst=ConstraintSetMover()
		cnst.constraint_file(cnstfile)
		cnst.add_constraints(True)
		cnst.apply(pose)		
		sf.set_weight(atom_pair_constraint,self.distConst)
		sf.set_weight(approximate_buried_unsat_penalty,1)
		sf.set_weight(fa_intra_rep,0.2)

	def runRelax(self,pose,sf,relBackBone = True, numRelax = 1,mMap = False,dumpTraj=False,dumpName='traj'):
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
			relax.apply(pose)
		# undo the dump_traj setting in the sfxn

		return pose

	def runFastDesign(self,pose, sf, designStr, designRun=1):

		for i in range(designRun):
			mmp = MoveMap()
			mmp.set_bb(True)
			mmp.set_chi(True)

			task_design = TaskFactory()
			read = ReadResfile(self.dRes)
			task_design.push_back(read)
			
			design = FastDesign()
			# If a custom ramping scheme is given, apply it here
			relax_script = vector_std_string(designStr)
			design.set_script_from_lines(relax_script)
			print('\n\n\n\nCustom relax script loaded:\n',designStr,'\n\n\n')
			design.set_scorefxn(sf)
			design.set_movemap(mmp)
			design.set_task_factory(task_design)
			design.apply(pose)

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
		mc.recover_low(pose)

		for i in range(Min_cycles-1):
			min_mover = MinMover(mMap, sf, 'dfpmin', 0.01, True)
			min_mover.apply(pose)
		min_mover = MinMover(mMap, sf, 'dfpmin', 0.00001, True)
		min_mover.apply(pose)

		return pose

	def RBTranslate(self, pose, step_size=500):
		# Count ligand No.
		ppdb2 = PandasPdb().read_pdb(self.pdbFile)
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

	def get_SASA(self, pose, resiSelector=False):
		sasa_metric = sm.PerResidueSasaMetric()
		if bool(resiSelector):
			sasa_metric.set_residue_selector(resiSelector)
		resi_sasa = sasa_metric.calculate(pose)
		resi_sasa = sorted(resi_sasa.items(), key=operator.itemgetter(1), reverse=False)
		relist, salist = zip(*resi_sasa)
		SASA = salist[0]

		return SASA

	def calRMSD(self, pose1, pose2, resiSelector=False):
		rmsd = RMSDMetric()
		if bool(resiSelector):
			rmsd.set_residue_selector(resiSelector) 
			rmsd.set_residue_selector_reference(resiSelector)
		rmsd.set_comparison_pose(pose2)
		RMSD = rmsd.calculate(pose1)

		return RMSD

	def get_sc_score(self, pose, ligand_resi_string, whole_protein_resi_string):
		sc = ShapeComplementarityFilter()
		sc.residues1(ligand_resi_string)
		sc.residues2(whole_protein_resi_string)
		sc_score = sc.score(pose)

		return sc_score

	def get_whole_protein_resi_string(self, pose):
		resi_string = ''
		for i in range(1,pose.total_residue()-4):
			resi_string += str(i)+','
		return resi_string

	def runProcesses(self):
		# Define function to load pickle files
		def load(fileName):
			fileObject2 = open(fileName, 'rb')
			modelInput = pkl.load(fileObject2)
			fileObject2.close()
			return modelInput

		# Load WT PDB
		pose = pose_from_pdb(self.pdbFile)

		# Set up the a default scorefunction
		#sf = get_fa_scorefxn()
		sf = create_score_function('beta_nov16')
		sf.set_weight(approximate_buried_unsat_penalty,1)
		sf.set_weight(fa_intra_rep,0.2)
		emo = EnergyMethodOptions()
		emo.hbond_options().decompose_bb_hb_into_pair_energies( True )
		sf.set_energy_method_options(emo)
		self.initial_score = sf(pose)
		print('\n\n\nInitiail score = ',self.initial_score,'\n\n\n')

		# Add initial bb constraint and initiating first flexible bb relax process
		self.addConst(pose, sf, self.cnst_init)
		print('\n\n\nInitiating first flexible bb relax process\n\n\n')
		self.runRelax(pose,sf,relBackBone=True)
		
		# Add des bb constraint and initiating hb and hp processes
		print('\n\n\nFlexible bb relax completed. Now initiating hbDesign process\n\n\n')
		self.addConst(pose, sf, self.cnst_des)
		designStr = load(self.hbRelaxScript)
		self.runFastDesign(pose, sf, designStr)
		pose.dump_pdb(self.output+'_hb'+self.designSuffix)

		print('\n\n\nhbDesign process completed. Now initiating hpDesign process\n\n\n')
		designStr = load(self.hpRelaxScript)
		self.runFastDesign(pose, sf, designStr)
		pose.dump_pdb(self.output+'_hp'+self.designSuffix)

		# Initiating post Design processes
		print('\n\n\nhpDesign process completed. Now initiating postDesign process\n\n\n')
		
		# Load in residue numbers indicated in resfile
		df = pd.read_csv(self.topPath+self.enzymeName+'/'+self.ResiList,delimiter=' ')
		resString = ''
		for i in range(len(df)): 
			resString +=(str(df.iloc[i]['resi'])+df.iloc[i]['chain'])+','
		resString = resString[:-1]
		sf(pose) # score before the selector for maximum performance

		# Set up residue selectors and movemap
		resiSelector = ResidueIndexSelector(resString) # Select only designed residues
		mmf = MoveMapFactory()
		mmf.all_bb(False) # Turn off all backbone movement
		mmf.all_chi(False) # Turn off all sidechain movement
		mmf.add_chi_action(mm_enable, resiSelector)

		mm  = mmf.create_movemap_from_pose(pose)

		# Initiate first Repack-Min to minimize BS residues
		pose_4_rlx = Pose()
		pose_4_rlx.assign(pose)
		Rlx_resiString = resString+',997Z'
		Rlx_resiSelector = ResidueIndexSelector(resString)
		self.runPac_Min(mm,sf,pose,set_PackerFocus=Rlx_resiSelector)
		print('\n\n\nFirst repack min=',sf(pose),'\n\n\n')
		Sub_selector = ResidueNameSelector('25W')
		
		# Initialte Relax test and calculate RMSD
		print('\n\n\ninitiating fix bb relax\n\n\n')
		'''
		rmmf = MoveMapFactory()
		rmmf.all_bb(False) # Turn off all backbone movement
		rmmf.all_chi(False) # Turn off all sidechain movement
		rmmf.add_chi_action(mm_enable, Rlx_resiSelector)
		rmmf.add_bb_action(mm_enable,Sub_selector)

		rmm  = rmmf.create_movemap_from_pose(pose)
		'''
		self.runRelax(pose_4_rlx,sf, relBackBone=False)
		print('\n\n\nscore after relax=',sf(pose_4_rlx),'\n\n\n')
		self.RMSD_Sub = self.calRMSD(pose, pose_4_rlx, resiSelector=Sub_selector)
		
		self.bound = sf(pose)
		pose.dump_pdb(self.output+'_min'+self.designSuffix)

		resi = pose.pdb_info().pdb2pose(self.tsrCode2,self.tsrPDBNum2)

		hSet = hbonds.HBondSet()
		pose.update_residue_neighbors()
		hbonds.fill_hbond_set(pose,False,hSet)

		resHSet = hSet.residue_hbonds(resi)
		self.total_hb = len(resHSet)

		
		# Checking hydrogen bonding energy for the TSR
		self.hEn = pose.energies().residue_total_energies(resi)[hbond_sc]+pose.energies().residue_total_energies(resi)[hbond_bb_sc]
		# Repulsive energy for TSR
		self.rep = pose.energies().residue_total_energies(resi)[fa_rep]
		# Intramolecular rep energy for TSR
		self.intra = pose.energies().residue_total_energies(resi)[fa_intra_rep]
		# Modified Leonard-Jones attractive energy
		self.atr = pose.energies().residue_total_energies(resi)[fa_atr]
		# Total TSR energy
		self.resiEn = pose.energies().residue_total_energy(resi)
		# Calculating SASA
		TSR_selector = ResidueNameSelector(self.lig_name)
		self.SASA = self.get_SASA(pose,resiSelector=TSR_selector)
		# approximate_buried_unsat_penalty
		self.Bunsat = pose.energies().residue_total_energies(resi)[approximate_buried_unsat_penalty]


		# Initiating translation
		pose2 = Pose()
		pose2.assign(pose)
		self.RBTranslate(pose2)
		print('\n\n\nLigand translated\n\n\n')
		pose2.dump_pdb(self.output+'_trans'+self.designSuffix)
		self.unbound2 =sf(pose2)
		print('\n\n\nMinimizing\n\n\n')
		self.runPac_Min(mm,sf,pose2,set_PackerFocus=resiSelector)

		self.unbound3 = sf(pose2)

		pose2.dump_pdb(self.output+'_pfc'+self.designSuffix)

		#Calculating IE and RMSD
		self.IE1_2 = self.bound - self.unbound2
		self.IE1_3 = self.bound - self.unbound3
		self.RMDS_BS = self.calRMSD(pose,pose2,resiSelector=resiSelector)

		# Get Shape Complementarity
		p_string = self.get_whole_protein_resi_string(pose)
		p_string = p_string+'503C,999X'
		lig_string = '997Z'
		self.SC = self.get_sc_score(pose,lig_string,p_string)


		# Set up output
		self.headStr = 'pdbFN\ttotalEn\tnumHB\thbEn\tSC\tSASA\tBunsat\tfa_rep\tintra\tfa_atr\tIE1_2\tIE1_3\tRMSD_BS\tRMSD_Sub\n'
		self.scoreStr = self.pdbFile+'\t'+str(self.resiEn)+'\t'+str(self.total_hb)+'\t'+str(self.hEn)+'\t'+str(self.SC)+'\t'+str(self.SASA)+'\t'+str(self.Bunsat)+'\t'+str(self.rep)+'\t'+str(self.intra)+'\t'+str(self.atr)+'\t'+str(self.IE1_2)+'\t'+str(self.IE1_3)+'\t'+str(self.RMDS_BS)+'\t'+str(self.RMSD_Sub)+'\n'

		f = open(self.saveFN, "w+")
		f.write(self.headStr)
		f.write(self.scoreStr)
		f.close()


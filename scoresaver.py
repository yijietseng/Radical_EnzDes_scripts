from pyrosetta import *
from pyrosetta.rosetta.core.scoring import fa_rep,hbond_sc,hbond_bb_sc,atom_pair_constraint,dihedral_constraint,hbonds,fa_atr, fa_intra_rep
from pyrosetta.rosetta.core.scoring.methods import EnergyMethodOptions
from pyrosetta.rosetta.core.select.residue_selector import ResidueNameSelector,ResidueIndexSelector
import pyrosetta.rosetta.core.simple_metrics.per_residue_metrics as sm
from pyrosetta.rosetta.core.pack.task.operation import OperateOnResidueSubset, RestrictToRepacking,PreventRepackingRLT,RestrictToRepackingRLT,ReadResfile
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, ResidueNameSelector
from pandas import read_csv
import operator, os, glob, json
from sys import argv
init('-load_PDB_components False','-beta_nov16 True')

# load parameter file
with open(argv[1]) as jsonfile:
	xParams = json.load(jsonfile)

stage = argv[2] # enter cls,hb,or hp
stage = stage.lower()

seqNo = int(argv[3])

# input a file_increment, or it will grab the folder with the largest increment to start with
if len(argv) > 4:
	file_incre = "{:02d}".format(int(argv[4]))
else:
	for i in itertools.count(start=1):
		if not os.path.exists('../'+xParams['PDBID']+'/CLUSTERX_'+"{:02d}".format(i)+'/'):
			file_incre = "{:02d}".format(i-1)
			break

# for cls the seqNo should be 1-375
if stage == 'cls':
	PDB = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/ClashChecked/'+xParams['PDBID']+'_PS'+"{:03d}".format(seqNo)+'_C.pdb'
	pose = pose_from_pdb(PDB)
	sf = create_score_function('beta_nov16')
	resi = pose.pdb_info().pdb2pose(xParams['tsrCode'],xParams['tsrPDBNum'])
	
	totalEnzEn = str(sf(pose))
	FASTA = pose.sequence()
	emo = EnergyMethodOptions()
	emo.hbond_options().decompose_bb_hb_into_pair_energies(True)
	sf.set_energy_method_options(emo)
	
	hSet = hbonds.HBondSet()
	pose.update_residue_neighbors()
	hbonds.fill_hbond_set(pose,False,hSet)

	resHSet = hSet.residue_hbonds(resi)
	# Loops through the hbonds and adds energy contribution if it involves one of the hydroxyls
	hydroHEn = 0
	hydroCount = len(resHSet)

	# Calculating the intra-hbonding score
	for i in range(1, len(resHSet)+1):
		if resHSet[i].don_res() != resHSet[i].acc_res():
			hydroHEn += resHSet[i].energy()*resHSet[i].weight()

	rep = str(pose.energies().residue_total_energies(resi)[fa_rep])
	intra = str(pose.energies().residue_total_energies(resi)[fa_intra_rep])
	atr = str(pose.energies().residue_total_energies(resi)[fa_atr])
	resTotal = str(pose.energies().residue_total_energy(resi))

	header = 'pdbFN\tFASTA\tRotamer\ttotalEnzEn\treshydroHEn\tresNumHB\tresFaRep\tresFaIntraRep\tresFaAtr\tresTotal\n'
	f = open('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/ClashChecked/logs/'+i[-16:-4]+'.txt', a)
	f.write(header)
	f.write(PDB+'\t'+FASTA+'\t'+PDB[-9:-6]+'\t'+totalEnzEn+'\t'+str(hydroHEn)+'\t'+str(hydroCount)+'\t'+rep+'\t'+intra+'\t'+atr+'\t'+resTotal+'\n')
	f.close()

# for hb the stage should be 1-375
if stage == 'hb':
	header = 'pdbFN\tFASTA\tRotamer\ttotalEnzEn\treshydroHEn\tresNumHB\tresFaRep\tresFaIntraRep\tresFaAtr\tresTotal\n'
	f = open('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hbDesigned/logs/'+xParams['PDBID']+'_PS'+"{:03d}".format(seqNo)+'_C.txt', 'a+')
	f.write(header)
	f.close()

	for j in range(1,11):
		PDB = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hbDesigned/'+xParams['PDBID']+'_PS'+"{:03d}".format(seqNo)+'_C_'+"{:02d}".format(j)+'.pdb'

		pose = pose_from_pdb(PDB)

		sf = create_score_function('beta_nov16')
		resi = pose.pdb_info().pdb2pose(xParams['tsrCode'],xParams['tsrPDBNum'])
		
				
		totalEnzEn = str(sf(pose))
		FASTA = pose.sequence()
		emo = EnergyMethodOptions()
		emo.hbond_options().decompose_bb_hb_into_pair_energies(True)
		sf.set_energy_method_options(emo)
		
		hSet = hbonds.HBondSet()
		pose.update_residue_neighbors()
		hbonds.fill_hbond_set(pose,False,hSet)

		resHSet = hSet.residue_hbonds(resi)
		# Loops through the hbonds and adds energy contribution if it involves one of the hydroxyls
		hydroHEn = 0
		hydroCount = len(resHSet)

		# Calculating the intra-hbonding score
		for i in range(1, len(resHSet)+1):
			if resHSet[i].don_res() != resHSet[i].acc_res():
				hydroHEn += resHSet[i].energy()*resHSet[i].weight()

		rep = str(pose.energies().residue_total_energies(resi)[fa_rep])
		intra = str(pose.energies().residue_total_energies(resi)[fa_intra_rep])
		atr = str(pose.energies().residue_total_energies(resi)[fa_atr])
		resTotal = str(pose.energies().residue_total_energy(resi))


		f = open('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hbDesigned/logs/'+xParams['PDBID']+'_PS'+"{:03d}".format(seqNo)+'_C.txt', 'a+')
		f.write(PDB+'\t'+FASTA+'\t'+PDB[-9:-6]+'\t'+totalEnzEn+'\t'+str(hydroHEn)+'\t'+str(hydroCount)+'\t'+rep+'\t'+intra+'\t'+atr+'\t'+resTotal+'\n')
		f.close()


# for hp the seqNo should be 1-375
if stage == 'hp':
	pdbList = glob.glob('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hpDesigned/*.pdb')

	for j in range(1,11):
		header = 'pdbFN\tFASTA\tRotamer\ttotalEnzEn\treshydroHEn\tresSASA\tresNumHB\tresFaRep\tresFaIntraRep\tresFaAtr\tresTotal\n'
		f = open('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hpDesigned/logs/'+xParams['PDBID']+'_PS'+"{:03d}".format(seqNo)+'_C_'+"{:02d}".format(j)+'.txt', 'a+')
		f.write(header)
		f.close()

		for k in range(1,11):
			PDB = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hpDesigned/'+xParams['PDBID']+'_PS'+"{:03d}".format(seqNo)+'_C_'+"{:02d}".format(j)+'_'+"{:02d}".format(k)+'.pdb'
			print(PDB)

			if PDB in pdbList:

				pose = pose_from_pdb(PDB)

				sf = create_score_function('beta_nov16')
				resi = pose.pdb_info().pdb2pose(xParams['tsrCode'],xParams['tsrPDBNum'])
				
						
				totalEnzEn = str(sf(pose))
				FASTA = pose.sequence()
				emo = EnergyMethodOptions()
				emo.hbond_options().decompose_bb_hb_into_pair_energies(True)
				sf.set_energy_method_options(emo)
				
				hSet = hbonds.HBondSet()
				pose.update_residue_neighbors()
				hbonds.fill_hbond_set(pose,False,hSet)

				resHSet = hSet.residue_hbonds(resi)
				# Loops through the hbonds and adds energy contribution if it involves one of the hydroxyls
				hydroHEn = 0
				hydroCount = len(resHSet)

				# Calculating the intra-hbonding score
				for i in range(1, len(resHSet)+1):
					if resHSet[i].don_res() != resHSet[i].acc_res():
						hydroHEn += resHSet[i].energy()*resHSet[i].weight()

				TSR_selector = ResidueNameSelector(xParams['fusligName']) 
				sasa_metric = sm.PerResidueSasaMetric()
				sasa_metric.set_residue_selector(TSR_selector)
				# Getting the SASA
				resi_sasa = sasa_metric.calculate(pose)
				# Formatting the SASA
				resi_sasa = sorted(resi_sasa.items(), key=operator.itemgetter(1), reverse=False)
				relist, salist = zip(*resi_sasa)
				SASA = salist[0]

				rep = str(pose.energies().residue_total_energies(resi)[fa_rep])
				intra = str(pose.energies().residue_total_energies(resi)[fa_intra_rep])
				atr = str(pose.energies().residue_total_energies(resi)[fa_atr])
				resTotal = str(pose.energies().residue_total_energy(resi))


				f = open('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hpDesigned/logs/'+xParams['PDBID']+'_PS'+"{:03d}".format(seqNo)+'_C_'+"{:02d}".format(j)+'.txt', 'a+')
				f.write(PDB+'\t'+FASTA+'\t'+PDB[-15:-12]+'\t'+totalEnzEn+'\t'+str(hydroHEn)+'\t'+str(SASA)+'\t'+str(hydroCount)+'\t'+rep+'\t'+intra+'\t'+atr+'\t'+resTotal+'\n')
				f.close()

			else:
				print('\n', PDB, 'not in pdblist')
				print('no structure found!!!!')
				
if stage == 'pfc':
	pdbList = glob.glob('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hpDesigned/*.pdb')

	for j in range(1,11):
		header = 'pdbFN\tFASTA\tRotamer\ttotalEnzEn\treshydroHEn\tresSASA\tresNumHB\tresFaRep\tresFaIntraRep\tresFaAtr\tresTotal\n'
		f = open('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hpDesigned/logs/'+xParams['PDBID']+'_PS'+"{:03d}".format(seqNo)+'_C_'+"{:02d}".format(j)+'.txt', 'a+')
		f.write(header)
		f.close()

		for k in range(1,11):
			PDB = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hpDesigned/'+xParams['PDBID']+'_PS'+"{:03d}".format(seqNo)+'_C_'+"{:02d}".format(j)+'_'+"{:02d}".format(k)+'.pdb'
			print(PDB)

			if PDB in pdbList:

				pose = pose_from_pdb(PDB)

				sf = create_score_function('beta_nov16')
				resi = pose.pdb_info().pdb2pose(xParams['tsrCode'],xParams['tsrPDBNum'])
				
						
				totalEnzEn = str(sf(pose))
				FASTA = pose.sequence()
				emo = EnergyMethodOptions()
				emo.hbond_options().decompose_bb_hb_into_pair_energies(True)
				sf.set_energy_method_options(emo)
				
				hSet = hbonds.HBondSet()
				pose.update_residue_neighbors()
				hbonds.fill_hbond_set(pose,False,hSet)

				resHSet = hSet.residue_hbonds(resi)
				# Loops through the hbonds and adds energy contribution if it involves one of the hydroxyls
				hydroHEn = 0
				hydroCount = len(resHSet)

				# Calculating the intra-hbonding score
				for i in range(1, len(resHSet)+1):
					if resHSet[i].don_res() != resHSet[i].acc_res():
						hydroHEn += resHSet[i].energy()*resHSet[i].weight()

				TSR_selector = ResidueNameSelector(xParams['fusligName']) 
				sasa_metric = sm.PerResidueSasaMetric()
				sasa_metric.set_residue_selector(TSR_selector)
				# Getting the SASA
				resi_sasa = sasa_metric.calculate(pose)
				# Formatting the SASA
				resi_sasa = sorted(resi_sasa.items(), key=operator.itemgetter(1), reverse=False)
				relist, salist = zip(*resi_sasa)
				SASA = salist[0]

				rep = str(pose.energies().residue_total_energies(resi)[fa_rep])
				intra = str(pose.energies().residue_total_energies(resi)[fa_intra_rep])
				atr = str(pose.energies().residue_total_energies(resi)[fa_atr])
				resTotal = str(pose.energies().residue_total_energy(resi))


				f = open('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/hpDesigned/logs/'+xParams['PDBID']+'_PS'+"{:03d}".format(seqNo)+'_C_'+"{:02d}".format(j)+'.txt', 'a+')
				f.write(PDB+'\t'+FASTA+'\t'+PDB[-15:-12]+'\t'+totalEnzEn+'\t'+str(hydroHEn)+'\t'+str(SASA)+'\t'+str(hydroCount)+'\t'+rep+'\t'+intra+'\t'+atr+'\t'+resTotal+'\n')
				f.close()

			else:
				print('\n', PDB, 'not in pdblist')
				print('no structure found!!!!')

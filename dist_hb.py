from pyrosetta import *
from pyrosetta.rosetta.core.scoring import atom_pair_constraint, fa_atr, fa_rep,hbonds, hbond_sc, hbond_bb_sc, cart_bonded,fa_intra_rep,approximate_buried_unsat_penalty
from pyrosetta.rosetta.core.scoring.methods import EnergyMethodOptions
from sys import argv
import pandas as pd
import os, json, itertools, glob
#init(extra_options='-dunbrack_prob_buried 1.0 -dunbrack_prob_nonburied 1.0 -dunbrack_prob_buried_semi 1.0')
init('-load_PDB_components False','-beta_nov16 True')


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

try:
	os.mkdir('../'+xParams['PDBID']+'/dist_hb')
except:
	print('Folder exist!!!!')

PDBList = glob.glob('../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/PFCheck/Rlx_test/*rlx.pdb')

df = pd.DataFrame(columns=['PDBFN','numHB','total_hEn',"dist_CZ_C5", 'FASTA'])

PDBFNlist = []
numHBlist = []
HBEnlist = []
distlst = []
seqlist = []

for i in PDBList:

	PDBFNlist.append(i)

	pose = pose_from_pdb(i)

	sf = create_score_function('beta_nov16')
	sf.set_weight(approximate_buried_unsat_penalty,1)
	sf.set_weight(fa_intra_rep, 0.2)

	emo = EnergyMethodOptions()
	emo.hbond_options().decompose_bb_hb_into_pair_energies( True )
	sf.set_energy_method_options(emo)

	totalEnzEn = sf(pose)
	FASTA = pose.sequence()
	seqlist.append(FASTA)

	hSet = hbonds.HBondSet()
	pose.update_residue_neighbors()
	hbonds.fill_hbond_set(pose,False,hSet)

	# Identify the number of the TSR
	resi = pose.pdb_info().pdb2pose(xParams['tsrCode2'],xParams['tsrPDBNum2'])
	resHSet = hSet.residue_hbonds(resi)
	hydroHEn = 0
	hydroCount = len(resHSet)
	for i in range(1, len(resHSet)+1):
		if resHSet[i].don_res() != resHSet[i].acc_res():
			hydroHEn += resHSet[i].energy()*resHSet[i].weight()    
	
	numHBlist.append(hydroCount)
	HBEnlist.append(hydroHEn)

	Lig_resi = pose.total_residue()
	AD_resi = Lig_resi - 1

	CZ_xyz = pose.residue(Lig_resi).xyz('CZ')
	C5_xyz = pose.residue(AD_resi).xyz("C5'")
	dist_vector = CZ_xyz - C5_xyz
	dist = round(dist_vector.norm(),2)

	distlst.append(dist)



df.PDBFN = PDBFNlist
df.numHB = numHBlist
df.total_hEn = HBEnlist
df.dist_CZ_C5 = distlst
df.FASTA = seqlist

result = '../'+xParams['PDBID']+'/dis_hb/'+file_incre+'_dist_hb.txt'

df.to_csv(result,sep='\t',index=False)


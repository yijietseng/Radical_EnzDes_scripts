'''
Make sure to change the ligand's residue number to 998 and chain to Z
'''

from pyrosetta import *
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.scoring import fa_rep,hbond_sc,hbond_bb_sc,hbonds,fa_atr, fa_intra_rep
from pyrosetta.rosetta.protocols.simple_filters import ShapeComplementarityFilter
from pyrosetta.rosetta.core.scoring.methods import EnergyMethodOptions
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, ResidueNameSelector
import pyrosetta.rosetta.core.simple_metrics.per_residue_metrics as sm
from pyro_scratch import bio_dist_select # This is a script that Talmage wrote
from sys import argv
import operator, itertools
init('-load_PDB_components False','-beta_nov16 True')



def get_whole_protein_resi_string(pose):
	resi_string = ''
	for i in range(1,pose.total_residue()-4):
		resi_string += str(i)+','
	return resi_string

def get_sc_score(pose, ligand_resi_string, enzyme_resi_string):
	sc = ShapeComplementarityFilter()
	sc.residues1(str(ligand_resi_string))
	sc.residues2(enzyme_resi_string)
	sc_score = sc.score(pose)

	return sc_score

def get_whole_protein_resi_string(pose):
	resi_string = ''
	for i in range(1,pose.total_residue()-4):
		resi_string += str(i)+','
	return resi_string

def get_pocket_res(pose, pdbFile):
	sc_pocket_str = ''
	pocket_res = bio_dist_select(pdbFile,resi=998)
	for i in range(len(pocket_res)):
		sc_pocket_str += str(pose.pdb_info().pdb2pose(pocket_res.iloc[i]['chain_id'],pocket_res.iloc[i]['residue_number']))+','
	sc_pocket_str = sc_pocket_str[:-4]
	return sc_pocket_str


#PDBID = argv[1]
Lig_name = input('\nPlease enter the 3-letter name of the ligand:')
Lig_name = Lig_name.upper()
#PDB = '../'+PDBID+'/PREDES/'+PDBID+'.pdb'
PDB = argv[1]
pose = pose_from_pdb(PDB)
sf = create_score_function('beta_nov16')
resi = pose.pdb_info().pdb2pose('Z',998)


movemap = MoveMap()
movemap.set_bb(True)
movemap.set_chi(True)

relax = FastRelax()
relax.set_movemap(movemap)
relax.set_scorefxn(sf)
relax.apply(pose)

'''
for i in itertools.count(start=1):
	if not os.path.exists('../'+PDBID+'/PREDES/'+PDBID+'_Rlx'+"{:02d}".format(i)+'pdb'):
		incre = "{:02d}".format(i)
		break
'''
#pose.dump_pdb('../'+PDBID+'/PREDES/'+PDBID+'_Rlx'+incre+'.pdb')
pose.dump_pdb('../4R34/CLUSTERX_01/top10/200_9_6_M363E_Rlx.pdb')

p_string = get_whole_protein_resi_string(pose)
lig_string = '998Z'

# Active Site SC
active_sc = round(get_sc_score(pose, ligand_resi_string=resi, enzyme_resi_string = get_pocket_res(pose, PDB)),2)

sf(pose)
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

TSR_selector = ResidueNameSelector(Lig_name) 
sasa_metric = sm.PerResidueSasaMetric()
sasa_metric.set_residue_selector(TSR_selector)
# Getting the SASA
resi_sasa = sasa_metric.calculate(pose)
# Formatting the SASA
resi_sasa = sorted(resi_sasa.items(), key=operator.itemgetter(1), reverse=False)
relist, salist = zip(*resi_sasa)
SASA = round(salist[0],2)

rep = str(round(pose.energies().residue_total_energies(resi)[fa_rep],2))
intra = str(round(pose.energies().residue_total_energies(resi)[fa_intra_rep],2))
atr = str(round(pose.energies().residue_total_energies(resi)[fa_atr],2))
resTotal = str(round(pose.energies().residue_total_energy(resi),2))

header = 'pdbFN\tFASTA\ttotalEnzEn\treshydroHEn\tresSASA\tresNumHB\tSC\tresFaRep\tresFaIntraRep\tresFaAtr\tresTotal\n'
#f = open('../'+PDBID+'/PREDES/'+PDBID+'Stats.txt', 'w+')
f = open('../4R34/CLUSTERX_01/top10/200_9_6_M363E_Rlx.txt','w+')
f.write(header)
f.write(PDB+'\t'+FASTA+'\t'+totalEnzEn+'\t'+str(round(hydroHEn,2))+'\t'+str(hydroCount)+'\t'+str(active_sc)+'\t'+rep+'\t'+intra+'\t'+atr+'\t'+resTotal+'\n')
f.close()
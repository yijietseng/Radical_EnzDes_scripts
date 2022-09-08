'''
Created and modified by Josh Tseng --- 2/24/2022

Description:
This script will take the structures that passed the fidility checks, prep them for autodock and execute autoduck. I have not yet figure out a way to automate the analyzing process. So far manual evaluation is required.  


Instructions:
This script will take the selected top structures (the fidelity tested structure) and prep them for autodock. The process flow is stated below:

1.Automatically identifying search space. The script will automatically identify search space by using the CA2 atom as the central coordinates and will box a 6x6x6 area by default as the search space. You can also run Box_search_area_visual.py in pymol to visuallize and modify the search area. 

2. Remove the ligand from the binding site and save it as RotXXX.pdb

3. Prep the RotXXX.pdbqt by running a script(prepare_ligand4.py, path stored in .bashrc or echo $prepLig) from mgltools, which assigns partial charges to all atoms.

4. Prep the receptor pdbqt by running a script (prepare_receptor4.py, path stored in .bashrc or echo $prepReceptor)from mgltools, which assigns partial charges to all atoms (this process can take a while based on protein size).

5. Modify SF4 partial charges. Because mgltools does not have a database of SF4. So after the run, the partial charges for SF4 is 0. So I mangaed a way to automatically replace the partial charges as shown below:
	S: -2
	Fe: +2.5
	There was a paper about it, but I forgot which one it was.

6. Run autodock.

7. Analyze data
'''

import pymol, os, json
from glob import glob
from sys import argv
from pymol import cmd



# loading parameter file
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



top_struct_folder = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/top_struct/'
AD_folder = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/top_struct/autodock/'
Receptor_folder = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/top_struct/autodock/receptor/'
Lig_folder = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/top_struct/autodock/Lig/'
output_folder = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/top_struct/autodock/output/'
log_foler = '../'+xParams['PDBID']+'/CLUSTERX_'+file_incre+'/top_struct/autodock/log/'


#---------------------------------------------------------------------------------------------------------------------------

def prep_input(PDB,lig_name,Receptor_folder=Receptor_folder,Lig_folder=Lig_folder):
	pythonsh = '/home/apps/autodock_vina/mgltools/bin/pythonsh '
	prepLig = '/home/apps/autodock_vina/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py '
	prepReceptor = '/home/apps/autodock_vina/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py '

	RotNo = 'Rot'+PDB[-19:-16]+'_'+PDB[-13:-10].strip('0')+PDB[-10:-8].strip('0')
	EnzShort = 'd'+file_incre.strip('0')+'_'+RotNo[3:]
	tempLigPDB = Lig_folder+RotNo+'.pdb'
	LigPDBQT = Lig_folder+RotNo+'.pdbqt'
	tempRecpPDB = Receptor_folder+EnzShort+'.pdb'
	RecpPDBQT = Receptor_folder+EnzShort+'.pdbqt'
	
	cmd.delete('*')
	cmd.load(PDB)
	#Take ligand and prep for ligand PDBQT file
	cmd.select(lig_name,'resn '+lig_name)
	cmd.save(tempLigPDB,lig_name)
	'''
	The prepare_ligand4.py cannot identify paths to ligand. So it is required to change the directory to the coresponding directory.
	'''
	os.chdir(Lig_folder)
	os.system(pythonsh+prepLig+'-l '+tempLigPDB)
	os.chdir('../../../../../scripts/')
	os.remove(tempLigPDB)
	print('\n'+LigPDBQT+' is generated!!!!')

	#Prep the protein, 5AD and SF4 into PDBQT
	cmd.select('Receptor','chain A+B+C+X')
	cmd.save(tempRecpPDB,'Receptor')
	os.system(pythonsh+prepReceptor+'-r '+tempRecpPDB+' -o '+RecpPDBQT+" -U waters")
	os.remove(tempRecpPDB)
	print('\n'+RecpPDBQT+' is generated!!!!')
	print('\nModifying partial charges for SF4!!!!')

	os.system('sed -i "s/     0.000 SA/    -2.000 SA/g" '+RecpPDBQT)
	os.system('sed -i "s/     0.000 Fe/     2.500 Fe/g" '+RecpPDBQT)

	print('\nPartial charges for S atoms in SF4 is set to: -2')
	print('Partial charges for Fe atoms in SF4 is set to +2.5\n\n')

	cmd.delete('*')

	return RecpPDBQT, LigPDBQT

def center_atom_coord(lig_name,atom='CA2'):
	'''
	By default it uses the CA atom of the ligand as the central atom. This can be varied according to different ligand
	'''
	cmd.select(lig_name,'resn '+lig_name+' and name '+atom)
	cen_coords = cmd.get_coords(lig_name,1)[0]
	print('Selected central coordinates is:', cen_coords)

	return cen_coords


def run_vina(RecpPDBQT, LigPDBQT, cen_coords,x_size=6,y_size=6,z_size=6, repeat=500):

	vina = '/home/apps/autodock_vina/bin/vina '
	cen_x = cen_coords[0]
	cen_y = cen_coords[1]
	cen_z = cen_coords[2]

	out = output_folder+LigPDBQT.rsplit('/',1)[-1]
	log = log_foler+LigPDBQT.rsplit('/',1)[-1][:-4]+'.txt'

	os.system(vina+'--receptor '+RecpPDBQT+' --ligand '+LigPDBQT+' --center_x '+str(cen_x)+' --center_y '+str(cen_y)+' --center_z '+str(cen_z)+' --size_x '+str(x_size)+' --size_y '+str(y_size)+' --size_z '+str(z_size)+' --out '+out+' --log '+log+' --exhaustiveness '+str(repeat))


# ------------------------------------------------------------------

try:
	os.mkdir(AD_folder)
	os.mkdir(Receptor_folder)
	os.mkdir(Lig_folder)
	os.mkdir(output_folder)
	os.mkdir(log_foler)
	print('Folders created')
except:
	print('Folders exist!!!!')

# Get PDB list
PDBlist = glob(top_struct_folder+'*rlx.pdb')

for i in PDBlist:
	cmd.delete('*')
	cmd.load(i)
	cen_coord = center_atom_coord(xParams['cutligName'])
	inputs = prep_input(i,xParams['cutligName'],Receptor_folder=Receptor_folder,Lig_folder=Lig_folder)
	run_vina(inputs[0], inputs[1], cen_coord,x_size=6,y_size=6,z_size=6)

# Need to figure out how to analyze the structures

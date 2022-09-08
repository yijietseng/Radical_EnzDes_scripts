import os

PDBIDLS = ["1R30","2A5H","1TV8","5FEW","3T7V",
           "4M7T","4NJK","4R34","4WCX","4U0P",
           "5EXK","5VSM","6FD2","6P78"]

# Make necessary folders
for i in PDBIDLS:
	if not os.path.exists('../'+i):
		os.mkdir('../'+i)

if not os.path.exists('./params4marylou/'):
	os.mkdir('./params4marylou/')

# Create param files
for i in PDBIDLS:
	with open('params.json',"rt") as fin:
		with open('./params4marylou/'+i+"_params.json","wt") as fout:
			lines = fin.readlines()
			lines[1] = lines[1].replace("Designs against CPF with proper charge", i+" scaffold check")
			lines[2] = lines[2].replace('4M7T',i)
			lines[3] = lines[3].replace('F5M_27.pdb', 'FCP.pdb')
			lines[6] = lines[6].replace('F5M','FCP')
			lines[9] = lines[9].replace('UNK','CPF')
			# change the torsion atom list
			lines[10] = lines[10].replace("O","N1'")
			lines[10] = lines[10].replace("C2",'CA')
			for ele in lines:
				fout.write(ele)

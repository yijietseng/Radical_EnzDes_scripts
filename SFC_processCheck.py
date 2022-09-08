import os, glob

PDBID = glob.glob('../????/')
PDBID.sort()

for i in PDBID:
    CLS = glob.glob(i+'CLUSTERX_01/Scaffold_Check/*.pdb')
    CLSPDB = len(CLS)
    print(i+':'+str(CLSPDB)+'('+str(round(CLSPDB/375*100,2))+')')

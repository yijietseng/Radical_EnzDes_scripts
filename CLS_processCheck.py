import os, glob

PDBID = glob.glob('../????/')

for i in PDBID:
    CLS = glob.glob(i+'CLUSTERX_01/ClashChecked/*.pdb')
    CLSPDB = len(CLS)
    print(i+':'+str(CLSPDB)+'('+str(round(CLSPDB/375*100,2))+')')

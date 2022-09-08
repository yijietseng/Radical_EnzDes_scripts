#!/bin/bash
:<<'comment'
There might be more artifacts in other structures, so manual addition of those
artifacts are sometimes required.
comment

# Identify how many arguments there are
# Input arguments are the 1. PDBID 2. original name for SAM or MET and 3. 5AD.
if [ "$#" -eq 2 ]; then
    PDBID=$1;
    MET=$2;
    AD=$2;
    #echo $PDBID ${2}
fi
if [ "$#" -eq 3 ]; then
    PDBID=$1;
    MET=$2;
    AD=$3;
    #echo $PDBID $MET $AD
fi

struct_PATH=../${PDBID}/PREDES/


sed -i "s/.$MET ...../ SAM X 999/g" "${struct_PATH}5AD_raw.pdb"
sed -i "s/.$AD ...../ SAM X 999/g" "${struct_PATH}5AD_raw.pdb"
#sed -i "s/SF4....../SF4 B 500/g" "${struct_PATH}5AD_raw.pdb"
grep -v ANISOU "${struct_PATH}5AD_raw.pdb" > "${struct_PATH}5AD.pdb"
sed -i "s/      D    /           /g" "${struct_PATH}5AD.pdb"
sed -i "s/      D   /          /g" "${struct_PATH}5AD.pdb"
sed -i "s/      E    /           /g" "${struct_PATH}5AD.pdb"
sed -i "s/      E   /          /g" "${struct_PATH}5AD.pdb"
sed -i "s/      C    /           /g" "${struct_PATH}5AD.pdb"
sed -i "s/      B    /           /g" "${struct_PATH}5AD.pdb"
sed -i "s/      B   /          /g" "${struct_PATH}5AD.pdb"
sed -i "s/      F    /           /g" "${struct_PATH}5AD.pdb"
sed -i "s/  G    /       /g" "${struct_PATH}5AD.pdb"
sed -i "s/      H    /           /g" "${struct_PATH}5AD.pdb"
sed -i "s/      I    /           /g" "${struct_PATH}5AD.pdb"
sed -i "s/      J    /           /g" "${struct_PATH}5AD.pdb"
sed -i "s/  J   /      /g"  "${struct_PATH}5AD.pdb"
sed -i "s/      K    /           /g" "${struct_PATH}5AD.pdb"
sed -i "s/SE / SD/ " "${struct_PATH}5AD.pdb"
sed -i "s/N   SAM X 999/N   MLF C 503/ " "${struct_PATH}5AD.pdb"
sed -i "s/CA  SAM X 999/CA  MLF C 503/ " "${struct_PATH}5AD.pdb"
sed -i "s/C   SAM X 999/C   MLF C 503/ " "${struct_PATH}5AD.pdb"
sed -i "s/O   SAM X 999/O   MLF C 503/ " "${struct_PATH}5AD.pdb"
sed -i "s/CB  SAM X 999/CB  MLF C 503/ " "${struct_PATH}5AD.pdb"
sed -i "s/CG  SAM X 999/CG  MLF C 503/ " "${struct_PATH}5AD.pdb"
sed -i "s/CE  SAM X 999/CE  MLF C 503/ " "${struct_PATH}5AD.pdb"
sed -i "s/OXT SAM X 999/OXT MLF C 503/ " "${struct_PATH}5AD.pdb"
sed -i "s/SD  SAM X 999/SD  MLF C 503/ " "${struct_PATH}5AD.pdb"
sed -i "s/C5' SAM X 999/C5' 5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/C4' SAM X 999/C4' 5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/O4' SAM X 999/O4' 5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/C3' SAM X 999/C3' 5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/O3' SAM X 999/O3' 5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/C2' SAM X 999/C2' 5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/O2' SAM X 999/O2' 5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/C1' SAM X 999/C1' 5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/N1  SAM X 999/N1  5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/C2  SAM X 999/C2A 5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/H2  SAM X 999/H2A 5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/N3  SAM X 999/N3  5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/C4  SAM X 999/C4  5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/C5  SAM X 999/C5  5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/C6  SAM X 999/C6  5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/N6  SAM X 999/N6  5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/N7  SAM X 999/N7  5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/C8  SAM X 999/C8  5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/N9  SAM X 999/N9  5AD X 999/ " "${struct_PATH}5AD.pdb"
sed -i "s/C4  SAM X 999/C4  5AD X 999/ " "${struct_PATH}5AD.pdb"
#sed -i "s/     D    /          / " "${struct_PATH}5AD.pdb"
#sed -i "s/     E  /        / " "${struct_PATH}5AD.pdb"




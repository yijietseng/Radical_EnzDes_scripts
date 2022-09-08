#!/bin/bash
:<<'comment'
This script will fix some of the name artifacts in the cluster center pdb files.
And will also create rotamer directories necessary for later process
comment

ID=${1} #raw pdb file
PDB=${ID:3:4} # extract PDBID

struct_PATH=../${PDB}/PREDES/

grep -v ANISOU "${ID}" |
grep -v BARG |
grep -v BVAL |
grep -v BALA |
grep -v BGLY |
grep -v BLEU |
grep -v BPRO |
grep -v BPHE |
grep -v BTRP |
grep -v BTYR |
grep -v BASP |
grep -v BGLU |
grep -v BARG |
grep -v BHIS |
grep -v BLYS |
grep -v BSER |
grep -v BTHR |
grep -v BCYS |
grep -v BMET |
grep -v BASN |
grep -v BGLN |
grep -v BILE > "${struct_PATH}${PDB}.pdb"

sed -i 's/AARG/ ARG/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/AVAL/ VAL/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/AALA/ ALA/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/AGLY/ GLY/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/ALEU/ LEU/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/APRO/ PRO/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/APHE/ PHE/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/ATRP/ TRP/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/ATYR/ TYR/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/AASP/ ASP/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/AGLU/ GLU/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/AARG/ ARG/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/AHIS/ HIS/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/ALYS/ LYS/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/ASER/ SER/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/ATHR/ THR/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/ACYS/ CYS/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/AMET/ MET/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/AASN/ ASN/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/AGLN/ GLN/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/AILE/ ILE/g' "${struct_PATH}${PDB}.pdb"
sed -i 's/  A    /       /g' "${struct_PATH}${PDB}.pdb"
sed -i 's/      A   /          /g' "${struct_PATH}${PDB}.pdb"
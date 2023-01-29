from spyrmsd import io, rmsd

import os
import glob
import csv


# Folder path and count the number of output ligands
path = r'/home/yuyang/temp02'
out_filenames = glob.glob(path + "/vina/Trial05/output" + "/*.pdb")
out_length = len(out_filenames)

print("Output ligands out: ", out_length)

# Load models
data, ref_filenames = [], []

for i in range(out_length):
	mols = io.loadallmols(out_filenames[i])
	ref_file = path + "/ligand/" + out_filenames[i][-12:-8] + "_ligand.pdb"
	ref = io.loadmol(ref_file)

	# Remove hydrogens
	ref.strip()

	for mol in mols:
		mol.strip()

	# Symmetry-corrected RMSD
	coords_ref = ref.coordinates
	anum_ref = ref.atomicnums
	adj_ref = ref.adjacency_matrix

	coords = [mol.coordinates for mol in mols]
	anum = mols[0].atomicnums
	adj = mols[0].adjacency_matrix

	print("Processing ", out_filenames[i][-12:-8])

	RMSD = rmsd.symmrmsd(coords_ref, coords, anum_ref, anum, adj_ref, adj)

	data.append(RMSD)

file = open('Trial05_rmsd_results.csv', 'a+', newline = '')

with file:
	write = csv.writer(file)
	write.writerows(data)


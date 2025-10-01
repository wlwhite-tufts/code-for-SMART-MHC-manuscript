import argparse
from Bio import PDB
from tqdm import tqdm

def get_args():

	parser = argparse.ArgumentParser(description="Align the input PDBs to the target by chain B and save the position on chain A as a new PDB.")
	parser.add_argument("-inputs", nargs='*', type=str, help="Filenames of pdbs to align and resave.")
	parser.add_argument("-target", type=str, help="Filenames of the target pdb.")
	parser.add_argument("-out_folder", type=str, help="Where to save the outputs.")

	return parser.parse_args()

if __name__ == '__main__':

	args = get_args()

	target_struct = PDB.PDBParser().get_structure('target', open(args.target,'r'))[0]
	target_struct = list(target_struct.get_chains())[0]

	for pdb in tqdm(args.inputs):

		input_struct = PDB.PDBParser().get_structure('target', open(pdb,'r'))[0]
		input_chains = list(input_struct.get_chains())

		target_atoms = []
		input_atoms = []
		for target_res, input_res in zip(target_struct, input_chains[1]) :
			assert target_res.resname == input_res.resname
			
			target_atoms.append(target_res['CA'])                
			input_atoms.append(input_res['CA'])

		super_imposer = PDB.Superimposer()
		super_imposer.set_atoms(target_atoms, input_atoms)
		super_imposer.apply(input_chains[0])

		io = PDB.PDBIO()
		io.set_structure(input_chains[0])
		io.save(args.out_folder + '/' + pdb.split('/')[-1])
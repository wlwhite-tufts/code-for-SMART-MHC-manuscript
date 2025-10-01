import argparse

def get_args():

	parser = argparse.ArgumentParser(description="Run a certain number of simulations with the given parameters.")
	parser.add_argument("-pdb", type=str, help="pdb file to modify")
	parser.add_argument("-length", type=int, help="offset length")
	
	return parser.parse_args()

if __name__ == '__main__':

	args = get_args()

	with open(args.pdb,'r') as f:
		lines = f.readlines()

	for line in lines[::-1]:
		if line.startswith('ATOM') or line.startswith('HETATM'):
			break_idx = int(line[23:27]) - args.length
			break

	with open(args.pdb,'w') as f:
		for line in lines:
			if  (line.startswith('ATOM') or line.startswith('HETATM') )and int(line[23:27]) > break_idx:
				line = line[:21] + 'B' + line[22:]
			f.write(line)
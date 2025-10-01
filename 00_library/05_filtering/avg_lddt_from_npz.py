import numpy as np
import pandas as pd
import argparse

def get_argparse():
	parser = argparse.ArgumentParser("Get summary stats from each provided .npz file in .csv format.")
	parser.add_argument('npzs', nargs='*', type=str, help='The files to aggregate.')
	parser.add_argument('-out_name',type=str,help='Name of output csv file.')

	return parser

if __name__ == '__main__':
	parser = get_argparse()
	args = parser.parse_args()

	scores = []
	for npz in args.npzs:
		name = npz.split('/')[-1].replace('.npz','')

		data = np.load(npz)
		lddts = data['lddt']
		data.close()

		design_mean = np.mean(lddts[:-177])
		MHC_mean = np.mean(lddts[-177:])
		design_min = np.min(lddts[:-177])
		MHC_min = np.min(lddts[-177:])

		scores.append([name,design_mean,design_min,MHC_mean,MHC_min])

	scores = pd.DataFrame(scores,columns=['name','design_mean','design_min','MHC_mean','MHC_min'])
	scores.to_csv(args.out_name)
from pyrosetta import *
import numpy as np
import argparse
import pandas as pd
import os

def get_args():
	parser = argparse.ArgumentParser(description="Evaluate loopability")
	parser.add_argument("-input", help="The input silent file", type=str)
	parser.add_argument("-tags", nargs='*', default=[], help="The name of the structure in the silent file to score.", type=str)
	parser.add_argument("-output", help="The output score file", type=str)
	parser.add_argument("-len", help="The length of the loop to test", type=int)
	parser.add_argument("-pattern", default='GGS' ,help="The linker pattern to repeat", type=str)
	parser.add_argument("-n_top", default=5 ,help="Loop checker will consider the top N outputs form Remodel", type=int)
	parser.add_argument('-save_pdbs', action='store_true', default=False)
	return parser.parse_args()

def get_poses_from_silent(file,tags):

	sfd_in = rosetta.core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
	sfd_in.read_file(file)

	if len(tags) == 0:
		tags = sfd_in.tags()
	
	for tag in tags:
		pose = Pose()
		sfd_in.get_structure(tag).fill_pose(pose)
		yield tag,pose

def make_blueprint(pose,length):

	#get ss info to use for residues on either side of the linker
	pose_aa = pose.sequence()

	#make blueprint for inserting the loop
	bp_list = []

	#chain A
	for i in range(pose.chain_end(1)-1):
		bp_list.append('{} {} .'.format(i+1,pose_aa[i]))

	#pre-loop residue has to be designable
	bp_list.append('{} {} L'.format(i+2,pose_aa[i+1]))

	#loop
	for ii in range(length):
		bp_list.append('0 x L')

	#post-loop residue has to be designable
	bp_list.append('{} {} L'.format(i+3,pose_aa[i+2]))

	#chain B
	for i in range(pose.chain_end(1)+1,pose.size()):
		bp_list.append('{} {} .'.format(i+1,pose_aa[i]))

	return '\n'.join(bp_list), pose.chain_end(1)

def relax_loop(pose,chA_end,length,pattern):
	full_name = {
	'A':'ALA',
	'C':'CYS',
	'D':'ASP',
	'E':'GLU',
	'F':'PHE',
	'G':'GLY',
	'H':'HIS',
	'I':'ILE',
	'K':'LYS',
	'L':'LEU',
	'M':'MET',
	'N':'ASN',
	'P':'PRO',
	'Q':'GLN',
	'R':'ARG',
	'S':'SER',
	'T':'THR',
	'V':'VAL',
	'W':'TRP',
	'Y':'TYR'
	}

	#select loopresidues
	loop_selector = rosetta.core.select.residue_selector.ResidueIndexSelector()
	loop_selector.set_index('{}-{}'.format(chA_end+1,chA_end+length))

	#select residues near loop
	near_loop_selector = rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
	near_loop_selector.set_distance(8)
	near_loop_selector.set_include_focus_in_subset(True)
	near_loop_selector.set_focus(loop_selector.apply(pose))

	#mutate residues to selected pattern
	selector = rosetta.core.select.residue_selector.ResidueIndexSelector()
	mutator = rosetta.protocols.simple_moves.MutateResidue()
	for i in range(chA_end+1,chA_end+length+1):
		selector.set_index(str(i))

		mutator.set_res_name(full_name[pattern[(i-chA_end-1)%len(pattern)]])
		mutator.set_selector(selector)

		mutator.apply(pose)

	#movemap to hold binder and MHC constant
	mm = rosetta.core.kinematics.MoveMap()
	mm.set_chi(near_loop_selector.apply(pose)) #allow loop and nearby to move sidechains
	mm.set_bb(loop_selector.apply(pose)) #allow only loop to move backbone

	#standard scorefunction, with chainbreak term turned on
	sfx = rosetta.core.scoring.ScoreFunction()
	sfx.add_weights_from_file("beta_nov16.wts")
	sfx.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.chainbreak,(1))

	#set foldtree so that only loop moves
	ft = rosetta.core.kinematics.FoldTree()
	ft.add_edge(1,chA_end,-1) #peptide edge for binder
	ft.add_edge(1,chA_end+1,1) #jump to attach linker to root of ft
	ft.add_edge(1,chA_end+length+1,2) #jump to attach MHC to root of ft
	ft.add_edge(chA_end+1,chA_end+length,-1) #peptide edge for linker
	ft.add_edge(chA_end+length+1,pose.size(),-1) #peptide edge for MHC
	pose.fold_tree(ft)

	# add chainbreak constraints to ends of linker to prevent ft from causing werid bond lengths/angles
	add_chainbreak = rosetta.protocols.protein_interface_design.movers.AddChainBreak()
	add_chainbreak.change_foldtree(False)
	add_chainbreak.resnum(str(chA_end))
	add_chainbreak.apply(pose)
	add_chainbreak.resnum(str(chA_end+length))
	add_chainbreak.apply(pose)

	pose.update_residue_neighbors()

	#standard minnmover
	min_mover = rosetta.protocols.minimization_packing.MinMover(mm,sfx,"lbfgs_armijo_nonmonotone",0.01,True,False,False)

	min_mover.apply(pose)
	return

def get_score(pose,chA_end,length):

	energy = pose.energies()
	residue_energies = np.array([energy.residue_total_energy(x) for x in range(chA_end-2,chA_end+length+4)])

	return max(residue_energies),np.mean(residue_energies)

if __name__ == '__main__':

	args = get_args()

	#flags modified from a script Hua gave me
	flags = ["-beta",
			 "-mute all",
			 "-in:file:silent_struct_type binary",

			 "-remodel::blueprint tmp.bp",
			 "-remodel::save_top {}".format(args.n_top),
			 "-generic_aa G",
			 "-remodel:quick_and_dirty",
			 "-remodel:design:no_design",
			 "-remodel:use_same_length_fragments false",
			 "-sample_over_loops"]
	init(" ".join(flags))

	with open(args.output,'w') as f:
		f.write('description,worst_res_linker,score_per_res_linker\n')

	for tag,in_pose in get_poses_from_silent(args.input,args.tags):

		#reset chains so remodel doesn't get confused
		in_pose.pdb_info().set_chains('A')

		#create real blueprint before remodel has to use it
		bp,chA_end = make_blueprint(in_pose,args.len)
		with open('tmp.bp','w') as f:
			f.write(bp)

		# #remodel the loop
		remodel = rosetta.protocols.forge.remodel.RemodelMover()
		remodel.apply(in_pose) #creates default-named output pdb files (cannot be changed with flags)

		#for each of the remodel outputs, mutate to the selected pattern, minimize and score
		scores = []
		for i in range(1,args.n_top+1):

			pdb_name = '{}.pdb'.format(i)

			if os.path.isfile(pdb_name): #might not exist if not that many solutions found
				pose = pose_from_pdb(pdb_name)
				relax_loop(pose,chA_end,args.len,args.pattern)

				if args.save_pdbs:
					pose.dump_pdb('{}_relax.pdb'.format(i))
				
				scores.append(get_score(pose,chA_end,args.len))

				#clean up pdb
				os.remove(pdb_name)

		if len(scores) == 0: #if no remodel outputs, store NaNs

			with open(args.output,'a') as f:
				f.write('{},{},{}\n'.format(tag,'NaN','NaN'))

			continue

		scores = pd.DataFrame(scores,columns=['worst_res_linker','score_per_res_linker'])
		scores = scores.sort_values(by='worst_res_linker',ascending=True)
		scores.index = range(len(scores))

		with open(args.output,'a') as f:
			f.write('{},{},{}\n'.format(tag,*list(scores.iloc[0,:])))

	#clean up tmp file
	os.remove('tmp.bp')



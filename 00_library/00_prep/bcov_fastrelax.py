#!/usr/bin/env python
import sys,os
import pyrosetta

def repack_rotamers_in_list(in_pose=None,
                            pack_scorefxn=None,
                            packable_residues_list=[] #Rosetta ndx
                           ):

    task = pyrosetta.standard_packer_task( in_pose )
    task.initialize_extra_rotamer_flags_from_command_line()
    for iaa in range( 1 , in_pose.total_residue() + 1 ):
        if iaa in packable_residues_list:
            task.nonconst_residue_task( iaa ).restrict_to_repacking()
        else:
            task.nonconst_residue_task( iaa ).prevent_repacking()

    pyrosetta.rosetta.core.pack.pack_rotamers(in_pose, pack_scorefxn, task)


def bcov_fasterlax(in_pose,
                   repack_residue_list=[],
                   relax_bb_residue_list=[],
                   frelax_rounds=1,
                   relax_jumps=[],
                   b_constraint_pose=True,
                   constraint_val=1.0
                  ):

    sfx_A=pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function("beta")
    sfx_A.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.fa_rep,(0.02))
    sfx_A.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.coordinate_constraint,constraint_val)

    sfx_B=pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function("beta")
    sfx_B.set_weight(pyrosetta.rosetta.core.scoring.ScoreType.coordinate_constraint,constraint_val)

    min_chi_vector=pyrosetta.rosetta.utility.vector1_bool(in_pose.total_residue())
    min_bb_vector=pyrosetta.rosetta.utility.vector1_bool(in_pose.total_residue())
    for ires in range(in_pose.total_residue()):
        if (ires in repack_residue_list): #flex_sc and
            min_chi_vector[ires+1]=True
        else:
            min_chi_vector[ires+1]=False
        if (ires in relax_bb_residue_list): #flex_bb and
            min_bb_vector[ires+1]=True
        else:
            min_bb_vector[ires+1]=False
    mm=pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_chi( min_chi_vector );
    mm.set_bb( min_bb_vector );
    if relax_jumps:
        for ijump in relax_jumps:
            mm.set_jump(ijump, True);

    min_mover_A=pyrosetta.rosetta.protocols.minimization_packing.MinMover( mm,
                                                            sfx_B, #sfx_A
                                                            "lbfgs_armijo_nonmonotone",
                                                            0.01,
                                                            True,
                                                            False,
                                                            False)
    min_mover_B=pyrosetta.rosetta.protocols.minimization_packing.MinMover( mm,
                                                                sfx_B,
                                                                "lbfgs_armijo_nonmonotone",
                                                                1e-6,
                                                                True,
                                                                False,
                                                                False)
    for x in range(frelax_rounds):
        if b_constraint_pose:
            for kres in range(1,in_pose.total_residue()+1):
                for jatom in range( 1, in_pose.residue(kres).natoms()+1): #In principle the first four atoms are the BB always
                    if not(in_pose.residue(kres).atom_is_hydrogen(jatom)):
                        in_pose.add_constraint( pyrosetta.rosetta.core.scoring.constraints.CoordinateConstraint(pyrosetta.rosetta.core.id.AtomID(jatom, kres),
                                                                                                                 pyrosetta.rosetta.core.id.AtomID(1, 1), #Fixed reference
                                                                                                                 in_pose.residue(kres).atom(jatom).xyz(), #Reference
                                                                                                                 pyrosetta.rosetta.core.scoring.func.HarmonicFunc(0.0, 1.0) #Functional form is: mean,        HarmonicFunc(mean, std)
                                                                                                                 )
                                                 )

        repack_rotamers_in_list(in_pose=in_pose,
                                pack_scorefxn=sfx_A,
                                packable_residues_list=repack_residue_list
                               )
        min_mover_A.apply( in_pose )

        repack_rotamers_in_list(in_pose=in_pose,
                                pack_scorefxn=sfx_B,
                                packable_residues_list=repack_residue_list
                               )

        min_mover_B.apply( in_pose )

        if b_constraint_pose:
            in_pose.remove_constraints()


if __name__ == "__main__":
    pyrosetta.init("-beta -mute all")
    pdb_file = sys.argv[1]

    print("Working on {}".format(pdb_file))

    input_pose = pyrosetta.pose_from_file(pdb_file)

    bcov_fasterlax(input_pose,
                   repack_residue_list=range(1, input_pose.size()+1),
                   relax_bb_residue_list=range(1, input_pose.size()+1),
                   frelax_rounds=1)

    pdb_id = os.path.basename(pdb_file).split(".")[0]
    pdb_dir = os.path.dirname(pdb_file)

    input_pose.dump_pdb("{}/{}_bcov_relaxed.pdb".format(pdb_dir,pdb_id))
    print("Result saved in {}/{}_bcov_relaxed.pdb".format(pdb_dir,pdb_id))

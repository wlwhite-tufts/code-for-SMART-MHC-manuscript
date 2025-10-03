IMPORTANT: flag files and jupyter notebook need to be edited to match where you put your files. Any text <that looks like this> (starting with "<" and ending with ">") should be replaced with file path or directory information that matches the locations of your files.

1. Create a new RIF, similar to the one created during the RIFdock stage of the 00_library protocol. Use rifgen as in part1:step2 of ../../00_library/01_RIFdock/README.txt
	Use the adjusted ./rifgen.flag file which references ./rif_seed.pdb which contains only W60 from B2m and no other residues
	Note: these flags use the same ../01_RIFdock/rifgen_residues.list as in the RIFdock stage
2. Align each of the final designs with the file from which the Trp position in ./rif_seed.pdb was extracted (for us this was ./target.pdb).
	Use ./align_keep_chain_A.py to do this
	Example command: python align_keep_chain_A.py -inputs <path/to/full/pdbs>/* -out_folder chain_A_aligned -target target.pdb
	Replace <path/to/full/pdbs> with the actual file path to a folder containing pdb files for each design to be scored
3. Run rifdock using a command of the following format:
	rif_dock @rifdock.flag -scaffolds <path/to/aligned/design.pdb> --seed_with_these_pdbs <path/to/aligned/design.pdb> -rif_dock:dokfile <path/to/output.dok> > <path/to/log> 2>&1

	Use the structures produced in step 2 (above) as the inputs (in place of <path/to/aligned/design.pdb>)
	Replace <path/to/output.dok> with the desired output file name, and <path/to/log> with the desired path to a log file
	Use the ./rifdock.flag file rather than the one in ../../00_library/01_RIFdock

4. Collect the dok scores and save them in one large csv file using the ./collect_trp_info.ipynb notebook
	The result of this can be downloaded at https://files.ipd.uw.edu/pub/SMART_MHC_2025/01_select_hits-01_ckeck_Trp-near_perfect_matches_with_rif_scores.csv

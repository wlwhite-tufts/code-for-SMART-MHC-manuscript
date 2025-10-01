IMPORTANT: jupyter notebooks need to be edited to match where you put your files. Any text <that looks like this> (starting with "<" and ending with ">") should be replaced with file path or directory information that matches the locations of your files.

1. (Optional) Select peptide barcodes (similar to those described in https://doi.org/10.1038/s41592-019-0389-8) to add to each design.
	Barcode sequences are listed in ./standard_barcodes_ms1_70000_iRT100_aa.fasta
	Note: we added these barcodes but did not use them in any experiments presented in the paper.
2. Reverse translate barcode sequences. Any reverse translation tool could be used at this step - we used a custom python script, set to avoid BamHI and XhoI cut sites. An example reverse translation command is provided below:
	python reverse_translate.py -aa_fasta standard_barcodes_ms1_70000_iRT100_aa.fasta -dna_fasta standard_barcodes_ms1_70000_iRT100_dna.fasta
	The resulting DNA sequences from our reverse translation are listed in ./standard_barcodes_ms1_70000_iRT100_dna.fasta
3. Create sequence-scrambled control "designs" using the ./make_seqs_for_300mer.ipynb notebook
	a. Randomly select the desired number of designs to use as the basis for controls. We selceted a number such that there would be 11000 total sequences.
	b. Randomly shuffle amino acids within each selected design sequence, preserving only the amino acid category:
		Polar: RHKDESTNQC
		Non-polar: AVILMFYW
		Special: GP
4. Combine the list of scrambled sequences with the list of control sequences (using the ./make_seqs_for_300mer.ipynb notebook).
	The resulting file used in this publication is ./designs_and_scrambles_aa.fasta
5. Reverse translate the design and control sequences. We used the same method as above (step 2), but any other reverse translation method is also acceptable.
	The resulting file used in this publication is ./designs_and_scrambles_dna.fasta
6. (Optional) Randomly assign peptide barcodes to each design and control sequence (using the ./make_seqs_for_300mer.ipynb notebook)
7. (Optional) Filter designs and scrambles to reduce library complexity. This step was performed because of complexity constraints in the DNA library, and can be skipped if the the DNA synthesis method you are using can accomodate more sequences. We sorted the designs by decreasing interface_lddt, randomly removed half of the control sequences, and took from the top of the list until no additional sequences could be accomodated in the DNA library. The sorting and removal of controls was performed with the sort_designs_for_limited_300mer_chip.ipynb notebook.
	The resulting sorted list is contained in ./full_seqs_dna_300mer_sorted.list
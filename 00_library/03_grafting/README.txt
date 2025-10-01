IMPORTANT: jupyter notebooks need to be edited to match where you put your files. Any text <that looks like this> (starting with "<" and ending with ">") should be replaced with file path or directory information that matches the locations of your files.

1. Follow steps 15 and 16 of the Cao et al guide to extract and cluster the motifs from the designs created in the previous step.
2. Use ./longxing_motif_selection.ipynb to select the best motif from each cluster, and then select the best of those to use in the grafting step.
3. Follow step 19 of the Cao et al guide, using the same scaffold set as in the RIFdock stage:
	http://files.ipd.uw.edu/pub/robust_de_novo_design_minibinders_2021/supplemental_files/scaffolds.tar.gz
	Use ./tutorial_grafting_v2.xml with the flags set in ./grafting.flag
4. Repeat the fastdesign process using exactly the same scripts (in ../02_fastdesign) but a new Jupyter notebook (./graft_predictor.ipynb) and using the new outputs from motif grafting (step 3, above) as inputs.
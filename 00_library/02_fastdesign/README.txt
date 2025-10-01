IMPORTANT: flag files, python scripts, and jupyter notebooks need to be edited to match where you put your files. Any text <that looks like this> (starting with "<" and ending with ">") should be replaced with file path or directory information that matches the locations of your files.

1. Collect RIFdock outputs
	To save space, convert the outputs from RIFdock to one large silent file. Make sure to adjust the name of each design such that they are unique and the RIF that was used (standard or native) can be determined. See Silent Tools documentation for details (https://github.com/bcov77/silent_tools).
2. Split the full silent file into reasonable sized batches following recommendations in step 11 of the Cao et al guide. Do not delete the full silent file.
3. Run the ./predictor_v11.xml RosettaScript using the ./predictor.flags file on each batch from step 2. The command should look like the one below:
	<path/to/your_rosetta>/rosetta/latest/bin/rosetta_scripts @predictor.flags -parser:protocol predictor_v11.xml -in:file:silent <path/to/your/batch.silent>
	where <path/to/your_rosetta> is the location of your Rosetta installation and <path/to/your/batch.silent> is the location of a small silent file created in the previous step.
	IMPORTANT: this format of command can be used in all other cases where Rosetta scripts are run. In some cases you may also need to specify an output file by setting the -out:file:silent flag to the desired location of the output.
4. Select a random sample of RIFdock outputs to run as a pilot test, as in step 12 of the Cao et al guide.
5. Run the ./tutorial_design_v5.xml RosettaScript using the pilot.flags file on the random sample from step 4
6. Calculate the distance between the C-terminus of the scaffold and the N-terminus of the target on all batches from step 2 using ./CtermA_NtermB_dist.py, using commands like the following:
	python CtermA_NtermB_dist.py -in:file:silent <path/to/batch.silent> -out_file <path/to/desired/output.csv>
	where <path/to/batch.silent> is replaced with the location of the silent file to use as input, and <path/to/desired/output.csv> is the desired location of the output file of distances
7. Use the ./dock_predictor.ipynb Jypyter Notebook to evaluate the likelihood that each design evaluated in step 3 will actually pass the full filter set. As suggested in step 13 of the Cao et al guide, the exact filter cutoffs may need to be adjusted to get reasonable predictive power. This notebook will also filter out any designs where the C-terminus of the design is too far from the N-terminus of the target (in order to minimize linker lengths).
8. Follow step 14 of the Cao et al guide to run ./tutorial_design_v5.xml with ./pilot.flags on the designs selected in step 7.
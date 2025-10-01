IMPORTANT: jupyter notebook needs to be edited to match where you put your files. Any text <that looks like this> (starting with "<" and ending with ">") should be replaced with file path or directory information that matches the locations of your files.

1. Assemble full-length sequences from paired-end read data.
	We used a custom python script to do this, but any method for assembling paired-end read data could be used instead. This script finds the best overlap between forward and reverse reads, creates an assembled DNA sequence, translates it into an amino acid sequence, and counts the number of times each amino acid sequence was observed.

	For each pair of files [sample_name]_S[n]_L[x]_R001.fastq.gz and [sample_name]_S[n]_L[x]_R002.fastq.gz (where [sample_name] is the name of the sample, [n] is the sample number, and [x] is the lane number), run the following command:

	python fastq_gz_assembly.py -lib_name [sample_name]_S[n]_L[x] -primer_fasta pETCON_adapters.fasta -in_dir <path/to/fastq/folder/> -out_dir <path/to/output/folder/>

	replacing <path/to/fastq/folder/> with the path to the folder where the fastq.gz files are stored, and <path/to/output/folder/> with the desired output folder.

	Note: the python script is located in ./fastq_gz_assembly.py, and required a file containing the sequences of the sequencing primers (./pETCON_adapters.fasta)

2. Merge data across all sequencing samples
	Use the ./combine_counts.ipynb notebook to read in the csvs created above and create a data table where each row represents an observed sequence, and each column represents a sample, with entries representing the number of times a sequence was observed in a sample. The result of this tabulation for our experiments is ./all_counts.csv
	Note: this notebook also creates a fasta file with all observed sequences that is used in the next step.

3. Create a BLAST database from the sequences of the design library (to query assembled sequences against)
	Use the following command to create a BLAST database:

	makeblastdb -in <path/to/library_aa.fasta> -input_type fasta -dbtype prot -out <path/to/db>

	replacing <path/to/library_aa.fasta> with a path to the fasta file containing the amino acid sequences of every design and control in the library (for us, this is ../../00_library/06_encoding/designs_and_scrambles_aa.fasta) and <path/to/db> with the desired prefix for the database files that will be created.

4. For each sequence observed in the sequencing data, search the library database for the closest match
	We used a custom python script to do this, but it could also be done with the BLAST command line tool,

	a. Split the observed sequences into batches of 500 sequences (each batch represented by a fasta file containing all sequences in the batch)

	b. For each batch, use the ./find_sequencing_matches.py script to identify matches:
	python find_sequencing_matches.py <path/to/batch.fasta> <path/to/db> <path/to/output.csv>
	replacing <path/to/batch.fasta> with the location of the batch fasta file, <path/to/db> with the same database prefix as you used above, and <path/to/output.csv> with the desired output file name

	c. concatenate all the resulting batch outputs into a large csv. The result of this for our dataset was ./all_matches.csv
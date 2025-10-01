IMPORTANT: jupyter notebook needs to be edited to match where you put your files. Any text <that looks like this> (starting with "<" and ending with ">") should be replaced with file path or directory information that matches the locations of your files.

This section is performed entirely within the ./analyze_ngs.ipynb notebook
Note: data from several additional sorting experiments not described in the paper are included but not used in the selection of hits.

1. Load sequence matching, design scoring, and RIF scoring, and NGS count data for all sequences (../01_check_Trp/near_perfect_matches_with_rif_scores.csv)

2. remove noisy sequences
	a. label sequences with confident count values in each experiment
		i. find the line of best fit (forced to have intercept = 0) to the relationship between abundance in one replicate and the other
		ii. calculate the distance of each sequence to this line, in log-transformed space
		iii. "confident" sequences have a log-distance greater than 0.5, and at least two counts in each replicate
	b. plot correlation between replicates of each sorting experiment to visualize the noisy and confident sequences
	c. remove sequences that are not labeled as confident in any experiment

3. Calculate enrichment values for each sequence in each sorting experiment
	a. calculate the fraction of counts representing each sequence in the pre- and post- sorted samples
	b. calculate the ratio of the post-sort fraction to the pre-sort fraction
	c. plot enrichment distributions to compare sequences that match controls, designs with W60-like Trp plpacement, and other designs

4. Select hits (manual set)
	a. remove sequences that are not positively enriched during the expression sort (enrichment > 1)
	b. remove sequences that are not well-enrched during the peptide binding sort (enrichment > 2)
	c. remove sequences with more than one mutation relative to the clisest design
	d. remove sequences with low counts after expression sorting (> 5)
	e. sort equences in order of descending peptide binding enrichment and manually select 6 sequences from the top 20 based on:
		i. enrichment in the peptide binding sort
    	ii. presence of multiple mutants of the same design in the top 20
    	iii. using the unmutated version if both are present in the top 20
    Note: the result of this selection for our dataset is saved in ./manual_hits.fasta

5. Select hits (W60-containing)
	a. apply filters a-d in step 4, above
	b. remove sequences that don't have a Trp residue placed in a similar manner to W60 from B2m
	c. select the top 24 most enriched designs (in the peptide binding sort) excluding the 6 selected in step 4, above
	Note: the result of this selection for our dataset is saved in ./trp_hits.fasta
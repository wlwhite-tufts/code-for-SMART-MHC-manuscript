IMPORTANT: jupyter notebook needs to be edited to match where you put your files. Any text <that looks like this> (starting with "<" and ending with ">") should be replaced with file path or directory information that matches the locations of your files.

This series of steps relies on DeepAccNet, which can be installed following the instructions prvided here: https://github.com/hiranumn/DeepAccNet

1. Split the 04_surface outputs into batches and run the ./extra_paper_metrics.xml script on each batch using ./metric.flags
2. Evaluate each of the batches from setp 1 with DeepAccNet using the following command as a template:
	python <path/to/your_deepaccnet>/DeepAccNet-SILENT.py <path/to/batch.silent> <name/of/score.csv> --binder
	replacing <path/to/your_deepaccnet> with the location where you intalled the DeepAccNet package and <path/to/batch.silent> and <name/of/score.csv> with the locations of batch silent file, and the score file to be created for the batch, respectively
3. Extract sequence information from the batch silent files using a command like the following:
	silentsequence <path/to/batch.silent> > <path/to/designs.seq>
	replacing <path/to/batch.silent> and <path/to/designs.seq> with the locations of batch silent file, and the sequence file to be created, respectively
4. Collect all the batches of scores into one large score file for each score type (Rosetta metrics and DeepAccNet predictions). You can do this with whatever tool suits you best.
5. Use the ./final_filtering.ipynb to select the best designs. Filter cutoff can be adjusted based on the recommendations in step 21 of the Cao et al guide. Note that the Cao et al guide does not have recommendationds for DeepAccNet scores, but that these should generally be maximized under the constraint of producing the desired number of passing designs. Finally, the length cutoff for designs was set to 57 due to practical limitations on DNA sequence length. In order to get roughly 10k designs that passed the filters, we took designs that were roughly in the top 50% of designs in all metrics. Structures of the final selected designs are listed in ./selected_designs_with_metrics_renamed.csv

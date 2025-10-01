import pandas as pd
import json
import pickle
from Bio.Blast.Applications import NcbiblastpCommandline
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Find matching sequences in library")
parser.add_argument('sequencing_fasta',type=str,help='The FASTA file containing the sequencing results for hits from a library screen.')
parser.add_argument('library_db_name',type=str,help='The BLAST db name containing the library sequences.')
parser.add_argument('output_file',type=str,help='The name of the pkl file to write the results to.')
args = parser.parse_args()

if __name__ == '__main__':

	blastp_cline = NcbiblastpCommandline(query=args.sequencing_fasta,
										db=args.library_db_name,
										outfmt=15,
										max_target_seqs=1)
	std_out,std_err = blastp_cline()
	result = json.loads(std_out)['BlastOutput2']

	out_df = pd.DataFrame(columns=['query_name','query_len','hit_len','match_name','mismatches_query','mismatches_alignment','algn_len','score','e_val'])

	#for each query in sequencing_fasta
	for q in range(len(result)):
		#just look at top hit because no danger of chimeras here
		query_len = result[q]['report']['results']['search']['query_len']
		if len(result[q]['report']['results']['search']['hits']) > 0:
			hit = result[q]['report']['results']['search']['hits'][0]
			n_match = hit['hsps'][0]['identity']

			out_df.loc[q,'query_name'] = result[q]['report']['results']['search']['query_title']
			out_df.loc[q,'query_len'] = query_len
			out_df.loc[q,'hit_len'] = hit['len']
			out_df.loc[q,'match_name'] = hit['description'][0]['title']
			out_df.loc[q,'mismatches_query'] = query_len - n_match
			out_df.loc[q,'mismatches_alignment'] = hit['hsps'][0]['align_len'] - n_match
			out_df.loc[q,'algn_len'] = hit['hsps'][0]['align_len']
			out_df.loc[q,'score'] = hit['hsps'][0]['score']
			out_df.loc[q,'e_val'] = hit['hsps'][0]['evalue']

		else:
			out_df.loc[q,'query_name'] = result[q]['report']['results']['search']['query_title']
			out_df.loc[q,'query_len'] = query_len
			out_df.loc[q,'hit_len'] = np.nan
			out_df.loc[q,'match_name'] = ''
			out_df.loc[q,'mismatches_query'] = np.nan
			out_df.loc[q,'mismatches_alignment'] = np.nan
			out_df.loc[q,'algn_len'] = np.nan
			out_df.loc[q,'score'] = np.nan
			out_df.loc[q,'e_val'] = np.nan

	out_df.to_csv(args.output_file)


#!/usr/bin/env python
import os,sys,re
from glob import glob
import subprocess
from multiprocessing import Pool
import time
import pandas as pd
from collections import defaultdict
import gzip
import numpy as np
from tqdm import tqdm

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.pairwise2 import align

import argparse

def get_args():
    parser = argparse.ArgumentParser(description="Find and assemble pairs of forward and reverse reads.")
    parser.add_argument("-lib_name", help="The string that uniquely identifies the library among all the files in the input directory.", type=str)
    parser.add_argument("-primer_fasta", help="The fasta file with the primer sequences to check for.", type=str)
    parser.add_argument("-in_dir", help="The input directory.", type=str)
    parser.add_argument("-out_dir", help="The output directory.", type=str)
    parser.add_argument("-min_overlap", help="The minimum number of basepairs of overlap between the paired sequences.", type=int, default=6)
    parser.add_argument("-max_bad_end", help="The maximum number of allowed bp on each end of the overlap that don't match.", type=int, default=3)
    parser.add_argument("-max_primer_mismatch", help="The maximum number of allowed bp in the primer alignment that don't match.", type=int, default=3)

    return parser.parse_args()

def find_overlap(seq_F, seq_R, max_bad_end):
    seq_F = np.array(list(seq_F))
    seq_R = np.array(list(seq_R))

    max_match = 0
    best_F_bad_end_count = np.inf
    best_R_bad_end_count = np.inf
    best_match = np.array([])
    for i in range(1,min(len(seq_F),len(seq_R))):
        matches = seq_F[-i:] == seq_R[:i]
        streak_starts = matches[:-1] != matches[1:]
        streak_starts = np.insert(streak_starts, 0, True)
        streak_ids = np.cumsum(streak_starts)*matches #zero out the streaks of false values
        unique_streaks,streak_lens = np.unique(streak_ids, return_counts=True)
        #remove 0 from the entries
        unique_streaks = unique_streaks[1:]
        streak_lens = streak_lens[1:]

        if len(streak_lens) == 0:
            continue

        long_streak_id = unique_streaks[streak_lens == np.max(streak_lens)][0] #take index zero in case of ties
        long_streak_locs = np.where(streak_ids == long_streak_id)[0]
        F_bad_end_count = long_streak_locs[0]
        R_bad_end_count = len(streak_ids) - long_streak_locs[-1] - 1

        if (np.max(streak_lens) > max_match) and (max(F_bad_end_count, R_bad_end_count) <= max_bad_end):
            max_match = np.max(streak_lens)
            best_F_bad_end_count = F_bad_end_count
            best_R_bad_end_count = R_bad_end_count
            best_match = matches

    return max_match, best_F_bad_end_count, best_R_bad_end_count, best_match

def look_for_primer(primer, read, direction):

    algn = align.localxd(primer,read, -np.inf, -np.inf, -0.5, -0.1, one_alignment_only=True)[0]
    algn_primer = np.array(list(algn[0]))
    algn_read = np.array(list(algn[1]))
    primer_locs = np.where(algn_primer != '-')[0]
    start = primer_locs[0]
    end = primer_locs[-1]

    n_mismatch = np.sum(algn_primer[start:end+1] != algn_read[start:end+1])

    if direction == 'fwd':
        read_start = np.sum(algn_read[:end+1] != '-') #how many bases of the read exist before the primer ends
    if direction == 'rev':
        read_start = np.sum(algn_read[start:] != '-') #how many bases of the read exist after the primer starts

    # print(''.join(algn_primer))
    # print(''.join(algn_read))
    # print(n_mismatch)
    # print(read[:-read_start])
    # print('*******************')

    return n_mismatch, read_start


if __name__ == "__main__":
    pwd = os.getcwd()
    debug = False

    args = get_args()

    # ## step0: collect the AA sequence ordered
    # order_fasta = f"{pwd}/designs_and_scrambles_aa.fasta"
    # order_dict = defaultdict(list)
    # for each_record in SeqIO.parse(order_fasta, "fasta"):
    #     print( each_record.id )
    #     print( each_record.seq )
    #     order_dict[each

    ## ================= step1: get adaptor information ================== ##
    adaptors_fasta = args.primer_fasta
    adaptors_dict = defaultdict(lambda: defaultdict(str))
    for each_record in SeqIO.parse(adaptors_fasta, "fasta"):
        base_id, term_type = each_record.id.split("_")
        adaptors_dict[base_id][term_type] = str(each_record.seq)


    print(adaptors_dict)

    ## ================= step2: gather NGS libraries information ================== ##
    fastq_dir = args.in_dir

    # collect all file ids and names
    gz_list = glob(f"{fastq_dir}/*{args.lib_name}*fastq.gz")
    print(gz_list)

    ## example name: Hua-Will-pETCON-library-pep1-reseq1_S9_L004_R1_001.fastq.gz
    id_pattern = re.compile("Hua-Will-(.+)-(library-.+)_S([0-9]+)_L([0-9]+)_R([12]+)_.+fastq.gz")

    FASTQ_files = defaultdict(lambda:defaultdict(lambda:defaultdict(str)))

    for each_gz in gz_list:
        id_info = id_pattern.findall(each_gz)
        if not id_info:
            print("no id match")
            continue

        vec_id   = str(id_info[0][0])
        if "pETCON" in vec_id:
            vec_id = "pETCON3"
        if "pET29" in vec_id:
            vec_id = "pET29b"

        lib_name = str(id_info[0][1])
        S_num    = int(id_info[0][2])
        L_num    = int(id_info[0][3])
        R_num    = int(id_info[0][4])

        lib_name = f"{vec_id}_{lib_name}_{S_num}"

        if R_num == 1:
            FASTQ_files[lib_name][L_num]["F"] = each_gz
        elif R_num == 2:
            FASTQ_files[lib_name][L_num]["RC"] = each_gz
        else:
            continue

    print(FASTQ_files)
    ## ================= step3: go through all libraries ================== ##
    print( FASTQ_files.keys() )
    for lib_name in FASTQ_files.keys():
        print(f"================ working on lib {lib_name} ===================")
        vec_id = lib_name.split("_")[0]

        adaptor_5prime = adaptors_dict[vec_id]["5prime"]
        adaptor_3prime = adaptors_dict[vec_id]["3prime"]

        AAseq_counter = defaultdict(int)

        ## each library will have a couple of channels, the L_num
        for L_num in FASTQ_files[lib_name].keys():
            fastq_pair = FASTQ_files[lib_name][L_num]
            print(fastq_pair)

            ## A_dict is the R1 sequences, the forward sequencing
            A_dict = {}

            with gzip.open(fastq_pair["F"], "rt") as gz_f:
                for record in SeqIO.parse(gz_f, "fastq"):

                    n_mismatch, read_start = look_for_primer(adaptor_5prime, str(record.seq), 'fwd')

                    if n_mismatch > args.max_primer_mismatch:
                        continue

                    DNA_str_A = str(record.seq)[read_start:]

                    ## remove bad sequencing?
                    if ("N" in DNA_str_A)  or (len(DNA_str_A) == 0):
                        continue

                    DNA_seq_A = Seq(DNA_str_A, IUPAC.unambiguous_dna)

                    A_dict[record.id] = DNA_seq_A

            ## B_dict is the R2 sequences, the reverse complement sequencing
            B_dict = {}

            with gzip.open(fastq_pair["RC"], "rt") as gz_f:
                for record in SeqIO.parse(gz_f, "fastq"):

                    n_mismatch, read_start = look_for_primer(adaptor_3prime, str(record.seq.reverse_complement()), 'rev')

                    if n_mismatch > args.max_primer_mismatch:
                        continue

                    DNA_str_B = str(record.seq.reverse_complement())[:-read_start]

                    ## remove bad sequencing?
                    if ("N" in DNA_str_B)  or (len(DNA_str_B) == 0):
                        continue

                    DNA_seq_B = Seq(DNA_str_B, IUPAC.unambiguous_dna)

                    B_dict[record.id] = DNA_seq_B

            ## and make a spliced sequence
            for each_key in A_dict.keys() & B_dict.keys():
                A_DNA = str( A_dict[each_key] )
                B_DNA = str( B_dict[each_key] )

                max_match, F_bad_end_count, R_bad_end_count, best_match = find_overlap(A_DNA, B_DNA, args.max_bad_end)
                # print(max_match, F_bad_end_count, R_bad_end_count, best_match)

                if max_match >= args.min_overlap and max(F_bad_end_count, R_bad_end_count) <= args.max_bad_end:
                    full_DNA = Seq(A_DNA[:len(A_DNA)-F_bad_end_count] + B_DNA[R_bad_end_count+max_match:], IUPAC.unambiguous_dna)

                    if len(full_DNA)%3:
                        continue

                    full_AA = full_DNA.translate()

                    # print(A_DNA)
                    # print(' '*(len(A_DNA) - F_bad_end_count - R_bad_end_count - max_match) + B_DNA)
                    # print(full_DNA)
                    # print(full_AA)
                    # print('---------------------')

                    # if len(full_DNA)%3:
                    #     print('bad length')

                    AAseq_counter[full_AA] += 1

        AAseq_df = pd.DataFrame( AAseq_counter.items(), columns=['AA', "count"] )
        AAseq_df.sort_values(by=['count'], ascending=False, inplace=True, ignore_index=True)
        out_csv_name = f"result_df_{args.lib_name}.csv"
        AAseq_df.to_csv(args.out_dir+out_csv_name, index=False)

import os
import numpy as np
import sys
from optparse import OptionParser
from Bio import AlignIO
import tempfile
import subprocess

sys.path.append('<path/to/your/>/02_cut_site/')
import pyrosetta
import pyrosetta.rosetta
from pyrosetta.rosetta.core.scoring import *
import scoring_utils

# Sliding window: Only keep PSSM vector for focal position -- but look 4 positions up and downstream
def get_window(aln_arr, aln_pos, n=4):
    """
    Input:  aln_arr: np.array of np.characters
            aln_pos: which position to construct the window around
            n:       size of the window up/downstream up aln_pos
    Output: Window indexes
    
    """
    assert aln_pos >= 0
    assert aln_pos < len(aln_arr)
    assert aln_arr[aln_pos] != b'-' # expecting this to be an amino acid...
    
    # First we need to find the idx at the upstream end of the seq window
    upstream_idx, aa_counter = aln_pos, 0
    for c in aln_arr[aln_pos+1:]:
        if aa_counter == n:
            break
        if c != b'-':
            aa_counter+=1
        upstream_idx += 1
        
    # Next we need the downstream sequence bit
    downstream_idx, aa_counter = aln_pos, 0
    for c in reversed(aln_arr[:aln_pos]):
        if aa_counter == n:
            break
        if c != b'-':
            aa_counter += 1
        downstream_idx -= 1
        
    return downstream_idx, upstream_idx

def write_charArr_2_fasta(charArr, filename):
    seqs = [''.join(s.astype(str)) for s in charArr]
    with open(filename, 'w') as f_out:
        for i,seq in enumerate(seqs):
            f_out.write(f'>s_{i}\n')
            f_out.write(f'{seq}\n')

def read_pssm(filepath):
    datalines = []
    with open(filepath, 'r') as f_open:
        for line in f_open:
            is_pssm_line = (len(line.split()) == 44)
            if is_pssm_line:
                datalines.append(line)
    assert len(datalines) > 4 # if this fails, check that the format of the PSSM is as expected (pssm vectors when split should be 44 long)
    return datalines

def save_pssm(pssm_vector_lines, out_f):
    out_str = ''
    out_str += '\n'
    out_str += 'Last position-specific scoring matrix computed, weighted observed percentages rounded down, information per position, and relative weight of gapless real matches to pseudocounts\n'
    out_str += '            A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V\n'
    for l in pssm_vector_lines:
        out_str += l   
    out_str += '\n'
    out_str += '                      K         Lambda\n'
    out_str += 'PSI Ungapped         0.1946     0.3179\n'
    out_str += 'PSI Gapped           0.0596     0.2670\n'
    with open(out_f,'w') as f_out:
        f_out.write(out_str)
    

def msa2pssm(msa):
    tmpdir = tempfile.TemporaryDirectory()
    tmp_dirname = tmpdir.name
    
    # write the msa
    msa_f = f'{tmp_dirname}/msa.fasta'
    write_charArr_2_fasta(msa, msa_f)

    # Write the query
    query_f = f'{tmp_dirname}/query.fasta'
    seq_no_gap = msa[0][msa[0]!=b'-']
    write_charArr_2_fasta([seq_no_gap], f'{tmp_dirname}/query.fasta')
    
    # Get the PSSM
    pssm_path = f'{tmp_dirname}/pssm.fasta'
    tmp_out = f'{tmp_dirname}/tmp.out'
    tmp_err = f'{tmp_dirname}/tmp.out'
    psiblast_exe = '<path/to/your/psiblast>/bin/psiblast'
    psiblast_cmd = f'{psiblast_exe} -subject {query_f} -in_msa {msa_f} -ignore_msa_master -out_ascii_pssm {pssm_path} > {tmp_out} 2>{tmp_err}'
    subprocess.run(psiblast_cmd, shell=True)
    #os.system(psiblast_cmd)
    
    pssm = read_pssm(pssm_path)
    
    return pssm

def fix_scorefxn(sfxn, allow_double_bb=False, repulsive=False):
    if not repulsive:
        sfxn.set_weight(fa_rep, 0.0)
    else:
        sfxn.set_weight(fa_rep, 1.0)
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)

def aln2pssm_window_and_contacts_filter(aln_arr, pssm_out, pdb, win_size=4, nbr_energy_cutoff=-1.5):
    """
    Input
        aln_f: Path to fasta alignment file
        win_size: Size of the sequence window used
                  for indel matching during alignment
                  filtering.
    
    Output: recombined_pssm
        A PSSM constructed from filtered subalignments
        that must all share the indel-pattern of the
        top sequence in the alignment.
    """
    # Load pose and find 2 body energies
    pyrosetta.init('-beta')
    pose_by_chain = pyrosetta.pose_from_pdb(pdb).split_by_chain()
    sf_wo_rep = pyrosetta.get_fa_scorefxn()
    fix_scorefxn(sf_wo_rep, repulsive=False)
    onebody_energies, twobody_energies = scoring_utils.get_one_and_twobody_energies(pose_by_chain[2], sf_wo_rep)
    
    # Load alignment and convert to something useful
    curr_seq = aln_arr[0,:]
    aln_depth = aln_arr.shape[0]
    
    # Setup aln to seq map
    aa2aln_idx = {}
    aa_seen = 0
    for i, c in enumerate(curr_seq):
        if c!=b'-': 
            aa2aln_idx[aa_seen] = i
            aa_seen += 1
    
    # Initialize various variables
    pssm_lines = []
    aa_counter = 0
    tgt_seq_len = sum(curr_seq!=b'-')
   
    # Iterate over the gapped alignment
    for i,c in enumerate(curr_seq):
        # Skip if focal is gap
        if c == b'-': 
            continue

        # Find nbrs
        contacts = np.where(twobody_energies[aa_counter,:]<nbr_energy_cutoff)[0]

        # Which positions do they correspond to in the alignment?
        aln_contacts = np.array([aa2aln_idx[resi0] for resi0 in contacts]).astype(int)
        
        # Find the sequence window and combine with aln_contacts
        assert win_size >= 1
        down_idx, up_idx = get_window(curr_seq, i, n=win_size)
        idx_range = np.arange(down_idx, up_idx+1)
        all_idxes = np.unique(np.concatenate([aln_contacts,idx_range]))
        
        # select submsa accordingly
        submsa = aln_arr[:,all_idxes]
        
        # Get rid of sequences that has non-matching indels
        indel_pattern = (submsa == b'-')
        query_indel_pattern = indel_pattern[0]
        non_indel_match_idxes = np.where((indel_pattern==query_indel_pattern).sum(axis=1)!=len(query_indel_pattern))[0]
        indel_match_idxes = np.where((indel_pattern==query_indel_pattern).sum(axis=1)==len(query_indel_pattern))[0]
        
        tmp_aln = aln_arr.copy()
        
        # Strategy 1: For the non-matching indels insert gaps, so they don't polute the PSSM information
        # This is apperently not a good way to do it. Very gapped
        # columns ends up with extremely low information content for some reason
        # tmp_aln[non_indel_match_idxes, i] = b'-'
        # print(f'For pos {aa_counter} censoring {len(non_indel_match_idxes)} seqs out of {aln_depth}')
        
        # Stategy 2: Let's just subselect all the sequences that match
        tmp_aln = tmp_aln[indel_match_idxes]
        print(f'For pos {aa_counter} subselecting {len(indel_match_idxes)} sequences out of {aln_depth}')
        log_handle.write(f'For pos {aa_counter} subselecting {len(indel_match_idxes)} sequences out of {aln_depth}\n')

        make_pssm = True
        if make_pssm:
            # Compute the PSSM
            pssm = msa2pssm(tmp_aln)
            
            # The pssm has the length of the tgt sequence
            tgt_pssm_vector = pssm[aa_counter]
            
            # Store the pssm vector for the cleaned positions
            pssm_lines.append(tgt_pssm_vector)
        aa_counter += 1
        
    native_seq = ''.join(curr_seq.astype(str)).replace('-','')
    ctrl_seq = ''.join([l.split()[1] for l in pssm_lines])

    # Check that the sequence from the combined PSSM
    # matches with the query sequence. (it must)
    assert native_seq == ctrl_seq
    
    # Finally save the pssm
    save_pssm(pssm_lines, pssm_out)


def aln2pssm_window_filter(aln_arr, pssm_out, win_size=4):
    """
    Input
        aln_f: Path to fasta alignment file
        win_size: Size of the sequence window used
                  for indel matching during alignment
                  filtering.
    
    Output: recombined_pssm
        A PSSM constructed from filtered subalignments
        that must all share the indel-pattern of the
        top sequence in the alignment.
    """
    
    # Load alignment and convert to something useful
    curr_seq = aln_arr[0,:]
    aln_depth = aln_arr.shape[0]
    
    # Initialize various variables
    pssm_lines = []
    tgt_aa_seen = 0
    tgt_seq_len = sum(curr_seq!=b'-')
    
    # Iterate over the gapped alignment
    for i,c in enumerate(curr_seq):
        # Skip if focal is gap
        if c == b'-': 
            continue
        
        # Find the sequence window and select submsa accordingly
        assert win_size >= 1
        down_idx, up_idx = get_window(curr_seq, i, n=win_size)
        submsa = aln_arr[:,down_idx:up_idx+1]
        
        # Get rid of sequences that has non-matching indels
        indel_pattern = (submsa == b'-')
        query_indel_pattern = indel_pattern[0]
        non_indel_match_idxes = np.where((indel_pattern==query_indel_pattern).sum(axis=1)!=len(query_indel_pattern))[0]
        indel_match_idxes = np.where((indel_pattern==query_indel_pattern).sum(axis=1)==len(query_indel_pattern))[0]
        
        tmp_aln = aln_arr.copy()
        
        # Strategy 1: For the non-matching indels insert gaps, so they don't polute the PSSM information
        # This is apperently not a good way to do it. Very gapped
        # columns ends up with extremely low information content for some reason
        # tmp_aln[non_indel_match_idxes, i] = b'-'
        # print(f'For pos {tgt_aa_seen} censoring {len(non_indel_match_idxes)} seqs out of {aln_depth}')
        
        # Stategy 2: Let's just subselect all the sequences that match
        tmp_aln = tmp_aln[indel_match_idxes]
        print(f'For pos {tgt_aa_seen} subselecting {len(indel_match_idxes)} sequences out of {aln_depth}')
        log_handle.write(f'For pos {tgt_aa_seen} subselecting {len(indel_match_idxes)} sequences out of {aln_depth}\n')
        
        # Compute the PSSM
        pssm = msa2pssm(tmp_aln)
        
        # The pssm has the length of the tgt sequence
        tgt_pssm_vector = pssm[tgt_aa_seen]
        
        # Store the pssm vector for the cleaned positions
        pssm_lines.append(tgt_pssm_vector)
        tgt_aa_seen += 1
        
    native_seq = ''.join(curr_seq.astype(str)).replace('-','')
    ctrl_seq = ''.join([l.split()[1] for l in pssm_lines])
    
    # Check that the sequence from the combined PSSM
    # matches with the query sequence. (it must)
    assert native_seq == ctrl_seq
    
    # Finally save the pssm
    save_pssm(pssm_lines, pssm_out)

def aln2pssm(aln_arr, pssm_out):
    """
    Input
        aln_f: Path to alignment file

    Output
        pssm_out: output PSSM
    """
    
    # Load alignment and convert to something useful
    curr_seq = aln_arr[0,:]
    pssm = msa2pssm(aln_arr)
    native_seq = ''.join(curr_seq.astype(str)).replace('-','')
    ctrl_seq = ''.join([l.split()[1] for l in pssm])
    
    # Check that the sequence from the combined PSSM
    # matches with the query sequence. (it must)
    assert native_seq == ctrl_seq
    
    # Finally save the pssm
    save_pssm(pssm, pssm_out)

if __name__ == "__main__":
    parser = OptionParser(usage="usage: %prog [options] FILE", version="0.1")
    parser.add_option("--in_fasta", type="string", dest="in_fasta", metavar="STR", help="Input MSA in fasta format")
    parser.add_option("--in_a3m", type="string", dest="in_a3m", metavar="STR", help="Input MSA in a3m format")
    
    parser.add_option("--pdb", type="string", dest="pdb", metavar="STR", help="Input PDB")
    parser.add_option("--out_pssm", type="string", dest="out_pssm", metavar="STR", help="Path to output PSSM")
    parser.add_option("--win_size", type="string", dest="win_size", metavar="STR", help="Filtering window size")
    parser.add_option("--filter_mode", type="string", dest="filter_mode", metavar="STR", help="0: normal psiblast, 1: indel filter by window, 2: indel filter by window and contacts")
    parser.add_option("--query_seq_min_id", type="string", dest="query_seq_min_id", metavar="STR", help="Minimum sequence identity to the query sequence")
    
    (opts, args) = parser.parse_args()
    
    assert opts.query_seq_min_id is not None
    assert opts.in_fasta is None or opts.in_a3m is None # Only provide the alignment in one way
    
    #if os.path.isfile(opts.out_pssm):
    #sys.exit()

    min_seqid = int(opts.query_seq_min_id)
    
    # Get the alignment
    if opts.in_a3m is not None:
        tmpdir = tempfile.TemporaryDirectory()
        tmp_dirname = tmpdir.name
        
        convert_script = '<path/to/your/>/02_cut_site/reformat.pl'
        filter_script = '<path/to/your/hh-suite>/build/bin/hhfilter'
        assert os.path.isfile(convert_script)
        assert os.path.isfile(filter_script)
        
        # We do some conservative redundancy reduction here, just to speed up psiblast
        # (psiblast is internally doing some sequence reweighing)
        max_seqid = 90
        n_hits = 1000000
        while n_hits>3000:
            cmd_filter = f'{filter_script} -qid {min_seqid} -id {max_seqid} -cov 75 -i {opts.in_a3m} -o {tmp_dirname}/filtered.a3m > {tmp_dirname}/filt.log 2>{tmp_dirname}/filt.err'
            subprocess.run(cmd_filter, shell=True)
            n_hits = sum([1 for l in open(f'{tmp_dirname}/filtered.a3m', 'r') if '>' in l])
            log_str = f'Filtering the alignment (id<{max_seqid} & cov>70) and seqid vs. query: {min_seqid}. Number seq hits: {n_hits}'
            print(log_str)
            max_seqid -= 1
            if n_hits<30: # This likely happens because the first iteration resulted in very few sequences! 
                           # This can happen because we are iterating over higher and higher min_seqid
                print("Found too few sequences... returning")
                sys.exit()
        
        print("Converting a3m to fasta")
        cmd_convert = f'{convert_script} a3m fas {tmp_dirname}/filtered.a3m {tmp_dirname}/tmp.fasta'
        subprocess.run(cmd_convert, shell=True)
        aln = AlignIO.read(f'{tmp_dirname}/tmp.fasta','fasta')
    else:
        aln = AlignIO.read(opts.in_fasta,'fasta')
    
    # Convert the alignment to some useful format
    aln_arr = np.array([list(rec) for rec in aln], np.character)
    
    pssm_f = opts.out_pssm
    win_size = 4 if opts.win_size is None else int(opts.win_size)
    pdb = opts.pdb
    log_handle = open(pssm_f+'.log', 'w')
    log_handle.write(log_str+'\n')
    
    if opts.filter_mode=='0':
        print("Making PSSM without alignment filtering")
        aln2pssm(aln_arr, pssm_f)
    elif opts.filter_mode=='1':
        print("Making PSSM using local sequence indel match filtering")
        aln2pssm_window_filter(aln_arr, pssm_f, win_size=win_size)
    elif opts.filter_mode=='2':
        assert opts.pdb is not None
        print("Making PSSM using local sequence + contacts indel match filtering")
        aln2pssm_window_and_contacts_filter(aln_arr, pssm_f, pdb, win_size=win_size, nbr_energy_cutoff=-1.5)
    else:
        sys.exit("You must provide a filter_mode.")
    





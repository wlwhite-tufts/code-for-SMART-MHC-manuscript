"""Modified from Robby Divine's original script
"""

import random
import sys

codon_dict = {
#AGG, AGA, CGG, CGA (R), GGA (G), ATA (I), CTA (L), CCC (P), ACA (T) associated with poor expression in E. coli. Indexed by frequency found in E. coli B
    'A': ['GCG', 'GCC', 'GCA', 'GCT'],
    'C': ['TGC', 'TGT'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['TTT', 'TTC'],
    'G': ['GGC', 'GGT', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC'],
    'K': ['AAA', 'AAG'],
    'L': ['CTG', 'TTG', 'TTA', 'CTC', 'CTT'],
    'M': ['ATG'],
    'N': ['AAT', 'AAC'],
    'P': ['CCG', 'CCA', 'CCT'],
    'Q': ['CAG', 'CAA'],
    'R': ['CGC', 'CGT'], 
    'S': ['AGC', 'TCG', 'AGT', 'TCT', 'TCC', 'TCA'],
    'T': ['ACC', 'ACG', 'ACT'],
    'V': ['GTG', 'GTT', 'GTC', 'GTA'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    '*': ['TAA', 'TAG', 'TGA']
}


codon_freq = {'TTT': 28.9, 'TTC': 18.8, 'TTA': 17.5, 'TTG': 18.6, 'CTT': 12.7, 'CTC': 14.1, 'CTA': 3.4, 'CTG': 54.9, 'ATT': 33.9, 'ATC': 31.0, 'ATA': 5.0, 'ATG': 37.4, 'GTT': 19.6, 'GTC': 14.3, 'GTA': 10.6, 'GTG': 33.9, 'TCT': 8.5, 'TCC': 8.0, 'TCA': 6.1, 'TCG': 11.4, 'CCT': 5.8, 'CCC': 2.4, 'CCA': 7.4, 'CCG': 24.9, 'ACT': 7.7, 'ACC': 25.2, 'ACA': 6.1, 'ACG': 14.6, 'GCT': 13.8, 'GCC': 25.5, 'GCA': 19.6, 'GCG': 32.6, 'TAT': 18.6, 'TAC': 8.5, 'TAA': 1.9, 'TAG': 0.3, 'CAT': 9.3, 'CAC': 7.2, 'CAA': 13.5, 'CAG': 24.7, 'AAT': 21.2, 'AAC': 15.9, 'AAA': 29.2, 'AAG': 8.8, 'GAT': 30.0, 'GAC': 15.1, 'GAA': 29.4, 'GAG': 18.0, 'TGT': 4.2, 'TGC': 5.8, 'TGA': 0.8, 'TGG': 12.7, 'CGT': 16.4, 'CGC': 18.8, 'CGA': 2.4, 'CGG': 5.0, 'AGT': 9.0, 'AGC': 14.3, 'AGA': 2.4, 'AGG': 2.1, 'GGT': 24.4, 'GGC': 33.1, 'GGA': 8.2, 'GGG': 14.3} 
#codon_freq amounts from https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=413997

genscript_codon_freq = {'TTT': 22.1, 'TTC': 16.0, 'TTA': 14.3, 'TTG': 13.0, 'CTT': 11.9, 'CTC': 10.2, 'CTA': 4.2, 'CTG': 48.4, 'ATT': 29.8, 'ATC': 23.7, 'ATA': 6.8, 'ATG': 26.4, 'GTT': 19.8, 'GTC': 14.3, 'GTA': 11.6, 'GTG': 24.4, 'TCT': 10.4, 'TCC': 9.1, 'TCA': 8.9, 'TCG': 8.5, 'CCT': 7.5, 'CCC': 5.4, 'CCA': 8.6, 'CCG': 20.9, 'ACT': 10.3, 'ACC': 22.0, 'ACA': 9.3, 'ACG': 13.7, 'GCT': 17.1, 'GCC': 24.2, 'GCA': 21.2, 'GCG': 30.1, 'TAT': 17.5, 'TAC': 12.2, 'TAA': 2.0, 'TAG': 0.3, 'CAT': 12.5, 'CAC': 9.3, 'CAA': 14.6, 'CAG': 28.4, 'AAT': 20.6, 'AAC': 21.4, 'AAA': 35.3, 'AAG': 12.4, 'GAT': 32.7, 'GAC': 19.2, 'GAA': 39.1, 'GAG': 18.7, 'TGT': 5.2, 'TGC': 6.1, 'TGA': 1.0, 'TGG': 13.9, 'CGT': 20.0, 'CGC': 19.7, 'CGA': 3.8, 'CGG': 5.9, 'AGT': 9.9, 'AGC': 15.2, 'AGA': 3.6, 'AGG': 2.1, 'GGT': 25.5, 'GGC': 27.1, 'GGA': 9.5, 'GGG': 11.3} 
#genscript_codon_freq amounts from https://www.genscript.com/tools/codon-frequency-table, can use by editing codon_frequency

restriction_enzymes = {
    'bamhi':  'GGATCC',
    'ndei': 'CATATG', 
    'pvuii': 'CAGCTG',
    'saci': 'GAGCTC',
    'xhoi': 'CTCGAG',
    'ECRBS': 'AGGAGG',
    'bsai': 'GGTCTC',
}

def first_pass(sequence):
    #creates the first reverse translation by randomly choosing among possible codons for each position
    rev_trans = []
    for aa in sequence:
        try:
            rev_trans.append(codon_dict[aa][0])
        except KeyError:
            sys.exit('Error: sequence contains non-canonical amino acids!')
    return rev_trans

def create_sliding_windows(codon_sequence, window_len):
    #takes codon list and returns list of n-numbered nucleotide sliding windows
    dna = ''.join(codon_sequence)
    return [dna[i:i+window_len] for i in range(len(dna) - window_len + 1)]

def get_codon_pos(nt_pos):
    #takes in a nucleotide position, outputs the codon position
    return int((nt_pos - 1) / 3) + 1

def get_aa_from_codon(codon):
    #returns the aa corresponding to the input codon
    for aa in codon_dict:
        if codon in codon_dict[aa]:
            return aa

def return_new_codon(starting_codon):
    #generates the next most frequent codon some percent of the time (currently set at ~80%), otherwise generates the previous codon
    if starting_codon == 'M':
        return starting_codon, 'M'
    elif starting_codon == 'W':
        return starting_codon, 'W'
    aa = get_aa_from_codon(starting_codon)
    codon_options = codon_dict[aa]
    codon_index = codon_options.index(starting_codon)
    random_int = random.randint(1,2)
    if random_int == 2:
        try:
            new_codon = codon_options[codon_index + 1]
        except:
            new_codon = codon_options[0]
    else:
        if codon_index != 0:
            new_codon = codon_options[codon_index - 1]
        else:
            new_codon = codon_options[codon_index + 1]
    return new_codon, aa

def replace_codon(starting_codon, low_gc=False):
    #replaces codon with the next-most codon chosen by return_new_codon
    if starting_codon == 'ATG' or starting_codon == 'TGG':
        return starting_codon
    new_codon, aa = return_new_codon(starting_codon)
    count = 0
    if low_gc == True:
        starting_codon_gc = (starting_codon.count('G') + starting_codon.count('C'))
        new_codon_gc = (new_codon.count('G') + new_codon.count('C'))
        if new_codon_gc < starting_codon_gc:
            return new_codon
        else:
            while new_codon_gc >= starting_codon_gc and count < 5:
                new_codon, aa = return_new_codon(starting_codon)
                new_codon_gc = (new_codon.count('G') + new_codon.count('C'))
                count += 1
            return new_codon
    return new_codon

def find_duplicates(lst):
    #returns the duplicate sequence with the window position
    first_time = set()
    duplicates = []
    for i, item in enumerate(lst):
        if item not in first_time:
            first_time.add(item)
        else:
            duplicates.append([item, i])
    return duplicates

def replace_duplicates(codon_sequence, duplicate_len = 8):
    #takes in a list of codons and outputs a new sequence (with any remaining duplicate positions) where the duplicate positions have been randomly changed
    sliding_windows = create_sliding_windows(codon_sequence, duplicate_len)
    duplicates = find_duplicates(sliding_windows)
    for i in duplicates:
        duplicate_codon_pos = get_codon_pos(i[1])
        old_codon = codon_sequence[duplicate_codon_pos]
        new_codon = replace_codon(old_codon)
        codon_sequence[duplicate_codon_pos] = new_codon
    new_duplicates = find_duplicates(create_sliding_windows(codon_sequence, duplicate_len))
    return codon_sequence, new_duplicates

def find_high_gc(sliding_windows, gc_cutoff = 0.65):
    #returns the windows with %gc above a specified cutoff
    high_gc_windows = []
    best_gc_window = []
    for count, window in enumerate(sliding_windows):
        gc_in_window = window.count('G') + window.count('C')
        gc_ratio = float(gc_in_window) / float(len(window))
        if gc_ratio >= gc_cutoff:
            high_gc_windows.append([window, count])
    return high_gc_windows

def find_least_gc(sequences):
    #finds the lowest gc in a list
    least_gc = sequences[0]
    for count, seq in enumerate(sequences):
        least_gc_so_far = (least_gc.count('G') + least_gc.count('C'))
        percent_least_gc_so_far = float(least_gc_so_far) / float(len(least_gc))
        try:
            next_seq = sequences[count]
            next_gc_in_seq = next_seq.count('G') + next_seq.count('C')
            percent_next_gc = float(next_gc_in_seq) / float(len(next_seq))
            if percent_next_gc <= percent_least_gc_so_far:
                least_gc = next_seq
        except:
            return least_gc

def reduce_gc(codon_sequence, window_length, gc_cutoff = 0.65):
    #attempts to reduce gc content by replacing higher gc codons with lower gc codons
    sliding_windows = create_sliding_windows(codon_sequence, window_length)
    high_gc_windows = find_high_gc(sliding_windows, gc_cutoff)
    count = 0
    for i in high_gc_windows:
        first_codon_pos = get_codon_pos(i[1])
        old = codon_sequence[first_codon_pos - 1]
        new = replace_codon(old, low_gc = True)
        if (old.count('G') + old.count('C')) > (new.count('G') + new.count('C')):
            codon_sequence[first_codon_pos - 1] = new
    return codon_sequence, high_gc_windows

def reduce_gc_last_pos(codon_sequence, window_length, gc_cutoff = 0.65):
    #attempts to reduce gc content by replacing higher gc codons with lower gc codons
    sliding_windows = create_sliding_windows(codon_sequence, window_length)
    high_gc_windows = find_high_gc(sliding_windows, gc_cutoff)
    count = 0
    already_mutated = []
    for i in high_gc_windows:
        #first_codon_pos = get_codon_pos(i[1]) - 1
        last_codon_pos = get_codon_pos(i[1] + window_length - 1) - 1
        old = codon_sequence[last_codon_pos]
        if last_codon_pos not in already_mutated:
            already_mutated.append(last_codon_pos)
            new = replace_codon(old, low_gc = True)
            if (old.count('G') + old.count('C')) > (new.count('G') + new.count('C')):
                codon_sequence[last_codon_pos] = new
    return codon_sequence, high_gc_windows

def codon_frequency(codon_sequence):
    #returns a 'score' based off how often that codon is found in e. coli
    total_score = 0
    for codon in codon_sequence:
        total_score += codon_freq[codon]
        #total_score += genscript_codon_freq[codon]
    return total_score

def find_best_expressing_seq(sequences):
    #determines best codon_frequency-scoring sequence. currently only really works if all sequences are the same length, which is fine for now.
    best_seq = sequences[0]
    best_score = 0

    count = 0

    for seq in sequences:
        current_score = codon_frequency(sequences[count])
        if current_score > best_score:
            best_score = current_score
            best_seq = sequences[count]
        count += 1
    return best_seq

def return_low_freq_codons(codon_sequence):
    output = []
    for codon in codon_sequence:
        if codon_freq[codon] <= 5.0:
            output.append('Low frequency codon ' + codon + ' at position ' + str(codon_sequence.index(codon) + 1))
    return output

def change_restriction_sites(codon_sequence, enzymes):
    sites = set()
    for enzyme in enzymes:
        site = restriction_enzymes[enzyme.lower()]
        sites |= {site, reverse_complement(site)}

    for site in sites:
        width = len(site)
        dna = ''.join(codon_sequence)
        if site not in dna:
            continue
        
        for i in range(len(dna)):
            # if we encounter a restriction site
            if dna[i:i + len(site)] == site:
                # change any of these codons
                overlapped_codons = sorted(set([int((i + offset)/ 3) for offset in range(width)]))
                # accept first change that removes restriction site
                for j in overlapped_codons:
                    # change this codon
                    new_codon = replace_codon(codon_sequence[j])
                    local_dna = ''.join([new_codon if k == j else codon_sequence[k] 
                                         for k in overlapped_codons])
                    # if codon removes this site, keep it
                    if site not in local_dna:
                        codon_sequence = codon_sequence[:j] + [new_codon] + codon_sequence[j+1:]
                        break
    dna = ''.join(codon_sequence)
    for site in sites:
        assert site not in dna
    return codon_sequence

def optimize_sequence(aa_sequence):
    """Loops through several times to remove gc, then removes short stretches of very high gc, 
    finally removes duplicates.
    """
    codon_sequence = first_pass(aa_sequence)
    high_gc_windows = find_high_gc(create_sliding_windows(codon_sequence, 20))
    lower_gc_seq = codon_sequence
    for _ in range(2):
        lower_gc_seq, high_gc_windows = reduce_gc(lower_gc_seq, 20, 0.85)
        if not high_gc_windows:
            break
    for _ in range(2):
        lower_gc_seq, high_gc_windows = reduce_gc_last_pos(lower_gc_seq, 100, 0.65)
        if not high_gc_windows:
            break
    lower_duplicate_seq, duplicates = replace_duplicates(lower_gc_seq, 8)
    for _ in range(5):
        lower_duplicate_seq, duplicates = replace_duplicates(lower_duplicate_seq, 8)
        if not duplicates:
            break
        
    return lower_duplicate_seq

watson_crick = {'A': 'T',
                'T': 'A',
                'C': 'G',
                'G': 'C',
                }

def reverse_complement(seq):
    return ''.join(watson_crick[x] for x in seq)[::-1]

def main(aa_sequence, num_times_to_loop=3, enzymes=('ndei', 'xhoi', 'bamhi', 'bsai')):
    #makes multiple (# depends on num_times variable) reverse translations and returns the most codon-optimized. Breaks duplicates one last time after that.
    rev_trans_sequences = []
    for i in range(num_times_to_loop):
        restriction_opt_seq = change_restriction_sites(
            optimize_sequence(aa_sequence), enzymes)
        rev_trans_sequences.append(restriction_opt_seq)
    optimized_seq = find_best_expressing_seq(rev_trans_sequences)
    optimized_seq, high_gc_windows = reduce_gc(optimized_seq, 100, 0.68)
    #optimized_seq, high_gc_windows = reduce_gc_last_pos(optimized_seq, 80, 0.68)
    for _ in range(3):
        optimized_seq, duplicates = replace_duplicates(optimized_seq, 8)
    removed_rs_seq = change_restriction_sites(optimized_seq, enzymes)
    return ''.join(removed_rs_seq)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("usage: reverse_translate.py AA_SEQUENCE")
        exit(0)

    seq_to_translate = sys.argv[1].upper()
    random.seed(hash(seq_to_translate))

    print(main(seq_to_translate, 50))
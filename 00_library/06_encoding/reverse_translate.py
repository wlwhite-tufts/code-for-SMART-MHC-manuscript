def parse_fasta(txt):
    entries = []
    txt = '\n' + txt.strip()
    for raw in txt.split('\n>'):
        name = raw.split('\n')[0].strip()
        seq = ''.join(raw.split('\n')[1:]).replace(' ', '')
        if name:
            entries += [(name, seq)]
    return entries

def fasta_to_table(filename, name='name', sequence='sequence'):
    """Turn a fasta file into a delimited table.

    :param filename: fasta file or "stdin"
    :param name: first column name
    :param sequence: second column name
    """
    import pandas as pd
    import sys
    
    if filename == 'stdin':
        text = sys.stdin.read()
    else:
        with open(filename, 'r') as fh:
            text = fh.read()

    return pd.DataFrame(parse_fasta(text), columns=(name, sequence))

def seq_df_to_fasta(df, filename, idx_col='name', seq_col='dna'):
    """write a fasta file using sequences from a dataframe.

    :param df: dataframe containing sequences to write
    :param filename: the name fo the fasta file to create
    :param idx_col: name of the column with a unique identifier for each sequence
    :param seq_col: column name containing animo acid sequences
    """
    with open(filename,'w') as f:
        for i,row in df.iterrows():
            name = row[idx_col]
            seq = row[seq_col]
            f.write(f'>{name}\n{seq}\n')

def reverse_translate(df, idx_col='name', seq_col='sequence', enzymes=('NdeI', 'XhoI', 'BamHI', 'BsaI'), repeats=1, seed=0):
    """Reverse translate amino acids to DNA using reverse_translate_robby.

    :param df: dataframe containing sequences to reverse translate
    :param idx_col: name of the column with a unique identifier for each sequence
    :param seq_col: column name containing animo acid sequences
    :param enzymes: comma-separated list of restriction sites to avoid
    :param repeats: number of attempts, one seems OK
    :param seed: random seed, incremented for each attempt
    """
    import sys
    from reverse_translate_robby import main
    import random
    import pandas as pd

    random.seed(seed)

    dna = [main(s, repeats, enzymes=enzymes) for s in df[seq_col]]

    return pd.DataFrame({idx_col: df[idx_col], 'dna': dna})

def get_args():
    import argparse

    parser = argparse.ArgumentParser(description="Reverse translate a fasta file")
    parser.add_argument('-aa_fasta',type=str,help='The FASTA file containing the amino acid sequences to be reverse translated.')
    parser.add_argument('-dna_fasta',type=str,help='The FASTA file to write the resulting DNA sequences to.')
    return parser.parse_args()

if __name__ == '__main__':

    args = get_args()

    aa_df = fasta_to_table(args.aa_fasta)
    dna_df = reverse_translate(aa_df,enzymes=['BamHI','XhoI'])
    seq_df_to_fasta(dna_df,args.dna_fasta)

#! /usr/bin/env python3
from Bio import Seq,SeqIO
import gzip
from collections import defaultdict
import pandas as pd
import click

def create_empty_dicts():
    '''Create empty dictionaries to populate with base content information'''
    D_tail= {"A":defaultdict(int),"T":defaultdict(int),"C":defaultdict(int),"G":defaultdict(int)}
    D_content= {"A":defaultdict(int),"T":defaultdict(int),"C":defaultdict(int),"G":defaultdict(int)}
    D_bpb = {"A":defaultdict(int),"T":defaultdict(int),"C":defaultdict(int),"G":defaultdict(int)}
    return((D_tail,D_content,D_bpb))

def process_seqs(file_list, D_tail, D_content, D_bpb):
    '''Iterate theough the given fastq files and adds the information to the given dictionaries'''
    for file in file_list:
        extension = file.split('.')[-1]
        if extension == "gz":
            with gzip.open(file, "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    analyze_sequence(record.seq, D_tail, D_content, D_bpb)
        else:
            with open(file, "r") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    analyze_sequence(record.seq, D_tail, D_content, D_bpb)
    return(None)

def dict_to_df(D_tail:dict, D_content:dict, D_bpb_compositional:dict):
    '''Convert base content dictionaries to a single output dataframe'''
    out_df = pd.DataFrame({"position":range(1,150)})
    outer_dicts = [D_tail, D_content, D_bpb_compositional ]
    for i,outer_D in enumerate(outer_dicts):
        for base,D in outer_D.items():
            temp_df = pd.DataFrame.from_dict(D, orient="index")
            temp_df=temp_df.reset_index()
            if i==0:
                temp_df.columns= ["position","{}_tail_sequences".format(base)]
            elif i==1:
                temp_df.columns= ["position","{}_total_sequences".format(base)]
            elif i==2:
                temp_df.columns= ["position","{}_pct_content".format(base)]
            out_df = out_df.merge(temp_df, right_on="position", left_on="position", how="outer")
    out_df = out_df.fillna(0)
    return(out_df)

def analyze_sequence(seq, D_tail, D_content, D_bpb):
    '''Analyze the base content of a single sequence and add the results to the given dictionaries.
    Modifies the given dictionaires, does not return anything.'''
    #Set up dictionaries
    tail_D = {}
    sequence_conent_D = {"A":0, "T":0, "C":0, "G":0}
    #Get sequence info
    seq_len = len(seq)
    reverse_seq = seq[::-1]
    #Get tail base
    tail_base = reverse_seq[0]
    tail_len=0
    continue_tail = True
    #Iterate through sequence
    for i,base in enumerate(reverse_seq):
        if base not in "ATCG":
            continue
        position = seq_len - i
        #Handle Tail
        if continue_tail and base == tail_base:
            tail_len += 1
        else:
            continue_tail = False
        #Count content for this sequence
        sequence_conent_D[base] += 1
        #Handle content by position
        D_bpb[base][position] += 1
    #Modify output dictionaries
    if tail_base in "ATCG":
        #Add one count to for the tail length.
        D_tail[tail_base][tail_len] += 1
        #Add one count to tail length 0 for all the non-tail bases
        other_bases = "ATCG".replace(tail_base,"")
        for b in other_bases:
            D_tail[b][0] += 1
    for base,base_count in sequence_conent_D.items():
        D_content[base][base_count] += 1
    return(None)

def normalize_per_position_base_columns(df):
    '''Normalize the base content by position columns so they add to 1'''
    columns_to_normalize = ["A_pct_content","T_pct_content","C_pct_content","G_pct_content"]
    for i in df.position:
        position_sum = df.loc[df.position==i, columns_to_normalize].values.sum()
        for col in columns_to_normalize:
            df.loc[df.position==i, col] = df.loc[df.position==i, col] / position_sum

@click.command()
@click.option('-i',"--input", required=True, help="Input: Text file list of input fastq files. 1 file per line.")
@click.option('-o',"--output", required=True, help='Output: Name of output tsv file.')

def seq_analysis(input, output):
    #Get file list
    file_list = [line.strip() for line in open(input,'r')]
    #Make empty dictionaries
    D_tail, D_content, D_bpb = create_empty_dicts()
    #Process sequences
    process_seqs(file_list, D_tail, D_content, D_bpb)
    #Make dataframe from dictionaries
    output_dataframe = dict_to_df(D_tail, D_content, D_bpb)
    #Normalize the base content by position columns
    normalize_per_position_base_columns(output_dataframe)
    #Save output dataframe
    output_dataframe.to_csv(output, sep='\t', index=None)

if __name__ == "__main__":
    seq_analysis()

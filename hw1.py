#!/usr/bin/python
__author__ = "Chang Lu"
__email__ = "chang.lu@yale.edu"
__copyright__ = "Copyright 2021"
__license__ = "GPL"
__version__ = "1.0.0"

### Usage: python hw1.py -i <input file> -s <score file> -out <output file>
### Example:  python hw1.py -i input.txt -s blosum62.txt -out output.txt
### Note: Smith-Waterman Algorithm

import argparse
import pandas as pd
import numpy as np

### This is one way to read in arguments in Python. 
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=False,default="sample-input1.txt")
parser.add_argument('-s', '--score', help='score file', required=False,default="blosum62.txt")
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
parser.add_argument('-out', '--output', help='output file', required=False, default="output.txt")
args = parser.parse_args()


### Implement your Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap):
    # read in the inputFile and scoreFile
    with open(inputFile, "r") as f:
        seq1_print, seq2_print = f.read().split("\n")
    seq1 = list(seq1_print)
    seq2 = list(seq2_print)
    score_df = pd.read_csv(scoreFile, delim_whitespace=True, index_col=0)

    # create the empty matrix and fill in the first column and row
    M_mat = np.zeros([(len(seq2)+1),(len(seq1)+1)], dtype="int64")
    X_mat = np.zeros([(len(seq2)+1),(len(seq1)+1)], dtype="int64")
    Y_mat = np.zeros([(len(seq2)+1),(len(seq1)+1)], dtype="int64")
    final_mat = np.zeros([(len(seq2)+1),(len(seq1)+1)], dtype="int64")

    # Fill the matrix
    for i in np.arange(1,len(seq2)+1):
        for j in np.arange(1,len(seq1)+1):
            X_mat[i,j] = max(X_mat[i-1, j] + int(extGap), final_mat[i-1, j] + int(openGap), 0)
            Y_mat[i,j] = max(Y_mat[i, j-1] + int(extGap), final_mat[i, j-1] + int(openGap), 0)
            match_score = score_df.loc[seq2[i-1], seq1[j-1]]
            M_mat[i,j] = max(final_mat[i-1, j-1] + match_score,X_mat[i-1, j-1]+match_score,Y_mat[i-1, j-1]+match_score,  0)
            final_mat[i,j] = max(M_mat[i,j],X_mat[i,j],Y_mat[i,j])


    # calculate the alignment score
    alignment_score = np.max(final_mat)

    # tracing back
    save_seq1 = []
    save_seq2 = []
    start_ind = np.where(final_mat == np.max(final_mat))
    ini_seq2_ind = row_ind = start_ind[0][0]
    ini_seq1_ind = col_ind = start_ind[1][0]
    ini_seq2_ind -= 1
    ini_seq1_ind -= 1
    while row_ind>0 and col_ind>0:
        if final_mat[row_ind,col_ind] == 0:
            break
        elif final_mat[row_ind,col_ind] == M_mat[row_ind,col_ind]:
            save_seq1.append(seq1[col_ind-1])
            save_seq2.append(seq2[row_ind-1])
            row_ind -= 1
            col_ind -= 1
        elif final_mat[row_ind,col_ind] == X_mat[row_ind,col_ind]:
            save_seq1.append("-")
            save_seq2.append(seq2[row_ind-1])
            row_ind -= 1
        elif final_mat[row_ind,col_ind] == Y_mat[row_ind,col_ind]:
            save_seq1.append(seq1[col_ind-1])
            save_seq2.append("-")
            col_ind -= 1
    final_seq1_ind = max(0,col_ind)
    final_seq2_ind = max(0,row_ind)

    save_seq1.reverse()
    save_seq2.reverse()

    # deal with the sequence for printing
    difference = abs(final_seq1_ind - final_seq2_ind)
    if final_seq1_ind > final_seq2_ind:
        seq1_out = seq1[:final_seq1_ind] + ["("] + save_seq1 + [")"] + seq1[ini_seq1_ind+1:]
        seq2_out = [" "]*difference + seq2[:final_seq2_ind] + ["("] + save_seq2 + [")"] + seq2[ini_seq2_ind+1:]
    else:
        seq1_out = [" "]*difference + seq1[:final_seq1_ind] + ["("] + save_seq1 + [")"] + seq1[ini_seq1_ind+1:]
        seq2_out = seq2[:final_seq2_ind] + ["("] + save_seq2 + [")"] + seq2[ini_seq2_ind+1:]
    middle = [" "]*(max(final_seq1_ind,final_seq2_ind)+1)
    for i in np.arange(0,len(save_seq1)):
        if save_seq1[i] == save_seq2[i]:
            middle.append("|")
        else:
            middle.append(" ")

    seq1_out = ''.join(seq1_out)
    seq2_out = ''.join(seq2_out)
    middle = ''.join(middle)

    # Generate output
    output = """-----------
|Sequences|
-----------
sequence1
{seq1}
sequence2
{seq2}
--------------
|Score Matrix|
--------------
    """.format(seq1=seq1_print, seq2=seq2_print)

    with open(args.output, "w") as outfile:
        outfile.write(output)

    final_df = pd.DataFrame(final_mat)
    final_df.insert(0, column='', value=[""] +seq2, allow_duplicates=True)
    final_df.columns = [""] + seq1 +[""]
    #final_df.insert(-1, column='', value=[""]*(len(seq2)+1), allow_duplicates=True)
    #final_df.insert(len(final_df.columns), "", "" * len(final_df),allow_duplicates=True)

    with open(args.output, 'a') as f:
        final_df.to_csv(f, sep='\t', index=False)



    output = """----------------------
|Best Local Alignment|
----------------------
Alignment Score:{score}
Alignment Results:
{seq1_out}
{middle}
{seq2_out}
""".format(score=alignment_score, seq1_out=seq1_out, middle=middle, seq2_out=seq2_out)


    with open(args.output, "a") as outfile:
        outfile.write(output)


### Run Smith-Waterman Algorithm
runSW(args.input, args.score, args.opengap, args.extgap)

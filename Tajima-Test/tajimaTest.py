#!/usr/bin/env python

import argparse
import math
from Bio import SeqIO


# Get sequences from the FASTA file
def getSeqs(f):
    sequences = []
    for seq_record in SeqIO.parse(f, "fasta"):
        sequences.append(str(seq_record.seq))
    return sequences


# Finds the pairwise difference between two sequences
def differences(s1, s2):
    count = 0
    for n in range(0, len(s1)):
        if s1[n] != s2[n]:
            count += 1
    return int(count)


def run(args):
    # Loading in the sequences
    seqs = getSeqs(args.input)
    # Calculating number and length of sequences
    N = len(seqs)
    nucleotides = len(seqs[0])
    # Calculate Segregating sites:
    Sn = nucleotides
    for nuc in range(0, nucleotides):
        if seqs[0][nuc] == seqs[1][nuc] == seqs[2][nuc] == seqs[3][nuc] == seqs[4][nuc]:
            Sn -= 1
    print("The number of segregation sites is =", Sn)
    # Calculating harmonic sum:
    harmonic_sum = sum([1 / i for i in range(1, N)])
    # Calculating Theta:
    theta = Sn / harmonic_sum
    # Calculates pi:
    N_choose_2 = N * (N - 1) / 2
    diffs = []
    for i in range(0, N - 1):  # 0-4
        for j in range(i + 1, N):  # i+1 - 4
            kij = differences(seqs[i], seqs[j])
            diffs.append(kij)
    numerator = sum(diffs)
    pi = numerator / N_choose_2
    # Calculates the normalization constant C
    a1 = sum([1 / i for i in range(1, N)])
    a2 = sum([1 / i ** 2 for i in range(1, N)])
    b1 = (N + 1) / (3 * (N - 1))
    b2 = 2 * (N ** 2 + N + 3) / (9 * (N - 1))
    c1 = b1 - (1 / a1)
    c2 = b2 - ((N + 2) / (a1 * N)) + (a2 / a1 ** 2)
    x = (c1 / a1 * Sn) + (c2 / ((a1 ** 2) + a2)) * Sn * (Sn - 1)
    C = math.sqrt(x)
    # Calculates the Tajima statistic D:
    Dt = (pi - theta) / C
    print("Tajima D statistic =", round(Dt, 5))


def main():
    parser = argparse.ArgumentParser(description="Performs a Tajima D Test",
                                     epilog="This tool was created by Tyler Aprati")
    parser.add_argument("-i", help="FASTA input file", dest="input", type=str, required=True)
    parser.set_defaults(func=run)
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()

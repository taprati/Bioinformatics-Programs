# Tyler Aprati
#
# Performs a tajima test on sequences
#

import math

def differences(s1,s2):     # Finds the pairwise difference between two sequences
    count = 0
    for n in range(0,len(s1)):
        if s1[n] != s2[n]:
            count += 1
    return int(count)

# Loading in the sequences

##s1 = "GATACGGACAAGGCTGTTGTTAAGGCCATCTGGGCTAAGATCAGCCCCAAGGCCGATGAA"
##s2 = "CCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAG"
##s3 = "GCCGCCGACAAGGGCAATGTCAAGGCCGCCTGGGGCAAGGTTGGCGGCCACGCTGCAGAG"
##s4 = "GGGGAAGACAAAAGCAACATCAAGGCTGCCTGGGGGAAGATTGGTGGCCATGGTGCTGAA"
##s5 = "GCAGATGACAAAACCAACATCAAGAACTGCTGGGGGAAGATTGGTGGCCATGGTGGTGAA"
##seqs = [s1,s2,s3,s4,s5]

with open("tajima_sequences.txt") as file:
    seqs = [line.rstrip() for line in file]

# Calculating number and length of sequences
N = len(seqs)
nucleotides = len(seqs[0])

#Calculate Segregating sites:
Sn = nucleotides
for nuc in range(0,nucleotides):
    if seqs[0][nuc] == seqs[1][nuc] == seqs[2][nuc] == seqs[3][nuc] == seqs[4][nuc]:
        Sn -= 1
print("The number of segregation sites is =",Sn)

# Calculating harmonic sum:
harmonic_sum = sum([ 1/i for i in range(1,N)])
print("Harmonic Sum =",harmonic_sum)

# Calculating Theta:
theta = Sn / harmonic_sum
print("Theta =",theta)

# Calculates pi:
N_choose_2 = N*(N-1)/2
diffs = []
for i in range(0,N-1): # 0-4
    for j in range(i+1,N): # i+1 - 4
        kij = 0
        kij = differences(seqs[i],seqs[j])
        diffs.append(kij)
print("Completed correctly =",(len(diffs) == N_choose_2))
numerator = sum(diffs)
pi = numerator / N_choose_2
print("Pi =",pi)

# Calculates the normalization constanst C
a1 = sum([ 1/i for i in range(1,N)])
a2 = sum([ 1/i**2 for i in range(1,N)])
b1 = (N+1) / (3*(N-1))
b2 = 2 * (N**2 + N + 3) / (9*(N-1))
c1 = b1 - (1/a1)
c2 = b2 - ((N+2)/(a1*N)) + (a2/a1**2)
x = (c1/a1 * Sn) + (c2 / ((a1**2) + a2)) * Sn * (Sn-1)

C = math.sqrt(x)
print("Normalization Constanct C =",C)

# Calculates the Tajima statistic D:
Dt = (pi - theta) / C
print("Tajima statistic D =",Dt)


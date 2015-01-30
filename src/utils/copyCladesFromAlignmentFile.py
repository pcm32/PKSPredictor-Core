import os
import sys

__author__ = 'pmoreno'

# read the alignment file, keep all seqs in a map, read the clades files, for each id, print it from
# the alignment file into a new clade file.

if __name__ == "__main__":
    # read the overall alignment
    from Bio import SeqIO
    handle = open(sys.argv[1], "rU")
    sequences = {}
    for record in SeqIO.parse(handle, "fasta") :
        sequences[record.id] = record;
    handle.close()

    listing = os.listdir(sys.argv[2])

    for item in listing:
        if item.endswith(".fas"):
            handle = open(sys.argv[2]+item)
            handleCladeOut = open(sys.argv[2]+item+".new","w")
            for record in SeqIO.parse(handle, "fasta") :
                if record.id in sequences :
                    SeqIO.write(sequences.get(record.id), handleCladeOut, "fasta")
            handleCladeOut.close()
            handle.close()



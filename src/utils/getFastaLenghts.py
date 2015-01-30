import sys

__author__ = 'pmoreno'


if __name__ == "__main__":
    from Bio import SeqIO
    handle = open(sys.argv[1], "rU")
    for record in SeqIO.parse(handle, "fasta") :
        print record.id+"\t"+str(len(record.seq))

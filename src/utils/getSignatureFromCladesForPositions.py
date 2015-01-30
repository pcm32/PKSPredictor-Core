#! /usr/bin/python
import sys
import os
# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="pmoreno"
__date__ ="$Jul 14, 2010 3:42:49 PM$"

if __name__ == "__main__":
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    os.chdir(sys.argv[1])
    positionStr = sys.argv[2]
    listOfPos = positionStr.split(",");
    listing = os.listdir(".")
    consensus = {}
    consensusThres = 0.7

    print "Clade",
    for position in listOfPos:
        print position,

    print

    for item in listing:
        if item.endswith(".fas"):
            print item,
            alignments = AlignIO.parse(item,"fasta")
            for alignment in alignments:
                summ = AlignInfo.SummaryInfo(alignment)
                consensus[item] = summ.gap_consensus(consensusThres)
                for position in listOfPos:
                    print consensus[item][int(position)],
                print

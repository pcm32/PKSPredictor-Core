#! /usr/bin/python
import sys
import os
# To change this template, choose Tools | Templates
# and open the template in the editor.
#
# This script requires as first argument the directory where all the individual clades are, and as second
# argument the path to the complete (all clade) alignments.

__author__="pmoreno"
__date__ ="$Mar 28, 2011 4:17:40 PM$"

if __name__ == "__main__":

    from Bio.Align import AlignInfo, AlignIO
    from Bio.Align import MultipleSeqAlignment
    from Bio.Alphabet import IUPAC, Gapped
    from Bio.Seq import SeqIO
    from Bio.SeqRecord import SeqRecord
    # Directory where files are
    os.chdir(sys.argv[1])
    listing = os.listdir(".")
    consensus = {}
    genConsensus = ''
    pssmGen = ''
    # this value should be read from the arguments or else use a default
    consensusThres = 0.7 
    # sys.argv[2] holds the path to the general alignment
    generalAlignment = AlignIO.parse(sys.argv[2],"fasta",alphabet=Gapped(IUPAC.ExtendedIUPACProtein(),"-"))
    lengthGenAl = 0
    positionsToMask = []
    for genAlignment in generalAlignment:
        sumGen = AlignInfo.SummaryInfo(genAlignment)
        genConsensus = sumGen.gap_consensus(consensusThres)
        for index, residue in enumerate(genConsensus):
            if genConsensus[index] == '-':
                continue
            if genConsensus[index] == 'X':
                continue
            positionsToMask.append(index)
        #pssmGen = sumGen.pos_specific_score_matrix(genConsensus,chars_to_ignore = ['-'])
        pssmGen = sumGen.pos_specific_score_matrix(genConsensus)
        lengthGenAl = len(genAlignment)

    print positionsToMask

    for item in listing:
        if item.endswith(".fas"):
            #alignments = AlignIO.parse(item,"fasta",alphabet=IUPAC.ExtendedIUPACProtein())
            alignments = AlignIO.parse(item, "fasta", alphabet=Gapped(IUPAC.ExtendedIUPACProtein(), "-"))
            for alignment in alignments:
                summ = AlignInfo.SummaryInfo(alignment)
                consensus[item] = summ.gap_consensus(consensusThres)
                for posToMask in positionsToMask:
                    if consensus[item][posToMask] == '-':
                        continue
                    for alignElement in alignment:
                        mutSeq = alignElement.seq.tomutable()
                        mutSeq[posToMask] = 'X'
                        alignElement.seq = mutSeq.toseq()
                SeqIO.write(alignment, item+"_noPKSsignal.faa", "fasta")
                summ = AlignInfo.SummaryInfo(alignment)
                consensus[item] = summ.gap_consensus(consensusThres)
                print item, consensus[item]

                

    

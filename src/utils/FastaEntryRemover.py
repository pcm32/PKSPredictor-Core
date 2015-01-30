#! /usr/bin/python
from Bio.Seq import Seq, SeqRecord
import sys
# This script should read the alignment file of a clade, extract one element
# write it in a separate fasta file (removing the gaps). Then it should also
# read the version of the clade that doesn't have the general PKS signal,
# extract the same element, and write to another fasta file the rest of the
# entries. The first fasta file written, the only only holding one sequence,
# should be the query fasta file. The second, bigger file, should be processed
# by the HMMMER tool to create a profile to search against.

__author__="pmoreno"
__date__ ="$May 29, 2011 5:16:28 PM$"

if __name__ == "__main__":
    #dirOfHMMModels = sys.argv[1]
    fastaFileCladeNoGeneralSignal = sys.argv[1]
    fastaFileClade = sys.argv[2]
    entryToTest = int(sys.argv[3])
    resultFolder = sys.argv[4]

    from Bio import AlignIO, SeqIO
    from Bio.Alphabet import IUPAC, Gapped
    from Bio.Align import MultipleSeqAlignment

    alignmentNoGenSignalIterator = AlignIO.parse(fastaFileCladeNoGeneralSignal,"fasta",alphabet=Gapped(IUPAC.ExtendedIUPACProtein(),"-"));
    alignmentIterator = AlignIO.parse(fastaFileClade,"fasta",alphabet=Gapped(IUPAC.ExtendedIUPACProtein(),"-"));

    noGenSignalAlignment = alignmentNoGenSignalIterator.next()
    queryFasta = resultFolder+"/"+"Query_%d.faa" % (entryToTest,)
    ownCladeProfile = resultFolder+"/"+"ForOwnCladeProfile_%d.faa" % (entryToTest,)
    #print testAlignment[entryToTest].id
    #print testAlignment[entryToTest].seq

    alignmentWithSignal = alignmentIterator.next()
    desiredSeqString = str(alignmentWithSignal[entryToTest-1].seq)
    desiredSeqString = desiredSeqString.replace("-", "")
    #print desiredSeqString
    seqNoGaps = Seq(desiredSeqString, alphabet=IUPAC.ExtendedIUPACProtein())
    #print seqNoGaps
    seqRecNoGaps = SeqRecord(seq=seqNoGaps, id=alignmentWithSignal[entryToTest-1].id)
    #print seqRecNoGaps.seq
    #print seqRecNoGaps.id
    SeqIO.write(seqRecNoGaps, queryFasta, "fasta")
    
    #print "Number of entries: ",len(testAlignment)
    # Here we remove the desired element
    newTestAlignment = []
    for i in range(len(noGenSignalAlignment)):
        if i != entryToTest-1 :
            newTestAlignment.append(noGenSignalAlignment[i])

    newAlignment = MultipleSeqAlignment(newTestAlignment)

    AlignIO.write(newAlignment, ownCladeProfile,"fasta")
    #print "Number of entries after: ",len(newTestAlignment)



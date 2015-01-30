import sys
import os
# To change this template, choose Tools | Templates
# and open the template in the editor.

__author__="pmoreno"
__date__ ="$May 9, 2010 4:21:18 AM$"

if __name__ == "__main__":
    from Bio import AlignIO
    from Bio.Align import AlignInfo
    from Bio.Align import MultipleSeqAlignment
    from Bio.Alphabet import IUPAC, Gapped
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    os.chdir(sys.argv[1])
    listing = os.listdir(".")
    consensus = {}
    genConsensus = '';
    pssmGen = '';
    consensusThres = 0.7

    #generalAlignment = AlignIO.parse(sys.argv[2],"fasta",alphabet=IUPAC.ExtendedIUPACProtein())
    generalAlignment = AlignIO.parse(sys.argv[2],"fasta",alphabet=Gapped(IUPAC.ExtendedIUPACProtein(),"-"))
    lengthGenAl = 0
    for genAlignment in generalAlignment:
        sumGen = AlignInfo.SummaryInfo(genAlignment)
        genConsensus = sumGen.gap_consensus(consensusThres)
        #pssmGen = sumGen.pos_specific_score_matrix(genConsensus,chars_to_ignore = ['-'])
        pssmGen = sumGen.pos_specific_score_matrix(genConsensus)
        lengthGenAl = len(genAlignment)

    for item in listing:
        if item.endswith(".fas"):
            #alignments = AlignIO.parse(item,"fasta",alphabet=IUPAC.ExtendedIUPACProtein())
            alignments = AlignIO.parse(item,"fasta",alphabet=Gapped(IUPAC.ExtendedIUPACProtein(),"-"))
            for alignment in alignments:
                summ = AlignInfo.SummaryInfo(alignment)
                consensus[item] = summ.gap_consensus(consensusThres)

    consensusAlignment = MultipleSeqAlignment([],alphabet=Gapped(IUPAC.ExtendedIUPACProtein(),"-"))
    for seqKey in consensus.keys():
        #consensusAlignment.append(SeqRecord(Seq(consensus[seqKey], Gapped(IUPAC.ExtendedIUPACProtein(),"-")), id=seqKey))
        consensusAlignment.append(SeqRecord(consensus[seqKey]))
    summConsensusAlign = AlignInfo.SummaryInfo(consensusAlignment)

    print "Clade","AlignedSeqs","Pos","NumConflicts","CladeCons","Support","GenConsensus","GenConsSupp","DiffLetterSeen","InfCont","InfConGen","InfConConsensus"

    for item in listing:
        if item.endswith(".fas"):
            #alignments = AlignIO.parse(item,"fasta",alphabet=IUPAC.ExtendedIUPACProtein())
            alignments = AlignIO.parse(item,"fasta",alphabet=Gapped(IUPAC.ExtendedIUPACProtein(),"-"))
            for alignment in alignments:
                summ = AlignInfo.SummaryInfo(alignment)
                cladeConsensus = summ.gap_consensus(consensusThres)
                #for otherClade in consensus.keys():
                #    if item != otherClade:
                #pssms = summ.pos_specific_score_matrix(genConsensus)#,
                                            #chars_to_ignore = ['N','X'])
                #pssms = summ.pos_specific_score_matrix(genConsensus,chars_to_ignore = ['-'])
                pssms = summ.pos_specific_score_matrix(genConsensus)
                #print item
                countImportantRes = 0
                relevantPositions = dict();
                setOfLettersSeen = dict();
                for otherClade in consensus.keys():
                    if otherClade == item:
                        continue
                    for index, residue in enumerate(consensus[otherClade]):
                        if cladeConsensus[index] == '-':
                                continue
                        if cladeConsensus[index] == 'X':
                                continue
                        if residue == 'X' or residue == '-':
#                            if (pssms[index][cladeConsensus[index]] / len(alignment)) > consensusThres:
#                                #print "Position ",index," specific for"," res:",cladeConsensus[index]," general res:",residue
#                                #countImportantRes+=1
#                                if not index in relevantPositions.keys():
#                                    relevantPositions[index]=1
#                                else:
#                                    relevantPositions[index]+=1
                            continue
                        if pssmGen[index][residue]/lengthGenAl > 0.05:
                            if not index in setOfLettersSeen.keys():
                                setOfLettersSeen[index] = list();
                            if not residue in setOfLettersSeen[index]:
                                setOfLettersSeen[index].append(residue)
                        score = pssms[index][residue] / len(alignment)
                        if score < consensusThres:
                            if (pssms[index][cladeConsensus[index]] / len(alignment)) > consensusThres:
                                #print "Position ",index," specific for"," res:",cladeConsensus[index]," general res:",residue
                                #countImportantRes+=1
                                if not index in relevantPositions.keys():
                                    relevantPositions[index]=1
                                else:
                                    relevantPositions[index]+=1

                    #print "Important residues: ",countImportantRes;

                

                for position in relevantPositions.keys():
                    if cladeConsensus[position] != genConsensus[position] or relevantPositions[position] > 10:
                        seenLettersInPos = 0
                        if position in setOfLettersSeen.keys():
                            seenLettersInPos = len(setOfLettersSeen[position])
                        genSupportAuxInPos = 0;
                        if genConsensus[position] != '-':
                            genSupportAuxInPos = pssmGen[position][genConsensus[position]]
                        print item,len(alignment),
                        print position,relevantPositions[position],cladeConsensus[position],pssms[position][cladeConsensus[position]]/len(alignment),
                        print genConsensus[position],genSupportAuxInPos/lengthGenAl,seenLettersInPos,
                        print summ.information_content(position, position+1, chars_to_ignore = ['-','X']),
                        print sumGen.information_content(position, position+1, chars_to_ignore = ['-','X']),
                        print summConsensusAlign.information_content(position,position+1, chars_to_ignore = ['-','X'])
                #print
                #print item
                #print pssms
                #print

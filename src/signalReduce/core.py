'''
Created on Sep 21, 2011

@author: pmoreno
'''
import sys
import os

class SignalReducer(object):
    '''
    classdocs
    '''


    def __init__(self,pathToCladesAlignments,generalAlignment,outputPath):
        '''
        Constructor
        '''
        self.pathToCladesAlignments = pathToCladesAlignments
        self.generalAlignment = generalAlignment
        self.outPutPath = outputPath
        
        
    def run(self,consensusThreshold):
        from Bio import AlignIO, SeqIO
        from Bio.Align import AlignInfo
        # from Bio.Align import MultipleSeqAlignment
        from Bio.Alphabet import IUPAC, Gapped
        # from Bio.Seq import Seq
        # from Bio.SeqRecord import SeqRecord
        # Directory where files are
        # os.chdir(sys.argv[1])
        # listing = os.listdir(".")
        listing = os.listdir(self.pathToCladesAlignments)
        consensus = {}
        genConsensus = '';
        pssmGen = '';
        # this value should be read from the arguments or else use a default
        consensusThres = consensusThreshold
        # sys.argv[2] holds the path to the general alignment
        generalAlignment = AlignIO.parse(self.generalAlignment,"fasta",alphabet=Gapped(IUPAC.ExtendedIUPACProtein(),"-"))
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
        print listing
    
        resultAlignFiles = []
        for item in listing:
            if item.endswith(".fas"):
                #alignments = AlignIO.parse(item,"fasta",alphabet=IUPAC.ExtendedIUPACProtein())
                alignments = AlignIO.parse(self.pathToCladesAlignments+item,"fasta",alphabet=Gapped(IUPAC.ExtendedIUPACProtein(),"-"))
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
                    SeqIO.write(alignment, self.outPutPath+item+"_noPKSsignal_Thres%d.faa" % (consensusThres*100,), "fasta")
                    resultAlignFiles.append(self.outPutPath+item+"_noPKSsignal_Thres%d.faa" % (consensusThres*100,))
                    summ = AlignInfo.SummaryInfo(alignment)
                    consensus[item] = summ.gap_consensus(consensusThres)
                    print item, consensus[item]
        return resultAlignFiles
                    
        
        
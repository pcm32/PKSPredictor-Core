'''
Created on Sep 21, 2011

@author: pmoreno
'''

class FastaGapRemover(object):
    '''
    classdocs
    '''
    def __init__(self,pathToFasta):
        '''
        Constructor
        '''
        self.pathToFasta = pathToFasta
    
    def writeNoGaps(self,outputfile):
        fastaFile = self.pathToFasta
        
        from Bio import AlignIO, SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import IUPAC, Gapped
        
        seqIterator = SeqIO.parse(fastaFile,"fasta",alphabet=Gapped(IUPAC.ExtendedIUPACProtein(),"-"));
        outFasta = outputfile
        records = list()
        for alignment in seqIterator:
            desiredSeqString = str(alignment.seq)
            desiredSeqString = desiredSeqString.replace("-", "")
            #print desiredSeqString
            seqNoGaps = Seq(desiredSeqString,alphabet=IUPAC.ExtendedIUPACProtein())
        
            seqRecNoGaps = SeqRecord(seq=seqNoGaps, id=alignment.id)
            #print seqRecNoGaps.seq
            #print seqRecNoGaps.id
            records.append(seqRecNoGaps)
        
        SeqIO.write(records, outFasta, "fasta")
        

class FastaEntryRemover(object):
    '''
    classdocs
    '''


    def __init__(self,pathToFasta):
        '''
        Constructor
        '''
        self.pathToFasta = pathToFasta
        
    def getFastaEntry(self,entryNum,resultFolder):
        fastaFile = self.pathToFasta
        
        from Bio import AlignIO, SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio.Alphabet import IUPAC, Gapped
        
        alignmentIterator = AlignIO.parse(fastaFile,"fasta",alphabet=Gapped(IUPAC.ExtendedIUPACProtein(),"-"));
        queryFasta = resultFolder+"/"+"Query_%d.faa" % (entryNum,)
        alignment = alignmentIterator.next()
        if entryNum > len(alignment):
            return None
        
        desiredSeqString = str(alignment[entryNum-1].seq)
        desiredSeqString = desiredSeqString.replace("-", "")
        #print desiredSeqString
        seqNoGaps = Seq(desiredSeqString,alphabet=IUPAC.ExtendedIUPACProtein())
        
        seqRecNoGaps = SeqRecord(seq=seqNoGaps, id=alignment[entryNum-1].id)
        #print seqRecNoGaps.seq
        #print seqRecNoGaps.id
        SeqIO.write(seqRecNoGaps, queryFasta, "fasta")
        return queryFasta
        
    def generateFastaWithOutEntry(self,entryNumToRemove,resultFolder):
        fastaFile = self.pathToFasta
        # fastaFileClade = sys.argv[2]
        # entryToTest = int(sys.argv[3])
        entryToTest = int(entryNumToRemove)

        from Bio import AlignIO
        from Bio.Alphabet import IUPAC, Gapped
        from Bio.Align import MultipleSeqAlignment
    
        alignmentIterator = AlignIO.parse(fastaFile,"fasta",alphabet=Gapped(IUPAC.ExtendedIUPACProtein(),"-"));
        
    
        alignment = alignmentIterator.next()
        pathToNewFile = resultFolder+"/"+"WithOutEntry_%d.faa" % (entryToTest,)
        #print testAlignment[entryToTest].id
        #print testAlignment[entryToTest].seq
        #print "Number of entries: ",len(testAlignment)
        # Here we remove the desired element
        newTestAlignment = []
        for i in range(len(alignment)):
            if i != entryToTest-1 :
                newTestAlignment.append(alignment[i])
    
        newAlignment = MultipleSeqAlignment(newTestAlignment)
    
        AlignIO.write(newAlignment, pathToNewFile,"fasta")
        #print "Number of entries after: ",len(newTestAlignment)
        return pathToNewFile

class GBKRankingBasedFeatureRemover(object):

    def __init__(self,maxRanking):
        self.maxRanking = maxRanking


    def processGBK(self,seqRecord):
        listToRemove = []
        for seqFeat in seqRecord.features:
            if seqFeat.type == "domain" and "ranking" in seqFeat.qualifiers:
                if int(seqFeat.qualifiers["ranking"]) > self.maxRanking:
                    listToRemove.append(seqFeat)

        for seqFeat in listToRemove:
            seqRecord.features.remove(seqFeat)

        return seqRecord


class SeqRecordAnnotator(object):

    def __init__(self,modelAnnotator):
        self.modelAnnotator = modelAnnotator

    def annotateSeqFeatures(self,seqRecord):
        for feature in seqRecord.features:
            clade = feature.qualifiers.get("name",None)
            if clade is not None:
                annot = self.modelAnnotator.getAnnotationForKey(clade)
                if annot is not None:
                    feature.qualifiers["note"] = annot



    
    
        
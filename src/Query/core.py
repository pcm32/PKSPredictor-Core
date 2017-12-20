from os.path import isfile

from Bio import SeqIO
from Bio.Alphabet import IUPAC
from hmmer.core import HMMERSuite, HMMERParse, HMMERModelBasedCutoffStore
import re
from Bio.SeqFeature import SeqFeature, FeatureLocation
from EMBOSS.core import FuzzProRunner
import tempfile



class QueryRunner(object):
    
    def __init__(self, fastaQuery, HMMModel, outPutPath,
                 HMMBinaryPath, fuzzProPath, NRPS2Path, allDomains, HMMModelNonClades):
        self.fastaQuery = fastaQuery
        self.HMMModel = HMMModel
        self.HMMModelNonClades = HMMModelNonClades
        self.outPutPath = outPutPath
        self.HMMBinaryPath = HMMBinaryPath
        self.fuzzProPath = fuzzProPath
        self.skipClade = list()
        self.useCutOff = False
        self.nrps2Path = NRPS2Path
        self.allDomains = allDomains

    def useModelCutOff(self):
        self.useCutOff = True
            
    def __setAlphabet(self,seq):
        if re.match("[ATCG]+$", str(seq)):
            seq.alphabet = IUPAC.unambiguous_dna
        else:
            seq.alphabet = IUPAC.protein
    
    def __generalLocationContainsCladeLoc(self,generalLoc,cladeLocation):
        if generalLoc.start <= cladeLocation.start and cladeLocation.end <= generalLoc.end:
            return True
        else:
            return False
                
    def setCladesToSkip(self,cladesToSkip):
        self.skipClade=cladesToSkip
    
    def runQuery(self):
        listOfAASeqs = list()
        for seq_record in SeqIO.parse(self.fastaQuery,"fasta"):
            self.__setAlphabet(seq_record.seq)
            # print "Processing "+seq_record.id+" alphabet "+str(seq_record.seq.alphabet)
            # if the sequence is a nucleotide, get 6 reading frames into AA.
            if seq_record.seq.alphabet == IUPAC.unambiguous_dna:
            #if seq_record.seq.alphabet != IUPAC.ExtendedIUPACProtein:
                try:
                    for i in range(3):
                        seqAA = seq_record[i::3].translate(table=11)
                        seqAA.id = seq_record.id+"_T%" % i
                        listOfAASeqs.append(seqAA)
                        seqAA = seq_record.reverse_complement()[i::3].translate(table=11)
                        seqAA.id = seq_record.id+"_RT%" % i
                        listOfAASeqs.append(seqAA)
                        # by default we use the bacterial traslation table
                except AttributeError:
                    print "Skipping "+seq_record.id+" due to wrongly detected alphabet "+str(seq_record.seq.alphabet)
            else:
                listOfAASeqs.append(seq_record)
                 
        # hmmerExec = HMMERSuite(path=self.HMMBinaryPath)
        sequenceMarkers = list()
        sequenceMarkers.append(CladeFeatureMarker(self.outPutPath, self.HMMModel, self.HMMBinaryPath, self.skipClade,
                                                  self.useCutOff, self.allDomains))
        sequenceMarkers.append(NRPSPred2FeatureMarker(self.nrps2Path))
        if self.HMMModelNonClades is not None:
            sequenceMarkers.append(HMMERHitFeatureMarker(self.HMMModelNonClades, self.HMMBinaryPath))

        if self.fuzzProPath:
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath,
                                                        "[DPACQRTELNHVSKYGF][IVLMAP][QDHNRSTGKLAEC][GAVT][VAIL][HFILVM][YHFNVQ][ILSATGMCFNV][VAGTPS][PGILMVATFR]X(3)D[DVNEGTLKRQASHPC][KTASPRGLFVENDYMIH][VILMFAC][QELATRHKSNGVFDMYIP][QNVKRGSEDATLHYM][HMTLIKSARNQF][KGTSPNDEARQH][EIPRQTWAFLMDVSKYHG][PAEVDTKRQNSHFG][LADQETRHSKGINM][KWFLYVIMARSGT][DRNHISQEVTACKYLMP]", "KR D-configured OH groups"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath,
                                                        "[DPACQRTELNHVSKYGF][IVLMAP][QDHNRSTGKLAEC][GAVT][VAIL][HFILVM][YHFNVQ][ILSATGMCFNV][VAGTPS][PGILMVATFR]X(3){D}[DVNEGTLKRQASHPC][KTASPRGLFVENDYMIH][VILMFAC][QELATRHKSNGVFDMYIP][QNVKRGSEDATLHYM][HMTLIKSARNQF][KGTSPNDEARQH][EIPRQTWAFLMDVSKYHG][PAEVDTKRQNSHFG][LADQETRHSKGINM][KWFLYVIMARSGT][DRNHISQEVTACKYLMP]", "KR L-configured OH groups"))

        self.listOfSeqRecords = list()
        print "Starting processing..."

            
        for seqAA in listOfAASeqs:
            # For each amino acid sequence, we run HMMER and mark
            print "Processing "+seqAA.id

            for seqMarker in sequenceMarkers:
                seqMarker.markFeaturesInSequenceRecord(seqAA)

            # Add additional motifs/domains marks, like the methytransferase trough fuzzpro

            # This can be written at the outside to GBK
            self.listOfSeqRecords.append(seqAA)
        print "Finished processing"
    
    def getSeqsWithResults(self):
        return self.listOfSeqRecords

from abc import ABCMeta, abstractmethod

class FeatureMarker:
    __metaclass__ = ABCMeta

    @abstractmethod
    def markFeaturesInSequenceRecord(self, seqRecord):
        pass

class PatternFeatureMarker(FeatureMarker):

    def __init__(self, fuzzProPath, pattern, patternName):
        self.runner = FuzzProRunner(fuzzProPath)
        self.runner.setPattern(pattern, patternName)

    def markFeaturesInSequenceRecord(self, seqRecord):
        self.runner.run(seqRecord)

class NRPSPred2FeatureMarker(FeatureMarker):

    def __init__(self, nrps2Path):
        from NRPS2Pred.core import NRPS2PredRunner
        self.runner = NRPS2PredRunner(nrps2Path)

    def markFeaturesInSequenceRecord(self, seqRecord):
        self.runner.run(seqRecord)

class HMMERHitFeatureMarker(FeatureMarker):

    def __init__(self, HMMModel, HMMBinaryPath):
        self.outPutPath = tempfile.gettempdir()+"/"
        self.HMMModel = HMMModel
        self.hmmerExec = HMMERSuite(path=HMMBinaryPath)

    def markFeaturesInSequenceRecord(self, seqAA):
        indfastaPath=self.outPutPath+"Query_"+seqAA.id+".faa"
        SeqIO.write(sequences=seqAA, handle=indfastaPath, format="fasta")
        scanOutputPath=self.outPutPath+seqAA.id+"_hmmerRes.txt"
        self.hmmerExec.runScan(query=indfastaPath, output=scanOutputPath, model=self.HMMModel)
        hmmResFH = open(scanOutputPath, 'r')
        hmmParse = HMMERParse(hmmerScanFileHandle=hmmResFH, maxBestModels=10)
        # check maxBestModels applicability
        model = hmmParse.goToNextModelAlignments()
        while model != None:
            domainHit = hmmParse.nextAlignmentResult()
            while domainHit != None:
                domain_loc = FeatureLocation(domainHit.getEnvStart(),domainHit.getEnvStop())
                qual = {"evalue" : domainHit.getCEvalue(),
                        "score" : domainHit.getScore(),
                        "name" : model,
                        "subtype": model}
                domain_feat = SeqFeature(domain_loc, type="domain", strand=1, id=model, qualifiers=qual)
                seqAA.features.append(domain_feat)
                domainHit = hmmParse.nextAlignmentResult()
            model = hmmParse.goToNextModelAlignments()





class CladeFeatureMarker(FeatureMarker):

    def __init__(self, outPutPath, HMMModel, HMMBinaryPath, skipClade, useCutOff, allDomains):
        self.outPutPath = outPutPath
        self.HMMModel = HMMModel
        if not isfile(self.HMMModel):
            print "Could not find HMMModel at "+self.HMMModel
            exit(-1)
        self.skipClade = skipClade
        self.useCutOff = useCutOff
        self.hmmerExec = HMMERSuite(path=HMMBinaryPath)
        self.allDomains = allDomains
        if self.useCutOff:
            self.__useModelCutOff()

    def __useModelCutOff(self):
        self.hmmerCutOffStore = HMMERModelBasedCutoffStore(self.HMMModel+".modelcutoffs")

    def markFeaturesInSequenceRecord(self, seqAA):

        indfastaPath=self.outPutPath+"Query_"+seqAA.id+".faa"
        SeqIO.write(sequences=seqAA, handle=indfastaPath, format="fasta")
        scanOutputPath=self.outPutPath+seqAA.id+"_hmmerRes.txt"
        self.hmmerExec.runScan(query=indfastaPath, output=scanOutputPath, model=self.HMMModel)
        hmmResFH = open(scanOutputPath,'r')
        hmmParse = HMMERParse(hmmerScanFileHandle=hmmResFH, maxBestModels=10)

        model = hmmParse.goToNextModelAlignments()
        while model is not None:
            domainHit = hmmParse.nextAlignmentResult()
            while domainHit is not None:
                domain_loc = FeatureLocation(domainHit.getEnvStart(), domainHit.getEnvStop())

                qual = {"evalue": domainHit.getCEvalue(),
                        "score": domainHit.getScore(),
                        "name": model,
                        "subtype": "KS" if model != "AllTransKS" else "General KS"}

                domain_feat = SeqFeature(domain_loc, type="domain",
                                         strand=1, id=model, qualifiers=qual)
                seqAA.features.append(domain_feat)
                domainHit = hmmParse.nextAlignmentResult()
            model = hmmParse.goToNextModelAlignments()

class SeqFeatureContainer(object):
    
    def __init__(self,generalSeqFeature):
        self.generalSeqFeature = generalSeqFeature
        self.listOfSubFeatures = list()
        
    def generalLocationContainsCladeLoc(self,cladeLocation):
        generalLoc = self.generalSeqFeature.location
        if generalLoc.start <= cladeLocation.start and cladeLocation.end <= generalLoc.end:
            return True
        else:
            return False
        
    def addCladeSeqFeature(self,cladeSeqFeature):
        self.listOfSubFeatures.append(cladeSeqFeature)

    def getAllCladesSorted(self):
        if len(self.listOfSubFeatures) == 0:
            return None
        sortedList = sorted(self.listOfSubFeatures, reverse=True,
                            key=lambda feature: feature.qualifiers["score"])
        for i in range(0, len(sortedList), 1):
            sortedList[i].qualifiers["ranking"] = i+1

        return sortedList

        
    def getBestCladeSeqFeature(self):
        sortedList = self.getAllCladesSorted()
        if len(sortedList) == 0:
            return None
        bestHit = sortedList[0]
        quals = bestHit.qualifiers
        if len(sortedList) > 1:
            nextEvalue = sortedList[1].qualifiers["evalue"]
            nextScore = sortedList[1].qualifiers["score"]
            nextClade = sortedList[1].qualifiers["name"]
            quals["note"] = "Next best result was "+nextClade+" Evalue: %0.E Score: %d " % (nextEvalue,nextScore)
        bestHit.qualifiers = quals
        return bestHit
        
            
            
            
    
    
    



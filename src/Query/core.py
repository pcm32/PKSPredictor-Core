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
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "GAGTGx(75,85)[LV]HAT", "methyltransferase_p"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "PVIA", "crotonase_p_1"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "FTPG", "crotonase_p_2"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "GG[ALVMTGS]G[GRYDAVTHSK][LIV]G", "KR motif stereospecificty S"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "[GDPSERVIHKNAT]x(3)[GVTA][IVAL][VHILF][HYFVQY][SIATGMLCFNV][AVTPS][GLIMRP]x(3)D","KR D-configured OH groups"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "[HN][GSA][TAVSIP][GAS][TMGS]","Elongating KS"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DLxx[MWL]x(38)[VTF]x(20)[LYMG]x[CDI]X(20)[TAS]x(7)Vx(185)K", "D-Ala"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "D[MLVA]xx[IFLG]x(38)[WFIG]x(20)[ICY]x[AG]x(20)[CGAM]x(7)[LPIV]x(185)K", "Ala"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "IDxxTx(38)Lx(20)IxSx(20)Sx(7)Gx(185)K", "beta-Ala"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DVxxAx(38)Dx(20)VxCx(20)Ax(7)Ix(185)K", "Arg"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DLxxTx(38)Kx(20)ILVxGx(20)Ex(7)Vx(185)K", "Asn"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DLxxTx(38)Kx(20)[IV]xGx(20)[AH]x(7)[VI]x(185)K", "Asp"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DLxx[YF]x(38)Nx(20)[ML]xSx(20)[MLP]x(7)Ix(185)K", "Cys"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DHxxEx(38)Sx(20)Dx[IV]x(20)Gx(7)Ix(185)K", "Cys"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DLxxEx(38)Wx(20)NxTx(20)Tx(7)Vx(185)K", "Dab"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DAxx[QW]x(38)[DQ]x(20)[LF]xGx(20)[VL]x(7)[VI]x(185)K", "Gln"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DAxxWx(38)Hx(20)FxGx(20)[SG]x(7)Vx(185)K", "Glu"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DAxxKx(38)Dx(20)[IL]xGx(20)Vx(7)Vx(185)K", "Glu"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DQxxGx(38)Kx(20)TxGx(20)Vx(7)Gx(185)K", "3Me-Glu"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "D[IL]xx[FLQ]x(38)[NQM]x(20)[NLF]x[AGV]x(20)[LA]x(7)[TIM]x(185)K", "Gly"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DSxxEx(38)Ax(20)TxAx(20)Ex(7)Vx(185)K", "His"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DGxxFx(38)Fx(20)[LF]xGx(20)[V]x(7)[V]x(185)K", "Ile"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DAxxFx(38)Fx(20)[Y]xGx(20)[I]x(7)Tx(185)K", "Ile"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DAxx[W]x(38)[F]x(20)[L]xGx(20)[N]x(7)Vx(185)K", "Leu"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "D[AG]xx[F]x(38)[FM]x(20)[LM]xGx(20)[MVC]x(7)[V]x(185)K", "Leu"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DGxx[A]x(38)Yx(20)[T]xGx(20)Ex(7)Vx(185)K", "Leu"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DAxxWx(38)Lx(20)YxGx(20)Ax(7)Vx(185)K", "Leu"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DAxxEx(38)Sx(20)IxGx(20)Sx(7)Vx(185)K", "Lys"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DVxxAx(38)Hx(20)PxGx(20)Fx(7)Ix(185)K", "6hLys"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DAxxQx(38)Dx(20)AxGx(20)Cx(7)Vx(185)K", "6haLys"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DVxxGx(38)Ex(20)Lx[MLV]x(20)Sx(7)[VIL]", "Orn"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "D[IM]xx[NE]x(38)[YN]x(20)[WL]xGx(20)[GL]x(7)Ix(185)K", "5hOrn"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DAxx[WPFLV]x(38)[ITVG]x(20)[MVI]x[GA]x(20)[GA]x(7)[VIT]x(185)K", "Phe"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DYxxQx(38)Yx(20)[CL]xGx(20)Hx(7)Lx(185)K", "Pip"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DVxxQx(38)[YFSV]x(20)[AI]xAx(20)Hx(7)Vx(185)K", "Pro"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DIxx[TA]x(38)[LV]x(20)[IV]x[AT]x(20)[GV]x(7)Lx(185)K", "Pro"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DVxxWx(38)Hx(20)FxSx(20)LxVx(185)K", "Ser"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "D[LV]xxWx(38)Hx(20)LxSx(20)LxIx(185)K", "Ser"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "D[IAV]xx[FVY]x(38)Hx(20)LxGx(20)Lx(7)L", "HPG"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DPxxYx(38)Lx(20)GxGx(20)Tx(7)Lx(185)K", "DHPG"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DIxxFx(38)Lx(20)WxGx(20)Lx(7)Lx(185)K", "PGly"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DFxxWx(38)Nx(20)[IVL]xGx(20)Mx(7)Vx(185)K", "Thr"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "D[GA]xxWx(38)[SA]x(20)Vx[GA]x(20)Sx(7)Vx(185)K", "Trp"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "D[TA]xxSx(38)[KT]x(20)[VL]x[GA]x(20)Ax(7)Ix(185)K", "3hTyr"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "D[IP]xx[LW]x(38)[QG]x(20)LxGx(20)Lx(7)[VI]x(185)K", "3h4mPhe"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DGxxTx(38)Ix(20)TxAx(20)Ex(7)Vx(185)K", "Tyr"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DAxxLx(38)Vx(20)TxGx(20)Ax(7)Vx(185)K", "Tyr"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DAxxSx(38)Tx(20)VxAx(20)Ax(7)Vx(185)K", "Tyr"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DAxxFx(38)Wx(20)IxGx(20)Gx(7)Tx(185)K", "Val"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DFxxEx(38)Sx(20)TxAx(20)Ax(7)Vx(185)K", "Val"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "DAxxWx(38)Mx(20)FxAx(20)Ax(7)Vx(185)K", "Val"))
            sequenceMarkers.append(PatternFeatureMarker(self.fuzzProPath, "EPxxRx(38)Nx(20)IxVxEx(20)Fx(7)Vx(185)K", "Aad"))



        self.listOfSeqRecords = list()
            
        for seqAA in listOfAASeqs:
            # For each amino acid sequence, we run HMMER and mark
            print "Processing "+seqAA.id

            for seqMarker in sequenceMarkers:
                seqMarker.markFeaturesInSequenceRecord(seqAA)
            # indfastaPath=self.outPutPath+"Query_"+seqAA.id+".faa"
            # SeqIO.write(sequences=seqAA, handle=indfastaPath, format="fasta")
            # scanOutputPath=self.outPutPath+seqAA.id+"_hmmerRes.txt"
            # hmmerExec.runScan(query=indfastaPath, output=scanOutputPath, model=self.HMMModel)
            # hmmResFH = open(scanOutputPath,'r')
            # hmmParse = HMMERParse(hmmerScanFileHandle=hmmResFH, maxBestModels=10)
            # # First model hit needs to be of the general PKS, get those coordinates
            # model = hmmParse.goToNextModelAlignments()
            # if not (model == "AllTransKS"):
            #     print "First hit is not the PKS General signal, skipping"
            #     continue
            #
            # listOfGeneralModelDomains = list()
            # # A list of SeqFeatureContainer
            # domainHit = hmmParse.nextAlignmentResult()
            # print "Processing domains for AllTransKS"
            # while domainHit != None:
            #     domain_loc = FeatureLocation(domainHit.getQueryStart(),domainHit.getQueryStop())
            #     qual = {"evalue" : domainHit.getCEvalue(),
            #                 "score" : domainHit.getScore(),
            #                 "name" : model}
            #     domain_feat = SeqFeature(domain_loc,type="domain",strand=1,id=model)
            #     # we could put scores and evalues in qualifiers field
            #     listOfGeneralModelDomains.append(SeqFeatureContainer(domain_feat))
            #     # seqAA.features.append(domain_feat)
            #     domainHit = hmmParse.nextAlignmentResult()
            #
            # print "Processing domains for Clades"
            # # For each subsequent model (the clades), if they fit within detected PKS general
            # # coordinates, store the clade specific hit as a SeqFeature of the SeqRecord
            # model = hmmParse.goToNextModelAlignments()
            # while model is not None:
            #     domainHit = hmmParse.nextAlignmentResult()
            #     if self.skipClade.count(model)>0:
            #         print "Skipping model "+model
            #         model = hmmParse.goToNextModelAlignments()
            #         continue
            #     print "Processing model "+model
            #     while domainHit != None:
            #         domain_loc = FeatureLocation(domainHit.getQueryStart(),domainHit.getQueryStop())
            #         found=False
            #         pickedSeqFeatCont = None
            #         for seqFeatCont in listOfGeneralModelDomains:
            #             if seqFeatCont.generalLocationContainsCladeLoc(domain_loc):
            #                 found=True
            #                 pickedSeqFeatCont = seqFeatCont
            #                 break
            #         if not found:
            #             domainHit = hmmParse.nextAlignmentResult()
            #             break
            #         if self.useCutOff and self.hmmerCutOffStore.hasEvalueCutoff(model):
            #             if self.hmmerCutOffStore.getCutEvalueCutOff(model) < float(domainHit.getCEvalue()):
            #                 print "Skipping model "+model+" due to evalue cutoff."
            #                 domainHit = hmmParse.nextAlignmentResult()
            #                 break
            #         qual = {"evalue" : domainHit.getCEvalue(),
            #                 "score" : domainHit.getScore(),
            #                 "name" : model}
            #         domain_feat = SeqFeature(domain_loc, type="domain",
            #                                  strand=1, id=model, qualifiers=qual)
            #         # seqAA.features.append(domain_feat)
            #         pickedSeqFeatCont.addCladeSeqFeature(domain_feat)
            #         domainHit = hmmParse.nextAlignmentResult()
            #
            #     model = hmmParse.goToNextModelAlignments()
            #
            # # Finally, add the seqRecord to the list of SeqRecords to be returned.
            # for seqFeatCont in listOfGeneralModelDomains:
            #     bestClade = seqFeatCont.getBestCladeSeqFeature()
            #     if bestClade != None:
            #         seqAA.features.append(bestClade)

            # Add additional motifs/domains marks, like the methytransferase trough fuzzpro

            # This can be written at the outside to GBK
            self.listOfSeqRecords.append(seqAA)
    
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
        self.outPutPath = tempfile.gettempdir()
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
                domain_loc = FeatureLocation(domainHit.getQueryStart(),domainHit.getQueryStop())
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
        # First model hit needs to be of the general PKS, get those coordinates
        model = hmmParse.goToNextModelAlignments()
        if not (model == "AllTransKS"):
            print "First hit is not the PKS General signal, skipping"
            return

        listOfGeneralModelDomains = list()
        # A list of SeqFeatureContainer
        domainHit = hmmParse.nextAlignmentResult()
        print "Processing domains for AllTransKS"
        while domainHit != None:
            domain_loc = FeatureLocation(domainHit.getQueryStart(),domainHit.getQueryStop())
            qual = {"evalue" : domainHit.getCEvalue(),
                        "score" : domainHit.getScore(),
                        "name" : model}
            domain_feat = SeqFeature(domain_loc, type="domain", strand=1, id=model)
            # we could put scores and evalues in qualifiers field
            listOfGeneralModelDomains.append(SeqFeatureContainer(domain_feat))
            # seqAA.features.append(domain_feat)
            domainHit = hmmParse.nextAlignmentResult()

        print "Processing domains for Clades"
        # For each subsequent model (the clades), if they fit within detected PKS general
        # coordinates, store the clade specific hit as a SeqFeature of the SeqRecord
        model = hmmParse.goToNextModelAlignments()
        while model is not None:
            domainHit = hmmParse.nextAlignmentResult()
            if self.skipClade.count(model) > 0:
                print "Skipping model "+model
                model = hmmParse.goToNextModelAlignments()
                continue
            print "Processing model "+model
            while domainHit is not None:
                domain_loc = FeatureLocation(domainHit.getQueryStart(),domainHit.getQueryStop())
                found = False
                pickedSeqFeatCont = None # The general clade picked, where things are added later.
                for seqFeatCont in listOfGeneralModelDomains:
                    if seqFeatCont.generalLocationContainsCladeLoc(domain_loc):
                        found = True
                        pickedSeqFeatCont = seqFeatCont
                        break
                if not found:
                    domainHit = hmmParse.nextAlignmentResult()
                    break
                if self.useCutOff and self.hmmerCutOffStore.hasEvalueCutoff(model):
                    if self.hmmerCutOffStore.getCutEvalueCutOff(model) < float(domainHit.getCEvalue()):
                        print "Skipping model "+model+" due to evalue cutoff."
                        domainHit = hmmParse.nextAlignmentResult()
                        break
                qual = {"evalue": domainHit.getCEvalue(),
                        "score": domainHit.getScore(),
                        "name": model,
                        "subtype": "KS"}
                domain_feat = SeqFeature(domain_loc, type="domain",
                                         strand=1, id=model, qualifiers=qual)
                # seqAA.features.append(domain_feat)
                pickedSeqFeatCont.addCladeSeqFeature(domain_feat)
                domainHit = hmmParse.nextAlignmentResult()

            model = hmmParse.goToNextModelAlignments()

        # The stack number refers to the lump of features marked under the same KS region (same hits for a single region).
        featureStackNumber = 1
        # Finally, add the seqRecord to the list of SeqRecords to be returned.
        for seqFeatCont in listOfGeneralModelDomains:
            if self.allDomains:
                if seqFeatCont is not None:
                    allClades = seqFeatCont.getAllCladesSorted()
                    if allClades is not None:
                        for cladeFeature in allClades:
                            cladeFeature.qualifiers["region"] = featureStackNumber
                            seqAA.features.append(cladeFeature)
            else: # if we don't want all the clades, we only add the best one.
                bestClade = seqFeatCont.getBestCladeSeqFeature()
                if bestClade is not None:
                    seqAA.features.append(bestClade)
            featureStackNumber += 1

    
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
        
            
            
            
    
    
    



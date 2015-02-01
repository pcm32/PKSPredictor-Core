'''
Created on Sep 24, 2011

@author: pmoreno
'''
import sys
from symbol import except_clause

if __name__ == '__main__':
    from Query.core import QueryRunner


    pathToFastaQuery=sys.argv[1]
    pathToHMMModel=sys.argv[2]
    pathToOutput=sys.argv[3]
    pathToHMMBinary=sys.argv[4]
    pathToFuzzPro=sys.argv[5]
    pathToNRPS2=sys.argv[6]
    indexForSkips = 7
    allDomains = False
    if len(sys.argv)>7 and sys.argv[7] == "allDomains":
        ++indexForSkips
        allDomains = True

    modelsToSkip=list()
    if len(sys.argv)>indexForSkips:
        if sys.argv[indexForSkips].count(",")>0:
            modelsToSkip=sys.argv[indexForSkips].rsplit(",")

    queryRunner = QueryRunner(fastaQuery=pathToFastaQuery,
                              HMMModel=pathToHMMModel,
                              outPutPath=pathToOutput,
                              HMMBinaryPath=pathToHMMBinary,
                              fuzzProPath=pathToFuzzPro,
                              NRPS2Path=pathToNRPS2,
                              allDomains=allDomains
                              )
    
    if sys.argv.count("useModelCutOff") == 1:
        queryRunner.useModelCutOff()
    
    if len(modelsToSkip)>0:
        queryRunner.setCladesToSkip(modelsToSkip)
    queryRunner.runQuery()
    
    listOfSeqRecords = queryRunner.getSeqsWithResults()
    
    from Bio import SeqIO
    from SimpleFeatWriter.core import SimpleFeatWriter
    from fastaModifier.core import GBKRankingBasedFeatureRemover, SeqRecordAnnotator
    from hmmer.core import ModelAnnotator

    gbkProcessor = GBKRankingBasedFeatureRemover(1)
    modelAnnotator = ModelAnnotator(pathToHMMModel)
    seqRecordAnnotator = SeqRecordAnnotator(modelAnnotator)
    
    for seq in listOfSeqRecords:
        gbkToWrite = pathToOutput+seq.id+".gbk"
        seqRecordAnnotator.annotateSeqFeatures(seq)
        featWriter = SimpleFeatWriter(pathToOutput, 5) # number is max ranking to print
        featWriter.writeFeatures(seq)
        if "|" in seq.id:
            seq.id = seq.id.split("|")[2]
        if len(seq.id) > 16:
            seq.id = seq.id[:16]
        seq.name = seq.id

        try:
            SeqIO.write(sequences=gbkProcessor.processGBK(seq), handle=gbkToWrite, format="genbank")
        except AttributeError as error:
            print "Error writing "+seq.id+" to genbank: "+str(error)

        # touch finished file.
        open(pathToOutput+seq.id+".finished", "w").close()
    

'''
Created on Sep 24, 2011

@author: pmoreno
'''
import sys
import getopt
from __future__ import print_function
from symbol import except_clause


def usage():
    print("runQueryAgainstHMMModels.py options:", file=sys.stderr)
    print("-q path to query fasta file", file=sys.stderr)
    print("-c path to HMMER model for clades", file=sys.stderr)
    print("-o output path", file=sys.stderr)
    print("-m path to HMMER models for other domains", file=sys.stderr)
    print("-a flag, if present, all clades are annotated for each "
          + "region, otherwise only the best clade is annotated for each region.", file=sys.stderr)
    print("-s optional, names of the clades to skip, comma separated if more than one.", file=sys.stderr)
    print("--HMMERPath the path to the HMMER executables.", file=sys.stderr)
    print("--FuzzProPath the path to the EMBOSS fuzzpro executable.", file=sys.stderr)
    print("--NRPS2Path the path to the NRPS2 executable.", file=sys.stderr)


def main():
    from Query.core import QueryRunner


    try:
        opts, args = getopt.getopt(sys.argv[1:], "q:c:o:as:m:", ["HMMERPath=", "FuzzProPath=", "NRPS2Path="])
    except getopt.GetoptError as err:
        print(str(err)) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    allDomains = False
    modelsToSkip=list()
    for opt, arg in opts:
        if opt == "-q":
            pathToFastaQuery = arg  # 1
        elif opt == "-c":
            pathToHMMModel = arg  # 2
        elif opt == "-o":
            pathToOutput = arg # 3
        elif opt == "--HMMERPath":
            pathToHMMBinary = arg  # 4
        elif opt == "--FuzzProPath":
            pathToFuzzPro = arg  # 5
        elif opt == "--NRPS2Path":
            pathToNRPS2 = arg  # 6
        elif opt == "-a":
            allDomains = True
        elif opt == "-s":
            modelsToSkip = arg.split(",")
        elif opt == "-m":
            pathToHMMModelOthers = arg


    queryRunner = QueryRunner(fastaQuery=pathToFastaQuery,
                              HMMModel=pathToHMMModel,
                              HMMModelNonClades=pathToHMMModelOthers,
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
            print("Error writing "+seq.id+" to genbank: "+str(error))

        # touch finished file.
        open(pathToOutput+seq.id+".finished", "w").close()

if __name__ == '__main__':
    main()

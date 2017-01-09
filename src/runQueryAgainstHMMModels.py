'''
Created on Sep 24, 2011

@author: pmoreno
'''
from __future__ import print_function

import os
import sys
import getopt
from symbol import except_clause

from Clades.core import CladificationAnnotationReader
from SeqFeatureSorter.core import SeqFeatureSorter


# TODO this is horrible, should be moved to argparse instead!!!
from verifier.core import KSDomainVerifier, DH_PS_DomainVerifier


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

    found_q = False
    found_c = False
    found_o = False

    pathToHMMModelOthers = None
    pathToHMMBinary = ""
    pathToFuzzPro = ""

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
            found_q = True
        elif opt == "-c":
            pathToHMMModel = arg  # 2
            found_c = True
        elif opt == "-o":
            pathToOutput = arg # 3
            found_o = True
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

    if not (found_c or found_q or found_o):
        usage()
        sys.exit(2)


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
    if os.path.isfile(CladificationAnnotationReader.get_expected_annot_path(pathToHMMModel)):
        # This is the case for the new annotation format, which is expected with a slightly different name
        # (different extension).
        reader = CladificationAnnotationReader()
        clade_annotation = reader.get_annotation(CladificationAnnotationReader.get_expected_annot_path(pathToHMMModel))
        modelAnnotator = ModelAnnotator(clade_annotation=clade_annotation)
    else:
        # For the legacy annotation format.
        modelAnnotator = ModelAnnotator(pathToHMMModel)
    seqRecordAnnotator = SeqRecordAnnotator(modelAnnotator)

    for seq in listOfSeqRecords:
        gbkToWrite = pathToOutput+seq.id+".gbk"
        # Sort features in the seq object geographically (by start).
        SeqFeatureSorter.geographic_sort(seq)
        seqRecordAnnotator.annotateSeqFeatures(seq)
        # Run verifiers and add verification field output to the SimpleFeatWriter
        if clade_annotation is not None:
            for domain_verifier in [KSDomainVerifier(cladififcation_annotation=clade_annotation, seqObj=seq),
                                    DH_PS_DomainVerifier(seq)]:
                domain_verifier.verify()

        featWriter = SimpleFeatWriter(pathToOutput, 5) # number is max ranking to print
        featWriter.writeFeatures(seq)
        seq.name = seq.id
        if "|" in seq.name:
            seq.name = seq.name.split("|")[2]
        if len(seq.name) > 16:
            seq.name = seq.name[:16]


        try:
            SeqIO.write(sequences=gbkProcessor.processGBK(seq), handle=gbkToWrite, format="genbank")
        except AttributeError as error:
            print("Error writing "+seq.id+" to genbank: "+str(error))

        # touch finished file.
        open(pathToOutput+seq.id+".finished", "w").close()

if __name__ == '__main__':
    main()

__author__ = 'pmoreno'

if __name__ == "__main__":
    """
    This script receives the clade definition file, the path to the complete fasta alignment
    and an output path to produce alignment fasta files for each clade.
    """
    import sys
    from Clades.core import CladeDefinitionReader

    if len(sys.argv) < 2:
        print "USAGE: getSeparateAlignmentsPerCladeFromDefinitionFile.py fasta_id-clade_id-assignment-file complete_fasta_alignment outputdir hmmname"
        exit(1)


    cladeDefReader = CladeDefinitionReader(sys.argv[1])
    pathForCompleteFastaAlign = sys.argv[2]
    outputDir = sys.argv[3]
    hmmName = sys.argv[4]

    print "Using assignment file in:" + sys.argv[0]

    # The reader contains now mappings and others.
    from Clades.core import FastaCladeSplitter
    splitter = FastaCladeSplitter(cladeDefReader.entry2Clade)
    splitter.writeFastas(pathForCompleteFastaAlign, outputDir)

    # We also write the description/annotation for clades
    # TODO Deprecate what follows somehow.
    fh = open(outputDir+"/"+hmmName+".annot", mode="w")

    for clade in cladeDefReader.cladeDesc.keys():
        fh.write(clade+"\t"+cladeDefReader.cladeDesc.get(clade)+"\n")

    fh.close()

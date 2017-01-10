__author__ = 'pmoreno'

if __name__ == "__main__":
    """
    This script receives the clade definition file, the path to the complete fasta alignment
    and an output path to produce alignment fasta files for each clade.
    """
    import sys
    from Clades.core import CladeDefinitionReader


    cladeDefReader = CladeDefinitionReader(sys.argv[0])
    pathForCompleteFastaAlign = sys.argv[1]
    outputDir = sys.argv[2]
    hmmName = sys.argv[3]

    # The reader contains now mappings and others.
    from Clades.core import FastaCladeSplitter
    splitter = FastaCladeSplitter(cladeDefReader.entry2Clade)
    splitter.writeFastas(pathForCompleteFastaAlign,outputDir)

    # We also write the description/annotation for clades
    # TODO Deprecate what follows somehow.
    fh = open(outputDir+"/"+hmmName+".annot", mode="w")

    for clade in cladeDefReader.cladeDesc.keys():
        fh.write(clade+"\t"+cladeDefReader.cladeDesc.get(clade)+"\n")

    fh.close()

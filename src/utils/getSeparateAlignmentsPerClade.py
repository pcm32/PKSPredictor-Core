__author__ = 'pmoreno'


if __name__ == "__main__":
    """
    This script receives the Newick tree file, the path to the complete fasta alignment
    and an output path to produce alignment fasta files for each clade.
    """
    import sys
    pathForNewickTree = sys.argv[0]
    pathForCompleteFastaAlign = sys.argv[1]
    outputDir = sys.argv[2]
    hmmName = sys.argv[3]

    # We first read the newick tree file
    from Clades.core import NewickCladeReader
    newickReader = NewickCladeReader(pathForNewickTree)

    # The reader contains now mappings and others.
    from Clades.core import FastaCladeSplitter
    splitter = FastaCladeSplitter(newickReader.entry2Clade)
    splitter.writeFastas(pathForCompleteFastaAlign,outputDir)

    # We also write the description/annotation for clades
    fh = open(outputDir+"/"+hmmName+".annot",mode="w")

    for clade in newickReader.cladeDesc.keys():
        fh.write(clade+"\t"+newickReader.cladeDesc.get(clade)+"\n")

    fh.close()
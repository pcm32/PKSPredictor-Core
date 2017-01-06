__author__ = 'pmoreno'

class CladeReader(object):
    """
    Base class for Clade definition reading
    """

    def __init__(self):
        self.cladeDesc = {}
        self.entry2Clade = {}

    def __fillIndexes(self):
        pass



class NewickCladeReader(CladeReader):
    '''
    This class reads a set of newick files to produce an index of
    clades, stored in the cladeDesc[clade] = description and entr2clade[entry]=clade
    '''

    def __init__(self, newickFilesPath):
        self.pathToNewickFiles = newickFilesPath
        self.__fillIndexes()

    def __fillIndexes(self):
        import glob
        from Bio import Phylo

        for newickFileName in glob.glob1(self.pathToNewickFiles, "*.nwk"):
            clade_desc = newickFileName
            cladeTokens = clade_desc.split("_")
            clade = "_".join((cladeTokens[0], cladeTokens[1]))
            desc = newickFileName.replace(clade + "_", "").replace("\.nwk", "")

            # print "# "+clade+"\t"+desc
            self.cladeDesc[clade] = desc

            tree = Phylo.read(self.pathToNewickFiles + "/" + newickFileName, "newick")
            for subclade in tree.find_clades():
                if subclade.is_terminal() and subclade.name:
                    #print clade+"\t"+subclade.name
                    self.entry2Clade[subclade.name] = clade

class CladeDefinitionReader(CladeReader):
    """
    This class reads a Clade definition file, which is an alternative route to Newick files to define clades
    within a larger fasta alignment. The Clade definition file follows the following format:

    # Clade_1	amino acids
    >PKSXKS1	Clade_1
    >BaeKS1	Clade_1
    >PksXKS11	Clade_1
    >BaeKS10	Clade_1
    # Clade_2	amino acids
    >OzmKS9	Clade_2
    >OzmKS12	Clade_2
    # Clade_3	non-elongating KS
    >TaiKS4	Clade_3
    >RhiKS1	Clade_3
    >SGKS3	Clade_3

    Commented lines are used for the annotation file.
    """

    def __init__(self, pathToDefFile):
        CladeReader.__init__(self)
        self.pathToDef = pathToDefFile
        self.__fillIndexes()

    def __fillIndexes(self):
        defFile = open(self.pathToDef)
        for line in defFile:
            tokens = line.split("\t")
            if line.startswith("#"):
                clade = tokens[0].replace("# ", "", count=1)
                if len(tokens) > 1:
                    self.cladeDesc[clade] = tokens[1]
            elif len(tokens) > 1:
                clade = tokens[1]
                entry = tokens[0].replace(">","")
                self.entry2Clade[entry] = clade


class FastaCladeSplitter(object):
    '''
    This class reads a fasta file and mappings of its entries to clades, to produce
    a fasta file for each clade.
    '''

    def __init__(self,mappings):
        self._mappings = mappings
        self._clade2fileHandle = {}

    def writeFastas(self,pathToCompleteFastaAlignment,outputPath):
        """
        Writes individual clades as fasta alignments, where the files are stored in the output path
        and the name of each file is composed of <clade_name>.fasta . Data is obtained from the complete
        fasta alignment provided (pathToCompleteFastaAlignment) and the mappings given to the constructor.
        The mappings are a dictionary of fasta ids to clade names. Elements in the alignment not present
        in the mappings won't end in any clade fasta file.

        :param pathToCompleteFastaAlignment:
        :param outputPath:
        """
        from Bio import SeqIO, AlignIO
        from Bio.Alphabet import IUPAC, Gapped
        # create file handles for each clade
        for cladeName in list(self._mappings.values()):
            self._clade2fileHandle[cladeName] = open(outputPath+"/"+cladeName+".fas")

        # read complete alignment
        alignment = AlignIO.read(pathToCompleteFastaAlignment, "fasta", alphabet=Gapped(IUPAC.ExtendedIUPACProtein(), "-"))
        for record in alignment:
            handle = self._clade2fileHandle[self._mappings[record.id]]
            SeqIO.write(record, handle, "fasta")

        for cladeName in self._clade2fileHandle.keys():
            self._clade2fileHandle[cladeName].close()
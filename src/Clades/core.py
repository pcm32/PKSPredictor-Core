from collections import defaultdict
from itertools import chain

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
        CladeReader.__init__(self)
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
        self.entry2Clade = defaultdict(list)
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
                clade = tokens[1].rstrip()
                entry = tokens[0].replace(">", "")
                if clade not in self.entry2Clade[entry]:
                    self.entry2Clade[entry].append(clade)


class FastaCladeSplitter(object):
    '''
    This class reads a fasta file and mappings of its entries to clades, to produce
    a fasta file for each clade.
    '''

    def __init__(self, mappings):
        self._mappings = mappings
        self._clade2fileHandle = {}

    def writeFastas(self, pathToCompleteFastaAlignment, outputPath):
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

        for cladeName in set(chain.from_iterable(self._mappings.values())):
            self._clade2fileHandle[cladeName] = open(outputPath+"/"+cladeName+".fas", mode="w")

        # read complete alignment
        alignment = AlignIO.read(pathToCompleteFastaAlignment, "fasta", alphabet=Gapped(IUPAC.ExtendedIUPACProtein(), "-"))
        for record in alignment:
            if record.id in self._mappings:
                for clade in self._mappings[record.id]:
                    handle = self._clade2fileHandle[clade]
                    SeqIO.write(record, handle, "fasta")
            else:
                print "Skipping "+record.id+" as it doesn't have a clade assignment"

        for cladeName in self._clade2fileHandle.keys():
            self._clade2fileHandle[cladeName].close()


class KSDomainRequirement(object):
    # TODO probably get rid of this object... is nothing but a list of domain identifiers.
    def __init__(self, domains):
        self.__requirements_list = list()
        self.__requirements_list.extend(filter(lambda a: a is not '', domains.split(";")))

    def get_requirements(self):
        return self.__requirements_list


class CladificationAnnotation(object):
    def __init__(self):
        self.__molfile = {}
        self.__desc = {}
        self.__desc_tool = {}
        self.__post_processors = {}
        self.__verification_domains = {}
        self.__termination_rule = {}
        self.__non_elongating = {}
        self.__verification_mandatory = {}

    def add_entry(self, clade_id, desc, desc_tool, molfile, termination_rule,
                  verification_domains=[], non_elongating=False, verification_mandatory=False, post_processors=[]):
        """
        Adds an entry of annotation for a clade id.
        """
        self.__desc[clade_id] = desc
        self.__desc_tool[clade_id] = desc_tool
        self.__molfile[clade_id] = molfile
        self.__post_processors[clade_id] = post_processors
        self.__verification_domains[clade_id] = verification_domains
        self.__termination_rule[clade_id] = termination_rule
        self.__non_elongating[clade_id] = non_elongating
        self.__verification_mandatory[clade_id] = verification_mandatory

    def get_domain_requirements(self, clade_id):
        """
        :returns list of domains required upstream for a clade to be present
        """
        return self.__verification_domains[clade_id].get_requirements()

    def get_description(self, clade_id):
        return self.__desc[clade_id]

    def get_all_clade_ids(self):
        return self.__desc.keys()


class CladificationAnnotationReader(object):
    """
    Reads the cladification annotation file
    """

    @staticmethod
    def get_expected_annot_path(hmmer_model_path):
        """
        The expected annotation file path for a given hmmer model.
        """
        return hmmer_model_path+".annot_v2"

    def get_annotation(self, annotation_path):
        """
        Loads the annotation file whose path the object received on init. Data is left in a dictionary
        of KS_requirement classes.
        :return: CladificationAnnotation
        """
        fh = open(annotation_path)
        annotation = CladificationAnnotation()
        line = fh.readline()  # header
        line = fh.readline()
        while line is not None and len(line) > 1:
            (clade_id, desc, desc_tool, molfile, post_procs, verification_domains,
             termination_rule, non_elongating, verification_mandatory) = line.split("\t")

            annotation.add_entry(clade_id, desc, desc_tool, molfile,
                                 post_processors=filter(None, post_procs.split(";")),
                                 verification_domains=KSDomainRequirement(verification_domains),
                                 termination_rule=termination_rule, non_elongating=(non_elongating == 'yes'),
                                 verification_mandatory=(verification_mandatory == 'yes')
                                 )
            line = fh.readline()
        fh.close()
        return annotation
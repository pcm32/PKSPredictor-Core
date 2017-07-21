from tempfile import mkstemp
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os.path

__author__ = 'pmoreno'


class NRPS2Parser(object):
    def __init__(self, pathToReport):
        self.file = open(pathToReport, "r")
        # get rid of header
        self.file.readline()

    def nextRegion(self):
        line = self.file.readline()
        if line is None or len(line) == 0:
            return None
        id, sig, stachelhausCode, thirdClassPred, largeClassPred, smallClassPred, singleClassPred, \
        nearesStachehausCode, nrps1PredLargeClass, nrps2PredLargeClass, outsideAppDomain, coords, score \
            = tuple(line.strip().split("\t"))

        start, stop = coords.split(":")
        domain_loc = FeatureLocation(int(start), int(stop))
        qual = {"score": score,
                "name": singleClassPred if singleClassPred != "N/A" else nearesStachehausCode,
                "subtype": "NRPS2",
                "tool": "NRPS2"}
        domain_feat = SeqFeature(domain_loc, type="domain",
                                 strand=1, id=id + "_" + singleClassPred, qualifiers=qual)
        return domain_feat


class NRPS2PredRunner(object):
    def __init__(self, path):
        self.path = path

    def run(self, seqRecord):
        """
        Runs NRPSPredictor2 for the seqRecord given, then parses the result and produces a
        and writes it back to the seqRecord as a sequence feature.
        Sequence is BioPython SeqRecord object.

        :param seqRecord:
        :return: seqRecord: same object with an annotated sequence feature for the result of
        NRPSPredictor2 run.
        """

        # First we write the seqRecord to a temporary fasta file
        fastaFile, fastaFilePath = mkstemp()
        from Bio import SeqIO
        # print(seqRecord.format("fasta"))
        SeqIO.write(seqRecord, fastaFilePath, "fasta")

        # Now we set up the call
        command = "java -Ddatadir=" + self.path + "/data -cp "
        for jar in "dom4j-1.6.1.jar", "jaxen-1.1-beta-6.jar", "java-getopt-1.0.13.jar", \
                   "Utilities.jar", "libsvm.jar":
            command = command + self.path + "/lib/" + jar + ":"

        command = command + self.path + "/build/NRPSpredictor2.jar "
        command = command + "org.roettig.NRPSpredictor2.NRPSpredictor2 "
        command = command + "-i " + fastaFilePath + " -s 0 -r " + fastaFilePath + ".nrps2Report"

        # And we execute it
        from subprocess import call
        execCode = call(command, shell=True)

        # Finally we parse the result and add sequence features to the seq record.
        resultFilePath = fastaFilePath + ".nrps2Report"
        if os.path.isfile(resultFilePath):
            parser = NRPS2Parser(resultFilePath)

            while True:
                seqFeat = parser.nextRegion()
                if seqFeat is None:
                    break
                seqRecord.features.append(seqFeat)
        else:
            print "NRPS2predictor failed probably, result file not generated."

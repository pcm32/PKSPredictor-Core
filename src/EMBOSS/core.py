import subprocess
import re
from Bio.SeqFeature import SeqFeature, FeatureLocation

__author__ = 'pmoreno'


class FuzzProParser(object):
    def __init__(self, file, patternName):
        self.file = file
        prog = re.compile('\s+Start\s+End\s+Score')
        # get to the header line
        self.regionPat = re.compile('^\s+(\d+)\s+(\d+)\s+(\d+)')
        self.patternName = patternName
        while True:
            line = self.file.readline()
            if len(line) == 0:
                break
            if prog.match(line):
                break


    def nextRegion(self):
        """ Returns the next appearance of the pattern, if any, as a SeqFeature
            None is returned if the pattern doesn't appear again """
        line = self.file.readline()
        if re.match('^#',line):
            return None
        res = self.regionPat.match(line)
        if res:
            domain_loc = FeatureLocation(int(res.group(1)), int(res.group(2)))
            qual = { "score" : res.group(3),
                 "name" : self.patternName}
            domain_feat = SeqFeature(domain_loc, type="pattern",
                        strand=1, id=self.patternName, qualifiers=qual)
            return domain_feat
        return None


class FuzzProRunner(object):

    def __init__(self, path2FuzzPro):
        self.path = path2FuzzPro
        self.executable = "fuzzpro"

    def setPattern(self, pattern, patternName):
        self.pattern = pattern
        self.patternName = patternName

    def run(self, seqRecord):
        """
        Runs fuzzpro on the seqRecord given the set pattern, parses the output, and writes it
        back to the seqRecord as a sequence feature. Sequence is BioPython SeqRecord object.
        """
        from subprocess import Popen, PIPE
        # print [self.path+"/"+self.executable, "-pattern", self.pattern,
        #                    "-sequence","asis::"+seqRecord.seq.tostring(),
        #                    "-auto", "-rformat", "table", "-stdout"]
        self.pScan = Popen([self.path+"/"+self.executable, "-pattern", self.pattern,
                            "-sequence","asis::"+seqRecord.seq.tostring(),
                            "-auto", "-rformat", "table", "-stdout"], stdout=PIPE)
        # get stdout somehow and then parse it
        self.pScan.wait()

        parser = FuzzProParser(self.pScan.stdout, self.patternName)

        while True:
            seqFeat = parser.nextRegion()
            if seqFeat is None:
                break
            seqRecord.features.append(seqFeat)



import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_protein

__author__ = 'pmoreno'

import sys
sys.path.append('../../')

from EMBOSS.core import FuzzProRunner


class TestFuzzProRunner(unittest.TestCase):
    def test_run(self):
        fuzzProRunner = FuzzProRunner("/Users/pmoreno/Downloads/EMBOSS-6.5.7/emboss")
        fuzzProRunner.setPattern("GGx(5)AA","testPattern")
        seq = Seq("MSTEGGGRRCQAAASRRISFSASHRGGSKFLSAAENLKLFGKCNNPNGHG",generic_protein)
        seqRec = SeqRecord(seq, id="test")
        fuzzProRunner.run(seqRec)
        self.assertTrue(len(seqRec.features)==3,"Should find 3 sequence features")



if __name__ == '__main__':
    unittest.main()

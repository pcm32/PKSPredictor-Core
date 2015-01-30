'''
Created on Sep 21, 2011

@author: pmoreno
'''
import unittest


class Test(unittest.TestCase):
    
    def test_hmmer(self):
        #from hmmer.core import HMMERScan
        from hmmer.core import HMMERSuite
        hmmerExec = HMMERSuite(path="/Applications/hmmer-3.0-macosx-intel/binaries/")
        hmmerExec.runScan(query='~/fileToScan.faa',model='~/pathToModel')
        
        
        
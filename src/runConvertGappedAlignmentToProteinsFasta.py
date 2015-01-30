'''
Created on Apr 17, 2012

@author: pmoreno
'''
import sys

if __name__ == '__main__':
    from fastaModifier.core import FastaGapRemover
    
    pathToFastaQuery=sys.argv[1]
    
    gapRemover = FastaGapRemover(pathToFastaQuery)
    gapRemover.writeNoGaps(pathToFastaQuery+"_noGaps.fas")
    
    
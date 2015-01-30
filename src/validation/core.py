import math
from hmmer.core import HMMERModelHit
class ValidationResult(object):
    
    def __init__(self,thisCladeName,listOfModelHits):
        self.thisCladeName = thisCladeName
        self.listOfModelHits = listOfModelHits
        
    def addModelHit(self,modelHit):
        self.listOfModelHits.append(modelHit)
        
    def isCladeMatch(self):
        if len(self.listOfModelHits) > 0:
            return (self.listOfModelHits[0].getModel() == self.thisCladeName)
        return False
    
    def getScoreDifferenceToSecondModel(self):
        if len(self.listOfModelHits) > 1:
            return self.listOfModelHits[0].getFullSeqScore() - self.listOfModelHits[1].getFullSeqScore()
        
    def getEvalueDifferenceToSecondModel(self):
        if len(self.listOfModelHits) > 1:
            return math.log10(self.listOfModelHits[0].getFullSeqEvalue()) - math.log10(self.listOfModelHits[1].getFullSeqEvalue())

    
    def getFirstModel(self):
        if len(self.listOfModelHits) > 0:
            return self.listOfModelHits[0]

    
    def getFirstModelName(self):
        if len(self.listOfModelHits) > 0:
            return self.listOfModelHits[0].getModel()

    
    def getFirstModelEvalue(self):
        if len(self.listOfModelHits) > 0:
            return self.listOfModelHits[0].getFullSeqEvalue()

    
    def getFirsModelScore(self):
        if len(self.listOfModelHits) > 0:
            return self.listOfModelHits[0].getFullSeqScore()
    
    
        
    
    
    
    
    
    
            
    



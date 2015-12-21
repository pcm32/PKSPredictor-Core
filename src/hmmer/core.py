import subprocess


class HMMERSuite(object):
    
    def __init__(self, path):
        self.path = path
        self.executableScan = 'hmmscan'
        self.executableBuild = 'hmmbuild'
        self.executablePress = 'hmmpress'

    def runScan(self, query, output, model):
        from subprocess import Popen
        fileOut=open(output,'w')
        self.pScan = Popen([self.path+"/"+self.executableScan,model,query],stdout=fileOut)
        self.pScan.wait()
        fileOut.close()
    
    def runBuild(self, fastaAlign, output, modelName):
        from subprocess import Popen
        self.pBuild = Popen([self.path+"/"+self.executableBuild,"-n",modelName,"--informat","afa",output,fastaAlign])
        self.pBuild.wait()
        
    def runPress(self,hmmModelCompendium):
        from subprocess import Popen
        self.pPress = Popen([self.path+"/"+self.executablePress,hmmModelCompendium])
        self.pPress.wait()
    #def getOutput(self):
     #   return self.pScan.stdout

    
class HMMERDomainHit(object):
    
    def __init__(self,valid,score,bias,cevalue,ievalue,hmmStart,hmmStop,hmmBound,qstart,qstop,qBound):
        self.valid = valid
        self.score = float(score)
        self.bias = float(bias)
        self.cevalue = float(cevalue)
        self.ievalue = float(ievalue)
        self.hmmStart = int(hmmStart)
        self.hmmStop = int(hmmStop)
        self.hmmBound = hmmBound
        self.qstart = int(qstart)
        self.qstop = int(qstop)
        self.qBound = qBound
        
    def getQueryStart(self):
        return self.qstart
    
    def getQueryStop(self):
        return self.qstop
    
    def getIsValid(self):
        if self.valid == "!":
            return True
        return False
            
    def getCEvalue(self):    
        return self.cevalue
    
    def getIEvalue(self):
        return self.ievalue
    
    def getScore(self):
        return self.score
    
    
class HMMERModelHit(object):
    
    def __init__(self,fsEvalue,fsScore,fsBias,bdEvalue,bdScore,bdBias,exp,numOfDomains,model,desc):
        self.fsEvalue = float(fsEvalue)
        self.fsScore = float(fsScore)
        self.fsBias = float(fsBias)
        self.bdEvalue = float(bdEvalue)
        self.bdScore = float(bdScore)
        self.bdBias = float(bdBias)
        self.exp = exp
        self.numOfDomains = int(numOfDomains)
        self.model = model
        self.desc = desc
    
    def getFullSeqEvalue(self):
        return self.fsEvalue
    
    def getFullSeqScore(self):
        return self.fsScore
    
    def getModel(self):
        return self.model


class HMMERModelBasedCutoffStore(object):
    
    def __init__(self, fileWithCutoffs):
        self.model2EvalueCutoff = dict()
        fwc = open(fileWithCutoffs,"r")
        for line in fwc:
            tokens = line.rstrip("\n\r").split(':')
            if tokens[1] == "evalue":
                self.model2EvalueCutoff[tokens[0]]=float(tokens[2])
        fwc.close()
    
    def hasEvalueCutoff(self,modelName):
        return self.model2EvalueCutoff.has_key(modelName)
    
    def getCutEvalueCutOff(self,modelName):
        return self.model2EvalueCutoff[modelName]


class HMMERParse(object):

    def __init__(self, hmmerScanFileHandle,maxBestModels=5):
        self.maxModels = maxBestModels
        self.fileHandle = hmmerScanFileHandle
        self.line = self.fileHandle.readline()

        self.__skipFirstLines()
        self.__parseModelHits()

    def __parseModelHits(self):
        self.modelHits = []
        import re
        while (not self.line.startswith("Domain annotation")) and len(self.modelHits) <= self.maxModels:
            self.line = self.fileHandle.readline()
            tokens = re.split('\s+', self.line)
            if len(tokens) < 4:
                break
            modelHit = HMMERModelHit(fsEvalue=tokens[1], 
                                     fsScore=tokens[2], 
                                     fsBias=tokens[3], 
                                     bdEvalue=tokens[4], 
                                     bdScore=tokens[5], 
                                     bdBias=tokens[6], 
                                     exp=tokens[7], 
                                     numOfDomains=tokens[8], 
                                     model=tokens[9], 
                                     desc=tokens[10])
            self.modelHits.append(modelHit)
       
    def __skipFirstLines(self):
        while not self.line.startswith("Scores for complete sequence"):
            self.line = self.fileHandle.readline()

        skippedLines = 1
        while not skippedLines > 3:
            self.line = self.fileHandle.readline()
            skippedLines += 1
            
    def nextModelHit(self):
        if len(self.modelHits) > 0:
            return self.modelHits.pop(0)
    
    def goToNextModelAlignments(self):
        while not self.line.startswith(">>"):
            self.line = self.fileHandle.readline()
            if len(self.line) == 0:
                return None
        #for l in self.fileHandle:
        #    self.line = l
        #    if self.line.startswith(">>"):
        #        break
        currentModel = self.line[3:]
        currentModel = currentModel.replace(" \n","")
        currentModel = currentModel.strip()
        self.fileHandle.readline()
        self.fileHandle.readline()
        return currentModel 
    
    def nextAlignmentResult(self):
        self.line = self.fileHandle.readline()
        import re
        tokens = re.split('\s+', self.line)
        if len(tokens) > 2:
            # number = tokens[1]
            valid = tokens[2]
            score = tokens[3]
            bias = tokens[4]
            cEvalue = tokens[5]
            iEvalue = tokens[6]
            hmmStart = tokens[7]
            hmmStop = tokens[8]
            hmmBoundaries = tokens[9]
            alignStart = tokens[10]
            alignStop = tokens[11]
            alignBoundaries = tokens[12]
            
            hmmDomainHit = HMMERDomainHit(valid, score, bias, cEvalue, iEvalue, hmmStart, hmmStop, hmmBoundaries, 
                                          qstart=alignStart, qstop=alignStop, qBound=alignBoundaries)
            return hmmDomainHit


class ModelAnnotator(object):

    def __init__(self, pathToModel):
        self.pathToModel = pathToModel
        annotFile = open(pathToModel+".annot","r")
        line = annotFile.readline()
        self.annotMap = {}
        while len(line) > 0:
            line = line.strip()
            annot = line.rsplit("\t")
            if len(annot) > 1:
                self.annotMap[annot[0]] = annot[1]
            line = annotFile.readline()

        annotFile.close()

    def getAnnotationForKey(self,key):
        return self.annotMap.get(key,None)


            
        
        


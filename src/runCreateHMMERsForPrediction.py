'''
Created on Sep 23, 2011

@author: pmoreno
'''
import sys
from hmmer.core import HMMERSuite
from fastaModifier.core import FastaEntryRemover
import os
from signalReduce.core import SignalReducer
from subprocess import Popen
import Query

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print "USAGE: runCreateHMMERsForPrediction.py path-to-dir-with-clade-fastas " \
              "path-to-complete-alignment consensus-threshold output-path-hmmer-model hmmer-binary-path"
        exit(1)


    pathToDirWithClades=sys.argv[1]
    pathToFastaOveralAlignment=sys.argv[2]
    consensusThres=float(sys.argv[3])
    outPutPathForHMMERModel=sys.argv[4]
    if len(sys.argv) > 5:
        hmmerBinaryPath=sys.argv[5]
    else:
        hmmerBinaryPath = None
    
    
    overallKSModelName=Query.allKSModelName()
    
    # Create the overall PKS HMM
    
    hmmerExec = HMMERSuite(path=hmmerBinaryPath)
    
    filesPathToHMMPress = list()
    filesPathToHMMPress.append("cat")
    allKSModelPath=outPutPathForHMMERModel+overallKSModelName+".hmm"
    filesPathToHMMPress.append(allKSModelPath)
    hmmerExec.runBuild(fastaAlign=pathToFastaOveralAlignment, 
                       output=allKSModelPath, 
                       modelName=overallKSModelName)
    # Remove signal from individual clades
    listing = os.listdir(pathToDirWithClades)
    for entry in listing:
        if entry.startswith("."):
            listing.remove(entry)
    
    sigReducer = SignalReducer(pathToCladesAlignments=pathToDirWithClades, 
                  generalAlignment=pathToFastaOveralAlignment, 
                  outputPath=outPutPathForHMMERModel)
    alignFilesNoSig = sigReducer.run(consensusThreshold=consensusThres)
    
    # Create HMMER model for individual clade no signal
    for i in range(len(alignFilesNoSig)):
        cladeName = listing[i].replace(".fas","")
        cladeFasta = alignFilesNoSig[i] 
        outPutPathHMMModelClade = outPutPathForHMMERModel+cladeName+".hmm"
        filesPathToHMMPress.append(outPutPathHMMModelClade)
        hmmerExec.runBuild(fastaAlign=cladeFasta, 
                           output=outPutPathHMMModelClade, modelName=cladeName)
        
    # Press the models
    outputPathAllModels=outPutPathForHMMERModel+"PKSAllPlusClades_ConsSignal%d.hmm" % int(consensusThres*100)
    fileWithAllHMMModels=open(outputPathAllModels, 'w')
    p1=Popen(args=filesPathToHMMPress,stdout=fileWithAllHMMModels)
    p1.wait()
    fileWithAllHMMModels.close()
    hmmerExec.runPress(hmmModelCompendium=outputPathAllModels)
    
    filesPathToHMMPress.remove("cat")
    # Delete files not necessary to prediction
    for fileName in filesPathToHMMPress:
        os.remove(fileName)
    for fileName in alignFilesNoSig:
        os.remove(fileName)
    
        
        
        
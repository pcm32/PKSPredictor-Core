'''
Created on Sep 21, 2011

@author: pmoreno
'''
from subprocess import Popen
from validation.core import ValidationResult

def catHmmModelsToDirectoryExceptFor(listOfHMMModelPaths,entryToSkip,destinationFile):
    
    filesToCat=list()
    for modelIndex in range(len(listOfHMMModelPaths)):
        if modelIndex == entryToSkip:
            continue
        filesToCat.append(listOfHMMModelPaths[modelIndex])
        # print "Added file to cat:"+listOfHMMModelPaths[modelIndex]
        # os.symlink(listOfHMMModelPaths[modelIndex], destinationDir)
    fileOut = open(destinationFile,'w')
    # print "Outputting to:"+destinationFile
    pCat = Popen(args=filesToCat,executable='cat',stdout=fileOut)
    pCat.wait()
    fileOut.close()
     

from hmmer.core import HMMERDomainHit

if __name__ == '__main__':
    
    # We need to create the models without general signal
    # for this we need to open the complete alignment (all clades),
    # obtain the positions that have above certain threshold.
    # With this positions, we need to open each clade alignment
    # remove one sequence and then alter all the other sequences
    # in the defined positions. Then we build the models and search the
    # removed sequences (no gaps) against all models.
    import sys
    import os
    
    pathToClades = sys.argv[1]
    pathToGeneralAlign = sys.argv[2]
    minConsensusToRemove = float(sys.argv[3])
    outputPath = sys.argv[4]
    hmmerPath = sys.argv[5]
    
    listing = os.listdir(pathToClades)
    for entry in listing:
        if entry.startswith("."):
            listing.remove(entry)
    
    
    # print listing
    
    from signalReduce.core import SignalReducer
    
    minConsAddendumForPath = "%d/" % int(minConsensusToRemove*100)
    dirForNoSignalFaaAlignments=outputPath+"NoSignalCompleteFaaAlignments/"+minConsAddendumForPath
    if not os.path.exists(dirForNoSignalFaaAlignments):
        os.makedirs(dirForNoSignalFaaAlignments)
    signalReducer = SignalReducer(pathToCladesAlignments=pathToClades, generalAlignment=pathToGeneralAlign,outputPath=dirForNoSignalFaaAlignments)
    # elements in the resulting list should be one-to-one to the original clades listing.
    resultingAlignmentFiles = signalReducer.run(consensusThreshold=minConsensusToRemove)
    
    if len(resultingAlignmentFiles)!=len(listing):
        print "Original clade listing and signal reduced listing have different lengths:"
        print listing
        print resultingAlignmentFiles
    
    # print resultingAlignmentFiles
    # these files have all the sequences, but with the required positions turned off
    # Now, for each clade we create the HMM model for the complete alignment that doesn't
    # have the general signal.
    
    from fastaModifier.core import FastaEntryRemover
    from hmmer.core import HMMERSuite, HMMERDomainHit, HMMERModelHit, HMMERParse
    
    dirForNoSignalHMMModels=outputPath+"NoSignalCompleteHMMModels/"+minConsAddendumForPath
    if not os.path.exists(dirForNoSignalHMMModels):
        os.makedirs(dirForNoSignalHMMModels)
        
    hmmer = HMMERSuite(path=hmmerPath)
    hmmerModelsNoGenSignal = list()    
    for clade in range(len(listing)):
        fastaNoSignalPath=resultingAlignmentFiles[clade]
        cladeName = listing[clade].replace(".fas","");
        
        hmmModelOutput = dirForNoSignalHMMModels+cladeName+"_NoGenSig.hmm"
        hmmerModelsNoGenSignal.append(hmmModelOutput)
        hmmer.runBuild(fastaAlign=fastaNoSignalPath, output=hmmModelOutput, modelName=cladeName)
        
    # Now, for each clade, we create a modified version of the alignment with no signal
    # that doesn't contain one entry (for each entry). We obtain that same entry X from the original listed files
    # which include all the signals. We press all models (excepting the original one we just modified)
    # + the newly modified one (comp Y). We search the entry X againts compendium Y
    dirForStatsOutputFiles = outputPath + "StatsFiles/"+minConsAddendumForPath
    if not os.path.exists(dirForStatsOutputFiles):
        os.makedirs(dirForStatsOutputFiles)
        
    filePerCladeEntry = open(dirForStatsOutputFiles+"perCladeEntry.txt",'w')
    filePerClade = open(dirForStatsOutputFiles+"perClade.txt",'w')
    filePerCladeEntry.write("Clade\tEntry\tConsRem\tBestModel\tEvalue\tScore\tEvalueDiff\tScoreDiff\n")
    filePerClade.write("Clade\tConsRem\tTP\tFP\tSens\n")
    
    for clade in range(len(listing)):
        cladeName = listing[clade].replace(".fas","")
        fastaNoSignalToMod=FastaEntryRemover(pathToFasta=resultingAlignmentFiles[clade])
        fastaOriginalToGetEntries=FastaEntryRemover(pathToFasta=pathToClades+listing[clade])
        # for each entry
        entryCounter=1
        positives=0
        negatives=0
        while True:
            print "Running for Clade "+cladeName+" entry %d" % entryCounter
            outputDirForCladeEntry = outputPath+cladeName+"/%d" % entryCounter+"/" + minConsAddendumForPath
            if not os.path.exists(outputDirForCladeEntry):
                os.makedirs(outputDirForCladeEntry)
            pathToSingleEntry=fastaOriginalToGetEntries.getFastaEntry(entryNum=entryCounter, resultFolder=outputDirForCladeEntry)
            if pathToSingleEntry is None:
                os.rmdir(outputDirForCladeEntry)
                break
            pathToModFile=fastaNoSignalToMod.generateFastaWithOutEntry(entryNumToRemove=entryCounter, resultFolder=outputDirForCladeEntry)
            # we take all HMM models and cat to a file, except for our entry
            pathToAllOtherCattedModels=outputDirForCladeEntry+"allOthers_NoGenSig.hmm"
            # print "cat files > "+pathToAllOtherCattedModels+"\n"
            catHmmModelsToDirectoryExceptFor(listOfHMMModelPaths=hmmerModelsNoGenSignal, 
                                                 entryToSkip=clade, 
                                                 destinationFile=pathToAllOtherCattedModels)
            # we create the fasta with the missing entry and build its HMM
            pathToModelWithThisMissingEntry=outputDirForCladeEntry+cladeName+"_NoGenSig.hmm"
            hmmer.runBuild(fastaAlign=pathToModFile, output=pathToModelWithThisMissingEntry, modelName=cladeName)
            # we cat the file with all the other models and the one for this entry
            pathToAllModelsCatted=outputDirForCladeEntry+cladeName+"AllModels_NoGenSig.hmm"
            fileAllModelsCatted=open(pathToAllModelsCatted,'w')
            # print "cat "+pathToAllOtherCattedModels+" "+pathToModelWithThisMissingEntry+" > "+pathToAllModelsCatted
            pCat = Popen(['cat',pathToAllOtherCattedModels,pathToModelWithThisMissingEntry],stdout=fileAllModelsCatted)
            pCat.wait()
            fileAllModelsCatted.close()
            # We make the binary file for the search
            hmmer.runPress(hmmModelCompendium=pathToAllModelsCatted)
            # Finally, we make the search
            pathToHmmScanResult=outputDirForCladeEntry+cladeName+"_result.hmmer"
            hmmer.runScan(query=pathToSingleEntry, output=pathToHmmScanResult, model=pathToAllModelsCatted)
            resultFile=open(pathToHmmScanResult)
            hmmParser = HMMERParse(hmmerScanFileHandle=resultFile, maxBestModels=4)
            # and delete the model files which are heavy, just leave main results for manual inspection.
            os.remove(pathToAllModelsCatted)
            os.remove(pathToAllModelsCatted+".h3f")
            os.remove(pathToAllModelsCatted+".h3i")
            os.remove(pathToAllModelsCatted+".h3m")
            os.remove(pathToAllModelsCatted+".h3p")
            os.remove(pathToAllOtherCattedModels)
            os.remove(pathToModelWithThisMissingEntry)
            
            validationRes = ValidationResult(thisCladeName=cladeName, listOfModelHits=list())
            modelCounter=0
            while modelCounter<3:
                modelHit = hmmParser.nextModelHit()
                if modelHit == None:
                    break
                validationRes.addModelHit(modelHit=modelHit)
                modelCounter += 1
            
            evalueDiff = validationRes.getEvalueDifferenceToSecondModel()
            if evalueDiff is None:
                evalueDiff = "NA"
            else:
                evalueDiff = "%d" % evalueDiff
            scoreDiff = validationRes.getScoreDifferenceToSecondModel()
            if scoreDiff is None:
                scoreDiff = "NA"
            else:
                scoreDiff = "%d" % scoreDiff
            firstModelName = validationRes.getFirstModelName()
            if firstModelName is None:
                filePerCladeEntry.write(cladeName+"\t%d" % entryCounter
                                    +"\t%d" % int(minConsensusToRemove*100)
                                    +"\tNA"
                                    +"\tNA"
                                    +"\tNA"
                                    +"\t"+evalueDiff
                                    +"\t"+scoreDiff
                                    +"\n")
            else:
                filePerCladeEntry.write(cladeName+"\t%d" % entryCounter
                                    +"\t%d" % int(minConsensusToRemove*100)
                                    +"\t"+validationRes.getFirstModelName()
                                    +"\t%E" % validationRes.getFirstModelEvalue()
                                    +"\t%d" % validationRes.getFirsModelScore()
                                    +"\t"+evalueDiff
                                    +"\t"+scoreDiff
                                    +"\n")
            
            if validationRes.isCladeMatch():
                positives+=1
            else:
                negatives+=1
                
            entryCounter+=1
        
        # end of the entry loop
        filePerClade.write(cladeName+
                           "\t%d" % int(minConsensusToRemove*100)+
                           "\t%d" % positives+
                           "\t%d" % negatives+
                           "\t%.2f" % (float(positives)/(float(positives+negatives)))+"\n")
        filePerClade.flush()
    # end of clade loop
    filePerClade.close()

    
    
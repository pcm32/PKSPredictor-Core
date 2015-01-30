__author__ = 'pmoreno'



class SimpleFeatWriter(object):
    '''
    This class receives a seqRecord object and write a simple, one line per feature.
    '''
    def __init__(self, pathForOutput, maxRanking):
        '''
        Constructor
        '''
        self.path = pathForOutput
        self.maxRanking = maxRanking


    def writeFeatures(self,seqObject):
        self.file2write = open(self.path+seqObject.id+".features", mode="w")

        for feature in seqObject.features:
            if "ranking" in feature.qualifiers and feature.qualifiers["ranking"] > self.maxRanking:
                continue
            start = feature.location.start
            end = feature.location.end
            evalue = "N/A"
            ranking = "N/A"
            region = "N/A"
            note = "N/A"
            if "evalue" in feature.qualifiers:
                evalue = feature.qualifiers["evalue"]
            if "ranking" in feature.qualifiers:
                ranking = feature.qualifiers["ranking"]
            if "region" in feature.qualifiers:
                region = feature.qualifiers["region"]
            if "note" in feature.qualifiers:
                note = feature.qualifiers["note"]
            score = feature.qualifiers["score"]
            name = feature.qualifiers["name"]
            type = feature.type

            self.file2write.write("\t".join([str(start), str(end), str(evalue),
                                             str(score), str(ranking), str(region),
                                             type, name, note]))
            self.file2write.write("\n")

        self.file2write.close()



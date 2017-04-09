from verifier.core import Domain_Verifier

__author__ = 'pmoreno'



class SimpleFeatWriter(object):
    '''
    This class receives a seqRecord object and write a simple, one line per feature. Include verification should
    only be used when the new annotation scheme is used.
    '''
    def __init__(self, pathForOutput, maxRanking, include_verification=False):
        '''
        Constructor
        '''
        self.path = pathForOutput
        self.maxRanking = maxRanking
        self.use_verification = include_verification


    def writeFeatures(self, seqObject):
        file2write = open(self.path+seqObject.id+".features", mode="w")

        for feature in seqObject.features:
            if "ranking" in feature.qualifiers and feature.qualifiers["ranking"] > self.maxRanking:
                continue
            start = feature.location.start
            end = feature.location.end
            evalue = "N/A"
            ranking = "N/A"
            region = "N/A"
            note = "N/A"
            subtype = "N/A"
            if "evalue" in feature.qualifiers:
                evalue = feature.qualifiers["evalue"]
            if "ranking" in feature.qualifiers:
                ranking = feature.qualifiers["ranking"]
            if "region" in feature.qualifiers:
                region = feature.qualifiers["region"]
            if "note" in feature.qualifiers:
                note = feature.qualifiers["note"]
            if "subtype" in feature.qualifiers:
                subtype = feature.qualifiers["subtype"]
            score = feature.qualifiers["score"]
            name = feature.qualifiers["name"]
            type = feature.type

            if self.use_verification:
                verification = "N/A"
                if Domain_Verifier.get_qualifier_key() in feature.qualifiers:
                    verification = str(feature.qualifiers[Domain_Verifier.get_qualifier_key()])

                file2write.write("\t".join([str(start), str(end), str(evalue),
                                             str(score), str(ranking), str(region),
                                             type, subtype, name, note, str(verification)]))
                # TODO add support for extra column on Java side.
            else:
                file2write.write("\t".join([str(start), str(end), str(evalue),
                                            str(score), str(ranking), str(region),
                                            type, subtype, name, note]))

            file2write.write("\n")

        file2write.close()



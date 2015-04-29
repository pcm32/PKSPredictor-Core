__author__ = 'pmoreno'


class SeqFeatureSorter(object):

    @staticmethod
    def geographic_sort(seqObject):
        seqObject.features = sorted(seqObject.features, key=lambda feature: feature.location.start)

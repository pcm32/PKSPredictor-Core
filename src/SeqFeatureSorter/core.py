__author__ = 'pmoreno'


class SeqFeatureSorter(object):

    @staticmethod
    def geographic_sort(seqObject):
        seqObject.features = sorted(seqObject.features, key=lambda feature: feature.location.start)

    @staticmethod
    def geographic_with_pile_sort(seqObject):

        def compare_features(f1, f2):
            if "region" in f1.qualifiers and "region" in f2.qualifiers:
                if f1.qualifiers["region"] == f2.qualifiers["region"]:
                    return f1.qualifiers["ranking"] - f2.qualifiers["ranking"]

            return f1.location.start - f2.location.start

        seqObject.features = sorted(seqObject.features, cmp=compare_features)


from abc import ABCMeta, abstractmethod
from Clades.core import CladificationAnnotation

__author__ = 'pmoreno'

class Domain_Verifier:
    __metaclass__ = ABCMeta

    @staticmethod
    def get_qualifier_key():
        return "pass_check"

    @abstractmethod
    def verify(self):
        pass


class KSDomainVerifier(Domain_Verifier):
    """
    This object verifies that a KS domain given is preceeded, in terms of annotated domains, of all the required
    domains that is supposed to be preceeded by. These domains ought to be between the current KS domain and the
    previous (geographically speaking, so upstream) annotated KS domain. The domains expected are defined in the
    annotation file for the cladification.
    """

    def __init__(self, cladififcation_annotation, seqObj):
        """

        :param cladififcation_annotation: A cladificaiton annotation object.
        :type CladificationAnnotation
        :param seqObj: The sequence object where KS preceeding domains should be verified.
        :type seqObj:
        """
        self.__cladification_annot = cladififcation_annotation
        self.__seq_obj = seqObj

    def verify(self):
        """
        Identifies KS seq features in the sequence object given on init, and then verifies that required preceding domains
        are found.
        :return:
        """
        sorted_features = sorted(self.__seq_obj.features, key=lambda feature: feature.location.start)
        previous_stack = 0
        previous_stack_end = -1
        current_preceding_module = None
        for feature in sorted_features:
            if feature.qualifiers["subtype"] is "KS":
                if previous_stack == 0:
                    '''
                    We have no previous stack recorded, we record this one and move on.
                    '''
                    previous_stack = feature.qualifiers["region"]
                    previous_stack_end = feature.location.end
                elif previous_stack == feature.qualifiers["region"]:
                    '''
                    We are on a KS of the same stack/region just recorded, move on if it is the first,
                    verify otherwise
                    '''
                    if previous_stack > 1:
                        self.__run_check(feature, current_preceding_module)

                elif previous_stack < feature.qualifiers["region"]:
                    current_preceding_module = self.__seq_obj[previous_stack_end+1:feature.location.start-1]
                    self.__run_check(feature, current_preceding_module)
                    previous_stack_end = feature.location.end
                    previous_stack = feature.qualifiers["region"]

    def __run_check(self, feature, current_preceding_module):
        domains_to_be_found = list(self.__cladification_annot.get_domain_requirements(feature.qualifiers["name"]))
        # Only work if the list of domains to be verified is not empty
        if domains_to_be_found:
            for feature_preceding_module in current_preceding_module.features:
                # afterwards we look at any remaining items in domains_to_be_found
                domains_to_be_found.remove(feature_preceding_module.qualifiers["name"])
        if len(domains_to_be_found)>0:
            feature.qualifiers[self.get_qualifier_key()] = False
        else:
            feature.qualifiers[self.get_qualifier_key()] = True


class DH_PS_DomainVerifier(Domain_Verifier):
    """
    Checks whether DH domains include geographically the fuzzpro pattern that
    either confirms it as a DH domain ("DH_conf"), or the pattern that leads to
    changing that domain to PS.
    """

    # TODO This functionality is quite hard coded and should be improved in future.

    dh_confirmation_pattern = "DH_conf"
    ps_confirmation_pattern = "PS_not_DH"

    def __init__(self, seq_obj):
        self.__seq_obj = seq_obj

    def verify(self):
        """
        Goes through the seq_obj features, to check each DH domain for the appearance
        of a particular signature.
        """
        for feature in self.__seq_obj.features:
            if feature.qualifiers["name"] in ["DH", "PS"]:
                slice = self.__seq_obj[feature.location.start:feature.location.end]
                found_confirmation = False
                found_ps_hint = False
                for sub_features in slice.features:
                    if sub_features.id is DH_PS_DomainVerifier.dh_confirmation_pattern and not found_confirmation:
                        found_confirmation = True
                    if sub_features.id is DH_PS_DomainVerifier.ps_confirmation_pattern and not found_ps_hint:
                        found_ps_hint = True
                if feature.qualifiers["name"] is "DH":
                    if found_confirmation and not found_ps_hint:
                        feature.qualifiers[self.get_qualifier_key()] = True
                    else:
                        # mark it as not verified, a PS should be given preference to this.
                        feature.qualifiers[self.get_qualifier_key()] = False
                else:
                    # PS domain case
                    if not found_confirmation and found_ps_hint:
                        feature.qualifiers[self.get_qualifier_key()] = True
                    else:
                        # mark it as not verified, a DH should be given preference to this.
                        feature.qualifiers[self.get_qualifier_key()] = False


















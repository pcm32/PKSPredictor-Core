from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from verifier.core import KSDomainVerifier
from Clades.core import CladificationAnnotationReader
import pkg_resources, os, random
import unittest

__author__ = 'pmoreno'

class test_ks_domain_verifier(unittest.TestCase):

    def test_run_Clade3_with_PS(self):
        """
        Checks that Clade 3 is correctly verified when finding a PS domain before it.
        """
        resource_package = "Clades.test"  ## Could be any module/package name.
        #resource_path = os.path.join('annotation_file_1.tab')
        resource_path = os.path.join('annotation_file_test1.tsv')
        test_annotation_file = pkg_resources.resource_filename(resource_package, resource_path)
        clade_annot = CladificationAnnotationReader().get_annotation(test_annotation_file)


        seqobj = SeqRecord(Seq(''.join(
            [random.choice('ACDEFGHIKLMNPQRSTVWY') for n in xrange(100)])))
        self.__add_domain_from_to(seqobj, 2, 20, "Clade_2", is_ks=True, region=1)
        self.__add_domain_from_to(seqobj, 25, 30, "PS")
        clade_3 = self.__add_domain_from_to(seqobj, 40, 50, "Clade_3", is_ks=True, region=2)
        verifier = KSDomainVerifier(seqObj=seqobj, cladififcation_annotation=clade_annot)

        changes = verifier.verify()
        self.assertTrue(clade_3.qualifiers['verified'])


        # print changes
        # apply changes?

    def __add_domain_from_to(self, seq_record, start, end, name, is_ks=False, region=None):
        qual = {"evalue": "1E-7",
                "score": 100,
                "name": name,
                "subtype": name if not is_ks else "KS"}
        if region is not None:
            qual["region"] = region
        sf = SeqFeature(
            FeatureLocation(start=start, end=end),
            type="domain", strand=1, id=name, qualifiers=qual
        )

        seq_record.features.append(sf)
        return sf

if __name__ == '__main__':
    unittest.main()



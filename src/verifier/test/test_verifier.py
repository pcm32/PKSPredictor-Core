from verifier.core import KS_Domain_Verifier
import pkg_resources, os
from _pytest import unittest

__author__ = 'pmoreno'

class test_ks_domain_verifier(unittest.TestCase):
    def test_run(self):
        resource_package = __name__  ## Could be any module/package name.
        resource_path = os.path.join('verifier', 'test', 'annotation_file_1')
        test_annotation_file = pkg_resources.resource_string(resource_package, resource_path)

        verifier = KS_Domain_Verifier(seqObj=seqAA, pathToAnnotationFile=test_annotation_file)

        changes = verifier.verify()

        # print changes
        # apply changes?



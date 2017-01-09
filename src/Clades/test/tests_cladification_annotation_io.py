import unittest
import os
import pkg_resources
__author__ = 'pmoreno'

import sys
sys.path.append('../../')

from Clades.core import CladificationAnnotationReader


class TestCladificationAnnotationReader(unittest.TestCase):
    def test_run(self):
        reader = CladificationAnnotationReader()
        resource_package = __name__  ## Could be any module/package name.
        resource_path = os.path.join('annotation_file_test1.tsv')
        test_annotation_file = pkg_resources.resource_filename(resource_package, resource_path)
        clad_annotation = reader.get_annotation(annotation_path=test_annotation_file)
        #self.assertTrue(len(seqRec.features)==3,"Should find 3 sequence features")



if __name__ == '__main__':
    unittest.main()

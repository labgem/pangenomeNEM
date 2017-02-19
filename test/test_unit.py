import unittest
from collections import defaultdict
import sys
sys.path.append('../src/')
from classes import *
import argparse

logger = logging.getLogger()
logger.level = logging.DEBUG
logger.addHandler(logging.StreamHandler(sys.stdout))

class TestPangenomeMethods(unittest.TestCase):

    def import_files_progenome(self): 
        annotation_file = open("../data/specI-specI_v2_Cluster335.gene_annotations.tsv","r")
        eggNOG_clusters_file = open("../data/specI_v2_Cluster335.eggNOG_groups.tsv","r")

        return(annotation_file,eggNOG_clusters_file)
    def test_import_progenome_with_singleton(self):
        (annotation_file,eggNOG_clusters_file) = self.import_files_progenome()
        pan = Pangenome("progenome", annotation_file, eggNOG_clusters_file, False)
        
    def test_import_progenome_without_singleton(self):
        (annotation_file,eggNOG_clusters_file) = self.import_files_progenome()
        pan = Pangenome("progenome", annotation_file, eggNOG_clusters_file, True)

        #self.assertEqual('foo'.upper(), 'FOO')
        #self.assertTrue('FOO'.isupper())
        #self.assertFalse('Foo'.isupper())
        #with self.assertRaises(TypeError):
        #    pass #test()
    def test_import_microscope(self):
        pass
    def test_sub_pangenome(self):
        pass

if __name__ == '__main__':
    unittest.main()

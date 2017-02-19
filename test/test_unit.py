import unittest
from collections import defaultdict
import sys
sys.path.append('../src/')
from classes import *
import argparse

class TestPangenomeMethods(unittest.TestCase):

    def test_import_progenome(self):
        
        parser = argparse.ArgumentParser()

	parser.add_argument("-u", "--use", nargs=1, help = "what source of data to import : 'progenome', 'microscope', 'prokka/roary' 'prokka/MMseqs2' 'MEG'", required=True)
	group_progenome.add_argument('-a', '--annotations', type=argparse.FileType('r'), nargs=1, help="The tsv file provided by progenome containing the gene annotations")
	group_progenome.add_argument('-c', '--eggNOG_clusters', type=argparse.FileType('r'), nargs=1, help="The tsv file provided by progenome containing eggNOG orthologous groups related to the annotated genomes")
	parser.add_argument('-d', '--outputdirectory', type=str, nargs=1, default="output.dir", help="The output directory", required=True)
        
        pan = Pangenome(options)


        self.assertEqual('foo'.upper(), 'FOO')
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())
        with self.assertRaises(TypeError):
            pass #test()
    def test_import_microscope(self):
        pass
    def test_sub_pangenome(self):
        pass

if __name__ == '__main__':
    unittest.main()

import unittest
from collections import defaultdict
import sys
sys.path.append('../src/')
from classes import *
import argparse

logging.basicConfig(level=logging.DEBUG, format='\n%(asctime)s %(filename)s:l%(lineno)d:t%(thread)d %(levelname)s\t%(message)s', datefmt='%H:%M:%S')

class TestPangenomeMethods(unittest.TestCase):

    def import_files_progenome(self): 
        annotation_file = open("../../../Téléchargements/specI-specI_v2_Cluster36.gene_annotations.tsv","r")#../data/specI-specI_v2_Cluster335.gene_annotations_subset.tsv
        eggNOG_clusters_file = open("../../../Téléchargements/specI_v2_Cluster36.eggNOG_groups.tsv","r")#../data/specI_v2_Cluster335.eggNOG_groups_subset.tsv

        return(annotation_file,eggNOG_clusters_file)
    def test_import_progenome_with_singleton(self):
        (annotation_file,eggNOG_clusters_file) = self.import_files_progenome()
        pan = Pangenome("progenome", annotation_file, eggNOG_clusters_file, False)
        self.assertEqual(len(Pangenome.annotations),799)
        self.assertEqual(Pangenome.annotations[-1],tuple(['Cther_0132', 'CDS', '492476.PRJNA28257', '492476.PRJNA28257.ABVG02000008', '0EZJM', 26975, 27403, '-']))
        self.assertEqual(Pangenome.annotations[-2],tuple(['Cther_0131', 'CDS', '492476.PRJNA28257', '492476.PRJNA28257.ABVG02000008', '0EZGW', 25713, 26951, '-']))
        self.assertEqual(Pangenome.annotations[-3],tuple(['Cther_0130', 'CDS', '492476.PRJNA28257', '492476.PRJNA28257.ABVG02000008', 'Cther_0130', 25367, 25504, '-']))
       
        for org in pan.organism_positions.keys():
            for i in pan.annotation_positions[org]:
                self.assertTrue(Pangenome.annotations[i][2] == org)
                #add a coverage test

        self.assertEqual(len(pan.familly_positions),595)
        self.assertEqual(len(pan.familly_positions),len(Pangenome.presences_absences))
        self.assertEqual(len(pan.organism_positions),5)   

        #test empty annoatation file
        #test empty eggNOGcluster file 
    def test_import_progenome_without_singleton(self):
        (annotation_file,eggNOG_clusters_file) = self.import_files_progenome()
        pan = Pangenome("progenome", annotation_file, eggNOG_clusters_file, True)

        #with self.assertRaises(TypeError):
        #    pass #test()
    def test_import_microscope(self):
        pass
    def test_sub_pangenome(self):
        (annotation_file,eggNOG_clusters_file) = self.import_files_progenome()
        pan = Pangenome("progenome", annotation_file, eggNOG_clusters_file, False) 
        sub_pan1 = pan.sub_pangenome(["492476.PRJNA28257"])
        sub_pan2 = pan.sub_pangenome(["492476.PRJNA28257","1094188.PRJNA73891"])
        sub_pan_all = pan.sub_pangenome(pan.organism_positions.keys())

if __name__ == '__main__':
    unittest.main()

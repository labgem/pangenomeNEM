import unittest
from collections import defaultdict
import sys
sys.path.append('../src/')
from classes import *
import argparse

logging.basicConfig(level=logging.DEBUG, format='\n%(asctime)s %(filename)s:l%(lineno)d:t%(thread)d %(levelname)s\t%(message)s', datefmt='%H:%M:%S')

class TestPangenomeMethods(unittest.TestCase):

    def import_files_progenome(self): 
        annotation_file = open("../data/specI-specI_v2_Cluster335.gene_annotations_subset.tsv","r")
        eggNOG_clusters_file = open("../data/specI_v2_Cluster335.eggNOG_groups_subset.tsv","r")

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
        self.assertEqual(pan.nb_organisms,5)
        self.assertEqual(pan.pan_size,595)
        self.assertEqual(pan.core_size,0)
        self.assertEqual(len(pan.familly_ratio),len(pan.familly_positions))
        self.assertEqual(len(pan.core_list),0)
        #test empty annoatation file
        #test empty eggNOGcluster file 
    def test_import_progenome_without_singleton(self):
        (annotation_file,eggNOG_clusters_file) = self.import_files_progenome()
        pan = Pangenome("progenome", annotation_file, eggNOG_clusters_file, True)

        self.assertEqual(Pangenome.annotations[-1],tuple(['Cther_0132', 'CDS', '492476.PRJNA28257', '492476.PRJNA28257.ABVG02000008', '0EZJM', 26975, 27403, '-']))
        self.assertEqual(Pangenome.annotations[-2],tuple(['Cther_0131', 'CDS', '492476.PRJNA28257', '492476.PRJNA28257.ABVG02000008', '0EZGW', 25713, 26951, '-']))
        self.assertEqual(Pangenome.annotations[-3],tuple(['Cther_0130', 'CDS', '492476.PRJNA28257', '492476.PRJNA28257.ABVG02000008', None, 25367, 25504, '-']))
        
        for org in pan.organism_positions.keys():
            for i in pan.annotation_positions[org]:
                self.assertTrue(Pangenome.annotations[i][2] == org)
                #add a coverage test
        self.assertEqual(len(pan.familly_positions),518)
        self.assertEqual(len(pan.familly_positions),len(Pangenome.presences_absences))
        self.assertEqual(len(pan.organism_positions),5)   
        self.assertEqual(pan.nb_organisms,5)
        self.assertEqual(pan.pan_size,518)
        self.assertEqual(pan.core_size,0)
        self.assertEqual(len(pan.familly_ratio),len(pan.familly_positions))
        self.assertEqual(len(pan.core_list),0)
        #with self.assertRaises(TypeError):
        #    pass #test()
    def test_import_microscope(self):
        pass
    def test_sub_pangenome(self):
        (annotation_file,eggNOG_clusters_file) = self.import_files_progenome()
        pan = Pangenome("progenome", annotation_file, eggNOG_clusters_file, False) 
        sub_pan1 = pan.sub_pangenome(["492476.PRJNA28257"])
        
        self.assertEqual(["492476.PRJNA28257"],sub_pan1.annotation_positions.keys())


        for org in sub_pan1.organism_positions.keys():
            for i in sub_pan1.annotation_positions[org]:
                self.assertTrue(Pangenome.annotations[i][2] == org)
                #add a coverage test

        self.assertEqual(len(sub_pan1.familly_positions),181)
        self.assertEqual(len(sub_pan1.organism_positions),1)   
        self.assertEqual(sub_pan1.nb_organisms,1)
        self.assertEqual(sub_pan1.pan_size,181)
        self.assertEqual(sub_pan1.core_size,181)
        self.assertEqual(len(sub_pan1.familly_ratio),len(sub_pan1.familly_positions))
        self.assertEqual(len(sub_pan1.core_list),181)

        sub_pan2 = pan.sub_pangenome(["492476.PRJNA28257","1094188.PRJNA73891"])
        
        self.assertEqual(["1094188.PRJNA73891","492476.PRJNA28257"],sub_pan2.annotation_positions.keys())


        for org in sub_pan2.organism_positions.keys():
            for i in sub_pan2.annotation_positions[org]:
                self.assertTrue(Pangenome.annotations[i][2] == org)
                #add a coverage test

        self.assertEqual(len(sub_pan2.familly_positions),208)
        self.assertEqual(len(sub_pan2.organism_positions),2)   
        self.assertEqual(sub_pan2.nb_organisms,2)
        self.assertEqual(sub_pan2.pan_size,208)
        self.assertEqual(sub_pan2.core_size,7)
        self.assertEqual(len(sub_pan2.familly_ratio),len(sub_pan2.familly_positions))
        self.assertEqual(len(sub_pan2.core_list),7)

        sub_pan_all = pan.sub_pangenome(pan.organism_positions.keys())
        print(pan.organism_positions.keys())
        print(sub_pan_all.annotation_positions.keys())

        for org in sub_pan_all.organism_positions.keys():
            for i in sub_pan_all.annotation_positions[org]:
                self.assertTrue(Pangenome.annotations[i][2] == org)
                #add a coverage test

        self.assertEqual(sub_pan_all.familly_positions,pan.familly_positions)
        self.assertEqual(sub_pan_all.organism_positions,pan.organism_positions) 
        self.assertEqual(sub_pan_all.nb_organisms,pan.nb_organisms)
        self.assertEqual(sub_pan_all.pan_size,pan.pan_size)
        self.assertEqual(sub_pan_all.core_size,pan.core_size)
        self.assertEqual(sub_pan_all.familly_ratio,pan.familly_ratio)
        self.assertEqual(sub_pan_all.core_list,pan.core_list)

        with self.assertRaises(ValueError):
            sub_bar = pan.sub_pangenome(["foo"])
if __name__ == '__main__':
    unittest.main()

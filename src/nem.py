# !/usr/bin/python
# -*- coding: iso-8859-1 -*-

from util import *
from classes import *

import argparse
import os
import sys
import errno
from collections import defaultdict
import warnings
import math
import random
import string
import subprocess

logging.basicConfig(level = logging.DEBUG, format = '\n%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s', datefmt='%H:%M:%S')

SUBPATH = "/NEMOUT"

OUTPUTDIR = ""
OUTPUTDIR_EVOL = ""

if __name__=='__main__':

	parser = argparse.ArgumentParser(description='Calculate and partition gene families of a pangenome from annotated genomes and gene families')
	parser.add_argument('--version', action='version', version='0.1')
	parser.add_argument('-o', '--organisms', type=argparse.FileType('r'), nargs=1, help="The tsv file provided by progenome containing the gene annotations", required=True)
	parser.add_argument('-g', '--gene_families', type=argparse.FileType('r'), nargs=1, help="The tsv file provided by progenome containing eggNOG orthologous groups related to the annotated genomes",  required=True)
	parser.add_argument('-d', '--outputdirectory', type=str, nargs=1, default="output.dir", help="The output directory")
	parser.add_argument('-r', '--remove_high_copy_number_families', type=int, nargs=1, default=[-1], help="Remove families having a number of copy of one families above or equal to this threshold in at least one organism (0 or negative value keep all families whatever their occurence)")
	parser.add_argument("-n", "--neighborcomputation", type=bool, default=1, help="Consider neighboors for the statiscal partionning")
	parser.add_argument("-k", "--classnumber", type=int, nargs=1, default=[3], help="Number of classes to consider for other results than evolution (default = 3)")
	parser.add_argument("-t", "--num_thread", default=1, nargs=1, help="The number of thread to use, 0 for autodetection")
	parser.add_argument("-w", "--verbose", default=False, action="store_true", help="verbose")

	options = parser.parse_args()

	OUTPUTDIR       = options.outputdirectory[0]
	NEMOUTPUTDIR    = OUTPUTDIR+"/NEM_results/"
	FIGUREOUTPUTDIR = OUTPUTDIR+"/"+"figures/"
	for directory in [OUTPUTDIR, NEMOUTPUTDIR, FIGUREOUTPUTDIR]:
		if not os.path.exists(directory):
			os.makedirs(directory)
		else:
			logging.getLogger().error(directory+" already exist")
			exit()
	EXACT_CORE_FILE            = OUTPUTDIR+"/"+"exact_core_cluster.txt"
	STATS_EXACT_FILE           = OUTPUTDIR+"/"+"stats_exact.txt"
	STATS_NEM_FILE             = OUTPUTDIR+"/"+"stats_nem_"+str(options.classnumber[0])+".txt" 
	ORGANISMS_FILE             = OUTPUTDIR+"/"+"organisms.txt"	

	pan = Pangenome("file", options.organisms[0],  options.gene_families[0], options.remove_high_copy_number_families[0])

	pan.partition(NEMOUTPUTDIR+"/nb"+str(pan.nb_organisms)+"_k"+str(options.classnumber[0])+"i_0", options.classnumber[0], write_graph = "gexf")	

	logging.getLogger().info(pan)

	logging.getLogger().info("> Creation of a exact_core_cluster.txt file with all core clusters")
	outfile = open(ORGANISMS_FILE,"w")
	outfile.writelines(["%s\n" % item  for item in (pan.organisms)])
	outfile.close()

	core_cluster_file = open(EXACT_CORE_FILE,"w")
	core_cluster_file.write("\n".join(pan.core_list))
	core_cluster_file.write("\n")
	core_cluster_file.close()

	stat_exact_file = open(STATS_EXACT_FILE,"w")
	stat_exact_file.write(str(pan.nb_organisms) + "\t"+ str(pan.core_size) + "\t"+str(pan.pan_size-pan.core_size) + "\t" + str(pan.pan_size)+ "\n")
	stat_nem_file = open(STATS_NEM_FILE,"w")
	stat_nem_file.write(str(pan.nb_organisms) + str(pan.partitions_size["P"]) + "\t"+str(pan.partitions_size["S"]) + "\t" + str(pan.partitions_size["C"]) + "\t" + str(pan.pan_size)+ "\t"+str(pan.BIC)+"\n")

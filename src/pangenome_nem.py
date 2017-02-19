# !/usr/bin/python
# -*- coding: iso-8859-1 -*-

from util import *
from classes import *

import argparse
import pandas as pd
import os
import sys
import errno
from collections import defaultdict
import numpy as np
import warnings
import math
from multiprocessing import Pool
import itertools
import shutil
import random
import string
from joblib import Parallel, delayed
import multiprocessing
import subprocess
from Bio import SeqIO
import json
import urllib2 
from goatools.obo_parser import GODag
from goatools.base import download_go_basic_obo
from goatools import obo_parser
from goatools import semantic
import requests
import re

warnings.filterwarnings("ignore", "Unknown table.*")

MASH_LOCATION = "/env/cns/proj/agc/home/ggautrea/ponderation/mash-Linux64-v1.1/"#ADD MASH
SUBPATH = "/NEMOUT"

OUTPUTDIR = ""
OUTPUTDIR_EVOL = ""

def run(i,prefix, pan, sub_organisms = None):
	if sub_organisms != None:
		
		pan = pan.sub_pangenome(sub_organisms)
		
	pan.classify(i,prefix)
	print(pan)
def findCOG(pangenome):

	print("HERE")	
	headers = {
	    'Origin': 'http://eggnogdb.embl.de',
	    'Accept-Encoding': 'gzip, deflate',
	    'Accept-Language': 'en-US,en;q=0.8,fr;q=0.6',
	    'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/53.0.2785.143 Safari/537.36',
	    'Content-Type': 'application/json;charset=UTF-8',
	    'Accept': 'application/json, text/plain, */*',
	    'Referer': 'http://eggnogdb.embl.de/',
	    'Connection': 'keep-alive',
	    'DNT': '1',
	}
		
	ortho_2_COG_funcat = pd.DataFrame("S", index = pangenome.families, columns = ["funcat"])
	for ortho in sorted(pangenome.families):
		if len(ortho) == 5:
			try:
				data = '{"desc":"","seqid":"","target_species":"","level":"","nognames":"'+ortho+'","page":0}'
				content = requests.post('http://eggnogapi.embl.de/meta_search', headers=headers, data=data)
			except :
				continue
			if content.status_code == 200:
				funcat = json.loads(content.text)["matches"][0]["funcat"]
				ortho_2_COG_funcat.loc[ortho,"funcat"] = funcat
			print(ortho+"  "+str(ortho_2_COG_funcat.loc[ortho,:]))
	return(ortho_2_COG_funcat)

def findGO(pangenome):

	obo_fname = download_go_basic_obo()
	go_parser = obo_parser.GODag("go-basic.obo")

	API= "http://eggnogapi.embl.de/nog_data/json/go_terms/"
	
	GO_roots = { "Biological Process": set(['GO:0008150']), "Cellular Component" : set(['GO:0005575']), "Molecular Function" : set(['GO:0003674'])}

	ortho_2_GO = pd.DataFrame("S", index = sorted(pangenome.families), columns = GO_roots.keys())
	for ortho in sorted(pangenome.families):
		link = API+ortho
		try:
			content = urllib2.urlopen(link)
		except urllib2.URLError:
			continue
		if content.getcode() == 200:
			go_terms = json.loads(content.read())["go_terms"]
			for namespace in GO_roots.keys():
				if namespace in go_terms.keys():

					go_terms_name_space = go_terms[namespace]
					result=semantic.common_parent_go_ids([go_id[0] for go_id in go_terms_name_space],go_parser)
					result -= GO_roots[namespace]
					ortho_2_GO.loc[ortho,namespace]=go_parser[result.pop()].name if len(result) > 0 else "Unknown"
		print(ortho+"  "+str(ortho_2_GO.loc[ortho,:]))
	return(ortho_2_GO)


def calc_mash_distance(fasta, OUTPUTDIR, num_thread):

	MASH_DIRECTORY = OUTPUTDIR+"/"+"mash/"
	if not os.path.exists(MASH_DIRECTORY):
		os.makedirs(MASH_DIRECTORY)
	
	mash_parameters = set()
	fasta_sequences = SeqIO.parse(fasta,'fasta')
	for fasta in fasta_sequences:
		elements = fasta.id.split(".")
		out_file_name = MASH_DIRECTORY+elements[0]+"."+elements[1]
		mash_parameters.add(out_file_name)

		with open(out_file_name+".fasta","a") as out_file:
			SeqIO.write(fasta,out_file,"fasta")

	command = MASH_LOCATION+'mash sketch -n -p ' + str(num_thread) + " -o "+MASH_DIRECTORY+"all_sketch.msh "+MASH_DIRECTORY+"*.fasta"
	print(command)
	proc = subprocess.Popen(command, shell=True)
	proc.communicate()

	command = MASH_LOCATION+'mash dist -t '+(MASH_DIRECTORY+'all_sketch.msh ') * 2 + " > "+OUTPUTDIR+"/mash_distance.csv"
	print(command)
	proc = subprocess.Popen(command, shell=True,stdout=subprocess.PIPE)
	proc.communicate()
	
	distances = pd.read_csv(OUTPUTDIR+"/mash_distance.csv", sep="\t", index_col =0)
	organisms_names = {i : os.path.splitext(os.path.basename(i))[0] for i in distances.index}
	distances.rename(index=organisms_names,columns=organisms_names, inplace=True)
	return(distances)

if __name__=='__main__':

	parser = argparse.ArgumentParser(description='Calculate and classify genes of a pangenome from annotated genomes and ortholog groups')
	parser.add_argument('--version', action='version', version='0.1')
		
	parser.add_argument("-u", "--use", nargs=1, help = "what source of data to import : 'progenome', 'microscope', 'prokka/roary' 'prokka/MMseqs2' 'MEG'", required=True)

	#TODO add_argument_group and .add_mutually_exclusive_group()
	group_progenome = parser.add_argument_group("progenome options")
	group_progenome.add_argument('-a', '--annotations', type=argparse.FileType('r'), nargs=1, help="The tsv file provided by progenome containing the gene annotations")
	group_progenome.add_argument('-c', '--eggNOG_clusters', type=argparse.FileType('r'), nargs=1, help="The tsv file provided by progenome containing eggNOG orthologous groups related to the annotated genomes")

	group_microscope = parser.add_argument_group("microscope options")
	group_microscope.add_argument('-o', '--organisms', type=str, nargs="*", default=False, help="Choose organisms to compare for the core/pan-genome analysis")
	group_microscope.add_argument("-i", "--micfamparam", type=int, nargs=1, default=1, help="MICFAM parameter : 1 (80/80), 2 (50/80)")	
	group_microscope.add_argument("-r", "--artefact", default=False, action="store_true",help="Include artefacts to the analysis (surpredicted genes)")	

	#group_roary = parser.add_argument_group("roary/prokka options")
	#group_roary.add_argument('-v',"--roary_csv_file", type=argparse.FileType('r'), nargs=1, help="roary 'gene_presence_absence.csv' file")
	#group_roary.add_argument('-g',"--gff", type=argparse.FileType('r'), nargs="*", help="gff files provided by prokka")
	
	group_mmseqs2 = parser.add_argument_group("mmseqs2/prokka options")
	group_mmseqs2.add_argument('-v',"--tsv_file", type=argparse.FileType('r'), nargs=1, help="MMseqs2 tsv file")
	group_mmseqs2.add_argument('-g',"--gff", type=argparse.FileType('r'), nargs="*", help=" 'gff files provided by prokka")

	group_MEG = parser.add_argument_group("MEG options") 

	parser.add_argument("-p", "--ponderation", default=False, action="store_true", help="use mash to calculate distance between genomes and based on these distances ponderate underepresented genomes via a MDS approach")
	ponderation = parser.add_mutually_exclusive_group()
	ponderation.add_argument("-x", "--mash_fasta",type =argparse.FileType('r'),nargs=1, help="a csv file containing a matrix of distances between organisms")
	ponderation.add_argument("-f", "--distances_file", type=argparse.FileType('r'), nargs=1, help="fasta file containing contig required to computed distance between genome")

	parser.add_argument("-s", "--remove_singleton", default=False, action="store_true", help="Remove singleton to the pan-genome")
	#parser.add_argument("-m", "--remove_ME", default=False, action="store_true", help="Remove the mobile elements (transposons, integrons, prophage gene) of the pan-genome")
	parser.add_argument('-d', '--outputdirectory', type=str, nargs=1, default="output.dir", help="The output directory", required=True)
	parser.add_argument("-n", "--neighborcomputation", type=int, default=1, help="Consider neighboors for the analysis with the max neighbor distance (integer) (0 = infinite distance)")
	parser.add_argument("-k", "--classnumber", type=int, nargs=1, default=[3], help="Number of classes to consider for other results than evolution (default = 3)")
	parser.add_argument("-e", "--evolution", default=False, action='store_true', help="compute several sample of organism to enable to build the evolution curves.")	
	parser.add_argument("-t", "--num_thread", default=1, nargs=1, help="The number of thread to use, 0 for autodetection")
	parser.add_argument("-w", "--verbose", default=False, action="store_true", help="verbose")
	parser.add_argument("-]", "--max.resampling", default = 30, nargs=1, help="Number max of subsamples in each combinaison of organisms")
	parser.add_argument("-[", "--min.resampling", default = 10, nargs=1, help="Number max of subsamples in each combinaison of organisms")
	
	options = parser.parse_args()

	OUTPUTDIR       = options.outputdirectory[0]
	NEMOUTPUTDIR    = OUTPUTDIR+"/NEM_results/"		
        FIGUREOUTPUTDIR = OUTPUTDIR+"/"+"figures/"
        for directory in [OUTPUTDIR, NEMOUTPUTDIR, FIGUREOUTPUTDIR]:
	        if not os.path.exists(directory):
 		        os.makedirs(directory)
		else:
			print(directory+" already exist")
			exit()
        EXACT_CORE_FILE = OUTPUTDIR+"/"+"exact_core_cluster.txt"
	EVOLUTION_STATS_NEM_FILE = OUTPUTDIR+"/"+"evolution_stats_nem_"+str(options.classnumber[0])+".txt" 
	EVOLUTION_STATS_EXACT_FILE = OUTPUTDIR+"/"+"evolution_stats_exact.txt"
	ORGANISMS_FILE       = OUTPUTDIR+"/"+"organisms.txt"	
	ONTOLOGY_FILE = OUTPUTDIR+"/"+"ontology.txt"
	COG_FILE = OUTPUTDIR+"/"+"COG.txt"

	max_neighbordistance = float("inf") if int(options.neighborcomputation)==0 else int(options.neighborcomputation)
	num_thread = multiprocessing.cpu_count() if options.num_thread==0 else int(options.num_thread)
	num_thread = 1 if num_thread==0 else num_thread 
	
	distances = None
	
	if options.ponderation:
		if options.mash_fasta[0]:
			distances = calc_mash_distance(options.mash_fasta[0],OUTPUTDIR,num_thread)
		#elif options.perdersen_distance:
		# inspired by Hierarchical Sets proposed by perdersen 2017
				
		else:
			distances = pd.read_csv(options.distancesfile, sep="\t", index_col = 0)
			
	
	pan = pangenome(options, distances)
	print(pan.organisms)
	
	outfile = open(ORGANISMS_FILE,"w")
	outfile.writelines(["%s\n" % item  for item in sorted(pan.organisms)])
	outfile.close()

	core_cluster_file = open(EXACT_CORE_FILE,"w")
	core_cluster_file.write("\n".join(pan.cluster_core.index))
	core_cluster_file.write("\n")
	core_cluster_file.close()
	print "> Creation of a exact_core_cluster.txt file with all core clusters"

	print(pan)
	
	if options.ponderation:	
		pan.ponderate()
	
	if options.evolution:
		arguments = list()       
		organisms = list(pan.organisms)	
		total_combinations = oidsCombinations(range(0,len(organisms)),10,1000,10)
		del total_combinations[len(organisms)]

		nb_total_combinations = 0	

		for c in total_combinations:
			nb_total_combinations += len(total_combinations[c])
		print "..... Preparing to classify 1 pangenome + "+str(nb_total_combinations)+" subsampled pangenomes ....."
		for comb_size in total_combinations:
			for combination in total_combinations[comb_size]:
				arguments.append([NEMOUTPUTDIR, pan, [organisms[i] for i in combination]])
		
		random.shuffle(arguments)
		arguments.insert(0, [NEMOUTPUTDIR, pan])
		Parallel(n_jobs=num_thread)(delayed(run)(i,*arguments[i]) for i in range(len(arguments)))

	else:
		pan.classify(0,NEMOUTPUTDIR)

	command = 'cat '+NEMOUTPUTDIR+"/*/nem.stat | sort -k1n > "+EVOLUTION_STATS_NEM_FILE
	print(command)
	proc = subprocess.Popen(command, shell=True)
	proc.communicate()	
	command = 'cat '+NEMOUTPUTDIR+"/*/exact.stat | sort -k1n > "+EVOLUTION_STATS_EXACT_FILE
	print(command)
	proc = subprocess.Popen(command, shell=True)
	proc.communicate()

	#ortho_2_GO = findGO(pan)
	#ortho_2_GO.to_csv(ONTOLOGY_FILE)
	ortho_2_COG_funcat = findCOG(pan)
	ortho_2_COG_funcat.to_csv(COG_FILE)

	with open(NEMOUTPUTDIR+"/nborg"+str(len(pan.organisms))+"_k"+str(options.classnumber[0])+"_i0/file.mf","r") as mf_file:
		for i,line in enumerate(mf_file):
			print(str(i)+line)
			if i==6:
				elements = line.split()
				print(line)
				print(elements)
				if options.ponderation:#meaning Normal model
					BIC = 2 * float(elements[3]) - (options.classnumber[0] * (len(pan.organisms) + 1) + 1 + options.classnumber[0] - 1) * math.log(len(pan.families)) 
				else:# meaning Bernoulli model
					BIC = 2 * float(elements[3]) - (options.classnumber[0] * len(pan.organisms) + 1 + options.classnumber[0] - 1) * math.log(len(pan.families))
	print(BIC)



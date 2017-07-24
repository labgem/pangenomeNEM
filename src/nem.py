# !/usr/bin/python
# -*- coding: iso-8859-1 -*-

from util import *
from classes import *

import argparse
import os
import sys
import errno
from collections import defaultdict
#import numpy as np
import warnings
import math
#from multiprocessing import Pool
#import itertools
#import shutil
import random
import string
#from joblib import Parallel, delayed
#import multiprocessing
import subprocess

#import json
#import requests
#import re
#from scipy.cluster.hierarchy import linkage
#from scipy.cluster.hierarchy import to_tree
#from sklearn import manifold
#from scipy.spatial import Voronoi, voronoi_plot_2d
#import matplotlib.pyplot as plt
#from scipy.spatial.distance import jaccard
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm

logging.basicConfig(level = logging.DEBUG, format = '\n%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s', datefmt='%H:%M:%S')

MASH_LOCATION = ""#ADD MASH
SUBPATH = "/NEMOUT"

OUTPUTDIR = ""
OUTPUTDIR_EVOL = ""

def run(i,prefix, pan, sub_organisms = None):
	if sub_organisms != None:
		
		pan = pan.sub_pangenome(sub_organisms)
		
	pan.partition(i,prefix)
	print(pan)

if __name__=='__main__':

	parser = argparse.ArgumentParser(description='Calculate and partition gene families of a pangenome from annotated genomes and gene families')
	parser.add_argument('--version', action='version', version='0.1')
		
	#TODO add_argument_group and .add_mutually_exclusive_group()
	group_progenome = parser.add_argument_group("progenome options")
	group_progenome.add_argument('-o', '--organisms', type=argparse.FileType('r'), nargs=1, help="The tsv file provided by progenome containing the gene annotations")
	group_progenome.add_argument('-g', '--gene_families', type=argparse.FileType('r'), nargs=1, help="The tsv file provided by progenome containing eggNOG orthologous groups related to the annotated genomes")

	# parser.add_argument("-p", "--ponderation", default=False, action="store_true", help="use mash to calculate distance between genomes and based on these distances ponderate underepresented genomes via a MDS approach")
	# ponderation = parser.add_mutually_exclusive_group()
	# ponderation.add_argument("-x", "--mash_fasta",type =argparse.FileType('r'),nargs=1, help="fasta file containing contig required to computed distance between genome")
	# ponderation.add_argument("-f", "--distances_file", type=argparse.FileType('r'), nargs=1, help="a csv file containing a matrix of distances between organisms")
	# ponderation.add_argument("-j", "--jacquard_distance", action="store_true", help="compute jacquard distances between organisms based on presence/absences of famillies among organisms")


	#parser.add_argument("-s", "--remove_singleton", default=False, action="store_true", help="Remove singleton to the pan-genome")
	#parser.add_argument("-m", "--remove_ME", default=False, action="store_true", help="Remove the mobile elements (transposons, integrons, prophage gene) of the pan-genome")
	parser.add_argument('-d', '--outputdirectory', type=str, nargs=1, default="output.dir", help="The output directory")
	parser.add_argument("-n", "--neighborcomputation", type=int, default=1, help="Consider neighboors for the analysis with the max neighbor distance (integer) (0 = infinite distance)")
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
	ONTOLOGY_FILE              = OUTPUTDIR+"/"+"ontology.txt"
	COG_FILE                   = OUTPUTDIR+"/"+"COG.txt"

	max_neighbordistance = float("inf") if int(options.neighborcomputation)==0 else int(options.neighborcomputation)
	# num_thread = multiprocessing.cpu_count() if options.num_thread==0 else int(options.num_thread)
	# num_thread = 1 if num_thread==0 else num_thread 
	
	distances = None
	
	pan = Pangenome("file", options.organisms[0],  options.gene_families[0])#, options.remove_singleton

	pan.partition(NEMOUTPUTDIR+"/nb"+str(pan.nb_organisms)+"_k"+str(options.classnumber[0])+"i_0", options.classnumber[0], write_graph = "gexf")	

	print(pan)

	print("> Creation of a exact_core_cluster.txt file with all core clusters")
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

	# if options.ponderation:
	# 	if options.mash_fasta:
	# 		distances = calc_mash_distance(options.mash_fasta[0],OUTPUTDIR,num_thread)
	# 	elif options.jacquard_distance:
	# 		distances = pd.DataFrame(0.0, index = pan.organism_positions.keys(), columns = pan.organism_positions.keys())		
	# 		for i in range(0,pan.nb_organisms):
	# 			i_vector = [1 if fam_vector[i] > 0 else 0 for fam_vector in Pangenome.presences_absences.values()]
	# 			for j in range(0,i):
	# 				j_vector = [1 if fam_vector[j] > 0 else 0 for fam_vector in Pangenome.presences_absences.values()]
	# 				dis = jaccard(i_vector, j_vector)
	# 				distances.iloc[i,j] = dis
	# 				distances.iloc[j,i] = dis
	# 		distances.to_csv(OUTPUTDIR+"/jacquard_distance.csv")
	# 	elif options.distances_file:
	# 		distances = pd.read_csv(options.distances_file[0], sep="\t", index_col = 0)

	# 	print(distances)

	# 	distances = distances.round(decimals=3)
	# 	triangle = np.triu(distances.values)
	# 	step = np.min(triangle[np.nonzero(triangle)])/10
	# 	step = 0.001 if step > 0.001 else step

	# 	mds     = manifold.MDS(n_components=pan.nb_organisms-1, dissimilarity="precomputed", random_state=1000, max_iter = 1000, n_init = 1000,n_jobs = num_thread)
	# 	results = mds.fit(distances.values)
	# 	print(results.stress_)
	# 	# xmin=np.min(results.embedding_[:,0])*1.5
	# 	# xmax=np.max(results.embedding_[:,0])*1.5
	# 	# ymin=np.min(results.embedding_[:,1])*1.5
	# 	# ymax=np.max(results.embedding_[:,1])*1.5
	# 	coords  = [tuple(coord) for coord in results.embedding_] #already centered around the centroid which is the origin of the space
	# 	logging.getLogger().debug(step)
	# 	logging.getLogger().debug(coords)
		
	# 	cpt=0
	# 	# cmap = cm.get_cmap(name='rainbow')
	# 	# centroid_x = sum([coord[0] for coord in coords])/len(coords)
	# 	# centroid_y = sum([coord[1] for coord in coords])/len(coords)
	# 	# centroid = (centroid_x,centroid_y,0) # 0 is the radius of a point

	# 	adresses = dict((i,set([i])) for i in range(0,len(coords)))
	# 	weights = [0] * len(coords)
	# 	while len(coords) > 1:
	# 		new_coords = []
	# 		fixed = set()
	# 		new_adresses = defaultdict(set)
	# 		new_weights = defaultdict(float)
	# 		while len(fixed)!=len(coords): 
	# 			# fig, ax = plt.subplots()
	# 			# ax.plot(centroid[0],centroid[1], 'ro')
	# 			# ax.text(centroid[0],centroid[1], s = "c")
	# 			# ax.axis((xmin,xmax,ymin,ymax))
	# 			for i in range(0,len(coords)):
	# 				if i not in fixed:
	# 					intersection = False
	# 					for j in range(0,len(coords)):
	# 						if j != i:
	# 							new_coord = hypersphere_intersection(coords[i],new_weights[i],coords[j],new_weights[j])
	# 							if type(new_coord) == tuple:
	# 								intersection = True
	# 								fixed.add(i)
	# 								fixed.add(j)
	# 								new_adresses[len(new_coords)].update(adresses[i],adresses[j])
	# 								new_coords.append(new_coord)
	# 								break

	# 					if not intersection:
	# 						new_weights[i] += step

	# 			# 	ax.add_artist(plt.Circle((coords[i][0],coords[i][1]), radius=new_weights[i], color=cmap(i), fill=True, alpha = 0.3))
	# 			# 	ax.text(coords[i][0],coords[i][1], s = "\n".join([str(distances.index[adress]) for adress in adresses[i]]))

	# 			# fig.savefig('circles_step'+str(cpt)+".png")
	# 			# cpt+=1
	# 		#weights are propagated to the root
	# 		for i, adress in adresses.items():
	# 			for elem in adress:
	# 				weights[elem] += len(coords)*new_weights[i]
	# 		logging.getLogger().debug(dict(zip(distances.index, weights)))
	# 		coords   = new_coords  
	# 		adresses = new_adresses

	# 	sum_weigths = sum(weights)
	# 	weights = [round(wei/sum_weigths,3) for wei in weights]
	# 	weights = dict(zip(distances.index, weights))

	# 	with open(OUTPUTDIR+"/weights.txt","w") as weight_file:
	# 		for org, wei in sorted(weights.items()):
	# 			weight_file.write(org+"\t"+str(wei)+"\n")

	# 	logging.getLogger().debug(weights)	
	# 	pan.assign_weights(weights)

	# if options.evolution:
	# 	arguments             = list()       
	# 	organisms             = list(pan.organisms)	
	# 	total_combinations    = oidsCombinations(range(0,len(organisms)),options.min_resampling,options.max_resampling,options.min_resampling)
	# 	del total_combinations[len(organisms)]
	# 	nb_total_combinations = 0
	# 	for c in total_combinations:
	# 		nb_total_combinations += len(total_combinations[c])
	# 	print("..... Preparing to classify 1 pangenome + "+str(nb_total_combinations)+" subsampled pangenomes .....")
	# 	for comb_size in total_combinations:
	# 		for combination in total_combinations[comb_size]:
	# 			arguments.append([NEMOUTPUTDIR, pan, [organisms[i] for i in combination]])
		
	# 	random.shuffle(arguments)
	# 	arguments.insert(0, [NEMOUTPUTDIR, pan])
	# 	Parallel(n_jobs=num_thread)(delayed(run)(i,*arguments[i]) for i in range(len(arguments)))
	# else:


	# command = 'cat '+NEMOUTPUTDIR+"/*/nem.stat | sort -k1n > "+EVOLUTION_STATS_NEM_FILE
	# print(command)
	# proc    = subprocess.Popen(command, shell=True)
	# proc.communicate()
	# command = 'cat '+NEMOUTPUTDIR+"/*/exact.stat | sort -k1n > "+EVOLUTION_STATS_EXACT_FILE
	# print(command)
	# proc    = subprocess.Popen(command, shell=True)
	# proc.communicate()
	# ortho_2_COG_funcat = findCOG(pan)
	# ortho_2_COG_funcat.to_csv(COG_FILE)
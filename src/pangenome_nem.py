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

import json
# from goatools.obo_parser import GODag
# from goatools.base import download_go_basic_obo
# from goatools import obo_parser
# from goatools import semantic
import requests
import re
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import to_tree
from sklearn import manifold
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
from scipy.spatial.distance import jaccard
import matplotlib.pyplot as plt
import matplotlib.cm as cm


logging.basicConfig(level = logging.DEBUG, format = '\n%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s', datefmt='%H:%M:%S')


MASH_LOCATION = ""#ADD MASH
SUBPATH = "/NEMOUT"

OUTPUTDIR = ""
OUTPUTDIR_EVOL = ""

def run(i,prefix, pan, sub_organisms = None):
	if sub_organisms != None:
		
		pan = pan.sub_pangenome(sub_organisms)
		
	pan.classify(i,prefix)
	print(pan)

# def findGO(pangenome):

# 	obo_fname = download_go_basic_obo()
# 	go_parser = obo_parser.GODag("go-basic.obo")

# 	API= "http://eggnogapi.embl.de/nog_data/json/go_terms/"
	
# 	GO_roots = { "Biological Process": set(['GO:0008150']), "Cellular Component" : set(['GO:0005575']), "Molecular Function" : set(['GO:0003674'])}

# 	ortho_2_GO = pd.DataFrame("S", index = sorted(pangenome.families), columns = GO_roots.keys())
# 	for ortho in sorted(pangenome.families):
# 		link = API+ortho
# 		try:
# 			content = urllib2.urlopen(link)
# 		except urllib2.URLError:
# 			continue
# 		if content.getcode() == 200:
# 			go_terms = json.loads(content.read())["go_terms"]
# 			for namespace in GO_roots.keys():
# 				if namespace in go_terms.keys():

# 					go_terms_name_space = go_terms[namespace]
# 					result=semantic.common_parent_go_ids([go_id[0] for go_id in go_terms_name_space],go_parser)
# 					result -= GO_roots[namespace]
# 					ortho_2_GO.loc[ortho,namespace]=go_parser[result.pop()].name if len(result) > 0 else "Unknown"
# 		print(ortho+"  "+str(ortho_2_GO.loc[ortho,:]))
# 	return(ortho_2_GO)


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
	ponderation.add_argument("-x", "--mash_fasta",type =argparse.FileType('r'),nargs=1, help="fasta file containing contig required to computed distance between genome")
	ponderation.add_argument("-f", "--distances_file", type=argparse.FileType('r'), nargs=1, help="a csv file containing a matrix of distances between organisms")
	ponderation.add_argument("-j", "--jacquard_distance", action="store_true", help="compute jacquard distances between organisms based on presence/absences of famillies among organisms")


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
			logging.getLogger().error(directory+" already exist")
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
	
	if options.use[0] == "progenome":
		pan = Pangenome(options.use[0], options.annotations[0],  options.eggNOG_clusters[0], options.remove_singleton)
	else:
		exit(1)

	print(pan)
	
	outfile = open(ORGANISMS_FILE,"w")
	outfile.writelines(["%s\n" % item  for item in (pan.organism_positions.keys())])
	outfile.close()

	core_cluster_file = open(EXACT_CORE_FILE,"w")
	core_cluster_file.write("\n".join(pan.core_list))
	core_cluster_file.write("\n")
	core_cluster_file.close()
	print("> Creation of a exact_core_cluster.txt file with all core clusters")

	if options.ponderation:
		if options.mash_fasta:
			distances = calc_mash_distance(options.mash_fasta[0],OUTPUTDIR,num_thread)
		elif options.jacquard_distance:
			distances = pd.DataFrame(0.0, index = pan.organism_positions.keys(), columns = pan.organism_positions.keys())		
			for i in range(0,pan.nb_organisms):
				i_vector = [1 if fam_vector[i] > 0 else 0 for fam_vector in Pangenome.presences_absences.values()]
				for j in range(0,i):
					j_vector = [1 if fam_vector[j] > 0 else 0 for fam_vector in Pangenome.presences_absences.values()]
					dis = jaccard(i_vector, j_vector)
					distances.iloc[i,j] = dis
					distances.iloc[j,i] = dis
			distances.to_csv(OUTPUTDIR+"/jacquard_distance.csv")
		elif options.distances_file:
			distances = pd.read_csv(options.distances_file[0], sep="\t", index_col = 0)

		print(distances)
		# cutoff = 0.025

		# adj_graph = nx.from_pandas_dataframe(distance_melted, "source", "target", 'weight')
		# logging.getLogger().debug(adj_graph)
		# adj_subgraph = nx.Graph( [(u,v,d) for u,v,d in adj_graph.edges(data=True) if d['weight']>cutoff])
		# logging.getLogger().debug(adj_subgraph)

		#same = dict()

		# identitical = dict()

		# distance_melted = pd.DataFrame.stack(distances, level=0).reset_index()
		# distance_melted.columns = ["x","y","distance"]
		# logging.getLogger().debug(len(distance_melted.index))

		# same = defaultdict(set)
		# insame = []
		# distance_melted.sort_values(["x","y","distance"], inplace=True)
		# drop = []
		# for i, row in distance_melted.iterrows():
		# 	if row['y'] != row['x'] and row['distance'] == 0:
		# 		same[row['y']].add(row['x'])
		# 		insame.append(row['x'])
		# 	if (row['y'] in same) or (row['x'] in same):
		# 		drop.append(i)
		# print(same)
		# exit()
		# distance_melted.drop(distance_melted.index[drop],inplace=True)
		# logging.getLogger().debug(len(distance_melted.index))
		# #distance_melted = distance_melted.loc[lambda row: row.distance > 0.1]
		# logging.getLogger().debug(len(distance_melted.index))
		# #TODO  add identical organism

		# size = len(distance_melted["x"])
		
		# logging.getLogger().debug(len(set(distance_melted["x"])))
		# logging.getLogger().debug(len(set(distance_melted["y"])))
		# condensed_distances = pd.DataFrame(0, index=set(distance_melted["x"]),columns=set(distance_melted["x"]))
		# logging.getLogger().debug(condensed_distances)
		# for i, row in distance_melted.iterrows():
		# 	logging.getLogger().debug(row['x'])
		# 	logging.getLogger().debug(row['y'])
		# 	logging.getLogger().debug(row['distance'])
		# 	condensed_distances.loc[row['x'],row['y']] = row['distance']
		# 	condensed_distances.loc[row['y'],row['x']] = row['distance']

		# logging.getLogger().debug(condensed_distances)
		# #while distances.values
		# # case of 0 distance

		distances = distances.round(decimals=3)
		triangle = np.triu(distances.values)
		step = np.min(triangle[np.nonzero(triangle)])
		step = 0.001 if step > 0.001 else step

		mds     = manifold.MDS(n_components=2, dissimilarity="precomputed", random_state=10)
		results = mds.fit(distances.values)
		xmin=np.min(results.embedding_[:,0])*1.5
		xmax=np.max(results.embedding_[:,0])*1.5
		ymin=np.min(results.embedding_[:,1])*1.5
		ymax=np.max(results.embedding_[:,1])*1.5
		coords  = [(coord[0],coord[1]) for coord in results.embedding_] #already centered around the centroid which is the origin of the space
		#weights = dict(zip(pan.organism_positions.keys(), [0] * len(pan.organism_positions)))
		logging.getLogger().debug(step)
		logging.getLogger().debug(coords)
		
		cpt=0
		cmap = cm.get_cmap(name='rainbow')
		centroid_x = sum([coord[0] for coord in coords])/len(coords)
		centroid_y = sum([coord[1] for coord in coords])/len(coords)
		centroid = (centroid_x,centroid_y,0) # 0 is the radius of a point

		adresses = dict((i,set([i])) for i in range(0,len(coords)))
		weights = [0] * len(coords)
		while len(coords) > 1:
			new_coords = []
			fixed = set()
			
			new_adresses = defaultdict(set)
			new_weights = defaultdict(float)
			while len(fixed)!=len(coords): 
				fig, ax = plt.subplots()
				ax.plot(centroid[0],centroid[1], 'ro')
				ax.text(centroid[0],centroid[1], s = "c")
				ax.axis((xmin,xmax,ymin,ymax))
				for i in range(0,len(coords)):
					if i not in fixed:
						intersection = False
						for j in range(0,len(coords)):
							if j != i:
								new_coord = circle_intersection((coords[i][0],coords[i][1],new_weights[i]),(coords[j][0],coords[j][1],new_weights[j]))
								if type(new_coord) == tuple:
									intersection = True
									fixed.add(i)
									fixed.add(j)
									new_adresses[len(new_coords)].update(adresses[i],adresses[j])
									new_coords.append(new_coord)
									break

						if not intersection:
							if i == 52:
								print(new_weights[i])
							new_weights[i] += step

					ax.add_artist(plt.Circle((coords[i][0],coords[i][1]), radius=new_weights[i], color=cmap(i), fill=True, alpha = 0.3))
					ax.text(coords[i][0],coords[i][1], s = "\n".join([str(distances.index[adress]) for adress in adresses[i]]))

				fig.savefig('circles_step'+str(cpt)+".png")
				cpt+=1
			#weights are propagated to the root
			for i, adress in adresses.items():
				for elem in adress:
					weights[elem] += len(coords)*new_weights[i]
			logging.getLogger().debug(dict(zip(distances.index, weights)))
			coords   = new_coords  
			adresses = new_adresses
		weights = dict(zip(distances.index, weights))
		with open(OUTPUTDIR+"/weights.txt","w") as weight_file:
			for org, wei in weights.items():
				weight_file.write(org+"\t"+str(wei)+"\n")

		logging.getLogger().debug(weights)	
		# vor         = Voronoi(coords)
		# lines       = [LineString(vor.vertices[line]) for line in vor.ridge_vertices if -1 not in line]
		# convex_hull = MultiPoint([Point(i) for i in coords]).convex_hull.buffer(2)
		# pts         = MultiPoint([Point(i) for i in coords])
		# mask        = pts.convex_hull.union(pts.buffer(0.001, resolution=50, cap_style=1))
		# polys       = MultiPolygon([poly.intersection(mask) for poly in polygonize(lines)])
		# for poly in polys:
		# 	plt.fill(*zip(*np.array(list(zip(poly.boundary.coords.xy[0][:-1], poly.boundary.coords.xy[1][:-1])))), alpha=0.4, color = colors.rgb2hex(np.random.rand(3)))
		# plt.plot(coords[:,0], coords[:,1], 'ko')
		# plt.show()
		# logging.getLogger().debug(polys)
		# regions, vertices = voronoi_finite_polygons_2d(vor_result)
		# pts = MultiPoint([Point(i) for i in coords])
		# mask = pts.convex_hull.union(pts.buffer(10, resolution=5, cap_style=3))
		# new_vertices = []
		# for region in regions:
		# 	polygon = vertices[region]
		# 	shape = list(polygon.shape)
		# 	shape[0] += 1
		# 	p = MultiPolygon([poly.intersection(mask) for poly in polygonize(lines)])
		# 	#p = Polygon(np.append(polygon, polygon[0]).reshape(*shape)).intersection(mask)
		# 	poly = np.array(list(zip(p.boundary.coords.xy[0][:-1], p.boundary.coords.xy[1][:-1])))
		# 	new_vertices.append(poly)
		# 	plt.fill(*zip(*poly), alpha=0.4)
		# plt.plot(coords[:,0], coords[:,1], 'ko')
		# plt.title("Clipped Voronois")
		# plt.show()

  #       plt.show()
		#adj_graph = nx.relabel_nodes(adj_graph, distances.columns.str)
		

		#hclust = scipy.cluster.hierarchy.linkage(distances.values, method='single', metric='euclidean')
		#tree   = scipy.cluster.hierarchy.to_tree(hclust, rd=True)
		#dendro = scipy.cluster.hierarchy.dendrogram(hclust)

	if options.evolution:
		arguments = list()       
		organisms = list(pan.organisms)	
		total_combinations = oidsCombinations(range(0,len(organisms)),10,1000,10)
		del total_combinations[len(organisms)]

		nb_total_combinations = 0	

		for c in total_combinations:
			nb_total_combinations += len(total_combinations[c])
		print("..... Preparing to classify 1 pangenome + "+str(nb_total_combinations)+" subsampled pangenomes .....")
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



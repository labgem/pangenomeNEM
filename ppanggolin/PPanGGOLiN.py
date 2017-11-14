#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

import argparse
from collections import defaultdict
from collections import OrderedDict
from collections import Counter
from ordered_set import OrderedSet
import networkx as nx
import os
import sys
import re
import math
import logging
import random 
import string
import shutil
import gzip
import operator
import numpy as np
import time
import community
import tempfile
from random import randrange
from random import shuffle
from multiprocessing import Pool
import subprocess
import progressbar
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from tqdm import tqdm
import mmap

#import forceatlas2 

NEM_LOCATION  = os.path.dirname(os.path.abspath(__file__))+"/../NEM/nem_exe"
(TYPE, FAMILY, START, END, STRAND, NAME, PRODUCT) = range(0, 7)#data index in annotation
(ORGANISM_ID, ORGANISM_GFF_FILE) = range(0, 2)#data index in the file listing organisms 

(GFF_seqname, GFF_source, GFF_feature, GFF_start, GFF_end, GFF_score, GFF_strand, GFF_frame, GFF_attribute) = range(0,9) 

RESERVED_WORDS = set(["id", "label", "name", "weight", "partition", "partition_exact", "length", "length_min", "length_max", "length_avg", "length_avg", "product", 'nb_gene', 'community'])

"""  
    :mod:`pangenome` -- Depict microbial diversity
===================================

.. module:: pangenome
   :platform: Unix
   :synopsis: Depict microbial diversity via a partionned pangenome graph.
    .. moduleauthor:: Guillaume GAUTREAU (LABGeM, Genoscope, France) ggautrea@genoscope.cns.fr

    Description
    -------------------
    Pangenomes are generally stored in a binary matrix denoting the presence or absence of each gene family across organisms. However, this structure does not handle the genomic organization of gene families in each organism. We propose a graph model where nodes represent families and edges chromosomal neighborhood information. Indeed, it is known that core gene families share conserved organizations whereas variable regions are rather randomly distributed along genomes.
Moreover, our method classifies gene families through an Expectation/Maximization algorithm based on Bernoulli mixture model. This approach splits pangenomes in three groups: (1) persistent genome, equivalent to a relaxed core genome (genes conserved in all but a few genomes); (2) shell genome, genes having intermediate frequencies corresponding to moderately conserved genes potentially associated to environmental adaptation capabilities; (3) cloud genome, genes found at very low frequency.

""" 

class PPanGGOLiN:

    """  
          The ``PPanGGOLiN`` class
        ======================
        .. class:: PPanGGOLiN

            Pangenomes are generally stored in a binary matrix denoting the presence or absence of each gene family across organisms. 
            However, this structure does not handle the genomic organization of gene families in each organism. 
            We propose a graph model where nodes represent families and edges chromosomal neighborhood information. 
            Indeed, it is known that core gene families share conserved organizations whereas variable regions are rather randomly distributed along genomes.
            The PPanGGOLiN class models the genemic diversity of a pangenome, this modelisation organize the genemic diveristy via a pangenome graph of chromosomal neigborhood.
            Moreover, our method classifies gene families through an Expectation/Maximization algorithm based on Bernoulli mixture model. 
            This approach splits pangenomes in three groups: 
                1. *persistent genome*, equivalent to a relaxed core genome (genes conserved in all but a few genomes); 
                2. *shell genome*, genes having intermediate frequencies corresponding to moderately conserved genes potentially associated to environmental adaptation capabilities; 
                3. *cloud genome*, genes found at very low frequency.

            .. attribute:: annotations

                multilevel dictionnaries containing a dictionary of contig for each organism, and a dictionary of lists containing annotations for each contig

            .. attribute:: neighbors_graph

                a networkx undirected graph. Node correspond to gene families and edges to chromosomal colocalization beween families. Organisms supporting each edge are stored in edge attribute as weel as the edge weight (number of organism coverinf each edge).
                Nodes attributes contains the gene identifiers of each organism supporting this node.

            .. attribute:: organisms

                an ordored-set contains the imported organisms 

            .. attribute:: nb_organisms

                a int giving the number of imported organisms 

            .. attribute:: circular_contig_size

                a dict containing the size of the contigs (contigs identifiers as keys) which are both well assembled and circular contigs (used to circularize the graph). The contigs wich are not circular are not in this dictionnaries.

            .. attribute:: families_repeted_th
   
                a int containing the threshold of the maximum number copy of each families. Families exceeding this threshold are removed and are listed in the next attribute.

            .. attribute:: families_repeted

                a set containing the family identifiers of ones having a maximun number of copy in at least one organism above the families_repeted_th threshold attribute.

            .. attribute:: pan_size

                The number of nodes into the graph.
                .. warning:: this number is not necesserally equal to the number of imported gene families. Indeed, if imported families are not present into organism, superfluous will be discard.

            .. attribute:: is_partitionned
            
                a boolean specifying if the pangenome graph has been partitionned or not

            .. attribute:: nem_intermediate_files

                a str storing the path to the nem nem intermediate_files

            .. attribute:: partitions

                a dict provoding the families present in each partion:
                    * partitions["core_exact"] contains the list of core families (core exact)
                    * partitions["accessory"] contains the list of families not in the core exact
                    * partitions["persistent"] contains the list of persistent families
                    * partitions["shell"] contains the list of shell families
                    * partitions["cloud"] contains the list of cloud families

            .. attribute:: BIC

                a float providing the Bayesian Information Criterion. This Criterion give an estimation of the quality of the partionning (a low value means a good one)
                . seealso:: https://en.wikipedia.org/wiki/Bayesian_information_criterion

            .. method:: partition(nem_dir_path)

                se

            .. method:: export_to_GEXF()

                bllla

            .. method:: import_from_GEXF()

                blbaal

            .. method:: neighborhood_computation()

                blbaal

            .. method:: delete_pangenome_graph()

                blbaal

            .. method:: delete_nem_intermediate_files()

                blbaal
    """ 
    def __init__(self, init_from = "args", *args):
        """ 
            :param init_from: specified the excepted input (can be "file", "args", "database")
            :param *args: depending on the previous paramter, args can take multiple forms
            :type init_from: str
            :type *args: list

            :Example:

            >>>pan = PPanGGOLiN("file", organisms, gene_families, remove_high_copy_number_families)
            >>>pan = PPanGGOLiN("args", annotations, organisms, circular_contig_size, families_repeted)# load direclty the main attributes
        """ 
        self.annotations              = dict()
        self.neighbors_graph          = None
        self.organisms                = OrderedSet()
        self.nb_organisms             = 0
        self.circular_contig_size     = dict()
        self.families_repeted_th      = 0
        self.families_repeted         = set()
        self.pan_size                 = 0
        self.is_partitionned          = False
        self.nem_intermediate_files   = None
        self.partitions               = {}
        self.partitions["undefined"]  = list()
        self.partitions["persistent"] = list()
        self.partitions["shell"]      = list()
        self.partitions["cloud"]      = list()
        self.partitions["core_exact"] = list()
        self.partitions["accessory"]  = list()
        self.BIC                      = None # Bayesian Index Criterion
        #self.partitions_by_organisms  = defaultdict(lambda: defaultdict(set))

        if init_from == "file":
            self.__initialize_from_files(*args)
        elif init_from == "args":
            (self.annotations,
             self.organisms,
             self.circular_contig_size,
             self.families_repeted) = args 
        elif init_from == "database":
            logging.getLogger().error("database is not yet implemented")
            pass
        else:
            raise ValueError("init_from parameter is required")
        self.nb_organisms = len(self.organisms)

    def __initialize_from_files(self, organisms_file, families_tsv_file, lim_occurence = 0, infere_singleton = False, already_sorted = False):
        """ 
            :param organisms_file: a file listing organims by compute, first column is organism name, second is path to gff file and optionnally other other to provide the name of circular contig
            :param families_tsv_file: a file listing families. The first element is the family identifier (by convention, we advice to use the identifier of the average gene of the family) and then the next elements are the identifiers of the genes belonging to this family.
            :param lim_occurence: a int containing the threshold of the maximum number copy of each families. Families exceeding this threshold are removed and are listed in the next attribute.
            :param infere_singleton: a bool specifying if singleton must be explicitely present in the families_tsv_file (False) or if single gene in gff files must be automatically infered as a singleton family (True)
            :type file: 
            :type file: 
            :type int: 
            :type bool: 
        """ 
        logging.getLogger().info("Reading "+families_tsv_file.name+" families file ...")
        families    = dict()
        first_iter  = True
        for line in families_tsv_file:
            elements = line.split()
            for gene in elements[1:]:
                families[gene]=elements[0]

        self.circular_contig_size = {}

        logging.getLogger().info("Reading "+organisms_file.name+" list of organism files ...")

        def get_num_lines(file_path):
            fp = open(file_path, "r+")
            buf = mmap.mmap(fp.fileno(), 0)
            lines = 0
            while buf.readline():
                lines += 1
            return lines

        bar = tqdm(organisms_file,total=get_num_lines(organisms_file.name), unit = "gff file", unit_scale = True)

        for line in bar:
            elements = line.split()
            bar.set_description("Processing "+elements[ORGANISM_GFF_FILE])
            bar.refresh()
            if len(elements)>2:
                self.circular_contig_size.update({contig_id: None for contig_id in elements[2:len(elements)]})# size of the circular contig is initialized to None (waiting to read the gff files to fill the dictionnaries with the correct values)
            self.annotations[elements[0]] = self.__load_gff(elements[ORGANISM_GFF_FILE], families, elements[ORGANISM_ID], lim_occurence, infere_singleton, already_sorted)

        check_circular_contigs = {contig: size for contig, size in self.circular_contig_size.items() if size == None }
        if len(check_circular_contigs) > 0:
            logging.getLogger().error("""
                The following identifiers of circular contigs in the file listing organisms have not been found in any region feature of the gff files: """+"\t".join(check_circular_contigs.keys()))
            exit()
    def __load_gff(self, gff_file_path, families, organism, lim_occurence = 0, infere_singleton = False, already_sorted = False):
        """
            Load the content of a gff file
            :param gff_file_path: a valid gff file path where only feature of the type 'CDS' will be imported as genes. Each 'CDS' feature must have a uniq ID as attribute (afterall called gene id).
            :param families: a dictionary having the gene as key and the identifier of the associated family as value. Depending on the infere_singleton attribute, singleton must be explicetly present on the dictionary or not
            :param organism: a str containing the organim name
            :param lim_occurence: a int containing the threshold of the maximum number copy of each families. Families exceeding this threshold are removed and are listed in the next attribute.
            :param lim_occurence: a bool specifying if singleton must be explicitely present in the families parameter (False) or if single gene automatically infered as a singleton family (True)
            :type str: 
            :type dict: 
            :type str: 
            :type int: 
            :type bool: 
            :return: annot: 
            :rtype: dict 
        """ 

        logging.getLogger().debug("Reading "+gff_file_path+" file ...")

        if organism not in self.organisms:
            self.organisms.add(organism)
            annot = defaultdict(OrderedDict)

            ctp_prev = 1
            cpt_fam_occ = defaultdict(int)

            gene_id_auto = False

            with open(gff_file_path,'r') as gff_file:
                for line in gff_file:
                    
                    if line.startswith('##',0,2):

                        if line.startswith('FASTA',2,7):
                            break
                        elif line.startswith('sequence-region',2,17):
                            fields = line.split(' ')
                            if fields[1] in self.circular_contig_size:
                                self.circular_contig_size[fields[1]] = int(3)
                        continue
                    gff_fields = line.split('\t')
                    if GFF_feature == 'region':
                        if GFF_seqname in self.circular_contig_size:
                            self.circular_contig_size = int(GFF_end)
                            continue
                    elif gff_fields[GFF_feature] == 'CDS':
                        attributes_feild = gff_fields[GFF_attribute].split(';')
                        attributes = {}
                        for att in attributes_feild:
                            (key, value) = att.strip().split('=')
                            attributes[key.upper()]=value
                        try:
                            protein = attributes["ID"]
                        except:
                            logging.getLogger().error("Each CDS feature of the gff files must own a unique ID attribute. Not the case for file: "+gff_file_path)
                            exit(1)
                        try:
                            family = families[protein]
                        except KeyError:
                            if infere_singleton:
                                families[protein] = protein
                                family           = families[protein]
                                logging.getLogger().info("infered singleton: "+protein)
                            else:
                                raise KeyError("Unknown families:"+protein, ", check your families file or run again the program using the option to infere singleton")

                        cpt_fam_occ[family]+=1
                        prev = families[protein]

                        try:
                            name = attributes.pop('NAME')
                        except KeyError:
                            try:
                                name = attributes.pop('GENE')
                            except KeyError:
                                name = ""

                        try:
                            product = attributes.pop('PRODUCT')
                        except KeyError:
                            product = ""

                        annot[gff_fields[GFF_seqname]][protein] = ["CDS",family,int(gff_fields[GFF_start]),int(gff_fields[GFF_end]),gff_fields[GFF_strand], name, product]

            if not already_sorted:
                for seq_id in list(annot):
                    annot[seq_id] = OrderedDict(sorted(annot[seq_id].items(), key= lambda item: item[1][START]))
                    
            if (lim_occurence > 0):
                fam_to_remove =[fam for fam, occ in cpt_fam_occ.items() if occ > lim_occurence]
                logging.getLogger().debug("highly repeted families found (>"+str(lim_occurence)+" in "+organism+"): "+" ".join(fam_to_remove))
                self.families_repeted = self.families_repeted.union(set(fam_to_remove))

            return(annot)
        else:
            raise KeyError("Redondant organism names was found ("+organism+")")

    def __str__(self):
        """ Return an overview of the statistics of the pangenome as a formated string """ 
        pan_str ="\n"
        pan_str += "----------- Statistics -----------\n"
        pan_str += "Number of organisms: "+str(self.nb_organisms)+"\n"

        if self.pan_size != 0:
            pan_str += "Pan-genome size:"+str(self.pan_size)+"\n"
            pan_str += "Exact core-genome size:"+str(len(self.partitions["core_exact"]))+"\n"
            pan_str += "Exact variable-genome size:"+str(self.pan_size-len(self.partitions["core_exact"]))+"\n"
            pan_str += "Persistent genome size:"+str(len(self.partitions["persistent"]))+"\n"
            pan_str += "Shell genome size:"+str(len(self.partitions["shell"]))+"\n"
            pan_str += "Cloud genome cloud:"+str(len(self.partitions["cloud"]))+"\n"
        else:
            pan_str += "No partitioning have been performed on this Pangenome instance\n"
            pan_str += "Run the partitioning function to obtain more detailled statistics...\n"
        pan_str += "----------------------------------"

        return(pan_str)    

    def __add_gene(self, fam_id, org, gene, name, length, product):#, multi_copy = None
        """
            Add gene to the pangenome graph
            :param fam_id: The family identifier
            :param org: The organism name
            :param gene : The gene identifier
            :param name: The biological name of the gene
            :param length: The number of nucleotide of the gene
            :param product: The name of the protein function
            :param length: The number of nucleotide of the gene
            :type str: 
            :type str:
            :type str: 
            :type str: 
            :type str: 
            :type str: 
            :type str: 
        """ 
        self.neighbors_graph.add_node(fam_id)

        try: 
            self.neighbors_graph.node[fam_id]["nb_gene"]+=1
        except KeyError:
            self.neighbors_graph.node[fam_id]["nb_gene"]=1
        try:
            self.neighbors_graph.node[fam_id][org].add(gene)
        except KeyError:
            self.neighbors_graph.node[fam_id][org] = set([gene])

        for attribute in ["name","length","product"]:
            try:
                self.neighbors_graph.node[fam_id][attribute].add(locals()[attribute])
            except KeyError:
                self.neighbors_graph.node[fam_id][attribute]=set([locals()[attribute]])

    def __add_link(self, fam_id, fam_id_nei, org, length):
        """
            Add line between families of a the pangenome graph
            :param fam_id: The family identifier the first node (need to be have at least one gene belonging to this family already added to the graph via the method __add_gene())
            :param fam_id_nei: The family identifier the second node (need to be have at least one gene belonging to this family already added to the graph via the method __add_gene())
            :param org : The identifier of the organism supporting this link
            :param length : The distance in number of base between the genes adding this link
            :type str: 
            :type str:
            :type str: 
        """ 

        if not self.neighbors_graph.has_edge(fam_id,fam_id_nei):
            self.neighbors_graph.add_edge(fam_id, fam_id_nei)
            # logging.getLogger().debug([str(i) for i in [fam_id, fam_id_nei, org]])
        try:
            self.neighbors_graph[fam_id][fam_id_nei][org]+=1
        except KeyError:
            self.neighbors_graph[fam_id][fam_id_nei][org]=1
            try:
                self.neighbors_graph[fam_id][fam_id_nei]["weight"]+=1.0
            except KeyError:
                self.neighbors_graph[fam_id][fam_id_nei]["weight"]=1.0
        try:
            self.neighbors_graph[fam_id][fam_id_nei]["length"].add(length)
        except KeyError:
            self.neighbors_graph[fam_id][fam_id_nei]["length"]=set([length])

    def neighborhood_computation(self, directed = True, light = False):#, untangle_multi_copy_families = False
        """ Use the information already loaded (annotation) to build the pangenome graph
            :param directed: a bool specifying if the graph is directed or undirected
            :param light: a bool specifying is the annotation attribute must be detroyed at each step to save memory
            :type bool: 
        """ 
        if self.neighbors_graph is None:
            self.neighbors_graph = nx.Graph()
        elif self.is_partitionned:
            raise Exception("The pangenome graph is already built and partionned, please use the function delete pangenome graph before build it again")

        for organism in list(self.annotations):
            for contig, contig_annot in self.annotations[organism].items():
                
                try:
                    (gene_start, gene_info_start) = contig_annot.popitem(last=False)
                    while (gene_info_start[FAMILY] in self.families_repeted):
                            (gene_start, gene_info_start) = contig_annot.popitem(last=False)
                except KeyError:
                    continue

                self.__add_gene(gene_info_start[FAMILY],
                                organism,
                                gene_start,
                                gene_info_start[NAME],
                                gene_info_start[END]-gene_info_start[START],
                                gene_info_start[PRODUCT])

                family_id_nei, end_family_nei  = gene_info_start[FAMILY], gene_info_start[END]
                logging.getLogger().debug(gene_info_start)
                for gene, gene_info in contig_annot.items():
                    logging.getLogger().debug(gene_info)
                    logging.getLogger().debug(gene)
                    if gene_info[FAMILY] not in self.families_repeted:

                        self.__add_gene(gene_info[FAMILY],
                                        organism,
                                        gene,
                                        gene_info[NAME],
                                        gene_info[END]-gene_info[START],
                                        gene_info[PRODUCT])
                        self.neighbors_graph.add_node(family_id_nei)
                        self.__add_link(gene_info[FAMILY],family_id_nei,organism, gene_info[START] - end_family_nei)
                        family_id_nei  = gene_info[FAMILY]
                        end_family_nei = gene_info[END]
                
                if contig in self.circular_contig_size:#circularization
                    self.__add_link(gene_info_start[FAMILY],family_id_nei,organism, (self.circular_contig_size[contig] - end_family_nei) + gene_info_start[START])

                contig_annot[gene_start]=gene_info_start

            if light:
                del self.annotations[organism]

        self.pan_size = nx.number_of_nodes(self.neighbors_graph)

    def partition(self, nem_dir_path = tempfile.mkdtemp(), beta = 1.00, free_dispersion = False):
        """
            Use the graph topology and the presence or absence of genes from each organism into families to partition the pangenome in three groups ('persistent', 'shell' and 'cloud')
            . seealso:: Read the Mo Dang's thesis to understand NEM and Bernouilli Mixture Model, a summary is available here : http://www.kybernetika.cz/content/1998/4/393/paper.pdf
            
            :param nem_dir_path: a str containing a path to store tempory file of the NEM program
            :param beta: a float containing the spatial coefficient of smothing of the clustering results using the weighted graph topology (0.00 turn off the spatial clustering)
            :param free_dispersion: a bool specyfing if the dispersion around the centroid vector of each paritition is the same for all the organisms or if the dispersion is free
            :type str: 
            :type float: 
            :type bool: 

            .. warning:: please use the function neighborhoodComputation before
        """ 
        if self.neighbors_graph is None:
            logging.getLogger().error("The neighbors_graph is not built, please use the function neighborhood_computation before")
        elif self.is_partitionned:
            logging.getLogger.error("The pangenome is already partionned")

        if not os.path.exists(nem_dir_path):
            #NEM requires 5 files: nem_file.index, nem_file.str, nem_file.dat, nem_file.m and nem_file.nei
            os.makedirs(nem_dir_path)
        self.nem_intermediate_files = nem_dir_path

        logging.getLogger().info("Writing nem_file.str nem_file.index nem_file.nei nem_file.dat and nem_file.m files")
        str_file   = open(nem_dir_path+"/nem_file.str", "w")
        index_file = open(nem_dir_path+"/nem_file.index", "w")
        org_file   = open(nem_dir_path+"/column_org_file", "w")
        nei_file   = open(nem_dir_path+"/nem_file.nei", "w")
        dat_file   = open(nem_dir_path+"/nem_file.dat", "w")
        m_file     = open(nem_dir_path+"/nem_file.m", "w")

        str_file.write("S\t"+str(self.pan_size)+"\t"+str(self.nb_organisms)+"\n")
        str_file.close()

        nei_file.write("1\n")#to enable weigthed partionning
        
        index = {node: index+1 for index, node in enumerate(self.neighbors_graph.nodes(data=False))}
        org_file.write(" ".join([org for org in self.organisms])+"\n")
        org_file.close()

        for node_name, node_organisms in self.neighbors_graph.nodes(data=True):

            index_file.write(str(index[node_name])+"\t"+str(node_name)+"\n")
            logging.getLogger().debug(node_organisms)
            logging.getLogger().debug(self.organisms)
            dat_file.write("\t".join(["1" if org in node_organisms else "0" for org in self.organisms])+"\n")

            row_fam         = []
            row_dist_score  = []
            neighbor_number = 0
            try:
                for neighbor in set(nx.all_neighbors(self.neighbors_graph, node_name)):
                    #nb_presences = sum([pre_abs for org, pre_abs in self.neighbors_graph[node_name][neighbor].items() if org not in RESERVED_WORDS])
                    #self.neighbors_graph[node_name][neighbor]["weight"]= nb_presences
                    w_1=0
                    w_2=0
                    coverage = 0
                    if self.neighbors_graph.is_directed():
                        try:
                            cov_sens = self.neighbors_graph[node_name][neighbor]["weight"]
                        except KeyError:
                            pass
                        try:
                            cov_antisens = self.neighbors_graph[neighbor][node_name]["weight"]
                        except KeyError:
                            pass
                        coverage = cov_sens + cov_antisens
                        
                    else:
                        distance_score = coverage/self.nb_organisms
                    row_fam.append(str(index[neighbor]))
                    row_dist_score.append(str(round(distance_score,4)))
                    neighbor_number += 1
                if neighbor_number>0:
                    nei_file.write("\t".join([str(item) for sublist in [[index[node_name]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
                    #logging.getLogger().debug("\t".join([str(item) for sublist in [[[index[node_name]]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
                else:
                    raise nx.exception.NetworkXError("no neighbors in selected organismss")
            except nx.exception.NetworkXError as nxe:
                logging.getLogger().debug("The family: "+node_name+" is an isolated family")
                nei_file.write(str(index[node_name])+"\t0\n")

        m_file.write("1 0.33333 0.33333 ") # 1 to initialize parameter, 0.333 and 0.333 for to give one third of initial portition to each class (last 0.33 is automaticaly determined by substraction)
        m_file.write(" ".join(["1"]*self.nb_organisms)+" ") # persistent binary vector
        m_file.write(" ".join(["1"]*self.nb_organisms)+" ") # shell binary vector (1 ou 0, whatever because dispersion will be of 0.5)
        m_file.write(" ".join(["0"]*self.nb_organisms)+" ") # cloud binary vector
        m_file.write(" ".join(["0.1"]*self.nb_organisms)+" ") # persistent dispersition vector (low)
        m_file.write(" ".join(["0.5"]*self.nb_organisms)+" ") # shell dispersition vector (high)
        m_file.write(" ".join(["0.1"]*self.nb_organisms)) # cloud dispersition vector (low)

        index_file.close()
        nei_file.close()
        dat_file.close()
        m_file.close()

        logging.getLogger().info("Running NEM...")
        # weighted_degree = sum(list(self.neighbors_graph.degree(weight="weight")).values())/nx.number_of_edges(self.neighbors_graph)
        # logging.getLogger().debug("weighted_degree: "+str(weighted_degree))
        # logging.getLogger().debug("org/weighted_degree: "+str(self.nb_organisms/weighted_degree))    
        #weighted_degree = sum(self.neighbors_graph.degree(weight="weight").values())/nx.number_of_edges(self.neighbors_graph)

        Q              = 3 # number of partitions
        ALGO           = "ncem" #fuzzy classification by mean field approximation
        ITERMAX        = 100 # number of iteration max 
        MODEL          = "bern" # multivariate Bernoulli mixture model
        PROPORTION     = "pk" #equal proportion :  "p_"     varying proportion : "pk"
        VARIANCE_MODEL = "skd" if free_dispersion else "sk_"#one variance per partition and organism : "sdk"      one variance per partition, same in all organisms : "sd_"   one variance per organism, same in all partion : "s_d"    same variance in organisms and partitions : "s__" 
        NEIGHBOUR_SPEC = "f"# "f" specify to use all neighbors, orther argument is "4" to specify to use only the 4 neighbors with the higher weight (4 because for historic reason due to the 4 pixel neighbors of each pixel)
        CONVERGENCE_TH = "clas "+str(0.000001)

        HEURISTIC      = "heu_d"# "psgrad" = pseudo-likelihood gradient ascent, "heu_d" = heuristic using drop of fuzzy within cluster inertia, "heu_l" = heuristic using drop of mixture likelihood
        STEP_HEURISTIC = 0.5 # step of beta increase
        BETA_MAX       = float(self.nb_organisms) #maximal value of beta to test,
        DDROP          = 0.8 #threshold of allowed D drop (higher = less detection)
        DLOSS          = 0.5 #threshold of allowed D loss (higher = less detection)
        LLOSS          = 0.02 #threshold of allowed L loss (higher = less detection)
        
        BETA           = ["-B",HEURISTIC,"-H",str(STEP_HEURISTIC),str(BETA_MAX),str(DDROP),str(DLOSS),str(LLOSS)] if beta == float('Inf') else ["-b "+str(beta)]

        command = " ".join([NEM_LOCATION, 
                            nem_dir_path+"/nem_file",
                            str(Q),
                            "-a", ALGO,
                            "-i", str(ITERMAX),
                            "-m", MODEL, PROPORTION, VARIANCE_MODEL,
                            "-s m "+ nem_dir_path+"/nem_file.m",
                            *BETA,
                            "-n", NEIGHBOUR_SPEC,
                            "-c", CONVERGENCE_TH,
                            "-f fuzzy",
                            "-l y"])
     
        logging.getLogger().info(command)

        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out,err) = proc.communicate()
        logging.getLogger().debug(out)
        logging.getLogger().debug(err)

        if os.path.isfile(nem_dir_path+"/nem_file.uf"):
            logging.getLogger().info("Reading NEM results...")
        else:
            logging.getLogger().error("No NEM output file found: "+ nem_dir_path+"/nem_file.uf")

        classification = ["undefined"] * self.pan_size
        try:
            with open(nem_dir_path+"/nem_file.uf","r") as classification_nem_file, open(nem_dir_path+"/nem_file.mf","r") as parameter_nem_file:

                sum_mu_k = []
                sum_epsilon_k = []
                proportion = []

                parameter = parameter_nem_file.readlines()
                M = float(parameter[6].split()[3]) # M is markov ps-like
                self.BIC = -2 * M - (Q * self.nb_organisms * 2 + Q - 1) * math.log(self.pan_size)
                logging.getLogger().info("The Bayesian Criterion Index of the partionning is "+str(self.BIC))

                for k, line in enumerate(parameter[-3:]):
                    logging.getLogger().debug(line)
                    vector = line.split() 
                    mu_k = [bool(float(mu_kj)) for mu_kj in vector[0:self.nb_organisms]]
                    logging.getLogger().debug(mu_k)
                    logging.getLogger().debug(len(mu_k))

                    epsilon_k = [float(epsilon_kj) for epsilon_kj in vector[self.nb_organisms+1:]]
                    logging.getLogger().debug(epsilon_k)
                    logging.getLogger().debug(len(epsilon_k))
                    proportion = float(vector[self.nb_organisms])
                    logging.getLogger().debug(proportion)
                    sum_mu_k.append(sum(mu_k))
                    logging.getLogger().debug(sum(mu_k))
                    sum_epsilon_k.append(sum(epsilon_k))
                    logging.getLogger().debug(sum(epsilon_k))

                #cloud is defined by a sum of mu_k near of 0 and a low sum of epsilon
                #shell is defined by an higher sum of epsilon_k
                #persistent is defined by a sum of mu near of nb_organism and a low sum of epsilon

                max_mu_k     = max(sum_mu_k)
                persistent_k = sum_mu_k.index(max_mu_k)

                max_epsilon_k = max(sum_epsilon_k)
                shell_k       = sum_epsilon_k.index(max_epsilon_k)

                cloud_k = set([0,1,2]) - set([persistent_k, shell_k])
                cloud_k = list(cloud_k)[0]

                logging.getLogger().debug(sum_mu_k)
                logging.getLogger().debug(sum_epsilon_k)

                logging.getLogger().debug(persistent_k)
                logging.getLogger().debug(shell_k)
                logging.getLogger().debug(cloud_k)

                partition                 = {}
                partition[persistent_k] = "persistent"
                partition[shell_k]      = "shell"
                partition[cloud_k]      = "cloud"

                for i, line in enumerate(classification_nem_file):
                    elements = [float(el) for el in line.split()]
                    max_prob = max([float(el) for el in elements])
                    positions_max_prob = [pos for pos, prob in enumerate(elements) if prob == max_prob]
                    logging.getLogger().debug(positions_max_prob)
                    logging.getLogger().debug(i)

                    if (len(positions_max_prob)>1):
                        classification[i]="shell"# in case of doubt (equiprobable partition), gene families is attributed to shell
                    else:
                        classification[i] = partition[positions_max_prob.pop()]

                logging.getLogger().debug(partition)
                #logging.getLogger().debug(index.keys())
        except FileNotFoundError:
            logging.getLogger().warning("The number of organisms is too low ("+str(self.nb_organisms)+" organisms in this pangenome) to partition the pangenome graph in persistent, shell, cloud partition, traditional partitions only (Core and Accessory genome) will be provided")

        for node, nem_class in zip(self.neighbors_graph.nodes(), classification):
            nb_orgs=0
            for key in list(self.neighbors_graph.node[node].keys()):
                try:
                    self.neighbors_graph.node[node][key]="|".join(self.neighbors_graph.node[node][key])
                except TypeError:
                    if key == "length":
                        l = list(self.neighbors_graph.node[node][key])
                        self.neighbors_graph.node[node]["length_avg"] = float(np.mean(l))
                        self.neighbors_graph.node[node]["length_med"] = float(np.median(l))
                        self.neighbors_graph.node[node]["length_min"] = min(l)
                        self.neighbors_graph.node[node]["length_max"] = max(l)
                        del self.neighbors_graph.node[node]["length"]

                if key not in RESERVED_WORDS:
                    last_org = key
                    #self.partitions_by_organisms[key][partition[int(nem_class)]].add(self.neighbors_graph.node[node][key])
                    nb_orgs+=1

            self.neighbors_graph.node[node]["partition"]=nem_class
            self.partitions[nem_class].append(node)

            if nb_orgs >= self.nb_organisms:
                self.partitions["core_exact"].append(node)
                self.neighbors_graph.node[node]["partition_exact"]="core_exact"
            else:
                self.partitions["accessory"].append(node)
                self.neighbors_graph.node[node]["partition_exact"]="accessory"

        for node_i, node_j, data in self.neighbors_graph.edges(data = True):
            l = list(data["length"])
            self.neighbors_graph[node_i][node_j]["length_avg"] = float(np.mean(l))
            self.neighbors_graph[node_i][node_j]["length_med"] = float(np.median(l))
            self.neighbors_graph[node_i][node_j]["length_min"] = min(l)
            self.neighbors_graph[node_i][node_j]["length_max"] = max(l)

            del self.neighbors_graph[node_i][node_j]["length"]

        if len(self.families_repeted)>0:
            logging.getLogger().info("Discarded families are:\t"+" ".join(self.families_repeted))
        else:
            logging.getLogger().info("No families have been Discarded")

        logging.getLogger().debug(nx.number_of_edges(self.neighbors_graph))


        self.is_partitionned=True

        # positions = forceatlas2.forceatlas2_networkx_layout(self.neighbors_graph, 
        #                                                    niter=10,
        #                                                    edgeWeightInfluence=0.8)
        # figure = plt.figure()
        # nx.draw_graphviz(self.neighbors_graph, ax=figure.add_subplot(111))#, positions
        # figure.savefig("graph2.png")
    
    # def export_to_DIMACS(self, DIMACS_output_path):
    #     cpt_n = 0
    #     n = ""
    #     e = ""
        
    #     index = {}

    #     for id_org, organism, annot_contigs in enumerate(self.annotations.items()):
    #         for contig, contig_annot in annot_contigs.items():
    #             n+="n label="+contig_annot[0][GENE]+"color="+id_org"\n"
    #             index[contig_annot[0][GENE]]=cpt_n
    #             cpt_n+=1
    #             for index, gene_row in enumerate(contig_annot[1:]):
    #                 n+="n label="+gene_row[GENE]+"color="+id_org"\n"
    #                 index[gene_row[GENE]]=cpt_n
    #                 e+="e "+str(cpt_n-1)+" "+str(cpt_n)+"\n"
    #                 cpt_n+=1

    #     for node_name, node_organisms in self.neighbors_graph.nodes(data=True):

    #     DIMACS_file = open(DIMACS_output_path,'w')
    #     DIMACS_file.write(n)
    #     DIMACS_file.write(p)
    #     DIMACS_file.close()



    def export_to_GEXF(self, graph_output_path, all_node_attributes = False, compressed=False):
        """
            Export the Partionned Pangenome Graph Of Linked Neighbors to a GEXF file  
            :param nem_dir_path: a str containing the path of the GEXF output file
            :param nem_dir_path: a bool specifying if organisms and genes attributes of each family node must be in the file or not. A GEXF file without this attribute can't be imported again after.
            :param nem_dir_path: a bool specifying if the file must be compressed in gzip or not
            :type str: 
            :type bool: 
            :type bool: 
        """ 

        if self.neighbors_graph is None:
            logging.getLogger().error("neighbors_graph is not built, please use the function neighborhoodComputation before")
        elif not self.is_partitionned:
            logging.getLogger().warnings("the pangenome graph is not partionned, please use the function partition before")
        else:
            # if compute_layout:
            #     from fa2l import force_atlas2_layout
            #     logging.getLogger().info("Compute Force Atlas layout")
            #     positions = force_atlas2_layout(self.neighbors_graph,None, iterations=1000, edge_weight_influence=0.7, jitter_tolerance = 20, barnes_hut_theta=1.2)
            #     print(positions)
            #     exit()
            logging.getLogger().info("Writing GEXF file")
            if compressed:
                graph_output_path = gzip.open(graph_output_path,"w")
            nx.write_gexf(self.neighbors_graph, graph_output_path)

    def import_from_GEXF(self, path_graph_to_update):
        """
            Import an already built Partionned Pangenome Graph Of Linked Neighbors from a GEXF file  
            :param nem_dir_path: a str containing the path of the GEXF input file (compressed or not)
            :type str: 

            .. warning:: please import a full GEXF input file having all the attributes
        """ 
        file = gzip.open(path_graph_to_update.name,"r")
        try:
            file.readline()
            file.seek(0)
            self.neighbors_graph = nx.read_gexf(file)
        except IOError as e:
            if e.message == "Not a gzipped file":
                self.neighbors_graph = nx.read_gexf(path_graph_to_update)
            else:
                logging.getLogger().error("Unable to open file "+path_graph_to_update)

        for node, data in self.neighbors_graph.nodes(data=True):
            for key, value in data.items():
                new_value = set(value.split('|'))
                self.neighbors_graph.node[node][key]=new_value
                if key not in RESERVED_WORDS:
                    self.organisms.add(key)
            if len(self.organisms) == 0:
                logging.getLogger().error("No attributes containing organisms names found on this node: "+str(node))
        logging.getLogger().debug(self.organisms)
        self.nb_organisms = len(self.organisms)

        for source, target, data in self.neighbors_graph.edges(data=True):
            try:
                del self.neighbors_graph[source][target]['id']
            except KeyError:
                logging.getLogger().warnings("No previous edge id found in gexf input file for edge: "+source+" <-> "+target)

    def csv_matrix(self, path, sep=",", header=True):
        """
            Exported the pangenome as a csv_matrix similar to the csv matrix exported by Roary (https://sanger-pathogens.github.io/Roary/)
            :param nem_dir_path: a str containing the path of the csv out file
            :param sep: a str containing the separator
            :param header: a bool specifying if the header must be added to the file or not
            :type str: 
            :type str: 
            :type bool: 
        """ 
        logging.getLogger().info("Writing csv matrix")
        with open(path,"w") as matrix:
            if header:
                matrix.write(sep.join(["family",
                                       "partition",
                                       "exact",
                                       "in_nb_org",
                                       "ratio_copy",
                                       "product",
                                       "length_avg",
                                       "length_med",
                                       "length_min",
                                       "length_max"]
                                       +list(self.organisms))+"\n")

            for node, data in self.neighbors_graph.nodes(data=True):
                genes  = [data[org] if org in data else "" for org in self.organisms]
                nb_org = len([gene for gene in genes if gene != ""])
                matrix.write(sep.join([node,
                                       data["partition"],
                                       data["partition_exact"],
                                       str(nb_org),
                                       str(round(data["nb_gene"]/nb_org,2)),
                                       data["product"],
                                       str(round(data["length_avg"],2)),
                                       str(round(data["length_med"],2)),
                                       str(data["length_min"]),
                                       str(data["length_max"])]
                                       +genes)+"\n")

    def delete_pangenome_graph(self, delete_NEM_files = False):
        """
            Delete all the pangenome graph of eventuelly the statistic of the partionning process (including the temporary file)
        """ 
        if delete_NEM_files:
            self.delete_nem_intermediate_files()

        self.nem_output               = None
        self.neighbors_graph          = None
        self.pan_size                 = 0
        self.nem_output               = None
        self.is_partitionned          = False
        self.partitions               = {}
        self.partitions["undefined"]  = list()
        self.partitions["persistent"] = list()
        self.partitions["shell"]      = list()
        self.partitions["cloud"]      = list()
        self.partitions["core_exact"] = list()
        self.partitions["accessory"]  = list()
        self.BIC                      = None
        #self.partitions_by_organisms  = defaultdict(lambda: defaultdict(set))

    def delete_nem_intermediate_files(self):
        """
            Delete all the tempory files used to partion the pangenome
        """ 
        
        if self.nem_intermediate_files is not None:
            logging.getLogger().info("delete "+self.nem_intermediate_files)
            shutil.rmtree(self.nem_intermediate_files)
            self.nem_intermediate_files = None

    def identify_communities_in_each_partition(self):
        """
            Use the Louvain's algorithm to label the nodes by their community in each partition
        """ 
        size_communities=defaultdict(lambda : defaultdict(int))
        for partition in ["persistent","shell", "cloud"]:
            subgraph = self.neighbors_graph.subgraph([nodes for nodes,data in self.neighbors_graph.nodes(data=True) if data['partition']==partition])
            comm = community.best_partition(subgraph)# = nx.algorithms.community.asyn_fluidc(subgraph, 100)
            for node, id_com in comm.items():
                self.neighbors_graph.node[node]['community'] = partition+"_"+str(id_com)
                size_communities[partition][id_com]+=1

    def identify_shell_subpaths(self, k_range = range(2,10),nem_dir_path = tempfile.mkdtemp()):
        
        subgraph = self.neighbors_graph.subgraph([nodes for nodes,data in self.neighbors_graph.nodes(data=True) if data['partition']=='Shell'])

        if not os.path.exists(nem_dir_path):
            #NEM requires 5 files: nem_file.index, nem_file.str, nem_file.dat, nem_file.m and nem_file.nei
            os.makedirs(nem_dir_path)
        self.nem_intermediate_files = nem_dir_path

        logging.getLogger().info("Writing nem_file.str nem_file.index nem_file.nei nem_file.dat and nem_file.m files")
        str_file   = open(nem_dir_path+"/nem_file.str", "w")
        index_file = open(nem_dir_path+"/nem_file.index", "w")
        org_file   = open(nem_dir_path+"/column_org_file", "w")
        nei_file   = open(nem_dir_path+"/nem_file.nei", "w")
        dat_file   = open(nem_dir_path+"/nem_file.dat", "w")

        str_file.write("S\t"+str(nx.number_of_nodes(subgraph))+"\t"+str(self.nb_organisms)+"\n")
        str_file.close()

        nei_file.write("1\n")#to enable weigthed partionning
        
        index = {node: index+1 for index, node in enumerate(subgraph.nodes(data=False))}
        index_inv = {i: node for node, i in index.items()}
        org_file.write(" ".join([org for org in self.organisms])+"\n")
        org_file.close()

        for node_name, node_organisms in subgraph.nodes(data=True):

            index_file.write(str(index[node_name])+"\t"+str(node_name)+"\n")
            logging.getLogger().debug(node_organisms)
            logging.getLogger().debug(self.organisms)
            dat_file.write("\t".join(["1" if org in node_organisms else "0" for org in self.organisms])+"\n")

            row_fam         = []
            row_dist_score  = []
            neighbor_number = 0
            try:
                for neighbor in nx.all_neighbors(subgraph, node_name):
                    #nb_presences = sum([pre_abs for org, pre_abs in self.neighbors_graph[node_name][neighbor].items() if org not in RESERVED_WORDS])
                    #self.neighbors_graph[node_name][neighbor]["weight"]= nb_presences
                    distance_score = subgraph[node_name][neighbor]["weight"]/self.nb_organisms
                    row_fam.append(str(index[neighbor]))
                    row_dist_score.append(str(round(distance_score,4)))
                    neighbor_number += 1
                if neighbor_number>0:
                    nei_file.write("\t".join([str(item) for sublist in [[index[node_name]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
                    #logging.getLogger().debug("\t".join([str(item) for sublist in [[[index[node_name]]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
                else:
                    raise nx.exception.NetworkXError("no neighbors in selected organismss")
            except nx.exception.NetworkXError as nxe:
                logging.getLogger().debug("The family: "+node_name+" is an isolated family")
                nei_file.write(str(index[node_name])+"\t0\n")

        index_file.close()
        nei_file.close()
        dat_file.close()

        for k in k_range:

            logging.getLogger().info("Running NEM uing "+str(k)+" class")
            # weighted_degree = sum(list(self.neighbors_graph.degree(weight="weight")).values())/nx.number_of_edges(self.neighbors_graph)
            # logging.getLogger().debug("weighted_degree: "+str(weighted_degree))
            # logging.getLogger().debug("org/weighted_degree: "+str(self.nb_organisms/weighted_degree))    
            #weighted_degree = sum(self.neighbors_graph.degree(weight="weight").values())/nx.number_of_edges(self.neighbors_graph)


            ALGO           = "ncem" #fuzzy classification by mean field approximation
            ITERMAX        = 100 # number of iteration max 
            MODEL          = "bern" # multivariate Bernoulli mixture model
            PROPORTION     = "pk" #equal proportion :  "p_"     varying proportion : "pk"
            VARIANCE_MODEL = "sk_" #one variance per partition and organism : "sdk"      one variance per partition, same in all organisms : "sd_"   one variance per organism, same in all partion : "s_d"    same variance in organisms and partitions : "s__" 
            NEIGHBOUR_SPEC = "f"# "f" specify to use all neighbors, orther argument is "4" to specify to use only the 4 neighbors with the higher weight (4 because for historic reason due to the 4 pixel neighbors of each pixel)
            CONVERGENCE_TH = "clas "+str(0.000001)

            HEURISTIC      = "heu_d"# "psgrad" = pseudo-likelihood gradient ascent, "heu_d" = heuristic using drop of fuzzy within cluster inertia, "heu_l" = heuristic using drop of mixture likelihood
            STEP_HEURISTIC = 1 # step of beta increase
            BETA_MAX       = float(self.nb_organisms) #maximal value of beta to test,
            DDROP          = 0.8 #threshold of allowed D drop (higher = less detection)
            DLOSS          = 0.5 #threshold of allowed D loss (higher = less detection)
            LLOSS          = 0.02 #threshold of allowed L loss (higher = less detection)
            
            BETA = ["-B",HEURISTIC,"-H",str(STEP_HEURISTIC),str(BETA_MAX),str(DDROP),str(DLOSS),str(LLOSS)]

            command = " ".join([NEM_LOCATION, 
                                nem_dir_path+"/nem_file",
                                str(k),
                                "-a", ALGO,
                                "-i", str(ITERMAX),
                                "-m", MODEL, PROPORTION, VARIANCE_MODEL,
                                "-s r 5",
                                *BETA,
                                "-n", NEIGHBOUR_SPEC,
                                "-c", CONVERGENCE_TH,
                                "-f fuzzy",
                                "-l y"])
         
            logging.getLogger().info(command)
            proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (out,err) = proc.communicate()  
            logging.getLogger().debug(out)
            logging.getLogger().debug(err)

            if os.path.isfile(nem_dir_path+"/nem_file.uf"):
                logging.getLogger().info("Reading NEM results...")
            else:
                logging.getLogger().error("No NEM output file found")
            
            with open(nem_dir_path+"/nem_file.uf","r") as classification_nem_file, open(nem_dir_path+"/nem_file.mf","r") as parameter_nem_file:
                classification = []
                
                parameter = parameter_nem_file.readlines()
                M = float(parameter[6].split()[3]) # M is markov ps-like
                BIC = -2 * M - (k * self.nb_organisms * 2 + k - 1) * math.log(len(self.partitions["Shell"]))

                logging.getLogger().info("The Bayesian Criterion Index of the partionning for "+str(k)+" is "+str(BIC))

                for i, line in enumerate(classification_nem_file):
                    elements = [float(el) for el in line.split()]
                    max_prob = max([float(el) for el in elements])
                    classes = [pos for pos, prob in enumerate(elements) if prob == max_prob]

                    self.neighbors_graph.node[index_inv[i+1]]["subshell"]=str(classes[0])


    def projection_polar_histogram(self, out_dir, organisms_to_project):
        
        for organism in organisms_to_project:
            with open(out_dir+"/"+organism+".csv","w") as out_file:
                out_file.write("gene\tcontig\tori\tfamily\tpartition\tpersistent\tshell\tcloud\n")
                for contig, contig_annot in self.annotations[organism].items():
                    for gene, gene_info in contig_annot.items():
                        if gene_info[FAMILY] not in self.families_repeted:
                            nei_partitions = [self.neighbors_graph.node[nei]["partition"] for nei in nx.all_neighbors(self.neighbors_graph,gene_info[FAMILY])]
                            out_file.write("\t".join([gene,
                                                      contig,
                                                      "T" if (gene_info[NAME].upper() == "DNAA" or gene_info[PRODUCT].upper() == "DNAA") else "F",
                                                      gene_info[FAMILY],
                                                      self.neighbors_graph.node[gene_info[FAMILY]]["partition"],
                                                      str(nei_partitions.count("persistent")),
                                                      str(nei_partitions.count("shell")),
                                                      str(nei_partitions.count("cloud"))])+"\n")

################ END OF CLASS PPanGGOLiN ################

# Calcul du nombre de combinaisons de k elements parmi n
def combinationNb(k,n):
        if (k == 0):
                return 1
        result = 1
        for i in range(0, k):
                result *= float(n - i)/(i + 1);
        return int(round(result))

# Calcul du nombre total de combinaisons uniques de n elements
def combinationTotalNb(size):
        return pow(2,size)-1

# Generation d'une sous-liste aléatoire de taille n
def randomSublist(items,n):
        item_array = np.array(items)
        index_array = np.arange(item_array.size)
        np.random.shuffle(index_array)
        ordered_index_array = sorted(index_array[:n])
        return list(item_array[ordered_index_array])

# Generation de toutes les combinaisons uniques (sans notion d'ordre) d'elements donnes
def exactCombinations(items):
        len_item  = len(items);
        combinations = defaultdict(list)
        for i in range(1, 1<<len_item):
                c = []
                for j in range(0, len_item):
                        if(i & (1 << j)):
                                c.append(items[j])
                combinations[len(c)].append(c)
        return combinations

# Echantillonage proportionnel d'un nombre donne de combinaisons (sans remise)
def samplingCombinations(items, sample_thr, sample_min):
        samplingCombinationList = defaultdict(list)
        item_size = len(items)
        combTotNb = combinationTotalNb(item_size)
        sample_coeff = (float(combTotNb)/sample_thr)
        for k in range(1,item_size+1):
                tmp_comb = []
                combNb = combinationNb(k,item_size)
                combNb_sample = math.ceil(float(combNb)/sample_coeff)
                # Plus petit echantillonage possible pour un k donné = sample_min
                if ((combNb_sample < sample_min) and k != item_size):
                        combNb_sample = sample_min
                i = 0;
                while (i < combNb_sample):
                        comb = randomSublist(items,k)
                        # Echantillonnage sans remise
                        if (comb not in tmp_comb):
                                tmp_comb.append(comb)
                                samplingCombinationList[len(comb)].append(comb)
                                i+=1
        return samplingCombinationList

# Generate list of combinations of organisms exaustively or following a binomial coeficient
def organismsCombinations(orgs, nbOrgThr, sample_thr, sample_min):
        if (len(orgs) <= nbOrgThr):
                comb_list = exactCombinations(orgs)
        else:
                comb_list = samplingCombinations(orgs, sample_thr, sample_min)
        return comb_list

def plot_Rscript(script_outfile, nem_dir, outpdf_Ushape, outpdf_matrix, evoltion_dir, outpdf_evolution, projection_dir, outputdir_pdf_projection, run_script = True):
    """
    run r script
    required the following package to be instaled : ggplot2, reshape2, data.table, ggrepel (last version)

    """

    rscript = """
#!/usr/bin/env R
options(show.error.locations = TRUE)

if(!require("ggplot2")) {{ install.packages("ggplot2", dep = TRUE,repos = "http://cran.us.r-project.org") }}
library("ggplot2")
if(!require("reshape2")) {{ install.packages("reshape2", dep = TRUE,repos = "http://cran.us.r-project.org") }}
library("reshape2")

color_chart = c(pangenome="black", "accessory"="#EB37ED", "core_exact" ="#FF2828", shell = "#00D860", persistent="#F7A507", cloud = "#79DEFF")

########################### START U SHAPED PLOT #################################

binary_matrix         <- read.table("{nem_dir}/nem_file.dat", header=FALSE)
occurences            <- rowSums(binary_matrix)
classification_vector <- apply (read.table("{nem_dir}/nem_file.uf", header=FALSE),1, FUN = function(x){{
ret = which(x==max(x))
if(length(ret)>1){{ret=2}}
return(ret)
}})

means <- data.frame(partition = c("1","2","3"), mean = rep(NA,3))

means[means$partition == "1","mean"] <- mean(occurences[classification_vector == "1"])
means[means$partition == "2","mean"] <- mean(occurences[classification_vector == "2"])
means[means$partition == "3","mean"] <- mean(occurences[classification_vector == "3"])

means <- means[order(means$mean),]

classification_vector[classification_vector == means[1,"partition"]] <- "cloud"
classification_vector[classification_vector == means[2,"partition"]] <- "shell"
classification_vector[classification_vector == means[3,"partition"]] <- "persistent"

c = data.frame(nb_org = occurences, partition = classification_vector)

plot <- ggplot(data = c) + 
    geom_bar(aes_string(x = "nb_org", fill = "partition")) +
    scale_fill_manual(name = "partition", values = color_chart, breaks=c("persistent","shell","cloud")) +
    scale_x_discrete(limits = seq(1, ncol(binary_matrix))) +
    xlab("# of organisms in which each familly is present")+
    ylab("# of families")

ggsave("{outpdf_Ushape}", device = "pdf", height= (par("din")[2]*1.5),plot)

########################### END U SHAPED PLOT #################################

########################### START RESENCE/ABSENCE MATRIX #################################

organism_names          <- unlist(strsplit(readLines("{nem_dir}/column_org_file")," "))
colnames(binary_matrix) <- organism_names
nb_org                  <- ncol(binary_matrix)

binary_matrix_hclust    <- hclust(dist(t(binary_matrix), method="binary"))
binary_matrix           <- data.frame(binary_matrix,"NEM partitions" = classification_vector, occurences = occurences, check.names=FALSE)

binary_matrix[occurences == nb_org, "Former partitions"] <- "core_exact"
binary_matrix[occurences != nb_org, "Former partitions"] <- "accessory"
binary_matrix = binary_matrix[order(match(binary_matrix$"NEM partitions",c("persistent", "shell", "cloud")),
                                    match(binary_matrix$"Former partitions",c("core_exact", "accessory")),
                                    -binary_matrix$occurences),
                              c(binary_matrix_hclust$label[binary_matrix_hclust$order],"NEM partitions","Former partitions")]

binary_matrix$familles <- seq(1,nrow(binary_matrix))
data = melt(binary_matrix, id.vars=c("familles"))

colnames(data) = c("fam","org","value")

data$value <- factor(data$value, levels = c(1,0,"persistent", "shell", "cloud", "core_exact", "accessory"), labels = c("presence","absence","persistent", "shell", "cloud", "core_exact", "accessory"))

plot <- ggplot(data = data)+
        geom_raster(aes_string(x="org",y="fam", fill="value"))+
        scale_fill_manual(values = c("presence"="green","absence"="grey80",color_chart)) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.border = element_blank(), panel.background = element_blank())

ggsave("{outpdf_matrix}", device = "pdf", plot)

########################### END PRESENCE/ABSENCE MATRIX #################################

########################### START EVOLUTION CURVE #################################

if(!require("ggrepel") || packageVersion("ggrepel") < "0.6.6") {{ install.packages("ggrepel", dep = TRUE, repos = "http://cran.us.r-project.org") }}
library("ggrepel")

if(!require("data.table")) {{ install.packages("data.table", dep = TRUE,repos = "http://cran.us.r-project.org") }}
library("data.table")

if (file.exists("{evoltion_dir}/stat_evol.txt")){{
    data <- read.table("{evoltion_dir}/stat_evol.txt", header = TRUE)


    data <- melt(data, id = "nb_org")
    colnames(data) <- c("nb_org","partition","value")

    final_state = data[data$nb_org == max(data$nb_org,na.rm=T),]
    final_state = final_state[!duplicated(final_state), ]
    final <- structure(names = as.character(final_state$partition), as.integer(final_state$value))

    #gamma and kappa are calculated according to the Tettelin et al. 2008 approach
    median_by_nb_org <- setDT(data)[,list(med=as.numeric(median(value))), by=c("nb_org","partition")]
    colnames(median_by_nb_org) <- c("nb_org_comb","partition","med")

    for (part in as.character(final_state$partition)){{
        regression  <- nls(med~kapa*(nb_org_comb^gama),median_by_nb_org[which(median_by_nb_org$partition == part),],start=list(kapa=1000,gama=1))
        coefficient <- coef(regression)
        final_state[final_state$partition == part,"formula" ] <- paste0("n == ", format(coefficient["kapa"],decimal.mark = ",",digits =2),"~N^{{",format(coefficient["gama"],digits =2),"}}")
    }}

    plot <- ggplot(data = data, aes_string(x="nb_org",y="value", colour = "partition"))+
            ggtitle(bquote(list("Rarefaction curve. Heap's law parameters based on Tettelin et al. 2008 approach", n == kappa~N^gamma)))+
            geom_smooth(data        = median_by_nb_org[median_by_nb_org$partition %in% c("pangenome","shell","cloud","accessory", "persistent", "core_exact") ,],# 
                        mapping     = aes_string(x="nb_org_comb",y="med",colour = "partition"),
                        method      = "nls",
                        formula     = y~kapa*(x^gama),method.args =list(start=c(kapa=1000,gama=1)),
                        linetype    ="twodash",
                        size        = 1.5,
                        se          = FALSE,
                        show.legend = FALSE)+
            stat_summary(fun.ymin = function(z) {{ quantile(z,0.25) }},  fun.ymax = function(z) {{ quantile(z,0.75) }}, geom="ribbon", alpha=0.1,size=0.1, linetype="dashed", show.legend = FALSE)+
            stat_summary(fun.y=median, geom="line",size=0.5)+
            stat_summary(fun.y=median, geom="point",shape=4,size=1, show.legend = FALSE)+
            stat_summary(fun.ymax=max,fun.ymin=min,geom="errorbar",linetype="dotted",size=0.1,width=0.2)+
            scale_x_continuous(breaks = as.numeric(unique(data$nb_org)))+
            scale_y_continuous(limits=c(0,max(data$value,na.rm=T)), breaks = seq(0,max(data$value,na.rm=T),1000))+
            scale_colour_manual(name = "NEM partitioning", values = color_chart, breaks=names(sort(final, decreasing = TRUE)))+
            geom_label_repel(data = final_state, aes_string(x="nb_org", y="value", colour = "partition", label = "value"), show.legend = FALSE,
                      fontface = 'bold', fill = 'white',
                      box.padding = unit(0.35, "lines"),
                      point.padding = unit(0.5, "lines"),
                      segment.color = 'grey50',
                      nudge_x = 45) +
            geom_label_repel(data = final_state, aes(x = nb_org*0.9, y = value, label = formula), size = 2, parse = TRUE, show.legend = FALSE, segment.color = NA) + 
            xlab("# of organisms")+
            ylab("# of families")+
            ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    
    ggsave("{outpdf_evolution}", device = "pdf", width = (par("din")[1]*2) ,plot)

}}
########################### END EVOLUTION CURVE #################################

########################### START PROJECTION #################################

for (org_csv in list.files(path = "{projection_dir}", pattern = "*.csv$")){{
    
    data <- read.table(paste0("{projection_dir}/",org_csv), header = T)
    data <- cbind(data, pos = seq(nrow(data)))

    max_degree_log2p1 <- max(apply(data,1,FUN = function(x){{
            sum(log2(as.numeric(x[6:8])+1))
        }}))

    ori <- which(data$ori == T, arr.ind=T)
    data$ori <- NULL

    duplicated_fam     <- unique(data[duplicated(data$family),"family"])
    data$family <- ifelse(data$family %in% duplicated_fam, data$family, NA)
    data$family = as.factor(data$family)
    colors_duplicated_fam <- rainbow(length(duplicated_fam))
    names(colors_duplicated_fam) <- duplicated_fam

    data_melted <- melt(data, id.var=c("contig", "gene","family","partition","pos"))
    data_melted$variable <- factor(data_melted$variable, levels = rev(c("persistent","shell","cloud")), ordered=TRUE)

    contig <- unique(data_melted$contig)
    contig_color <-  rainbow(length(contig))
    names(contig_color) <- contig

    data_melted$value <- log2(data_melted$value+1)

    plot = ggplot(data = data_melted)+
    ggtitle(paste0("plot corresponding to the file", org_csv))+
    geom_bar(aes_string(x = "gene", y = "value", fill = "variable"),stat="identity", show.legend = FALSE)+
    scale_y_continuous(limits = c(-30, max_degree_log2p1), breaks = seq(0,ceiling(max_degree_log2p1)))+
    geom_hline(yintercept = 0)+
    geom_rect(aes_string(xmin ="pos-1/2", xmax = "pos+1/2", fill = "partition"), ymin = -10, ymax=-1, color = NA, show.legend = FALSE)+
    geom_hline(yintercept = -10)+
    geom_rect(aes_string(xmin ="pos-1/2", xmax = "pos+1/2", fill = "family"), ymin = -20, ymax=-11,  color = NA, show.legend = FALSE)+
    geom_hline(yintercept = -20)+
    geom_rect(aes_string(xmin ="pos-1/2", xmax = "pos+1/2", fill = "contig"), ymin = -30, ymax=-21,  color = NA)+
    geom_vline(xintercept = ori)+
    scale_fill_manual(values = c(color_chart,colors_duplicated_fam, contig_color), na.value = "grey80")+
    coord_polar()+
    ylab("log2(degree+1) of the families in wich each gene is")+
    theme(axis.line        = ggplot2::element_blank(),
                        axis.text.x      = ggplot2::element_blank(),
                        axis.ticks.x       = ggplot2::element_blank(),
                        axis.title.x     = ggplot2::element_blank(),
                        panel.background = ggplot2::element_blank(),
                        panel.border     = ggplot2::element_blank(),
                        panel.grid.major.x = ggplot2::element_blank(),
                        panel.grid.minor.x = ggplot2::element_blank(),
                        plot.background  = ggplot2::element_blank(),
                        plot.margin      = grid::unit(c(0,0,0,0), "cm"),
                        panel.spacing    = grid::unit(c(0,0,0,0), "cm"))

    ggsave(paste0("{outputdir_pdf_projection}/projection_",org_csv,".pdf"), device = "pdf", height= 40, width = 49, plot)

}}

########################### END PROJECTION #################################

    """.format(nem_dir                  = nem_dir,
               outpdf_Ushape            = outpdf_Ushape,
               outpdf_matrix            = outpdf_matrix,
               evoltion_dir             = evoltion_dir,
               outpdf_evolution         = outpdf_evolution,
               projection_dir           = projection_dir,
               outputdir_pdf_projection = outputdir_pdf_projection)
    logging.getLogger().info("Writing R script generating plot")
    with open(script_outfile,"w") as script_file:
        script_file.write(rscript)
    logging.getLogger().info("Running R script generating plot")
    if run_script:
        logging.getLogger().info("Rscript "+script_outfile)
        proc = subprocess.Popen("Rscript "+script_outfile, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        (out,err) = proc.communicate()
        logging.getLogger().debug(out)
        logging.getLogger().debug(err)


if __name__=='__main__':
    """
    --organims is a tab delimited files containg at least 2 mandatory fields by row and as many optional field as circular contig. Each row correspond to an organism to be added to the pangenome.
    Reserved words are : "id", "label", "name", "weight", "partition", "partition_exact"
    The first field is the organinsm name (id should be unique not contain any spaces, " or ' and reserved words).
    The second field is the gff file associated to the organism. This path can be abolute or relative. The gff file must contain an id feature for each CDS (id should be unique not contain any spaces, " or ' and reserved words).
    The next fields contain the name of perfectly assemble circular contigs (contig name must should be unique and not contain any spaces, " or ' and reserved words).
    example:

    """
    parser = argparse.ArgumentParser(prog = "", description='Build a partitioned pangenome graph from annotated genomes and gene families')
    parser.add_argument('-?', '--version', action='version', version='0.1.1')
    parser.add_argument('-o', '--organisms', type=argparse.FileType('r'), nargs=1, help="""
A tab delimited file containg at least 2 mandatory fields by row and as many optional fields as the number of well assembled circular contigs. Each row correspond to an organism to be added to the pangenome.
Reserved words are : "id", "label", "name", "weight", "partition", "partition_exact"
The first field is the organinsm name (id should be unique not contain any spaces, " or ' and reserved words).
The second field is the gff file associated to the organism. This path can be abolute or relative. The gff file must contain an id feature for each CDS (id should be unique across all pangenome and not contain any spaces, " or ' and reserved words).
The next fields contain the name of perfectly assembled circular contigs.
Contig names should be unique and not contain any spaces, quote, double quote and reserved words.
example:""", required=True)
    parser.add_argument('-gf', '--gene_families', type=argparse.FileType('r'), nargs=1, help="""
A tab delimited file containg the gene families. Each row contain 2 fields.
Reserved words are : "id", "label", "name", "weight", "partition", "partition_exact"
The first field is the family name.
The second field is the gene name.
families are intended to be grouped by chuncks of row.
the families name can be any string but must should be unique and not contain any spaces, " or ' and reserved words
As a conventation, it is recommanded to use the name of the most reprensative gene of the families as the family name.
Gene name can be any string corresponding to the if feature in the gff files. they should be unique and not contain any spaces, " or ' and reserved words.
example:""",  required=True)
    parser.add_argument('-od', '--output_directory', type=str, nargs=1, default=["output.dir"], help="""
The output directory""")
    parser.add_argument('-f', '--force', action="store_true", help="""
Force overwriting existing output directory""")
    parser.add_argument('-r', '--remove_high_copy_number_families', type=int, nargs=1, default=[0], help="""
Remove families having a number of copy of one families above or equal to this threshold in at least one organism (0 or negative value keep all families whatever their occurence). 
When -u is set, only work on new organisms added""")
    parser.add_argument('-s', '--infere_singleton', default=False, action="store_true", help="""
If a gene id found in a gff file is absent of the gene families file, the singleton will be automatically infered as a gene families having a single element. 
if this argument is not set, the program will raise KeyError exception if a gene id found in a gff file is absent of the gene families file.""")
#    parser.add_argument("-u", "--update", default = None, type=argparse.FileType('r'), nargs=1, help="""
# Pangenome Graph to be updated (in gexf format)""")
    parser.add_argument("-b", "--beta_smoothing", default = [-1.00], type=float, nargs=1, help = """
Coeficient of smoothing all the partionning based on the Markov Random Feild leveraging the weigthed pangenome graph. A positive float, 0.0 means to discard spatial smoothing and 'inf' to find beta automaticaly (increase the computation time) 
""")
    parser.add_argument("-fd", "--free_dispersion", default = False, action="store_true", help = """
Specify if the dispersion around the centroid vector of each paritition is the same for all the organisms or if the dispersion is free
""")
    parser.add_argument("-df", "--delete_nem_intermediate_files", default=False, action="store_true", help="""
Delete intermediate files used by NEM. Do not delete these files if you can the gerate plot latter using the generated Rscript""")
    parser.add_argument("-c", "--compress_graph", default=False, action="store_true", help="""
Compress (using gzip) the file containing the partionned pangenome graph""")
    parser.add_argument("-ss", "--subpartition_shell", default = 0, type=int, nargs=1, help = """
Subpartition the shell genome in n subpartition, n can be ajusted automatically if n = -1, 0 desactivate shell genome subpartitioning""")
    parser.add_argument("-v", "--verbose", default=True, action="store_true", help="""
Show information message, otherwise only errors and warnings will be displayed""")
    parser.add_argument("-vv", "--verbose_debug", default=False, action="store_true", help="""
Show all messages including debug ones""")
    parser.add_argument("-as", "--already_sorted", default=False, action="store_true", help="""
Accelerate loading of gff files if there are sorted by the coordinate of gene annotations (starting point) for each contig""")
    parser.add_argument("-l", "--ligth", default=False, action="store_true", help="""
Free the memory elements which are no longer used""")
    parser.add_argument("-p", "--plots", default=False, action="store_true", help="""
Generate Rscript able to draw plots and run it. (required R in the path and the packages ggplot2, ggrepel, data.table and reshape2 to be installed)""")
    parser.add_argument("-di", "--directed", default=True, action="store_true", help="""
directed or not directed graph
""")
    parser.add_argument("-e", "--evolution", default=False, action="store_true", help="""
Relaunch the script using less and less organism in order to obtain a curve of the evolution of the pangenome metrics
""")
    parser.add_argument("-ep", "--evolution_resampling_param", type=int, nargs=4, default=[4,10,30,1], help="""
1st argument is the number of threads (int) to use to resemple in parallel the pangenome
2nd argument is the minimun number of resampling for each number of organisms
3nd argument is the maximun number of resampling for each number of organisms
4th argument is the step between each number of organisms
""")
    parser.add_argument("-pr", "--projection", type = int, nargs = "+", help="""
Project the graph as a circos plot on each organism.
Expected parameters are the line number (1 based) of each organism on which the graph will be projected providing a circos plot (well assembled representative organisms must be prefered).
0 means all organisms (it is discouraged to use -p and -pr 0 in the same time because the projection of the graph on all the organisms can take a long time).
""")

    options = parser.parse_args()

    level = logging.WARNING
    if options.verbose:
        level = logging.INFO
    if options.verbose_debug:
        level = logging.DEBUG

    logging.basicConfig(stream=sys.stdout, level = level, format = '\n%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s', datefmt='%H:%M:%S')

    

    logging.getLogger().info("Command: "+" ".join([arg for arg in sys.argv]))
    logging.getLogger().info("Python version: "+sys.version)
    logging.getLogger().info("Networkx version: "+nx.__version__)

    OUTPUTDIR       = options.output_directory[0]
    NEMOUTPUTDIR    = OUTPUTDIR+"/NEM_results/"
    FIGUREOUTPUTDIR = OUTPUTDIR+"/figures/"
    PROJECTION      = OUTPUTDIR+"/projection/"
    EVOLUTION       = OUTPUTDIR+"/evolution/"
    list_dir        = [OUTPUTDIR,NEMOUTPUTDIR,FIGUREOUTPUTDIR]
    if options.projection:
        list_dir.append(PROJECTION)
    if options.evolution:
        list_dir.append(EVOLUTION)
        EVOLUTION_STAT_FILE = EVOLUTION+"/"+"stat_evol.txt"
        (NB_THREADS, RESAMPLING_MIN, RESAMPLING_MAX, STEP) = options.evolution_resampling_param
    for directory in list_dir:
        if not os.path.exists(directory):
            os.makedirs(directory)
        elif not options.force:
            logging.getLogger().error(directory+" already exist")
            exit(1)

    GEXF_GRAPH_FILE     = OUTPUTDIR+"/"+"graph.gexf"
    MATRIX_CSV_FILE     = OUTPUTDIR+"/matrix.csv"

    #-------------
    start_loading_file = time.time()
    pan = PPanGGOLiN("file",
                     options.organisms[0],
                     options.gene_families[0],
                     options.remove_high_copy_number_families[0],
                     options.infere_singleton,
                     options.already_sorted)

    # if options.update is not None:
    #     pan.import_from_GEXF(options.update[0])
    end_loading_file = time.time()
    #-------------

    #-------------
    logging.getLogger().info("Neighborhood Computation...")
    start_neighborhood_computation = time.time()
    pan.neighborhood_computation(options.directed, options.ligth)
    end_neighborhood_computation = time.time()
    #-------------

    #-------------
    logging.getLogger().info("Partionning...")
    start_partitioning = time.time()
    pan.partition(nem_dir_path = NEMOUTPUTDIR, beta = options.beta_smoothing[0], free_dispersion = options.free_dispersion)
    end_partitioning = time.time()
    #-------------

    #-------------
    # start_identify_communities = time.time()
    # pan.identify_communities_in_each_partition()
    # end_identify_communities = time.time()
    #pan.identify_shell_subpaths()
    #-------------

    #-------------
    start_writing_output_file = time.time()
    if options.compress_graph:
        pan.export_to_GEXF(GEXF_GRAPH_FILE+".gz", compressed=True)
    else:
        pan.export_to_GEXF(GEXF_GRAPH_FILE, compressed=False)
    for filename, families in pan.partitions.items(): 
        file = open(OUTPUTDIR+"/"+filename+".txt","w")
        file.write("\n".join(families))
        file.close()
    pan.csv_matrix(MATRIX_CSV_FILE)
    end_writing_output_file = time.time()
    #-------------

    logging.getLogger().info(pan)

    #-------------
    if options.projection:
        logging.getLogger().info("Projection...")
        start_projection = time.time()
        pan.projection_polar_histogram(PROJECTION, [pan.organisms.__getitem__(index-1) for index in options.projection] if options.projection[0] >   0 else list(pan.organisms))
        end_projection = time.time()
    #-------------


    # print(pan.partitions_by_organisms)
    # partitions_by_organisms_file = open(OUTPUTDIR+"/partitions_by_organisms.txt","w")
    # exact_by_organisms_file = open(OUTPUTDIR+"/exacte_by_organisms.txt","w")
    # for org, partitions in pan.partitions_by_organisms.items(): 
    #     partitions_by_organisms_file.write(org+"\t"+str(len(partitions["persistent"]))+
    #                                            "\t"+str(len(partitions["shell"]))+
    #                                            "\t"+str(len(partitions["cloud"]))+"\n")
    #     exact_by_organisms_file.write(org+"\t"+str(len(partitions["core_exact"]))+
    #                                       "\t"+str(len(partitions["accessory"]))+"\n")
    # partitions_by_organisms_file.close()
    # exact_by_organisms_file.close()


    #-------------
    if options.evolution:
        logging.getLogger().info("Evolution...")
        start_evolution = time.time()
        logging.disable(logging.INFO)# disable message info

        combinations = organismsCombinations(list(pan.organisms), nbOrgThr=1, sample_thr=RESAMPLING_MAX, sample_min=RESAMPLING_MIN)
        del combinations[pan.nb_organisms]
        del combinations[1]
        shuffled_comb = [OrderedSet(comb) for nb_org, combs in combinations.items() for comb in combs if nb_org%STEP == 0]
        random.shuffle(shuffled_comb)

        with open(EVOLUTION_STAT_FILE,"w") as evol:

            evol.write("nb_org\tpersistent\tshell\tcloud\tcore_exact\taccessory\tpangenome\n")
            def write_stat(a_pan):#internal function
                evol.write("\t".join([str(a_pan.nb_organisms),
                                       str(len(a_pan.partitions["persistent"])),
                                       str(len(a_pan.partitions["shell"])),
                                       str(len(a_pan.partitions["cloud"])),
                                       str(len(a_pan.partitions["core_exact"])),
                                       str(len(a_pan.partitions["accessory"])),
                                       str(a_pan.pan_size)])+"\n")
                evol.flush()
            
            def pan_sample(index):#internal function
                
                pan_sample = PPanGGOLiN("args",
                                       {org:pan.annotations[org] for org in shuffled_comb[index]},
                                       shuffled_comb[index],
                                       pan.circular_contig_size,
                                       pan.families_repeted)

                pan_sample.neighborhood_computation(options.directed, light=True)
                pan_sample.partition(EVOLUTION+"/nborg"+str(len(shuffled_comb[index]))+"_"+str(index), options.beta_smoothing[0], options.free_dispersion)

                write_stat(pan_sample)
                evol.flush()
                
                pan_sample.delete_pangenome_graph(delete_NEM_files = options.delete_nem_intermediate_files)
                del pan_sample
                

            write_stat(pan)

            with ProcessPoolExecutor(NB_THREADS) as executor:
                futures = [executor.submit(pan_sample,i) for i in range(len(shuffled_comb))]

                for f in tqdm(as_completed(futures), total = len(shuffled_comb), unit = 'pangenome resampled',  unit_scale = True):
                    if f.exception() is not None:
                        logging.getLogger().error(f.exception())
                        executor.shutdown()
                        exit(1)
            end_evolution = time.time()
            logging.disable(logging.NOTSET)#restaure message info
    #-------------

    logging.getLogger().info("\n"+
    "Execution time of file loading: """ +str(round(end_loading_file-start_loading_file, 2))+" s\n"+
    "Execution time of neighborhood computation: " +str(round(end_neighborhood_computation-start_neighborhood_computation, 2))+" s\n"+
    "Execution time of partitioning: " +str(round(end_partitioning-start_partitioning, 2))+" s\n"+
    #"Execution time of community identification: " +str(round(end_identify_communities-start_identify_communities, 2))+" s\n"+
    "Execution time of writing output files: " +str(round(end_writing_output_file-start_writing_output_file, 2))+" s\n"+
    (("Execution time of projection: " +str(round(end_projection-start_projection, 2))+" s\n") if options.projection else "")+
    (("Execution time of evolution: " +str(round(end_evolution-start_evolution, 2))+" s\n") if options.evolution else "")+

    "Total execution time: " +str(round(time.time()-start_loading_file, 2))+" s\n")

    logging.getLogger().info("""
The pangenome computation is complete. 
Plots will be generated using R (in the directory: """+FIGUREOUTPUTDIR+""").
If R and the required package (ggplot2, reshape2, ggrepel(>0.6.6), data.table) are not instaled don't worry the R script will be saved in the directory allowing to generate the figures latter""")

    plot_Rscript(script_outfile           = OUTPUTDIR+"/generate_plots.R",
                 nem_dir                  = NEMOUTPUTDIR,
                 outpdf_Ushape            = FIGUREOUTPUTDIR+"/Ushaped_plot.pdf",
                 outpdf_matrix            = FIGUREOUTPUTDIR+"/Presence_absence_matrix_plot.pdf",
                 evoltion_dir             = EVOLUTION,
                 outpdf_evolution         = FIGUREOUTPUTDIR+"/evolution.pdf",
                 projection_dir           = PROJECTION,
                 outputdir_pdf_projection = FIGUREOUTPUTDIR,
                 run_script               = options.plots)

    if options.delete_nem_intermediate_files:
            pan.delete_nem_intermediate_files()   

    logging.getLogger().info("Finished !")
    exit(0)
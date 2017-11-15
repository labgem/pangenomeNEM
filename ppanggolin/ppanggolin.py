#!/usr/bin/python3
# -*- coding: iso-8859-1 -*-

from collections import defaultdict, OrderedDict
from ordered_set import OrderedSet
import networkx as nx
import os
import sys
import math
import logging
import shutil
import gzip
import numpy as np
#import community
import tempfile
import subprocess
from tqdm import tqdm
import mmap

#import forceatlas2 

NEM_LOCATION  = os.path.dirname(os.path.abspath(__file__))+"/NEM/nem_exe"

(TYPE, FAMILY, START, END, STRAND, NAME, PRODUCT) = range(0, 7)#data index in annotation
(ORGANISM_ID, ORGANISM_GFF_FILE) = range(0, 2)#data index in the file listing organisms 

(GFF_seqname, GFF_source, GFF_feature, GFF_start, GFF_end, GFF_score, GFF_strand, GFF_frame, GFF_attribute) = range(0,9) 

RESERVED_WORDS = set(["id", "label", "name", "weight", "partition", "partition_exact", "length", "length_min", "length_max", "length_avg", "length_avg", "product", 'nb_gene', 'community'])

"""  
    :mod:`ppanggolin` -- Depict microbial diversity
===================================

.. module:: ppanggolin
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

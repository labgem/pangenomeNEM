#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

import argparse
from collections import defaultdict
from collections import OrderedDict
from collections import Counter
from ordered_set import OrderedSet
import networkx as nx
import subprocess
import os
import sys
import re
import math
import logging
import gffutils
import random 
import string
import shutil
import gzip
import operator
import numpy as np
import time
import community
import tempfile

#import forceatlas2 

NEM_LOCATION  = os.path.dirname(os.path.abspath(__file__))+"/../NEM/nem_exe"
(GENE, TYPE, FAMILY, START, END, STRAND, NAME, PRODUCT) = range(0, 8)#data index in annotation
(ORGANISM_ID, ORGANISM_GFF_FILE) = range(0, 2)#data index in the file listing organisms 
#(ORGANISM, SEQUENCE, INDEX) = range(0, 3)# index for annotation in gene_location

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
                    * partitions["Core_Exact"] contains the list of core families (core exact)
                    * partitions["Accessory"] contains the list of families not in the core exact
                    * partitions["Persistent"] contains the list of persistent families
                    * partitions["Shell"] contains the list of shell families
                    * partitions["Cloud"] contains the list of cloud families

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
        #self.gene_location            = OrderedDict()
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
        self.partitions["Persistent"] = list()
        self.partitions["Shell"]      = list()
        self.partitions["Cloud"]      = list()
        self.partitions["Core_Exact"] = list()
        self.partitions["Accessory"]  = list()
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

    def __initialize_from_files(self, organisms_file, families_tsv_file, lim_occurence = 0, infere_singleton = False):
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

        for line in organisms_file:
            elements = line.split()
            if len(elements)>2:
                self.circular_contig_size.update({contig_id: 0 for contig_id in elements[2:len(elements)]})# size of the circular contig is initialized to zero (waiting to read the gff files to fill the dictionnaries with the correct values)
            self.annotations[elements[0]] = self.__load_gff(elements[ORGANISM_GFF_FILE], families, elements[ORGANISM_ID], lim_occurence, infere_singleton)

        check_circular_contigs = {contig: size for contig, size in self.circular_contig_size.items() if size == 0 }
        if len(check_circular_contigs) > 0:
            logging.getLogger().error("""
                The following identifiers of circular contigs in the file listing organisms have not been found in any region feature of the gff files: """+"\t".join(check_circular_contigs.keys()))
            exit()
    def __load_gff(self, gff_file, families, organism, lim_occurence = 0, infere_singleton = False):
        """
            Load the content of a gff file
            :param gff_file: a valid gff file where only feature of the type 'CDS' will be imported as genes. Each 'CDS' feature must have a uniq ID as attribute (afterall called gene id).
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

        logging.getLogger().info("Reading "+gff_file+" file ...")

        if organism not in self.organisms:
            self.organisms.add(organism)
            db_gff = gffutils.create_db(gff_file, ':memory:')
            annot = defaultdict(list)

            # prev = None
            ctp_prev = 1
            cpt_fam_occ = defaultdict(int)

            for row in db_gff.all_features(featuretype=('region','CDS'), order_by=('seqid','start')):
                if row.featuretype == 'region':
                    if row.seqid in self.circular_contig_size:
                        self.circular_contig_size = row.end
                else:
                    # logging.getLogger().debug(row)
                    protein = row.id
                    try:
                        family = families[protein]
                    except KeyError:
                        if infere_singleton:
                            families[protein] = protein
                            family           = families[protein]
                            logging.getLogger().info("infered singleton: "+protein)
                        else:
                            raise KeyError("Unknown families:"+protein, ", check your families file or run again the program using the option to infere singleton")

                    # if family == prev:
                    #     family = family+"-"+str(ctp_prev)

                    cpt_fam_occ[family]+=1
                    prev = families[protein]

                    try:
                        name = row.attributes['Name'].pop()
                    except:
                        try:
                            name = row.attributes['gene'].pop()
                        except:
                            name = ""

                    try:
                        product = row.attributes['product'].pop()
                    except:
                        product = ""

                    annot_row = [protein,"CDS",family,row.start,row.end,row.strand, name, product]
                    #self.gene_location[protein] = tuple([organism, row.seqid, len(annot[row.seqid])])
                    annot[row.seqid].append(annot_row)

                    # logging.getLogger().debug(annot[self.gene_location[protein][SEQUENCE]][self.gene_location[protein][INDEX]][GENE])

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
            pan_str += "Exact core-genome size:"+str(len(self.partitions["Core_Exact"]))+"\n"
            pan_str += "Exact variable-genome size:"+str(self.pan_size-len(self.partitions["Core_Exact"]))+"\n"
            pan_str += "Persistent genome size:"+str(len(self.partitions["Persistent"]))+"\n"
            pan_str += "Shell genome size:"+str(len(self.partitions["Shell"]))+"\n"
            pan_str += "Cloud genome cloud:"+str(len(self.partitions["Cloud"]))+"\n"
        else:
            pan_str += "No partitioning have been performed on this Pangenome instance\n"
            pan_str += "Run the partitioning function to obtain more detailled statistics...\n"
        pan_str += "----------------------------------\n"

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




        # try:
        #     self.neighbors_graph.node[fam_id][org].add(gene)
        #     # if multi_copy is not None:
        #     #     multi_copy[fam_id] = max(multi_copy[fam_id],len(self.neighbors_graph.node[fam_id][org]))#count max number of multicopy
        # except KeyError:# the exception is the rule in case of mono copy gene
        #     self.neighbors_graph.node[fam_id][org]=set([gene])
        # try:
        #     self.neighbors_graph.node[fam_id]["name"].add(name)
        # except KeyError:
        #     self.neighbors_graph.node[fam_id]["name"]=set([name])
        # try:
        #     self.neighbors_graph.node[fam_id]["length"].add(length)
        # except KeyError:
        #     self.neighbors_graph.node[fam_id]["length"]=set([length])
        # try:
        #     self.neighbors_graph.node[fam_id]["product"].add(product)
        # except KeyError:
        #     self.neighbors_graph.node[fam_id]["product"]=set([product])
    def __add_link(self, fam_id, fam_id_nei, org, length):
        """
            Add line between families of a the pangenome graph
            :param fam_id: The family identifier the first node (need to be have at least one gene belonging to this family already added to the graph via the method __add_gene())
            :param fam_id_nei: The family identifier the second node (need to be have at least one gene belonging to this family already added to the graph via the method __add_gene())
            :param gene : The organism name supporting this link
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
            self.neighbors_graph.node[fam_id]["length"].add(length)
        except KeyError:
            self.neighbors_graph.node[fam_id]["length"]=set([length])

    # def __extract_families_of_genomic_context(self, gene, context_size = 5):
    #     ret = list()
    #     try:
    #         (organism, sequence, index) = self.gene_location[gene]
    #     except KeyError:
    #         logging.getLogger().debug("not a gene id: "+str(gene))

    #     if sequence not in self.circular_contig: # linear contig
    #         min_boundary = index-context_size if (index-context_size)>0 else 0
    #         for annot in self.annotations[organism][sequence][min_boundary:(index+context_size+1)]:#python find automatically the max boudary
    #             ret.append(annot[FAMILY])
    #     else:# circular contig 
    #         seq_size      = len(self.annotations[organism][sequence])
    #         window_length = (context_size*2+1)
    #         if seq_size > window_length:# circular contig larger than window
    #             for i in range(index-context_size,index+context_size+1):
    #                 ret.append(self.annotations[organism][sequence][i%seq_size][FAMILY])
    #         else:# circular contig smaler than window
    #             for i in range(seq_size):
    #                 ret.append(self.annotations[organism][sequence][i][FAMILY])
    #     return tuple(ret)
    # 
    # def __untangle_multi_copy_families(self, multi_copy, context_size = 5):
    #     import editdistance
    #     from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
    #     from scipy.spatial.distance import squareform, pdist
    #     # from matplotlib import pyplot as plt

    #     result = defaultdict(set)

    #     def __group_similar_contexts(contexts, context_size, distance_max):
    #         def compute_minimal_edit_distance(contexts1, contexts2):
    #             reverse_distance = 0
    #             distance         = editdistance.eval(contexts1, contexts2)
    #             if distance != 0:
    #                 reverse_distance = editdistance.eval(contexts1, tuple(reversed(contexts2)))
    #             return min(distance, reverse_distance)
    #         contexts_keys = [context for context in contexts.keys()]
    #         # pdist(np.array(contexts_keys))
            
    #         if len(contexts_keys)>1:
    #             distance_matrix = []
    #             for i in range(len(contexts_keys)):
    #                 for j in range(i):
    #                     distance_matrix.append(compute_minimal_edit_distance(contexts_keys[j], contexts_keys[i]))
    #             logging.getLogger().debug(distance_matrix)
    #             logging.getLogger().debug(squareform(distance_matrix))
    #             hierarchical_linkage = linkage(squareform(distance_matrix), 'single')
    #             logging.getLogger().debug(hierarchical_linkage)
    #             groups = fcluster(hierarchical_linkage, t = distance_max, criterion = 'distance')

    #             #logging.getLogger().debug("|".join([str(p) for p in partition]))
    #             # plt.figure(figsize=(25, 10))
    #             # plt.title('Hierarchical Clustering Dendrogram (len), '+str(nb_copy))
    #             # plt.xlabel('sample index')
    #             # plt.ylabel('distance')
    #             # dendrogram(
    #             # hierarchical_linkage,
    #             # leaf_rotation=90., leaf_font_size=8.,)
    #             # plt.show()
                
    #             for index, group in enumerate(groups):
    #                 result[group] |= contexts[contexts_keys[index]]
    #         else:
    #             result[0] = contexts[contexts_keys[0]]


    #         return result
    #             # complete_contexts = [context for context in contexts.keys()]

    #             # computation_vector     = np.vectorize(compute_minimal_edit_distance)
    #             # logging.getLogger().debug(contexts)
    #             # logging.getLogger().debug(list(enumerate(complete_contexts)))
    #             # distance_matrix        = np.fromfunction(computation_vector, shape=(len(complete_contexts),len(complete_contexts)), dtype=int)
                
    #             # logging.getLogger().debug(distance_matrix)
    #             # hierarchical_linkage2   = linkage(squareform(distance_matrix), 'single')
    #             # all_c = cut_tree(hierarchical_linkage2, height = distance_max)
    #             # if len(partition) != len(all_c):
    #             #     logging.getLogger().warning(str(len(partition))+" "+str(len(all_c)))
    #                 # import time
    #                 # plt.figure(figsize=(25, 10))
    #                 # plt.title('Hierarchical Clustering Dendrogram (part), '+str(nb_copy))
    #                 # plt.xlabel('sample index')
    #                 # plt.ylabel('distance')
    #                 # txt = str(time.time())
    #                 # dendrogram(
    #                 # hierarchical_linkage,
    #                 # leaf_rotation=90., leaf_font_size=8.,)
    #                 # plt.savefig(txt+"part.png")
    #                 # plt.close()

    #                 # plt.figure(figsize=(25, 10))
    #                 # plt.title('Hierarchical Clustering Dendrogram (all), '+str(nb_copy))
    #                 # plt.xlabel('sample index')
    #                 # plt.ylabel('distance')
    #                 # dendrogram(
    #                 # hierarchical_linkage2,
    #                 # leaf_rotation=90., leaf_font_size=8.,)
    #                 # plt.savefig(txt+"part-all.png")
    #                 # plt.close()

    #     # def merge_overlapping_context(data):
    #     #     if len(data)>0:
    #     #         sets = (set(e) for e in data if e)
    #     #         results = [next(sets)]
    #     #         for e_set in sets:
    #     #             to_update = []
    #     #             for i,res in enumerate(results):
    #     #                 if not e_set.isdisjoint(res):
    #     #                     to_update.insert(0,i)
    #     #             if not to_update:
    #     #                 results.append(e_set)
    #     #             else:
    #     #                 last = results[to_update.pop(-1)]
    #     #                 for i in to_update:
    #     #                     last |= results[i]
    #     #                     del results[i]
    #     #                 last |= e_set
    #     #         return results
    #     #     else:
    #     #         return []

    #     logging.getLogger().debug(max(multi_copy.values()))
    #     logging.getLogger().debug(sum(multi_copy.values()))
    #     logging.getLogger().debug(len(multi_copy.values()))
    #     logging.getLogger().debug(multi_copy)
    #     #while len(resolved)!=len(multi_copy):
    #     for fam, nb_copy in sorted(multi_copy.items(), key=operator.itemgetter(1), reverse=True):
    #             contexts = defaultdict(set)
    #             for org, genes in self.neighbors_graph.node[fam].items():
    #                 for gene in genes:
    #                     if org not in ["name","length_min","length_max","length_avg","label"]:
    #                         contexts[self.__extract_families_of_genomic_context(gene)].add(gene)

    #             isolated_genes = None
    #             try:
    #                 isolated_genes=contexts.pop(tuple())
    #             except:
    #                 pass

    #             groups = __group_similar_contexts(contexts, context_size, 3)

    #             for id, (group, genes) in enumerate(groups.items()):
    #                 new_fam = fam+"#"+str(id)
    #                 for gene in genes:
    #                     (org, seq, index) = self.gene_location[gene]
                        
    #                     name  = self.annotations[org][seq][index][NAME]
    #                     self.__add_gene(new_fam, org, gene, name,  None)
    #                     if seq in self.circular_contig:
    #                         next_index = (index+1)%len(self.annotations[org][seq])
    #                         prev_index = (index-1)%len(self.annotations[org][seq])
    #                         if next_index != index:
    #                             new_nei1 = self.annotations[org][seq][(index+1)%len(self.annotations[org][seq])][FAMILY]
    #                             self.__add_link(new_fam, new_nei1, org)
    #                             new_nei2 = self.annotations[org][seq][(index-1)%len(self.annotations[org][seq])][FAMILY]
    #                             self.__add_link(new_fam, new_nei2, org)
    #                     else:
    #                         new_nei1=None
    #                         new_nei2=None
    #                         logging.getLogger().debug(str(org)+"|"+str(seq)+"|"+str(index))
    #                         try:
    #                             new_nei1 = self.annotations[org][seq][(index+1)][FAMILY]
    #                             # logging.getLogger().debug("1(len(seq)="+str(len(seq))+")")
    #                         except IndexError:
    #                              logging.getLogger().warning("No next neighbors for gene: "+gene+" (new family: "+new_fam+") (len(seq)="+str(len(self.annotations[org][seq]))+")")
    #                         if new_nei1 is not None:
    #                             self.__add_link(new_fam, new_nei1, org)
    #                         try:
    #                             new_nei2 = self.annotations[org][seq][(index-1)][FAMILY]
    #                             # logging.getLogger().debug("2(len(seq)="+str(len(seq))+")")
    #                         except IndexError:
    #                             logging.getLogger().warning("No previous neighbors for gene: "+gene+" (new family: "+new_fam+") (len(seq)="+str(len(self.annotations[org][seq]))+")")
    #                         if new_nei2 is not None:
    #                             logging.getLogger().debug("works "+str(len(self.annotations[org][seq]))+")")
    #                             self.__add_link(new_fam, new_nei1, org)
    #                     self.annotations[org][seq][index][FAMILY] = new_fam

    #             if isolated_genes is not None:
    #                 new_fam = fam+"#"+str(id+1)                    

    #             self.neighbors_graph.remove_node(fam)
    #             # contexts_size[fam]=len(contexts)
    #             # logging.getLogger().debug(fam)
    #             # logging.getLogger().debug(len(contexts))
    #             # logging.getLogger().debug(contexts)
                
    #             # ovl_contexts = merge_overlapping_context(contexts)
    #             # print(ovl_contexts)

    #             # ovl_contexts_size[fam]=len(ovl_contexts)
    #             # logging.getLogger().debug(len(ovl_contexts))
    #             # logging.getLogger().debug(ovl_contexts)

    #             # cpt = 0
    #             # for con in ovl_contexts:
    #             #     new_fam = fam+"#"+str(cpt)
    #             #     for context, genes in contexts.items():
    #             #         if len(context & con)>0:
    #             #             logging.getLogger().debug(cpt)
    #             #             logging.getLogger().debug(context)
    #             #             logging.getLogger().debug(genes)
                            
    #             #             for gene in genes:
    #             #                 logging.getLogger().debug(gene)
    #             #                 (org, seq, index) = self.gene_location[gene]
    #             #                 name  = self.annotations[org][seq][index][NAME]
    #             #                 self.__add_gene(new_fam, org, gene, name,  None)
    #             #                 for neighbor in context: 
    #             #                     self.__add_link(new_fam, neighbor, org)
    #             #                 self.annotations[org][seq][index][FAMILY] = new_fam
    #             #     cpt+=1
    #             # #resolved.add(context)
    #             # logging.getLogger().debug(nx.number_of_nodes(self.neighbors_graph)) 

    #     #         # self.neighbors_graph.remove_node(fam)
    #     # logging.getLogger().debug(max(contexts_size.values()))
    #     # logging.getLogger().debug(sum(contexts_size.values()))
    #     # logging.getLogger().debug(len(contexts_size.values()))
    #     # logging.getLogger().debug(max(ovl_contexts_size.values()))
    #     # logging.getLogger().debug(sum(ovl_contexts_size.values()))
    #     # logging.getLogger().debug(len(ovl_contexts_size.values()))
    #     # logging.getLogger().debug(ovl_contexts_size)
    #     # logging.getLogger().debug(contexts_size)
    
    def neighborhood_computation(self):#, untangle_multi_copy_families = False
        """ Use the information already loaded (annotation) to build the pangenome graph """
        if self.neighbors_graph is None:
            self.neighbors_graph = nx.Graph()
        elif self.is_partitionned:
            raise Exception("The pangenome graph is already built and partionned, please use the function delete pangenome graph before build it again")

        #multi_copy = defaultdict(int) if untangle_multi_copy_families else None

        for organism, annot_contigs in self.annotations.items():
            for contig, contig_annot in annot_contigs.items():
                at_least_2_families = False
                start = 0
                while (start < len(contig_annot) and contig_annot[start][FAMILY] in self.families_repeted):
                    start += 1
                if start == len(contig_annot):
                    continue

                family_id_nei  = contig_annot[start][FAMILY]
                end_family_nei = contig_annot[start][END]
                logging.getLogger().debug(contig_annot[start])
                for index, gene_row in enumerate(contig_annot[start+1:]):
                    logging.getLogger().debug(gene_row)
                    if gene_row[FAMILY] not in self.families_repeted:

                        self.__add_gene(gene_row[FAMILY],
                                        organism,
                                        gene_row[GENE],
                                        gene_row[NAME],
                                        gene_row[END]-gene_row[START],
                                        gene_row[PRODUCT])
                                        #multi_copy)
                        self.neighbors_graph.add_node(family_id_nei)
                        self.__add_link(gene_row[FAMILY],family_id_nei,organism, end_family_nei-gene_row[STRAT])
                        family_id_nei  = gene_row[FAMILY]
                        end_family_nei = gene_row[END]
                        at_least_2_families = True
                
                if contig in self.circular_contig_size and at_least_2_families:#circularization
                    self.__add_gene(contig_annot[start][FAMILY],
                                    organism,
                                    contig_annot[start][GENE],
                                    contig_annot[start][NAME],
                                    contig_annot[start][END]-contig_annot[start][START],
                                    contig_annot[start][PRODUCT])
                                    #multi_copy)
                    self.neighbors_graph.add_node(family_id_nei)
                    self.__add_link(contig_annot[start][FAMILY],family_id_nei,organism, (self.circular_contig_size - end_family_nei) + contig_annot[start][START])
                else:#no circularization
                    self.__add_gene(contig_annot[start][FAMILY],
                                    organism,
                                    contig_annot[start][GENE],
                                    contig_annot[start][NAME],
                                    contig_annot[start][END]-contig_annot[start][START],
                                    contig_annot[start][PRODUCT])
                                    #multi_copy)

        # if multi_copy is not None and len(multi_copy)>0:
        #     logging.getLogger().info("Untangleling ...")
        #     self.__untangle_multi_copy_families(multi_copy)

        self.pan_size = nx.number_of_nodes(self.neighbors_graph)

    def partition(self, nem_dir_path = tempfile.mkdtemp(), beta = 1.00):
        """
            Use the graph topology and the presence or absence of genes from each organism into families to partition the pangenome in three groups ('persistent', 'shell' and 'cloud')
            . seealso:: Read the Mo Dang's thesis to understand NEM and Bernouilli Mixture Model, a summary is available here : http://www.kybernetika.cz/content/1998/4/393/paper.pdf
            
            :param nem_dir_path: a str containing a path to store tempory file of the NEM program
            :param beta: a float containing the spatial coefficient of smothing of the clustering results using the weighted graph topology (0.00 turn off the spatial clustering)
            :type str: 
            :type float: 

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
                for neighbor in nx.all_neighbors(self.neighbors_graph, node_name):
                    #nb_presences = sum([pre_abs for org, pre_abs in self.neighbors_graph[node_name][neighbor].items() if org not in RESERVED_WORDS])
                    #self.neighbors_graph[node_name][neighbor]["weight"]= nb_presences
                    distance_score = self.neighbors_graph[node_name][neighbor]["weight"]/self.nb_organisms
                    row_fam.append(str(index[neighbor]))
                    row_dist_score.append(str(round(distance_score,4)))
                    neighbor_number += 1
                if neighbor_number>0:
                    nei_file.write("\t".join([str(item) for sublist in [[index[node_name]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
                    #logging.getLogger().debug("\t".join([str(item) for sublist in [[[index[node_name]]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
                else:
                    raise nx.exception.NetworkXError("no neighbors in selected organismss")
            except nx.exception.NetworkXError as nxe:
                logging.getLogger().warning("The family: "+node_name+" is an isolated family")
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

        K              = 3 # number of partitions
        ALGO           = "ncem" #fuzzy classification by mean field approximation
        ITERMAX        = 100 # number of iteration max 
        MODEL          = "bern" # multivariate Bernoulli mixture model
        PROPORTION     = "pk" #equal proportion :  "p_"     varying proportion : "pk"
        VARIANCE_MODEL = "sk_" #one variance per partition and organism : "sdk"      one variance per partition, same in all organisms : "sd_"   one variance per organism, same in all partion : "s_d"    same variance in organisms and partitions : "s__" 
        NEIGHBOUR_SPEC = "f"# "f" specify to use all neighbors, orther argument is "4" to specify to use only the 4 neighbors with the higher weight (4 because for historic reason due to the 4 pixel neighbors of each pixel)
        CONVERGENCE_TH = "clas "+str(0.000001)

        HEURISTIC      = "heu_d"# "psgrad" = pseudo-likelihood gradient ascent, "heu_d" = heuristic using drop of fuzzy within cluster inertia, "heu_l" = heuristic using drop of mixture likelihood
        STEP_HEURISTIC = 0.25 # step of beta increase
        BETA_MAX       = float(self.nb_organisms) #maximal value of beta to test,
        DDROP          = 0.8 #threshold of allowed D drop (higher = less detection)
        DLOSS          = 0.5 #threshold of allowed D loss (higher = less detection)
        LLOSS          = 0.02 #threshold of allowed L loss (higher = less detection)
        
        BETA           = ["-B",HEURISTIC,"-H",str(STEP_HEURISTIC),str(BETA_MAX),str(DDROP),str(DLOSS),str(LLOSS)] if beta < 0.00 else ["-b "+str(beta)]

        command = " ".join([NEM_LOCATION, 
                            nem_dir_path+"/nem_file",
                            str(K),
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
            logging.getLogger().error("No NEM output file found")

        classification = "Undefined" * self.pan_size
        try:
            with open(nem_dir_path+"/nem_file.uf","r") as classification_nem_file, open(nem_dir_path+"/nem_file.mf","r") as parameter_nem_file:

                sum_mu_k = []
                sum_epsilon_k = []
                proportion = []
                
                for k, line in enumerate(parameter[-3:]):
                    logging.getLogger().debug(line)
                    vector = line.split()
                    mu_k = [bool(int(mu_kj)) for mu_kj in vector[0:self.nb_organisms]]
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
                partition[persistent_k] = "Persistent"
                partition[shell_k]      = "Shell"
                partition[cloud_k]      = "Cloud"

                for i, line in enumerate(classification_nem_file):
                    elements = [float(el) for el in line.split()]
                    max_prob = max([float(el) for el in elements])
                    positions_max_prob = [pos for pos, prob in enumerate(elements) if prob == max_prob]
                    logging.getLogger().debug(classes)
                    logging.getLogger().debug(i)

                    if (len(positions_max_prob)>1):
                        classification[i]="Shell"# in case of doubt (equiprobable partition), gene families is attributed to shell
                    else:
                        classification[i] = partition[positions_max_prob.pop()]

                parameter = parameter_nem_file.readlines()

                M = float(parameter[6].split()[3]) # M is markov ps-like
                self.BIC = -2 * M - (K * self.nb_organisms * 2 + K - 1) * math.log(self.pan_size)

                logging.getLogger().info("The Bayesian Criterion Index of the partionning is "+str(self.BIC))
                logging.getLogger().debug(partition)
                #logging.getLogger().debug(index.keys())
        except FileNotFoundError:
            logging.getLogger().warning("The number of organisms is to low to partition the pangenome graph in persistent, shell, cloud partition, traditional partitions only (Core and Accessory genome) will be provided")

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
                self.partitions["Core_Exact"].append(node)
                self.neighbors_graph.node[node]["partition_exact"]="Core_Exact"
            else:
                self.partitions["Accessory"].append(node)
                self.neighbors_graph.node[node]["partition_exact"]="Accessory"

        for node_i, node_j, data in self.neighbors_graph.edges(data = True):
            self.neighbors_graph[node_i][node_j]["length_avg"] = float(np.mean(data["length"]))
            self.neighbors_graph[node_i][node_j]["length_med"] = float(np.median(data["length"]))
            self.neighbors_graph[node_i][node_j]["length_min"] = min(data["length"])
            self.neighbors_graph[node_i][node_j]["length_max"] = max(data["length"])

        if len(self.families_repeted)>0:
            logging.getLogger().info("Discarded families are:\n"+"\n".join(self.families_repeted))
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

    def delete_pangenome_graph(self):
        """
            Delete all the pangenome graph of eventuelly the statistic of the partionning process (including the temporary file)
        """ 
        self.delete_nem_intermediate_files()
        self.nem_output               = None
        self.neighbors_graph          = None
        self.pan_size                 = 0
        self.nem_output               = None
        self.is_partitionned          = False
        self.partitions               = {}
        self.partitions["Persistent"] = list()
        self.partitions["Shell"]      = list()
        self.partitions["Cloud"]      = list()
        self.partitions["Core_Exact"] = list()
        self.partitions["Accessory"]  = list()
        self.BIC                      = None
        #self.partitions_by_organisms  = defaultdict(lambda: defaultdict(set))

    def delete_nem_intermediate_files(self):
        """
            Delete all the tempory files used to partion the pangenome
        """ 
        logging.getLogger().info("delete "+self.nem_intermediate_files)
        if self.nem_intermediate_files is not None:
            shutil.rmtree(self.nem_intermediate_files)
            self.nem_intermediate_files = None

    def identify_communities_in_each_partition(self):
        """
            Use the Louvain's algorithm to label the nodes by their community in each partition
        """ 
        size_communities=defaultdict(lambda : defaultdict(int))
        for partition in ["Persistent","Shell", "Cloud"]:
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
                logging.getLogger().warning("The family: "+node_name+" is an isolated family")
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


if __name__=='__main__':
    """
    --organims is a tab delimited files containg at least 2 mandatory fields by row and as many optional field as circular contig. Each row correspond to an organism to be added to the pangenome.
    Reserved words are : "id", "label", "name", "weight", "partition", "partition_exact"
    The first field is the organinsm name (id should be unique not contain any spaces, " or ' and reserved words).
    The second field is the gff file associated to the organism. This path can be abolute or relative. The gff file must contain an id feature for each CDS (id should be unique not contain any spaces, " or ' and reserved words).
    The next fields contain the name of perfectly assemble circular contigs (contig name must should be unique and not contain any spaces, " or ' and reserved words).
    example:

    """
    parser = argparse.ArgumentParser(description='Build a partitioned pangenome graph from annotated genomes and gene families')
    parser.add_argument('-?', '--version', action='version', version='0.1')
    parser.add_argument('-o', '--organisms', type=argparse.FileType('r'), nargs=1, help="""
A tab delimited file containg at least 2 mandatory fields by row and as many optional fields as the number of well assembled circular contigs. Each row correspond to an organism to be added to the pangenome.
Reserved words are : "id", "label", "name", "weight", "partition", "partition_exact"
The first field is the organinsm name (id should be unique not contain any spaces, " or ' and reserved words).
The second field is the gff file associated to the organism. This path can be abolute or relative. The gff file must contain an id feature for each CDS (id should be unique across all pangenome and not contain any spaces, " or ' and reserved words).
The next fields contain the name of perfectly assembled circular contigs.
Contig names should be unique and not contain any spaces, quote, double quote and reserved words.
example:""", required=True)
    parser.add_argument('-f', '--gene_families', type=argparse.FileType('r'), nargs=1, help="""
A tab delimited file containg the gene families. Each row contain 2 fields.
Reserved words are : "id", "label", "name", "weight", "partition", "partition_exact"
The first field is the family name.
The second field is the gene name.
families are intended to be grouped by chuncks of row.
the families name can be any string but must should be unique and not contain any spaces, " or ' and reserved words
As a conventation, it is recommanded to use the name of the most reprensative gene of the families as the family name.
Gene name can be any string corresponding to the if feature in the gff files. they should be unique and not contain any spaces, " or ' and reserved words.
example:""",  required=True)
    parser.add_argument('-d', '--output_directory', type=str, nargs=1, default=["output.dir"], help="""
The output directory""")
    parser.add_argument('-ff', '--force', action="store_true", help="""
Force overwriting existing output directory""")
    parser.add_argument('-r', '--remove_high_copy_number_families', type=int, nargs=1, default=[0], help="""
Remove families having a number of copy of one families above or equal to this threshold in at least one organism (0 or negative value keep all families whatever their occurence). 
When -u is set, only work on new organisms added""")
    parser.add_argument('-s', '--infere_singleton', default=False, action="store_true", help="""
If a gene id found in a gff file is absent of the gene families file, the singleton will be automatically infered as a gene families having a single element. 
if this argument is not set, the program will raise KeyError exception if a gene id found in a gff file is absent of the gene families file.""")
    parser.add_argument("-u", "--update", default = None, type=argparse.FileType('r'), nargs=1, help="""
Pangenome Graph to be updated (in gexf format)""")
    parser.add_argument("-b", "--beta_smoothing", default = [-1.00], type=float, nargs=1, help = """
Coeficient of smoothing all the partionning based on the Markov Random Feild leveraging the weigthed pangenome graph. A positive float, 0.0 means to discard spatial smoothing and -1 to find beta automaticaly (increase the computation time) 
""")
    parser.add_argument("-i", "--delete_nem_intermediate_files", default=False, action="store_true", help="""
Delete intermediate files used by NEM""")
    parser.add_argument("-c", "--compress_graph", default=False, action="store_true", help="""
Compress (using gzip) the file containing the partionned pangenome graph""")
    parser.add_argument("-v", "--verbose", default=True, action="store_true", help="""
Show information message, otherwise only errors and warnings will be displayed""")
    parser.add_argument("-vv", "--verbose_debug", default=False, action="store_true", help="""
Show all messages including debug ones""")

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
    for directory in [OUTPUTDIR, NEMOUTPUTDIR, FIGUREOUTPUTDIR]:
        if not os.path.exists(directory):
            os.makedirs(directory)
        elif not options.force:
            logging.getLogger().error(directory+" already exist")
            exit()

    start_loading_file = time.time()

    pan = PPanGGOLiN("file",
                     options.organisms[0],
                     options.gene_families[0],
                     options.remove_high_copy_number_families[0],
                     options.infere_singleton)

    if options.update is not None:
        pan.import_from_GEXF(options.update[0])

    start_neighborhood_computation = time.time()
    pan.neighborhood_computation()
    start_partitioning = time.time()
    
    pan.partition(NEMOUTPUTDIR, beta = options.beta_smoothing[0])
    start_identify_communities = time.time()
    pan.identify_communities_in_each_partition()
    pan.identify_shell_subpaths()


    time_of_writing_output_file = time.time()
    GEXF_GRAPH_FILE  = OUTPUTDIR+"/"+"graph.gexf"
    if options.compress_graph:
        pan.export_to_GEXF(GEXF_GRAPH_FILE+".gz", compressed=True)
    else:
        pan.export_to_GEXF(GEXF_GRAPH_FILE, compressed=False)

    for filename, families in pan.partitions.items(): 
        file = open(OUTPUTDIR+"/"+filename+".txt","w")
        file.write("\n".join(families))
        file.close()

    pan.csv_matrix(OUTPUTDIR+"/matrix.csv")

    logging.getLogger().info(pan)
    logging.getLogger().info("\n\
    Execution time of file loading: " +str(round(start_neighborhood_computation-start_loading_file, 2))+" s\n"+
    "Execution time of neighborhood computation: " +str(round(start_partitioning-start_neighborhood_computation, 2))+" s\n"+
    "Execution time of partitioning: " +str(round(start_identify_communities-start_partitioning, 2))+" s\n"+
    "Execution time of community identification: " +str(round(time_of_writing_output_file-start_identify_communities, 2))+" s\n"+
    "Execution time of writing output files: " +str(round(time.time()-time_of_writing_output_file, 2))+" s\n"+
    "Total execution time: " +str(round(time.time()-start_loading_file, 2))+" s\n")

    # print(pan.partitions_by_organisms)
    # partitions_by_organisms_file = open(OUTPUTDIR+"/partitions_by_organisms.txt","w")
    # exact_by_organisms_file = open(OUTPUTDIR+"/exacte_by_organisms.txt","w")
    # for org, partitions in pan.partitions_by_organisms.items(): 
    #     partitions_by_organisms_file.write(org+"\t"+str(len(partitions["Persistent"]))+
    #                                            "\t"+str(len(partitions["Shell"]))+
    #                                            "\t"+str(len(partitions["Cloud"]))+"\n")
    #     exact_by_organisms_file.write(org+"\t"+str(len(partitions["Core_Exact"]))+
    #                                       "\t"+str(len(partitions["Accessory"]))+"\n")
    # partitions_by_organisms_file.close()
    # exact_by_organisms_file.close()

    if options.delete_nem_intermediate_files:
        pan.delete_pangenome_graph()

    logging.getLogger().info("Finished !")

    #####################################

    # start_size = pan.nb_organisms

    # with open(OUTPUTDIR+"/stat_evol.txt","w") as evol, open(OUTPUTDIR+"/stat_evol_exact.txt","w") as evol_exact:
    #     evol.write("\t".join([str(pan.nb_organisms),
    #                               str(len(pan.partitions["Persistent"])),
    #                               str(len(pan.partitions["Shell"])),
    #                               str(len(pan.partitions["Cloud"])),
    #                               str(pan.pan_size)])+"\n")
    #     evol_exact.write("\t".join([str(pan.nb_organisms),
    #                                   str(len(pan.partitions["Core_Exact"])),
    #                                   str(len(pan.partitions["Accessory"])),
    #                                   str(pan.pan_size)])+"\n")
    #     pan.delete_pangenome_graph()
    #     while pan.nb_organisms>2:
    #         if ((pan.nb_organisms%10)==0):
    #             pan.neighborhood_computation()
    #             pan.partition(OUTPUTDIR+"/"+str(pan.nb_organisms))
    #             evol.write("\t".join([str(pan.nb_organisms),
    #                               str(len(pan.partitions["Persistent"])),
    #                               str(len(pan.partitions["Shell"])),
    #                               str(len(pan.partitions["Cloud"])),
    #                               str(pan.pan_size)])+"\n")
    #             evol_exact.write("\t".join([str(pan.nb_organisms),
    #                                   str(len(pan.partitions["Core_Exact"])),
    #                                   str(len(pan.partitions["Accessory"])),
    #                                   str(pan.pan_size)])+"\n")
    #             evol.flush()
    #             evol_exact.flush()
    #             pan.nem_intermediate_files = None
    #         removed = pan.organisms.pop()
    #         pan.annotations.pop(removed, None)
    #         pan.nb_organisms-=1           
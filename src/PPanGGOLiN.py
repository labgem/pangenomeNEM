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
#import forceatlas2 

NEM_LOCATION  = "../NEM/nem_exe"
(GENE, TYPE, FAMILLY, START, END, STRAND, NAME) = range(0, 7)#data index in annotation
(ORGANISM, SEQUENCE, INDEX) = range(0, 3)# index for annotation in gene_location

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
                Nodes attributes contains the genes name in each organism belong the gene families of the node.

            .. attribute:: organisms

                an ordored-set contains the imported organisms 

            .. attribute:: nb_organisms

                a int giving the number of imported organisms 

            .. attribute:: circular_contig

                a set containing the name of both well assembled and circular contigs (used to circularize the graph) 

            .. attribute:: pan_size

                The number of nodes into the graph.
                .. warning:: this number is not necesserally equal to the number of imported gene families. Indeed, if imported famillies are not present into organism, superfluous will be discard.

            .. attribute:: is_partitionned
            
                a boolean specifying if the pangenome graph has been partitionned or not

            .. attribute:: nem_intermediate_files

                a str storing the path to the nem nem intermediate_files

            .. attribute:: partitions

                a dict provoding the families present in each partion:
                    * partitions["Core_Exact"] contains the list of core families (core exact)
                    * partitions["Persistent"] contains the list of persistent families
                    * partitions["Shell"] contains the list of shell families
                    * partitions["Cloud"] contains the list of cloud families

            .. attribute:: BIC

                a float providing the Bayesian Information Criterion. This Criterion give an estimation of the quality of the partionning (a low value means a good one)
                . seealso:: https://en.wikipedia.org/wiki/Bayesian_information_criterion

            .. attribute:: families_repeted_th

                a int specifying the threshold of presence in a single organism above which families are considered as repeted

            .. attribute:: families_repeted

                a list of str storing the repeted families

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
            :type *args: list of argument 

            :Example:

            >>>pan = PPanGGOLiN("file", organisms, gene_families, remove_high_copy_number_families)
            >>>pan = PPanGGOLiN("args", annotations, organisms, circular_contig, families_repeted)
        """ 
        self.annotations              = dict()
        self.gene_location            = OrderedDict()
        self.neighbors_graph          = None
        self.organisms                = OrderedSet()
        self.nb_organisms             = 0
        self.circular_contig          = set()
        self.families_repeted_th      = 0
        self.families_repeted         = set()
        self.pan_size                 = 0
        self.nem_intermediate_files   = None
        self.is_partitionned          = False
        self.familly_frequencies      = {}
        self.partitions               = {}
        self.partitions["Persistent"] = list()
        self.partitions["Shell"]      = list()
        self.partitions["Cloud"]      = list()
        self.partitions["Core_Exact"] = list()
        self.partitions["Accessory"]  = list()
        self.BIC                      = None # Bayesian Index Criterion
        self.partitions_by_organisms  = defaultdict(lambda: defaultdict(int))

        if init_from == "file":
            self.__initialize_from_files(*args)
        elif init_from == "args":
            (self.annotations,
             self.organisms,
             self.circular_contig,
             self.families_repeted) = args 
        elif init_from == "database":
            logging.getLogger().error("database is not yet implemented")
            pass
        else:
            raise ValueError("init_from parameter is required")
        self.nb_organisms = len(self.organisms)

    def __initialize_from_files(self, organisms_file, families_tsv_file, lim_occurence, infere_singleton = False):
        """ 
            :param organisms_file: 
            :param families_tsv_file: 
            :param lim_occurence: 
            :type file: 
            :type file: 
            :type int: 

            :Example:

            >>>pan.__initialize_from_files(organisms_file, families_tsv_file, lim_occurence)

The tsv file provided by progenome containing the gene annotations. --organims is a tab delimited files containg at least 2 mandatory fields by row and as many optional field as circular contig. Each row correspond to an organism to be added to the pangenome.
Reserved words are : "id", "label", "name", "weight", "partition", "partition_exact"
The first field is the organinsm name (id should be unique not contain any spaces, " or ' and reserved words).
The seconde field is the gff file associated to the organism. This path can be abolute or relative. The gff file must contain an id feature for each CDS (id should be unique not contain any spaces, " or ' and reserved words).
The next fields contain the name of perfectly assemble circular contigs (contig name must should be unique and not contain any spaces, " or ' and reserved words).
example:

        """ 
        logging.getLogger().info("Reading "+families_tsv_file.name+" families file ...")
        families    = dict()
        first_iter  = True
        for line in families_tsv_file:
            elements = line.split()
            for gene in elements[1:]:
                families[gene]=elements[0]

        self.circular_contig = []

        logging.getLogger().info("Reading "+organisms_file.name+" families file ...")

        for line in organisms_file:
            elements = line.split()
            self.organisms.add(elements[0])
            self.annotations[elements[0]] = self.__load_gff(elements[1], families, elements[0], lim_occurence, infere_singleton)
            if len(elements)>2:
                self.circular_contig += elements[2:len(elements)]

        self.circular_contig = set(self.circular_contig)

    def __load_gff(self, gff_file, families, organism, lim_occurence, infere_singleton = False):
        """
            Load the content of a gff file
            :param gff_file: 
            :param families: 
            :param organism: 
            :param lim_occurence: 
            :type str: 
            :type dict: 
            :type str: 
            :type int: 
            :return: annot: 
            :rtype: dict 
        """ 
        logging.getLogger().info("Reading "+gff_file+" file ...")
        db_gff = gffutils.create_db(gff_file, ':memory:')
        annot = defaultdict(list)

        # prev = None
        ctp_prev = 1

        cpt_fam_occ = defaultdict(int)

        for row in db_gff.all_features(featuretype='CDS', order_by=('seqid','start')):
            logging.getLogger().debug(row)
            protein = row.id
            try:
                familly = families[protein]
            except KeyError:
                if infere_singleton:
                    families[protein] = protein
                    familly           = families[protein]
                    logging.getLogger().info("infered singleton: "+protein)
                else:
                    logging.getLogger().error("unknown families:"+protein, ", check your families file or run again the program using the option to infere singleton")
                    raise KeyError("Unknown families while infere singleton is false")

            # if familly == prev:
            #     familly = familly+"-"+str(ctp_prev)

            cpt_fam_occ[familly]+=1
            prev = families[protein]

            try:
                name = row.attributes['Name'].pop()
            except:
                try:
                    name = row.attributes['gene'].pop()
                except:
                    name = ""

            annot_row = [protein,"CDS",familly,row.start,row.end,row.strand, name]
            self.gene_location[protein] = tuple([organism, row.seqid, len(annot[row.seqid])])
            annot[row.seqid].append(annot_row)

            logging.getLogger().debug(annot[self.gene_location[protein][SEQUENCE]][self.gene_location[protein][INDEX]][GENE])

        if (lim_occurence >0):
            fam_to_remove =[fam for fam, occ in cpt_fam_occ.items() if occ > lim_occurence]
            self.families_repeted = self.families_repeted.union(set(fam_to_remove))

        return(annot)

    def __str__(self):
        """ Return an overview of the object as a string """ 
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

    def __add_gene(self, fam_id, org, gene, name, multi_copy = None):
        """ """
        self.neighbors_graph.add_node(fam_id)
        try:
            self.neighbors_graph.node[fam_id][org].add(gene)
            if multi_copy is not None:
                multi_copy[fam_id] = max(multi_copy[fam_id],len(self.neighbors_graph.node[fam_id][org]))#count max number of multicopy
        except KeyError:# the exception is the rule in case of mono copy gene
            self.neighbors_graph.node[fam_id][org]=set([gene])
        try:
            self.neighbors_graph.node[fam_id]['name'].add(name)
        except KeyError:
            self.neighbors_graph.node[fam_id]['name']=set([name])
    def __add_link(self, fam_id, fam_nei, org):
        """  """
        if not self.neighbors_graph.has_edge(fam_id,fam_nei):
            self.neighbors_graph.add_edge(fam_id, fam_nei)
            logging.getLogger().debug([str(i) for i in [fam_id, fam_nei, org]])
        try:
            self.neighbors_graph[fam_id][fam_nei][org]+=1
        except KeyError:
            self.neighbors_graph[fam_id][fam_nei][org]=1
            try:
                self.neighbors_graph[fam_id][fam_nei]["weight"]+=1.0
            except KeyError:
                self.neighbors_graph[fam_id][fam_nei]["weight"]=1.0

    def __extract_families_of_genomic_context(self, gene, context_size = 1):
        ret = set()
        try:
            (organism, sequence, index) = self.gene_location[gene]
            min_boundary = index-context_size if index-context_size else 0
            for annot in self.annotations[organism][sequence][min_boundary:(index+context_size+1)]:
                if annot[GENE]!=gene:
                    ret.add(annot[FAMILLY])            
        except KeyError:
            logging.getLogger().debug("not a gene id: "+str(gene))
        return frozenset(ret)

    def __untangle_multi_copy_families(self, multi_copy):

        def merge_overlapping_context(data):
            if len(data)>0:
                sets = (set(e) for e in data if e)
                results = [next(sets)]
                for e_set in sets:
                    to_update = []
                    for i,res in enumerate(results):
                        if not e_set.isdisjoint(res):
                            to_update.insert(0,i)
                    if not to_update:
                        results.append(e_set)
                    else:
                        last = results[to_update.pop(-1)]
                        for i in to_update:
                            last |= results[i]
                            del results[i]
                        last |= e_set
                return results
            else:
                return []

        resolved=set()
        logging.getLogger().debug(max(multi_copy.values()))
        logging.getLogger().debug(sum(multi_copy.values()))
        logging.getLogger().debug(len(multi_copy.values()))
        logging.getLogger().debug(multi_copy)
        #while len(resolved)!=len(multi_copy):
        contexts_size={}
        ovl_contexts_size={}
        for fam, nb_copy in sorted(multi_copy.items(), key=operator.itemgetter(1)):
            if fam not in resolved:
                contexts = defaultdict(set)
                for org, genes in self.neighbors_graph.node[fam].items():
                    for gene in genes:
                        if org != "name":
                            contexts[self.__extract_families_of_genomic_context(gene)].add(gene)

                #all_neigh = nx.all_neighbors(self.neighbors_graph,fam)

                contexts_size[fam]=len(contexts)
                logging.getLogger().debug(fam)
                logging.getLogger().debug(len(contexts))
                logging.getLogger().debug(contexts)
                
                ovl_contexts = merge_overlapping_context(contexts)
                print(ovl_contexts)

                ovl_contexts_size[fam]=len(ovl_contexts)
                logging.getLogger().debug(len(ovl_contexts))
                logging.getLogger().debug(ovl_contexts)

                cpt = 0
                for con in ovl_contexts:
                    new_fam = fam+"#"+str(cpt)
                    for context, genes in contexts.items():
                        if len(context & con)>0:
                            logging.getLogger().debug(cpt)
                            logging.getLogger().debug(context)
                            logging.getLogger().debug(genes)
                            
                            for gene in genes:
                                logging.getLogger().debug(gene)
                                (org, seq, index) = self.gene_location[gene]
                                name  = self.annotations[org][seq][index][NAME]
                                self.__add_gene(new_fam, org, gene, name,  None)
                                for neighbor in context: 
                                    self.__add_link(new_fam, neighbor, org)
                                self.annotations[org][seq][index][FAMILLY] = new_fam
                    cpt+=1
                #resolved.add(context)
                logging.getLogger().debug(nx.number_of_nodes(self.neighbors_graph)) 

                self.neighbors_graph.remove_node(fam)
        logging.getLogger().debug(max(contexts_size.values()))
        logging.getLogger().debug(sum(contexts_size.values()))
        logging.getLogger().debug(len(contexts_size.values()))
        logging.getLogger().debug(max(ovl_contexts_size.values()))
        logging.getLogger().debug(sum(ovl_contexts_size.values()))
        logging.getLogger().debug(len(ovl_contexts_size.values()))
        logging.getLogger().debug(ovl_contexts_size)
        logging.getLogger().debug(contexts_size)
        
    def neighborhood_computation(self, untangle_multi_copy_families=True):
        """ """
        if self.neighbors_graph is None:
            self.neighbors_graph = nx.Graph()
        elif self.is_partitionned:
            logging.getLogger().error("The pangenome graph is already built and partionned, please use the function delete pangenome graph before build it again")

        multi_copy = defaultdict(int) if untangle_multi_copy_families else None

        for organism, annot_contigs in self.annotations.items():
            for contig, contig_annot in annot_contigs.items():
                at_least_2_families = False
                start = 0
                while (start < len(contig_annot) and contig_annot[start][FAMILLY] in self.families_repeted):
                    start += 1
                if start == len(contig_annot):
                    continue

                familly_id_nei = contig_annot[start][FAMILLY]
                logging.getLogger().debug(contig_annot[start])
                for index, gene_row in enumerate(contig_annot[start+1:]):
                    logging.getLogger().debug(gene_row)
                    if gene_row[FAMILLY] not in self.families_repeted:
                        self.__add_gene(gene_row[FAMILLY],
                                        organism,
                                        gene_row[GENE],
                                        gene_row[NAME],
                                        multi_copy)
                        self.neighbors_graph.add_node(familly_id_nei)
                        self.__add_link(gene_row[FAMILLY],familly_id_nei,organism)
                        familly_id_nei = gene_row[FAMILLY]
                        at_least_2_families = True
                
                if contig in self.circular_contig and at_least_2_families:#circularization
                    self.__add_gene(contig_annot[start][FAMILLY],
                                    organism,
                                    contig_annot[start][GENE],
                                    contig_annot[start][NAME],
                                    multi_copy)
                    self.neighbors_graph.add_node(familly_id_nei)
                    self.__add_link(contig_annot[start][FAMILLY],familly_id_nei,organism)
                else:#no circularization
                    self.__add_gene(contig_annot[start][FAMILLY],
                                    organism,
                                    contig_annot[start][GENE],
                                    contig_annot[start][NAME],
                                    multi_copy)

        if multi_copy is not None and len(multi_copy)>0:
            logging.getLogger().info("Untangleling ...")
            self.__untangle_multi_copy_families(multi_copy)

        self.pan_size = nx.number_of_nodes(self.neighbors_graph)

    def partition(self, nem_dir_path):
        """
            . seealso:: Read the Mo Dang's thesis to understand NEM and Bernouilli Mixture Model, a summary is available here : http://www.kybernetika.cz/content/1998/4/393/paper.pdf
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

        #if use_neighborhood:
        #    nei_file.write("1\n")
        #else:
        nei_file.write("0\n")
        
        index = {node: index+1 for index, node in enumerate(self.neighbors_graph.nodes(data=False))}

        org_file.write("\t".join([org for org in self.organisms])+"\n")
        org_file.close()

        for node in self.neighbors_graph.nodes(data=True):
            node_name, node_organisms = node

            index_file.write(str(index[node_name])+"\t"+str(node_name)+"\n")
            logging.getLogger().debug(node_organisms)
            logging.getLogger().debug(self.organisms)
            dat_file.write("\t".join(["1" if org in node_organisms else "0" for org in self.organisms])+"\n")

            # row_fam         = []
            # row_dist_score  = []
            # neighbor_number = 0
            # try:
            #     for neighbor in nx.all_neighbors(self.neighbors_graph, node_name):
            #         nb_presences = sum([pre_abs for org, pre_abs in self.neighbors_graph[node_name][neighbor].items() if org != 'weight'])
            #         self.neighbors_graph[node_name][neighbor]["weight"]= nb_presences
            #         distance_score = nb_presences/self.nb_organisms
            #         row_fam.append(str(index[neighbor]))
            #         row_dist_score.append(str(round(distance_score,4)))
            #         neighbor_number += 1
            #     if neighbor_number>0:
            #         nei_file.write("\t".join([str(item) for sublist in [[index[node_name]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
            #         #logging.getLogger().debug("\t".join([str(item) for sublist in [[[index[node_name]]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
            #     else:
            #         raise nx.exception.NetworkXError("no all_neighbors in selected organismss")
            # except nx.exception.NetworkXError as nxe:
            #     logging.getLogger().warning("The familly: "+node_name+" is an isolated familly")
            #     nei_file.write(str(index[node_name])+"\t0\n")

        m_file.write("1 0.333 0.333 ") # 1 to initialize parameter, 0.333 and 0.333 for to give one third of initial portition to each class (last 0.33 is automaticaly determined by substraction)
        m_file.write(" ".join(["1"]*self.nb_organisms)+" ") # persistent binary vector
        m_file.write(" ".join(["1"]*self.nb_organisms)+" ") # shell binary vector (1 ou 0, whatever because dispersion will be of 0.5)
        m_file.write(" ".join(["0"]*self.nb_organisms)+" ") # cloud binary vector
        m_file.write(" ".join(["0.01"]*self.nb_organisms)+" ") # persistent dispersition vector (low)
        m_file.write(" ".join(["0.5"]*self.nb_organisms)+" ") # shell dispersition vector (high)
        m_file.write(" ".join(["0.01"]*self.nb_organisms)) # cloud dispersition vector (low)

        index_file.close()
        nei_file.close()
        dat_file.close()
        m_file.close()

        logging.getLogger().info("Running NEM...")

        K              = 3 # number of partitions
        ALGO           = "nem" #fuzzy classification by mean field approximation
        BETA           = 0 # coeficient of spatial smoothing to apply, 0 is equivalent to EM for mixture model
        ITERMAX        = 100 # number of iteration max 
        MODEL          = "bern" # multivariate Bernoulli mixture model
        PROPORTION     = "pk" #equal proportion :  "p_"     varying proportion : "pk"
        VARIANCE_MODEL = "skd" #one variance per partition and organism : "sdk"      one variance per partition, same in all organisms : "sd_"   one variance per organism, same in all partion : "s_d"    same variance in organisms and partitions : "s__" 

        command = " ".join([NEM_LOCATION, 
                            nem_dir_path+"/nem_file",
                            str(K),
                            "-a", ALGO,
                            "-i", str(ITERMAX),
                            "-m", MODEL, PROPORTION, VARIANCE_MODEL,
                            "-s", "m", nem_dir_path+"/nem_file.m",
                            "-b", str(BETA),
                            "-f fuzzy",
                            "-l y -T" if logging.getLogger().getEffectiveLevel() < 20 else ""])
     
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
            for i, line in enumerate(classification_nem_file):
                elements = [float(el) for el in line.split()]
                max_prob = max([float(el) for el in elements])
                classes = [pos for pos, prob in enumerate(elements) if prob == max_prob]
                logging.getLogger().debug(classes)
                logging.getLogger().debug(i)

                if (len(classes)>1):
                    classification.append(2)#shell
                else:
                    if classes[0] == 0:
                        classification.append(1)#persistent
                    elif classes[0] == 2:
                        classification.append(3)#cloud
                    else:
                        classification.append(2)#shell

            parameter = parameter_nem_file.readlines()
            M = float(parameter[6].split()[3]) # M is markov ps-like
            self.BIC = -2 * M - (K * self.nb_organisms * 2 + K - 1) * math.log(self.pan_size)

            logging.getLogger().info("The Bayesian Criterion Index of the partionning is "+str(self.BIC))

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

            max_mu_k = max(sum_mu_k)
            persistent_k = sum_mu_k.index(max_mu_k)

            max_epsilon_k = max(sum_epsilon_k)
            shell_k = sum_epsilon_k.index(max_epsilon_k)

            cloud_k = set([0,1,2]) - set([persistent_k, shell_k])
            cloud_k = list(cloud_k)[0]

            logging.getLogger().debug(sum_mu_k)
            logging.getLogger().debug(sum_epsilon_k)

            logging.getLogger().debug(persistent_k)
            logging.getLogger().debug(shell_k)
            logging.getLogger().debug(cloud_k)

            partition                 = {}
            partition[persistent_k+1] = "Persistent"
            partition[shell_k+1]      = "Shell"
            partition[cloud_k+1]      = "Cloud"

            logging.getLogger().debug(partition)
            #logging.getLogger().debug(index.keys())
            for node, nem_class in zip(self.neighbors_graph.nodes(), classification):
                nb_orgs=0
                for key, items in self.neighbors_graph.node[node].items():

                    self.neighbors_graph.node[node][key]="|".join(items)    
                    if key!="name":
                        self.partitions_by_organisms[key][partition[int(nem_class)]]=+1
                        nb_orgs+=1

                self.neighbors_graph.node[node]["partition"]=partition[int(nem_class)]
                self.partitions[partition[int(nem_class)]].append(node)

                if nb_orgs >= self.nb_organisms:
                    self.partitions["Core_Exact"].append(node)
                    self.neighbors_graph.node[node]["partition_exact"]="Core_Exact"
                else:
                    self.partitions["Accessory"].append(node)
                    self.neighbors_graph.node[node]["partition_exact"]="Accessory"

        if len(self.families_repeted)>0:
            logging.getLogger().info("Discarded families are:\n"+"\n".join(self.families_repeted))
        else:
            logging.getLogger().info("No families have been Discarded")

        logging.getLogger().debug(nx.number_of_edges(self.neighbors_graph))

        self.is_partitionned=True

        #positions = forceatlas2.forceatlas2_networkx_layout(self.neighbors_graph, 
        #                                                    niter=10,
        #                                                    edgeWeightInfluence=0.8)
        # figure = plt.figure()
        # nx.draw_graphviz(self.neighbors_graph, ax=figure.add_subplot(111))#, positions
        # figure.savefig("graph2.png")
        
    def export_to_GEXF(self, graph_output_path, compressed=False):
        """ 

        
            .. warning:: please use the function neighborhoodComputation before please use the function partition before

        """
        if self.neighbors_graph is None:
            logging.getLogger().error("neighbors_graph is not built, please use the function neighborhoodComputation before")
        elif not self.is_partitionned:
            logging.getLogger().warnings("the pangenome graph is not partionned, please use the function partition before")
        else:
            logging.getLogger().info("Writing GEXF file")
            if compressed:
                graph_output_path = gzip.open(graph_output_path,"w")
            nx.write_gexf(self.neighbors_graph, graph_output_path)

    def import_from_GEXF(self, path_graph_to_update):
        """ """
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
                if key not in ["id", "name", "label", "weight", "partition", "partition_exact"]:
                    self.organisms.add(key)
        logging.getLogger().debug(self.organisms)
        self.nb_organisms = self.nb_organisms = len(self.organisms)
        for source, target, data in self.neighbors_graph.edges(data=True):
            try:
                del self.neighbors_graph[source][target]['id']
            except KeyError:
                logging.getLogger().warnings("No previous edge id found in gexf input file for edge: "+source+" <-> "+target)

    def delete_pangenome_graph(self):
        self.delete_nem_intermediate_files()
        self.nem_output               = None
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

    def delete_nem_intermediate_files(self):
        logging.getLogger().info("delete "+self.nem_intermediate_files)
        shutil.rmtree(self.nem_intermediate_files)



if __name__=='__main__':
#blabla-1
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
The first field is the familly name.
The second field is the gene name.
famillies are intended to be grouped by chuncks of row.
the families name can be any string but must should be unique and not contain any spaces, " or ' and reserved words
As a conventation, it is recommanded to use the name of the most reprensative gene of the famillies as the familly name.
Gene name can be any string corresponding to the if feature in the gff files. they should be unique and not contain any spaces, " or ' and reserved words.
example:""",  required=True)
    parser.add_argument('-d', '--output_directory', type=str, nargs=1, default="output.dir", help="""
The output directory""")
    parser.add_argument('-r', '--remove_high_copy_number_families', type=int, nargs=1, default=[-1], help="""
Remove families having a number of copy of one families above or equal to this threshold in at least one organism (0 or negative value keep all families whatever their occurence). 
When -u is set, only work on new organisms added""")
    parser.add_argument('-s', '--infere_singleton', default=False, action="store_true", help="""
If a gene id found in a gff file is absent of the gene families file, the singleton will be automatically infered as a gene families having a single element. 
if this argument is not set, the program will raise KeyError exception if a gene id found in a gff file is absent of the gene families file.""")
    parser.add_argument("-u", "--update", default = None, type=argparse.FileType('r'), nargs=1, help="""
Pangenome Graph to be updated (in gexf format)""")
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

    OUTPUTDIR       = options.output_directory[0]
    NEMOUTPUTDIR    = OUTPUTDIR+"/NEM_results/"
    FIGUREOUTPUTDIR = OUTPUTDIR+"/figures/"
    for directory in [OUTPUTDIR, NEMOUTPUTDIR, FIGUREOUTPUTDIR]:
        if not os.path.exists(directory):
            os.makedirs(directory)
        else:
            logging.getLogger().error(directory+" already exist")
            exit()

    pan = PPanGGOLiN("file",
                     options.organisms[0],
                     options.gene_families[0],
                     options.remove_high_copy_number_families[0],
                     options.infere_singleton)

    if options.update is not None:
        pan.import_from_GEXF(options.update[0])

    pan.neighborhood_computation()
    pan.partition(NEMOUTPUTDIR)
    logging.getLogger().info(pan)

    GEXF_GRAPH_FILE  = OUTPUTDIR+"/"+"graph.gexf"

    if options.compress_graph:
        pan.export_to_GEXF(GEXF_GRAPH_FILE+".gz", compressed=True)
    else:
        pan.export_to_GEXF(GEXF_GRAPH_FILE, compressed=False)

    for filename, families in pan.partitions.items(): 
        file = open(OUTPUTDIR+"/"+filename+".txt","w")
        file.write("\n".join(families))
        file.close()

    if options.delete_nem_intermediate_files:
        pan.delete_pangenome_graph()

    logging.getLogger().info("Finished !")
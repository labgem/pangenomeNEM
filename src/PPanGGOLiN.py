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
#import forceatlas2 

NEM_LOCATION  = "../NEM/nem_exe"
(GENE, TYPE, ORGANISM, FAMILLY_CODE, START, END, STRAND, NAME) = range(0, 8)


class PPanGGOLiN:
    """  """ 
    def __init__(self, init_from = "args", *args):
        """  """ 
        self.annotations          = dict()
        self.neighbors_graph      = None
        self.organisms            = OrderedSet()
        self.nb_organisms         = 0
        self.circular_contig      = set()
        self.pan_size             = 0
        self.intermediate_file    = None
        self.partitions           = {}
        self.partitions["P"]      = list()
        self.partitions["S"]      = list()
        self.partitions["C"]      = list()
        self.partitions["CE"]     = list()
        self.BIC                  = None # Bayesian Index Criterion
        #self.gene_location        = OrderedDict()
        self.families_repeted     = set()

        if init_from == "file":
            self.__initialize_from_files(*args)
        elif init_from == "args":
            (self.annotations,
             self.organisms,
             self.circular_contig,
             self.families_repeted) = args 
        elif init_from == "database":
            pass
        else:
            raise ValueError("init_from parameter is required")
        self.nb_organisms = len(self.organisms)

    def __initialize_from_files(self, organisms_file, families_tsv_file, lim_occurence):
        """ """
        logging.getLogger().info("Reading "+families_tsv_file.name+" families file ...")
        families    = dict()
        nb_families = 0
        first_iter  = True
        for line in families_tsv_file:
            elements = line.split()
            if first_iter:
                prec = elements[0]
                first_iter = False
            if elements[0] == prec:
                families[elements[1]]=elements[0]
            else :
                prec = elements[0]
                nb_families +=1
                families[elements[1]]=elements[0]
        organisms = []
        self.circular_contig = []

        logging.getLogger().info("Reading "+organisms_file.name+" families file ...")

        for line in organisms_file:
            elements = line.split()
            self.organisms.add(elements[0])
            self.annotations[elements[0]] = self.__load_gff(elements[1], families, elements[0], lim_occurence)
            if len(elements)>2:
                self.circular_contig += elements[2:len(elements)]

        self.circular_contig = set(self.circular_contig)

        return nb_families

    def __load_gff(self, gff_file, families, organism, lim_occurence):
        """ Load the content of a gff file """ 
        logging.getLogger().info("Reading "+gff_file+" file ...")
        db_gff = gffutils.create_db(gff_file, ':memory:')
        annot = defaultdict(list)

        tandem_repeat = {}
        prev = None
        ctp_prev = 1

        cpt_fam_occ = defaultdict(int)

        for row in db_gff.all_features(featuretype='CDS', order_by=('seqid','start')):
            logging.getLogger().debug(row)
            protein = row.id
            familly = families[protein]
            if familly == prev:
                familly = familly+"-"+str(ctp_prev)
                cpt+=1
            else:
                cpt=1

            cpt_fam_occ[familly]+=1
            prev = families[protein]

            try:
                name = row.attributes['Name'].pop()
            except:
                try:
                    name = row.attributes['gene'].pop()
                except:
                    name = ""

            annot_row               = [protein,"CDS",organism,familly,row.start,row.end,row.strand, name]
            info                    = [protein,"CDS",organism,familly,row.start,row.end,row.strand, name]
            #self.gene_location[protein] = [organism, row.seqid, len(annot[row.seqid])]
            annot[row.seqid].append(annot_row)
        if (lim_occurence >0):
            fam_to_remove =[fam for fam,occ in cpt_fam_occ.items() if occ>lim_occurence]
            self.families_repeted = self.families_repeted.union(set(fam_to_remove))

        return(annot)

    def __str__(self):
        """ Return an overview of the object as a string """ 
        pan_str ="\n"
        pan_str += "----------- Statistics -----------\n"
        pan_str += "Number of organisms: "+str(self.nb_organisms)+"\n"

        if self.pan_size != 0:
            pan_str += "Pan-genome size:"+str(self.pan_size)+"\n"
            pan_str += "Exact core-genome size:"+str(len(self.partitions["CE"]))+"\n"
            pan_str += "Exact variable-genome size:"+str(self.pan_size-len(self.partitions["CE"]))+"\n"
            pan_str += "Persistent genome size:"+str(len(self.partitions["P"]))+"\n"
            pan_str += "Shell genome size:"+str(len(self.partitions["S"]))+"\n"
            pan_str += "Cloud genome cloud:"+str(len(self.partitions["C"]))+"\n"
        else:
            pan_str += "No partitioning have been performed on this Pangenome instance\n"
            pan_str += "Run the partitioning function to obtain more detailled statistics...\n"
        pan_str += "----------------------------------\n"

        return(pan_str)    

    def partition(self, nem_dir_path):
        """ Read the Mo Dang's thesis to understand NEM and Bernouilli Mixture Model, a summary is available here : http://www.kybernetika.cz/content/1998/4/393/paper.pdf """ 

        if not os.path.exists(nem_dir_path):
            #NEM requires 5 files: nem_file.index, nem_file.str, nem_file.dat, nem_file.m and nem_file.nei
            os.makedirs(nem_dir_path)
        self.intermediate_file = nem_dir_path

        logging.getLogger().info("Writing nem_file.str nem_file.index nem_file.nei nem_file.dat and nem_file.m files")
        str_file = open(nem_dir_path+"/nem_file.str", "w")
        index_file = open(nem_dir_path+"/nem_file.index", "w")
        org_file = open(nem_dir_path+"/column_org_file", "w")
        nei_file = open(nem_dir_path+"/nem_file.nei", "w")
        dat_file = open(nem_dir_path+"/nem_file.dat", "w")
        m_file = open(nem_dir_path+"/nem_file.m", "w")

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

        m_file.write("1 0.33 0.33 ") # 1 to initialize parameter, 0.33 and 0.33 for to give one third of initial portition to each class (last 0.33 is automaticaly determined by substraction)
        m_file.write(" ".join(["1"]*self.nb_organisms)+" ") # persistent binary vector
        m_file.write(" ".join(["1"]*self.nb_organisms)+" ") # shell binary vector
        m_file.write(" ".join(["0"]*self.nb_organisms)+" ") # cloud binary vector
        m_file.write(" ".join(["0.01"]*self.nb_organisms)+" ") # persistent dispersition vector
        m_file.write(" ".join(["0.5"]*self.nb_organisms)+" ") # shell dispersition vector
        m_file.write(" ".join(["0.01"]*self.nb_organisms)) # cloud dispersition vector

        index_file.close()
        nei_file.close()
        dat_file.close()
        m_file.close()

        logging.getLogger().info("Running NEM...")

        K       = 3 # number of partitions
        ALGO    = "nem" #fuzzy classification by mean field approximation
        BETA    = 0 # coeficient of spatial smoothing to apply, 0 is equivalent to EM for mixture model
        ITERMAX = 100 # number of iteration max 
        MODEL   = "bern" # multivariate Bernoulli mixture model
        PROPORTION = "pk" #equal proportion :  "p_"     varying proportion : "pk"
        VARIANCE_MODEL = "skd" #one variance per partition and organism : "sdk"      one variance per partition, same in all organisms : "sd_"   one variance per organism, same in all partion : "s_d"    same variance in organisms and partions : "s__" 

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
            self.BIC = 2 * M - (K * self.nb_organisms * 2 + K - 1) * math.log(self.pan_size)

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

            cloud_k = set([0,1,2]) - set([persistent_k,shell_k])
            cloud_k = list(cloud_k)[0]

            logging.getLogger().debug(sum_mu_k)
            logging.getLogger().debug(sum_epsilon_k)

            logging.getLogger().debug(persistent_k)
            logging.getLogger().debug(shell_k)
            logging.getLogger().debug(cloud_k)

            partition                 = {}
            partition[persistent_k+1] = "P"
            partition[shell_k+1]      = "S"
            partition[cloud_k+1]      = "C"

            logging.getLogger().debug(partition)
            #logging.getLogger().debug(index.keys())
            for node, nem_class in zip(self.neighbors_graph.nodes(), classification):
                nb_orgs=0
                for key, items in self.neighbors_graph.node[node].items():

                    self.neighbors_graph.node[node][key]=" ".join(items)    
                    if key!="name":
                        nb_orgs+=1

                self.neighbors_graph.node[node]["partition"]=partition[int(nem_class)]
                self.partitions[partition[int(nem_class)]].append(node)

                if nb_orgs >= self.nb_organisms:
                    self.partitions["CE"].append(node)
                    self.neighbors_graph.node[node]["partition_exact"]="C"
                else:
                    self.neighbors_graph.node[node]["partition_exact"]="A"

        if len(self.families_repeted)>0:
            logging.getLogger().info("Discarded families are:\n"+"\n".join(self.families_repeted))
        else:
            logging.getLogger().info("No families have been Discarded")

        logging.getLogger().debug(nx.number_of_edges(self.neighbors_graph))
        #positions = forceatlas2.forceatlas2_networkx_layout(self.neighbors_graph, 
        #                                                    niter=10,
        #                                                    edgeWeightInfluence=0.8)
        # figure = plt.figure()
        # nx.draw_graphviz(self.neighbors_graph, ax=figure.add_subplot(111))#, positions
        # figure.savefig("graph2.png")
        
    def export_to_GEXF(self, graph_output_path, compressed=False):
        """ """
        if self.neighbors_graph is None:
            logging.getLogger().error("neighbors_graph is not built, please use the function neighborhoodComputation before")
        elif self.intermediate_file is None:
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
                new_value = set(value.split(' '))
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

    def __add_families(self, fam_id, fam_nei, org, gene, name, edge = True):
        """ """
        self.neighbors_graph.add_node(fam_id)
        self.neighbors_graph.add_node(fam_nei)
        try:
            self.neighbors_graph.node[fam_id][org].add(gene)
        except KeyError:
            self.neighbors_graph.node[fam_id][org]=set([gene])
        try:
            self.neighbors_graph.node[fam_id]['name'].add(name)
        except KeyError:
            self.neighbors_graph.node[fam_id]['name']=set([name])

        if edge == True:
            if not self.neighbors_graph.has_edge(fam_id,fam_nei):
                self.neighbors_graph.add_edge(fam_id, fam_nei)
                logging.getLogger().debug([str(i) for i in [fam_id, fam_nei, org, gene, name, edge]])
            try:
                self.neighbors_graph[fam_id][fam_nei][org]+=1
            except KeyError:
                self.neighbors_graph[fam_id][fam_nei][org]=1
                try:
                    self.neighbors_graph[fam_id][fam_nei]["weight"]+=1.0
                except KeyError:
                    self.neighbors_graph[fam_id][fam_nei]["weight"]=1.0

    def neighborhood_computation(self):
        """ """
        if self.neighbors_graph is None:
            self.neighbors_graph = nx.Graph()
        elif self.intermediate_file is not None:
            logging.getLogger().error("The pangenome graph is already built and partionned, please use the function delete pangenome graph before build it again")

        for organism, annot_contigs in self.annotations.items():
            for contig, contig_annot in annot_contigs.items():
                at_least_2_families = False
                start = 0
                while (start < len(contig_annot) and contig_annot[start][FAMILLY_CODE] in self.families_repeted):
                    start += 1
                if start == len(contig_annot):
                    continue

                familly_id_nei = contig_annot[start][FAMILLY_CODE]
                logging.getLogger().debug(contig_annot[start])
                for index, gene_row in enumerate(contig_annot[start+1:]):
                    logging.getLogger().debug(gene_row)
                    if gene_row[FAMILLY_CODE] not in self.families_repeted:
                        self.__add_families(gene_row[FAMILLY_CODE],
                                          familly_id_nei,
                                          organism,
                                          gene_row[GENE],
                                          gene_row[NAME],
                                          True)
                        familly_id_nei = gene_row[FAMILLY_CODE]
                        at_least_2_families = True
                
                if contig in self.circular_contig and at_least_2_families:
                    self.__add_families(contig_annot[start][FAMILLY_CODE],
                                      familly_id_nei,
                                      organism,
                                      contig_annot[start][GENE],
                                      contig_annot[start][NAME],
                                      True)
                else:
                    self.__add_families(contig_annot[start][FAMILLY_CODE],
                                      contig_annot[start][FAMILLY_CODE],
                                      organism,
                                      contig_annot[start][GENE],
                                      contig_annot[start][NAME],
                                      False)
        self.pan_size = nx.number_of_nodes(self.neighbors_graph) 

    def delete_pangenome_graph(self):

        self.delete_intermediate_file()
        self.nem_output           = None
        self.pan_size             = 0
        self.nem_output           = None
        self.partitions           = {}
        self.partitions["P"]      = list()
        self.partitions["S"]      = list()
        self.partitions["C"]      = list()
        self.partitions["CE"]     = list()
        self.BIC                  = None


    def delete_intermediate_file(self):
        logging.getLogger().warning("delete "+self.intermediate_file)
        shutil.rmtree(self.intermediate_file)

if __name__=='__main__':
#blabla-1
    parser = argparse.ArgumentParser(description='Build a partitioned pangenome graph from annotated genomes and gene families')
    parser.add_argument('-?', '--version', action='version', version='0.1')
    parser.add_argument('-o', '--organisms', type=argparse.FileType('r'), nargs=1, help="""
The tsv file provided by progenome containing the gene annotations. --organims is a tab delimited files containg at least 2 mandatory fields by row and as many optional field as circular contig. Each row correspond to an organism to be added to the pangenome.
Reserved words are : "id", "label", "name", "weight", "partition", "partition_exact"
The first field is the organinsm name (id should be unique not contain any spaces, " or ' and reserved words).
The seconde field is the gff file associated to the organism. This path can be abolute or relative. The gff file must contain an id feature for each CDS (id should be unique not contain any spaces, " or ' and reserved words).
The next fields contain the name of perfectly assemble circular contigs (contig name must should be unique and not contain any spaces, " or ' and reserved words).
example:""", required=True)
    parser.add_argument('-f', '--gene_families', type=argparse.FileType('r'), nargs=1, help="""
is a file containg contain the families of gene. Each row contain 2 fields.
Reserved words are : "id", "label", "name", "weight", "partition", "partition_exact"
first field is the familly name.
second field is the gene name.
famillies are intended to be grouped by chuncks of row.
the families name can be any string but must should be unique and not contain any spaces, " or ' and reserved words
As a conventation, it is recommanded to use the name of the most reprensative gene of the famillies as the familly name.
Gene name can be any string corresponding to the if feature in the gff files. they should be unique and not contain any spaces, " or ' and reserved words.
example:""",  required=True)
    parser.add_argument('-d', '--output_directory', type=str, nargs=1, default="output.dir", help="""
The output directory""")
    parser.add_argument('-r', '--remove_high_copy_number_families', type=int, nargs=1, default=[-1], help=""""
Remove families having a number of copy of one families above or equal to this threshold in at least one organism (0 or negative value keep all families whatever their occurence). 
When -u is set, only work on new organisms added""")
    parser.add_argument("-u", "--update", default = None, type=argparse.FileType('r'), nargs=1, help="""
Pangenome Graph to be updated (in gexf format)""")
    parser.add_argument("-i", "--delete_intermediate_files", default=False, action="store_false", help="""
Delete intermediate files used by NEM""")
    parser.add_argument("-c", "--compress_graph", default=False, action="store_true", help="""
Compress (using gzip) the file containing the partionned pangenome graph""")
    parser.add_argument("-v", "--verbose", default=True, action="store_true", help="""
Show information message, otherwise only errors and warnings will be displayed""")
    parser.add_argument("-vv", "--verbose_debug", default=False, action="store_true", help="""
Show all messages including debug ones""")

    options = parser.parse_args()

    level = logging.WARNING
    if options.verbose is not None:
        level = logging.INFO
    if options.verbose_debug is not None:
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

    GEXF_GRAPH_FILE  = OUTPUTDIR+"/"+"graph.gexf"
    EXACT_CORE_FILE  = OUTPUTDIR+"/"+"exact_core.txt"
    PERSISTENT_FILE  = OUTPUTDIR+"/"+"persistant.txt"
    SHELL_FILE       = OUTPUTDIR+"/"+"shell.txt"
    CLOUD_FILE       = OUTPUTDIR+"/"+"cloud.txt"
    STATS_EXACT_FILE = OUTPUTDIR+"/"+"stats_exact.txt"
    STATS_NEM_FILE   = OUTPUTDIR+"/"+"stats_nem.txt" 

    pan = PPanGGOLiN("file", options.organisms[0],  options.gene_families[0], options.remove_high_copy_number_families[0])

    if options.update is not None:
        pan.import_from_GEXF(options.update[0])

    pan.neighborhood_computation()
    pan.partition(NEMOUTPUTDIR)
    logging.getLogger().info(pan)
    if options.compress_graph:
        pan.export_to_GEXF(GEXF_GRAPH_FILE, compressed=True)
    else:
        pan.export_to_GEXF(GEXF_GRAPH_FILE, compressed=False)

    core_list_file = open(EXACT_CORE_FILE,"w")
    core_list_file.write("\n".join(pan.partitions["CE"]))
    core_list_file.write("\n")
    core_list_file.close()

    persistent_list_file = open(PERSISTENT_FILE,"w")
    persistent_list_file.write("\n".join(pan.partitions["P"]))
    persistent_list_file.write("\n")
    persistent_list_file.close()

    shell_list_file = open(SHELL_FILE,"w")
    shell_list_file.write("\n".join(pan.partitions["S"]))
    shell_list_file.write("\n")
    shell_list_file.close()

    cloud_list_file = open(CLOUD_FILE,"w")
    cloud_list_file.write("\n".join(pan.partitions["C"]))
    cloud_list_file.write("\n")
    cloud_list_file.close()

    stat_exact_file = open(STATS_EXACT_FILE,"w")
    stat_exact_file.write("\t".join([str(pan.nb_organisms),
                                     str(len(pan.partitions["CE"])),
                                     str(pan.pan_size-len(pan.partitions["CE"])),
                                     str(pan.pan_size)])+ "\n")

    stat_nem_file = open(STATS_NEM_FILE,"w")
    stat_nem_file.write("\t".join([str(pan.nb_organisms),
                                   str(len(pan.partitions["P"])),
                                   str(len(pan.partitions["S"])),
                                   str(len(pan.partitions["C"])),
                                   str(pan.pan_size),
                                   str(pan.BIC)])+ "\n")

    pan.delete_pangenome_graph()

    logging.getLogger().debug("Finished !")
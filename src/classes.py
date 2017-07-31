#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

from collections import defaultdict
from collections import OrderedDict
from collections import Counter
from ordered_set import OrderedSet
import networkx as nx
import subprocess
import os
import re
import math
import logging
import gffutils
#import matplotlib.pyplot as plt  
import random 
#import forceatlas2 

NEM_LOCATION  = "../NEM/"
(GENE, TYPE, ORGANISM, FAMILLY_CODE, START, END, STRAND, NAME) = range(0, 8)

class Pangenome:

    def __init__(self, init_from = "args", *args):
        self.annotations          = dict()
        self.neighbors_graph      = None
        self.organisms            = OrderedSet()
        self.nb_organisms         = 0
        self.circular_contig      = set()
        self.core_list            = list()
        self.pan_size             = 0
        self.core_size            = 0
        self.partitions_size      = {}
        self.partitions_size["P"] = 0
        self.partitions_size["S"] = 0
        self.partitions_size["C"] = 0
        self.k                    = None
        self.BIC                  = None # Bayesian Index Criterion
        self.gene_location        = OrderedDict()
        self.families_repeted     = set()

        if init_from == "file":
            self.__initialize_from_files(*args)
        elif init_from == "args":
            (self.nb_organisms,
             self.organism_positions,
             self.familly_positions,
             self.annotation_positions) = args 
        else:
            raise ValueError("init_from parameter is required")
 
        self.neighbors_graph = self.__neighborhoodComputation()
        self.pan_size        = nx.number_of_nodes(self.neighbors_graph) 
        
    def __initialize_from_files(self, organisms_file, families_tsv_file, lim_occurence):

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
            self.nb_organisms +=1
            self.organisms.add(elements[0])
            self.annotations[elements[0]] = self.__load_gff(elements[1], families, elements[0], lim_occurence)
            if len(elements)>2:
                self.circular_contig += elements[2:len(elements)]

        self.circular_contig = set(self.circular_contig)

        return nb_families

    def __load_gff(self, gff_file, families, organism, lim_occurence):
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
            self.gene_location[protein] = [organism, row.seqid, len(annot[row.seqid])]
            annot[row.seqid].append(annot_row)
        if (lim_occurence >0):
            fam_to_remove =[fam for fam,occ in cpt_fam_occ.items() if occ>lim_occurence]
            self.families_repeted = self.families_repeted.union(set(fam_to_remove))

        return(annot)

    def __str__(self):
        pan_str =""
        pan_str += "----------- Statistics -----------\n"
        pan_str += "Number of organisms: "+str(self.nb_organisms)+"\n"

        if self.pan_size != 0:
            pan_str = ""
            pan_str += "Pan-genome size:"+str(self.pan_size)+"\n"
            pan_str += "Exact core-genome size:"+str(self.core_size)+"\n"
            pan_str += "Exact variable-genome size:"+str(self.pan_size-self.core_size)+"\n"
            pan_str += "Persistent genome size:"+str(self.partitions_size["P"])+"\n"
            pan_str += "Shell genome size:"+str(self.partitions_size["S"])+"\n"
            pan_str += "Cloud genome cloud:"+str(self.partitions_size["C"])+"\n"
        else:
            pan_str += "No partitioning have been performed on this Pangenome instance\n"
            pan_str += "Run the partitioning function to obtain more detailled statistics...\n"
        pan_str += "----------------------------------\n"

        return(pan_str)    

    def partition(self, result_path, k = 3, use_neighborhood = True, write_graph = None, neighbor_jumps = 1):
        """ """ 
        if not k>1:
            raise ValueError("k must be at leat equal to 2")

        if not neighbor_jumps > 0 and not neighbor_jumps == float("Inf"):
            raise ValueError("neighbor_jumps must be at leat equal to 1 or equal to floar('Inf')")

        if write_graph is not None:
            accessed_graph_output = ["gexf","gml","graphml","gpickle"]
            if write_graph not in accessed_graph_output:
                raise ValueError("write_graph must contain a format in the following list:'"+"', ".join(accessed_graph_output)+"'")

        if not os.path.exists(result_path):
            #NEM requires 5 files: nem_file.index, nem_file.str, nem_file.dat, nem_file.m (optional) and nem_file.nei
            os.makedirs(result_path)

        else:
            raise ValueError("result_path already exist")

        logging.getLogger().info("Writing nem_file.str nem_file.index nem_file.nei nem_file.dat and nem_file.m nem_files")
        index_file = open(result_path+"/nem_file.index", "w")
        nei_file = open(result_path+"/nem_file.nei", "w")
        dat_file = open(result_path+"/nem_file.dat", "w")
        m_file = open(result_path+"/nem_file.m", "w")

        if use_neighborhood:
            nei_file.write("1\n")
        else:
            nei_file.write("0\n")
        
        index = {node: index+1 for index, node in enumerate(self.neighbors_graph.nodes(data=False))}

        str_file = open(result_path+"/nem_file.str", "w")
        str_file.write("S\t"+str(self.pan_size)+"\t"+str(self.nb_organisms)+"\n")
        str_file.close()

        for node in self.neighbors_graph.nodes(data=True):
            node_name, node_organisms = node

            index_file.write(str(index[node_name])+"\t"+str(node_name)+"\n")
            logging.getLogger().debug(node_organisms)
            logging.getLogger().debug(self.organisms)
            dat_file.write("\t".join(["1" if org in node_organisms else "0" for org in self.organisms])+"\n")

            if use_neighborhood:
                row_fam         = []
                row_dist_score  = []
                neighbor_number = 0
                try:
                    for neighbor in nx.all_neighbors(self.neighbors_graph, node_name):
                        nb_presences = sum([pre_abs for org, pre_abs in self.neighbors_graph[node_name][neighbor].items() if org != 'weight'])
                        self.neighbors_graph[node_name][neighbor]["weight"]= nb_presences
                        distance_score = nb_presences/self.nb_organisms
                        row_fam.append(str(index[neighbor]))
                        row_dist_score.append(str(round(distance_score,4)))
                        neighbor_number += 1
                    if neighbor_number>0:
                        nei_file.write("\t".join([str(item) for sublist in [[index[node_name]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
                        #logging.getLogger().debug("\t".join([str(item) for sublist in [[[index[node_name]]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
                    else:
                        raise nx.exception.NetworkXError("no all_neighbors in selected organismss")
                except nx.exception.NetworkXError as nxe:
                    logging.getLogger().warning("The familly: "+node_name+" is an isolated familly")
                    nei_file.write(str(index[node_name])+"\t0\n")

                

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

        #bernouli -> no weight or normal -> weight
        model = "bern" #if self.weights is None else "norm"
        print_log = " -l y" if logging.getLogger().getEffectiveLevel() < 20 else "" 
        #command = NEM_LOCATION+"nem_exe "+result_path+"/nem_file "+str(k)+" -a nem -i 2000 -m "+model+" pk skd -s r 10 -n f -B fix -b "+("1" if use_neighborhood else "0")+" -T -O random"+print_log
        command = NEM_LOCATION+"nem_exe "+result_path+"/nem_file "+str(k)+" -a nem -i 2000 -m "+model+" pk skd -s m "+result_path+"/nem_file.m -B fix -b 0 -T -O random -f fuzzy"+print_log
     
        logging.getLogger().info(command)
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (out,err) = proc.communicate()
        logging.getLogger().debug(out)
        logging.getLogger().debug(err)
        # M = float(re.search("\(M = (.+?)\)",output).group(1))# M=
        
        # logging.getLogger().info("Based on "+str(k)+" classes, BIC = "+str(round(self.BIC,4)))
        if os.path.isfile(result_path+"/nem_file.uf"):
            logging.getLogger().info("Reading NEM results")
        else:
            logging.getLogger().error("No NEM output found in nem_file.uf")
        
        with open(result_path+"/nem_file.uf","r") as classification_nem_file, open(result_path+"/nem_file.mf","r") as parameter_nem_file:
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
            if model == "norm":
                self.BIC = 2 * M - (k * (self.nb_organisms + 1) + (1 if use_neighborhood else 0) + k- 1) * math.log(nx.number_of_nodes(self.neighbors_graph)) 
            elif model == "bern":
                self.BIC = 2 * M - (k * self.nb_organisms * 2 + (1 if use_neighborhood else 0) + k - 1) * math.log(nx.number_of_nodes(self.neighbors_graph))

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

            partition               = {}
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

                self.partitions_size[partition[int(nem_class)]] += 1

                if nb_orgs == self.nb_organisms:
                    self.core_size+=1
                    self.neighbors_graph.node[node]["partition_exact"]="C"
                else:
                    self.neighbors_graph.node[node]["partition_exact"]="A"

            if write_graph is not None:
                logging.getLogger().info("Writing graphML file")
                getattr(nx,'write_'+write_graph)(self.neighbors_graph,result_path+"/graph."+write_graph)
        
        logging.getLogger().info("Discarded families are:\n"+"\n".join(self.families_repeted))
        #positions = forceatlas2.forceatlas2_networkx_layout(self.neighbors_graph, 
        #                                                    niter=10,
        #                                                    edgeWeightInfluence=0.8)
        # figure = plt.figure()
        # nx.draw_graphviz(self.neighbors_graph, ax=figure.add_subplot(111))#, positions
        # figure.savefig("graph2.png")
        
    def __neighborhoodComputation(self):#initNeighborDistance, maxNeighborDistance,

        neighbors_graph = nx.Graph()
        
        def add_families(fam_id, fam_nei,org, gene, name,edge = True):#, prec , gene_nei
            neighbors_graph.add_node(fam_id)
            neighbors_graph.add_node(fam_nei)


            try:
                neighbors_graph.node[fam_id][org].add(gene)
            except KeyError:
                neighbors_graph.node[fam_id][org]=set([gene])
            try:
                neighbors_graph.node[fam_id]['name'].add(name)
            except KeyError:
                neighbors_graph.node[fam_id]['name']=set([name])

            if edge == True and not neighbors_graph.has_edge(fam_id,fam_nei):
                neighbors_graph.add_edge(fam_id, fam_nei)
                try:
                    neighbors_graph[fam_id][fam_nei]["weight"]+=1
                except KeyError:
                    neighbors_graph[fam_id][fam_nei]["weight"]=1
                try:
                    neighbors_graph[fam_id][fam_nei][org]+=1
                except KeyError:
                    neighbors_graph[fam_id][fam_nei][org]=1

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
                        add_families(gene_row[FAMILLY_CODE],
                                     familly_id_nei,
                                     organism,
                                     gene_row[GENE],
                                     gene_row[NAME],
                                     True)
                        familly_id_nei = gene_row[FAMILLY_CODE]
                        at_least_2_families = True
                
                if contig in self.circular_contig and at_least_2_families:
                    add_families(contig_annot[start][FAMILLY_CODE],
                                 familly_id_nei,
                                 organism,
                                 contig_annot[start][GENE],
                                 contig_annot[start][NAME],
                                 True)
                else:
                    add_families(contig_annot[start][FAMILLY_CODE],
                                 contig_annot[start][FAMILLY_CODE],
                                 organism,
                                 contig_annot[start][GENE],
                                 contig_annot[start][NAME],
                                 False)
        
        return neighbors_graph
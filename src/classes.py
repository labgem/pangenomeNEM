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
(GENE, TYPE, ORGANISM, FAMILLY_CODE, START, END, STRAND) = range(0, 7)

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
        # th = self.nb_organisms
        # highly_connected_nodes = self.jump_highly_connected_nodes(th)
        # if len(highly_connected_nodes)>0:
        #     logging.getLogger().warning("This nodes was removed because there are hyperconnected (at least "+str(th)+" neighbors)")
        #     logging.getLogger().warning(highly_connected_nodes)
        self.pan_size        = nx.number_of_nodes(self.neighbors_graph) 
        
        # for fam in self.familly_positions.keys():
        #     conservation = (np.sum(np.count_nonzero(Pangenome.presences_absences[self.familly_positions[fam],list(self.organism_positions.values())]))/float(self.nb_organisms))
        #     if conservation == 1.00:
        #         self.core_size += 1
        #         self.core_list.append(fam)
        #     self.familly_ratio[fam] = conservation

    def __initialize_from_files(self, organisms_file, families_tsv_file):

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
                families[elements[1]]=str(nb_families)
            else :
                prec = elements[0]
                nb_families +=1
                families[elements[1]]=str(nb_families)
        organisms = []
        self.circular_contig = []

        logging.getLogger().info("Reading "+organisms_file.name+" families file ...")

        for line in organisms_file:
            elements = line.split()
            self.nb_organisms +=1
            self.organisms.add(elements[0])
            self.annotations[elements[0]] = self.__load_gff(elements[1], families, elements[0])
            if len(elements)>2:
                self.circular_contig += elements[2:len(elements)]

        self.circular_contig = set(self.circular_contig)

        return nb_families
        #Pangenome.presences_absences = coo_matrix((nb_families+1,self.nb_organisms), dtype=np.dtype('uint8')).todense()
        # cpt_org = 0
        # cpt_fam = 0
        # cur_org = self.annotations[0][ORGANISM]
        # self.organism_positions[cur_org]=cpt_org
        # self.annotation_positions[cur_org]="0-"
        # for i in range(0, len(self.annotations)):
        #     if self.annotations[i][ORGANISM] != cur_org:
        #         self.annotation_positions[cur_org] = intspan(self.annotation_positions[cur_org]+str(i-1))
        #         cur_org = self.annotations[i][ORGANISM]
        #         self.annotation_positions[cur_org] = str(i)+"-"
        #         cpt_org += 1
        #         self.organism_positions[cur_org]=cpt_org
        #     fam_id = self.annotations[i][FAMILLY_CODE]
        #     #self.familly_positions[fam_id] = int(fam_id)
        #     #Pangenome.presences_absences[int(fam_id),cpt_org] += 1   
        # self.annotation_positions[cur_org] = intspan(self.annotation_positions[cur_org]+str(i))
        # self.familly_positions             = OrderedDict(sorted(self.familly_positions.items(), key=lambda t: t[0]))

    def __load_gff(self, gff_file, families, organism):
        logging.getLogger().info("Reading "+gff_file+" file ...")
        db_gff = gffutils.create_db(gff_file, ':memory:')
        annot = defaultdict(list)

        tandem_repeat = {}
        prev = None
        ctp_prev = 1
        for row in db_gff.all_features(featuretype='CDS', order_by=('seqid','start')):
            logging.getLogger().debug(row)
            protein = row.id
            familly = families[protein]
            if familly == prev:
                familly = familly+"-"+str(ctp_prev)
                cpt+=1
            else:
                cpt=1
            prev = families[protein]

            annot_row               = [protein,"CDS",organism,familly,row.start,row.end,row.strand]
            info                    = [protein,"CDS",organism,familly,row.start,row.end,row.strand]
            self.gene_location[protein] = [organism, row.seqid, len(annot[row.seqid])]
            annot[row.seqid].append(annot_row)
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
            #NEM requires 5 files: nem_file.index, nem_file.str, nem_file.dat, nem_file.ck (optional) and nem_file.nei
            os.makedirs(result_path)

        else:
            raise ValueError("result_path already exist")

        logging.getLogger().info("Writing nem_file.str nem_file.index nem_file.nei nem_file.dat and nem_file.ck nem_files")
        index_file = open(result_path+"/nem_file.index", "w")
        #ck_file  = open(result_path+"/nem_file.ck", "w")
        nei_file = open(result_path+"/nem_file.nei", "w")
        dat_file = open(result_path+"/nem_file.dat", "w")

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

        index_file.close()
        #ck_file.close()
        nei_file.close()
        dat_file.close()             

        logging.getLogger().info("Running NEM...")

        #bernouli -> no weight or normal -> weight
        model = "bern" #if self.weights is None else "norm"
        print_log = " -l y" if logging.getLogger().getEffectiveLevel() < 20 else "" 
        command = NEM_LOCATION+"nem_exe "+result_path+"/nem_file "+str(k)+" -a nem -i 2000 -m "+model+" pk skd -s r 10 -B fix -b "+("1" if use_neighborhood else "0")+" -T -O random"+print_log
        logging.getLogger().info(command)
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.communicate()
        # M = float(re.search("\(M = (.+?)\)",output).group(1))# M=
        
        # logging.getLogger().info("Based on "+str(k)+" classes, BIC = "+str(round(self.BIC,4)))
        if os.path.isfile(result_path+"/nem_file.cf"):
            logging.getLogger().info("Reading NEM results")
        else:
            logging.getLogger().error("No NEM output found in nem_file.cf")
        
        with open(result_path+"/nem_file.cf","r") as classification_nem_file, open(result_path+"/nem_file.mf","r") as parameter_nem_file:
            classification = classification_nem_file.readline().split()
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
                for org, genes in self.neighbors_graph.node[node].items():
                    self.neighbors_graph.node[node][org]=" ".join(genes)
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
        
        #positions = forceatlas2.forceatlas2_networkx_layout(self.neighbors_graph, 
        #                                                    niter=10,
        #                                                    edgeWeightInfluence=0.8)
        # figure = plt.figure()
        # nx.draw_graphviz(self.neighbors_graph, ax=figure.add_subplot(111))#, positions
        # figure.savefig("graph2.png")
        
    def __neighborhoodComputation(self):#initNeighborDistance, maxNeighborDistance,

        neighbors_graph = nx.Graph()
        all_paths = defaultdict(lambda : defaultdict(list))  

        def add_neighbors(fam_id, fam_nei, prec,org, gene, gene_nei):
            neighbors_graph.add_node(fam_id)
            neighbors_graph.add_node(fam_nei)

            if prec is None:
                logging.getLogger().debug([fam_id, fam_nei, prec,org,gene])
            if prec is not None:
                all_paths[fam_nei][frozenset([prec,fam_id])].append(tuple([org,gene_nei]))
            else:
                all_paths[fam_nei][None].append(fam_id)
            try:
                neighbors_graph.node[fam_id][org]+=(" "+gene)
            except KeyError:
                neighbors_graph.node[fam_id][org]=gene
            if not neighbors_graph.has_edge(fam_id,fam_nei):
                neighbors_graph.add_edge(fam_id, fam_nei)
            # try:
            # neighbors_graph[fam_id][fam_nei]["weight"]=1
            # except KeyError:
            #     neighbors_graph[fam_id][fam_nei]["weight"]=1
            #try:
            #    neighbors_graph[fam_id][fam_nei][org]+=1
            #except KeyError:
            neighbors_graph[fam_id][fam_nei][org]=1

        for organism, annot_contigs in self.annotations.items():
            for contig, contig_annot in annot_contigs.items():
                if len(contig_annot)>1:
                    prec  = None
                    familly_id_nei = contig_annot[0][FAMILLY_CODE]
                    gene_nei       = contig_annot[0][GENE]
                    logging.getLogger().debug(contig_annot[0])
                    for index, row in enumerate(contig_annot[1:]):
                        logging.getLogger().debug(row)
                        gene       = row[GENE]
                        familly_id = row[FAMILLY_CODE]     
                        add_neighbors(familly_id, familly_id_nei, prec, organism, gene, gene_nei)
                        prec           = familly_id_nei
                        familly_id_nei = familly_id
                        gene_nei       = gene
                    
                    row = contig_annot[0]
                    familly_id = row[FAMILLY_CODE]
                    gene = row[GENE]
                    
                    if contig in self.circular_contig:
                        add_neighbors(familly_id, familly_id_nei, prec, organism, gene, gene_nei)
                        logging.getLogger().debug("first2 "+familly_id)
                        logging.getLogger().debug("prec "+str(all_paths[familly_id]))
                        prec = all_paths[familly_id].pop(None)[0]
                        all_paths[familly_id][frozenset([familly_id_nei,prec])].append(tuple([organism,gene]))
                    else:
                        logging.getLogger().debug(contig)
                        all_paths[familly_id].pop(None)[0]
                else:
                    gene       = contig_annot[0][GENE]
                    familly_id = contig_annot[0][FAMILLY_CODE]
                    neighbors_graph.add_node(familly_id)
                    try:
                        neighbors_graph.node[familly_id][organism]+=(" "+gene)
                    except KeyError:
                        neighbors_graph.node[familly_id][organism]=gene
        
        #inspired from http://stackoverflow.com/a/9114443/7500030
        def merge_overlapping_path(data):
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
 
        #untangling stage:
        for node in list(neighbors_graph.node):
            logging.getLogger().debug(node)
            logging.getLogger().debug(all_paths[node])
            
            path_groups = merge_overlapping_path(all_paths[node])
            logging.getLogger().debug(path_groups)
            logging.getLogger().debug(len(path_groups))

            if len(path_groups)>1:#if several path are non overlaping to browse this node
                logging.getLogger().debug("split "+str(len(path_groups)))
                for suffix, path_group in enumerate(path_groups):
                    new_node = node+"_"+str(suffix)
                    #all_paths[new_node][frozenset(path_group)]= -1
                    logging.getLogger().debug("new_node:"+new_node)
                    neibors_new_node = [neighbor for neighbor in nx.all_neighbors(neighbors_graph, node) if neighbor in path_group]
                    for new_neibor in neibors_new_node:
                        renamed_paths = dict()
                        while len(all_paths[new_neibor])>0:
                            to_rename = all_paths[new_neibor].popitem()
                            set_to_rename = set(to_rename[0])
                            logging.getLogger().debug("to_rename:"+str(set_to_rename))
                            if node in set_to_rename:
                                set_to_rename.remove(node)
                                set_to_rename.add(new_node)
                            renamed_paths[frozenset(set_to_rename)]=to_rename[1]
                            logging.getLogger().debug("renamed:"+str(set_to_rename))
                        all_paths[new_neibor].update(renamed_paths)
                    neighbors_graph.add_node(new_node)
                    
                    already_added_gene = set()
                    for neighbor in path_group:
                        neighbors_graph.add_edge(new_node,neighbor)
                        logging.getLogger().debug(neighbor)
                        for path, genes in all_paths[node].items():
                            logging.getLogger().debug(path)
                            if neighbor in path:
                                for org, gene in genes:
                                    if gene not in already_added_gene:
                                        all_paths[new_node][path].append(tuple([org,gene]))
                                        org, contig, pos = self.gene_location[gene]
                                        self.annotations[org][contig][pos][FAMILLY_CODE]=new_node
                                        try:
                                            neighbors_graph.node[new_node][org].add(gene)
                                        except KeyError:
                                            neighbors_graph.node[new_node][org]=set([gene])
                                        already_added_gene.add(gene)
                                    #try:
                                    #    neighbors_graph[new_node][neighbor][org]+=1
                                    #except:
                                    neighbors_graph[new_node][neighbor][org]=1
                                    # try:
                                    # neighbors_graph[new_node][neighbor]["weight"]=1
                                    # except:
                                        # neighbors_graph[new_node][neighbor]["weight"]=1
                
                del all_paths[node]
                neighbors_graph.remove_node(node)

        # for node_i, node_j, orgs in neighbors_graph.edges(data=True):
        #     if (neighbors_graph.has_edge(node_i, node_j)):
        #     new_node = "("+node_i+"#"+node_j+")"
        #     for path, genes in all_paths[node_i].items():
        #         for org, gene in genes:
        #             if org in orgs:
        #                 org, contig, pos = self.gene_location[gene]
        #                 try:
        #                     direction = 1
        #                     if contig in self.circular_contig:
        #                         nbgene = len(self.annotations[org][contig])
        #                         before_i = self.annotations[org][contig][(pos-(1*direction))%nbgene][FAMILLY_CODE]
        #                         after_i  = self.annotations[org][contig][(pos+(1*direction))%nbgene][FAMILLY_CODE]

        #                         if node_j == before_i:
        #                             before_i = after_i
        #                             after_i = node_j
        #                             direction = -1

        #                         logging.getLogger().debug(node_i)
        #                         logging.getLogger().debug(node_j)

        #                         logging.getLogger().debug([pos-(1*direction),self.annotations[org][contig][(pos-(1*direction))%nbgene]])
        #                         logging.getLogger().debug([pos,self.annotations[org][contig][pos%nbgene]])
        #                         logging.getLogger().debug([pos+(1*direction),self.annotations[org][contig][(pos+(1*direction))%nbgene]])
        #                         logging.getLogger().debug([pos+(2*direction),self.annotations[org][contig][(pos+(2*direction))%nbgene]])

        #                         after_j = self.annotations[org][contig][(pos+(2*direction))%nbgene][FAMILLY_CODE]
        #                         if node_j == after_i:
        #                             logging.getLogger().debug("ici")
        #                             all_paths[new_node][frozenset([before_i,after_j])].append(tuple([org,gene]))
        #                         else:
        #                             logging.getLogger().debug("la")
        #                         #exit()
        #                     else:
        #                         try:
        #                             before = self.annotations[org][contig][(pos-(1*direction))%len(self.annotations[org][contig])][FAMILLY_CODE]    
        #                             after  = self.annotations[org][contig][(pos+(1*direction))%len(self.annotations[org][contig])][FAMILLY_CODE]
        #                         except:
        #                             pass
        #                 except IndexError:
        #                     logging.getLogger().debug("Que faire ?")
        #                     logging.getLogger().debug([pos,self.annotations[org][contig][pos]])
        #                     logging.getLogger().debug([pos+(1*direction),self.annotations[org][contig][pos-(1*direction)]])
                            
        #                     logging.getLogger().debug(self.annotations[org][contig][pos+(1*direction)])
        #                     logging.getLogger().debug(self.annotations[org][contig][pos+(2*direction)])
        #                     #racrocher les gènes à la famille la pllus problable
        #                     pass
        #     logging.getLogger().debug(all_paths[new_node])
        #     path_groups = merge_overlapping_path(all_paths[new_node])
        #     if len(path_groups)>1:#if several path are non overlaping to browse this (meta)node
        #         logging.getLogger().debug("split "+str(len(path_groups)))
        #         logging.getLogger().debug("split "+str(path_groups))
        #         for suffix, path_group in enumerate(path_groups):
        #             new_node_i = node_i+"#"+str(suffix)
        #             new_node_j = node_j+"#"+str(suffix)
        #             neighbors_graph.add_node(new_node_i)
        #             neighbors_graph.add_node(new_node_j)
        #             neighbors_graph.add_edge(new_node_i, new_node_j)
        #             neighbors_graph[new_node_i][new_node_j]["weight"]=5
        #             logging.getLogger().debug("new_node:"+new_node)
        #             logging.getLogger().debug("new_node_i:"+new_node_i)
        #             logging.getLogger().debug("new_node_j:"+new_node_j)
        #             neibors_new_node_i = [neighbor for neighbor in nx.all_neighbors(neighbors_graph, node_i) if neighbor in path_group]
        #             logging.getLogger().debug("neibors_new_node_i:"+str(neibors_new_node_i))
        #             logging.getLogger().debug([neighbor for neighbor in nx.all_neighbors(neighbors_graph, node_i)])
        #             for new_neibor in neibors_new_node_i:
        #                 neighbors_graph.add_edge(new_node_i,new_neibor)
        #                 neighbors_graph.remove_edge(node_i,new_neibor)
        #                 neighbors_graph[new_node_i][new_neibor]["weight"]=5
        #                 logging.getLogger().debug("add_edge i:"+str([new_node_i,new_neibor]))
        #                 renamed_paths = dict()
        #                 while len(all_paths[new_neibor])>0:
        #                     to_rename = all_paths[new_neibor].popitem()
        #                     set_to_rename = set(to_rename[0])
        #                     logging.getLogger().debug("to_rename:"+str(set_to_rename))
        #                     if node_i in set_to_rename:
        #                         set_to_rename.remove(node_i)
        #                         set_to_rename.add(new_node_i)
        #                     renamed_paths[frozenset(set_to_rename)]=to_rename[1]
        #                     logging.getLogger().debug("renamed:"+str(set_to_rename))
        #                 all_paths[new_neibor].update(renamed_paths)
        #             neibors_new_node_j = [neighbor for neighbor in nx.all_neighbors(neighbors_graph, node_j) if neighbor in path_group]
        #             logging.getLogger().debug("neibors_new_node_j:"+str(neibors_new_node_j))
        #             logging.getLogger().debug([neighbor for neighbor in nx.all_neighbors(neighbors_graph, node_j)])
        #             for new_neibor in neibors_new_node_j:
        #                 neighbors_graph.add_edge(new_node_j,new_neibor)
        #                 neighbors_graph.remove_edge(node_j,new_neibor)
        #                 neighbors_graph[new_node_j][new_neibor]["weight"]=5
        #                 logging.getLogger().debug("add_edge j:"+str([new_node_j,new_neibor]))
        #                 renamed_paths = dict()
        #                 while len(all_paths[new_neibor])>0:
        #                     to_rename = all_paths[new_neibor].popitem()
        #                     set_to_rename = set(to_rename[0])
        #                     logging.getLogger().debug("to_rename:"+str(set_to_rename))
        #                     if node_j in set_to_rename:
        #                         set_to_rename.remove(node_j)
        #                         setsQ_to_rename.add(new_node_j)
        #                     renamed_paths[frozenset(set_to_rename)]=to_rename[1]
        #                     logging.getLogger().debug("renamed:"+str(set_to_rename))
        #                 all_paths[new_neibor].update(renamed_paths)
        #             already_added_gene = set()
        #             for neighbor in path_group:
        #                 logging.getLogger().debug(neighbor)
        #                 for path, genes in all_paths[node_i].items():
        #                     logging.getLogger().debug(path)
        #                     if neighbor in path:
        #                         logging.getLogger().debug(path)
        #                         for org, gene in genes:
        #                             if gene not in already_added_gene:
        #                                 all_paths[new_node_i][path].append(tuple([org,gene]))
        #                                 org, contig, pos = self.gene_location[gene]
        #                                 self.annotations[org][contig][pos][FAMILLY_CODE]=new_node_i
        #                                 try:
        #                                     neighbors_graph.node[new_node_i][org].add(gene)
        #                                 except KeyError:
        #                                     neighbors_graph.node[new_node_i][org]=set([gene])
        #                                 already_added_gene.add(gene)
        #                                 try:
        #                                     neighbors_graph[new_node_i][neighbor][org]+=1
        #                                 except KeyError:
        #                                     logging.getLogger().debug(new_node_i)
        #                                     logging.getLogger().debug(neighbor)
        #                                     logging.getLogger().debug(neighbors_graph[new_node_i][neighbor])
        #                                     neighbors_graph[new_node_i][neighbor][org]=1
        #                                 try:
        #                                     neighbors_graph[new_node_i][neighbor]["weight"]+=1
        #                                 except:
        #                                     neighbors_graph[new_node_i][neighbor]["weight"]=1
        #                         #neighbors_graph.node[new_node][org]=" ".join(neighbors_graph.node[new_node][org])
        #                 logging.getLogger().debug(neighbor)
        #                 for path, genes in all_paths[node_j].items():
        #                     logging.getLogger().debug(path)
        #                     if neighbor in path:
        #                         logging.getLogger().debug(path)
        #                         for org, gene in genes:
        #                             if gene not in already_added_gene:
        #                                 logging.getLogger().debug(gene)
        #                                 all_paths[new_node_j][path].append(tuple([org,gene]))
        #                                 org, contig, pos = self.gene_location[gene]
        #                                 self.annotations[org][contig][pos][FAMILLY_CODE]=new_node_j
        #                                 try:
        #                                     neighbors_graph.node[new_node_j][org].add(gene)
        #                                 except KeyError:
        #                                     neighbors_graph.node[new_node_j][org]=set([gene])
        #                                 already_added_gene.add(gene)
        #                                 try:
        #                                     neighbors_graph[new_node_j][neighbor][org]+=1
        #                                 except KeyError:
        #                                     logging.getLogger().debug(new_node_j)
        #                                     logging.getLogger().debug(neighbor)
        #                                     logging.getLogger().debug(neighbors_graph[new_node_j][neighbor])
        #                                     neighbors_graph[new_node_j][neighbor][org]=1
        #                                 try:
        #                                     neighbors_graph[new_node_j][neighbor]["weight"]+=1
        #                                 except:
        #                                     neighbors_graph[new_node_j][neighbor]["weight"]=1

                #del all_paths[node_i]
                #del all_paths[node_j]
                #neighbors_graph.remove_node(node_i)
                #neighbors_graph.remove_node(node_j)
            else:
                logging.getLogger().debug("no split")
        # for node, path_group in all_paths.items():
        #     for neighbors in path_group:
                
        #         res = [(node,node_bis) for path_group_neighbors in all_paths[neighbors] for node_bis in path_group_neighbors if node_bis == node]
        #         if len(res)>0:
        #             logging.getLogger().debug(res)
        #             exit()

        # for organism, annot_contigs in self.annotations.items():
        #     for contig, contig_annot in annot_contigs.items():
        #         if contig not in self.circular_contig:
        #                 continue
        #         order = 0
        #         for index, row in enumerate(contig_annot):
        #             fam = row[FAMILLY_CODE]
        #             try:
        #                 neighbors_graph.node[fam]["order_"+contig.replace(".", "_")]=order
        #             except:
        #                 try:
        #                     logging.getLogger().debug(row[FAMILLY_CODE])
        #                     neighbors_graph.node[row[FAMILLY_CODE]]["order_"+contig.replace(".", "_")]=order
        #                 except:
        #                     logging.getLogger().warning(row)
        #             order+=1
        return neighbors_graph

    def jump_highly_connected_nodes(self, max_degree):

        logging.getLogger().info("Jump over nodes having a degree higher or equal to "+str(max_degree))

        highly_connected_nodes = set()
        if self.neighbors_graph is not None:
            new_highly_connected_nodes = set([node for node, degree in self.neighbors_graph.degree_iter() if degree>=max_degree])
            while len(new_highly_connected_nodes - highly_connected_nodes) > 0:
                highly_connected_nodes = highly_connected_nodes.union(new_highly_connected_nodes)
                for hcn in new_highly_connected_nodes:
                    cpt = 0
                    for org, genes in self.neighbors_graph.node[hcn].items():
                        logging.getLogger().debug([org, genes])
                        for gene in genes:
                            org, contig, pos = self.gene_location[gene]
                            curr_gene = self.annotations[org][contig][pos]
                            try:
                                i=1
                                prev_gene = self.annotations[org][contig][pos-i]
                                while prev_gene[FAMILLY_CODE] in highly_connected_nodes:
                                    curr_gene = prev_gene
                                    prev_gene = self.annotations[org][contig][pos-i]
                                    i+=1
                                prev = True
                                try:
                                    self.neighbors_graph.remove_edge(prev_gene[FAMILLY_CODE], curr_gene[FAMILLY_CODE])
                                except nx.exception.NetworkXError:
                                    pass
                            except KeyError:
                                prev = False
                            try:
                                i=1
                                next_gene = self.annotations[org][contig][pos+i]
                                while next_gene[FAMILLY_CODE] in highly_connected_nodes:
                                    curr_gene = next_gene
                                    next_gene = self.annotations[org][contig][pos+i]
                                    i+=1
                                nex = True
                                try:
                                    self.neighbors_graph.remove_edge(next_gene[FAMILLY_CODE], curr_gene[FAMILLY_CODE])
                                except nx.exception.NetworkXError:
                                    pass
                            except KeyError:
                                nex  = False
                            logging.getLogger().debug(prev_gene)
                            logging.getLogger().debug(curr_gene)
                            logging.getLogger().debug(next_gene)

                            if prev == True and nex == True:
                                if not self.neighbors_graph.has_edge(prev_gene[FAMILLY_CODE],next_gene[FAMILLY_CODE]):
                                    self.neighbors_graph.add_edge(prev_gene[FAMILLY_CODE], next_gene[FAMILLY_CODE])
                                # try:
                                # self.neighbors_graph[prev_gene[FAMILLY_CODE]][next_gene[FAMILLY_CODE]]["weight"]=1
                                # except KeyError:
                                #     self.neighbors_graph[prev_gene[FAMILLY_CODE]][next_gene[FAMILLY_CODE]]["weight"]=1
                                #try:
                                #    self.neighbors_graph[prev_gene[FAMILLY_CODE]][next_gene[FAMILLY_CODE]][org]+=1
                                #except KeyError:
                                self.neighbors_graph[prev_gene[FAMILLY_CODE]][next_gene[FAMILLY_CODE]][org]=1
                
                new_highly_connected_nodes = set([node for node, degree in self.neighbors_graph.degree_iter() if degree>max_degree])
        else:
            raise ValueError("neighbors_graph is None")
        return(highly_connected_nodes)

    # def __neighborhoodComputation(self, max_degree = float("Inf"), high_degree_node = set()):#initNeighborDistance, maxNeighborDistance,
    #     neighbors_graph = nx.Graph()
    #     all_paths = defaultdict(lambda : defaultdict(list))  

    #     def add_neighbors(fam_id, fam_nei, prec,org,gene):
    #         neighbors_graph.add_node(fam_id)
    #         neighbors_graph.add_node(fam_nei)

    #         if prec is None:
    #             logging.getLogger().debug([fam_id, fam_nei, prec,org,gene])

    #         if prec is not None:
    #             all_paths[fam_nei][frozenset([prec,fam_id])].append(tuple([org,gene]))
    #         else:
    #             all_paths[fam_nei][None].append(fam_id)
    #         try:
    #             neighbors_graph.node[fam_id][org]+=" "+gene
    #         except KeyError:
    #             neighbors_graph.node[fam_id][org]=gene
    #         if not neighbors_graph.has_edge(fam_id,fam_nei):
    #             neighbors_graph.add_edge(fam_id, fam_nei)
    #         try:
    #             neighbors_graph[fam_id][fam_nei]["weight"]+=1
    #         except KeyError:
    #             neighbors_graph[fam_id][fam_nei]["weight"]=1
    #         try:
    #             neighbors_graph[fam_id][fam_nei][org]+=1
    #         except KeyError:
    #             neighbors_graph[fam_id][fam_nei][org]=1

    #     for organism, annot_contigs in self.annotations.items():
    #         for contig, contig_annot in annot_contigs.items():
    #             prec  = None
    #             first = None
    #             shift = 0
    #             while True:
    #                 row = contig_annot[shift]
    #                 familly_id_nei = row[FAMILLY_CODE]
    #                 gene_nei = row[GENE]
    #                 if familly_id_nei not in high_degree_node:
    #                    first = shift
    #                    break
    #                 else:
    #                     shift+=1

    #             logging.getLogger().debug(row)
                
    #             for index, row in enumerate(contig_annot[first+1:]):
    #                 logging.getLogger().debug(row)        
    #                 gene       = row[GENE]
    #                 familly_id = row[FAMILLY_CODE]
    #                 skipped_node = set()
    #                 if familly_id not in high_degree_node:
    #                     add_neighbors(familly_id, familly_id_nei, prec, organism, gene_nei)
    #                     prec           = familly_id_nei
    #                     familly_id_nei = familly_id
    #                     gene_nei       = gene
    #                 else:
    #                     logging.getLogger().debug("skipped"+gene)
    #                     skipped_node.add(index)
                
    #             row = contig_annot[first]
    #             familly_id = row[FAMILLY_CODE]
    #             gene = row[GENE]
                
    #             if contig in self.circular_contig:
    #                 add_neighbors(familly_id, familly_id_nei, prec, organism, gene_nei)
    #                 logging.getLogger().debug("first2 "+familly_id)
    #                 logging.getLogger().debug(shift)
    #                 logging.getLogger().debug("prec "+str(all_paths[familly_id]))
    #                 prec = all_paths[familly_id].pop(None)[0]
    #                 all_paths[familly_id][frozenset([familly_id_nei,prec])].append(tuple([organism,gene]))
    #             else:
    #                 logging.getLogger().debug(contig)
    #                 all_paths[familly_id].pop(None)[0]

    #     new_high_degree_node = set([node for node, degree in neighbors_graph.degree_iter() if degree>max_degree])
    #     if len(new_high_degree_node)>max_degree:
    #         logging.getLogger().debug(len(new_high_degree_node))
    #         neighbors_graph = self.__neighborhoodComputation(max_degree, high_degree_node.union(new_high_degree_node))
    #         return neighbors_graph
    #     else:
    #         untangeled_neighbors_graph = neighbors_graph.copy()
    #         #inspired from http://stackoverflow.com/a/9114443/7500030
    #         def merge_overlapping_path(data):
    #             sets = (set(e) for e in data if e)
    #             results = [next(sets)]
    #             for e_set in sets:
    #                 to_update = []
    #                 for i,res in enumerate(results):
    #                     if not e_set.isdisjoint(res):
    #                         to_update.insert(0,i)
    #                 if not to_update:
    #                     results.append(e_set)
    #                 else:
    #                     last = results[to_update.pop(-1)]
    #                     for i in to_update:
    #                         last |= results[i]
    #                         del results[i]
    #                     last |= e_set
    #             return results

    #         inversed_dict = dict()
    #         logging.getLogger().debug(all_paths)
    #         #untangling stage:
    #         for node in list(untangeled_neighbors_graph.node):
    #             logging.getLogger().debug(node)
    #             logging.getLogger().debug(all_paths[node])
    #             path_groups = merge_overlapping_path(all_paths[node])

    #             logging.getLogger().debug(path_groups)
    #             logging.getLogger().debug(len(path_groups))
    #             if len(path_groups)>1:
    #                 logging.getLogger().debug("split "+str(len(path_groups)))
    #                 for suffix, path_group in enumerate(path_groups):
    #                     new_node = node+"_"+str(suffix)
    #                     all_paths[new_node][frozenset(path_group)]= -1
    #                     logging.getLogger().debug("new_node:"+new_node)
    #                     neibors_new_node = [neighbor for neighbor in nx.all_neighbors(untangeled_neighbors_graph, node) if neighbor in path_group]
    #                     for new_neibor in neibors_new_node:
    #                         renamed_paths = dict()
    #                         while len(all_paths[new_neibor])>0:
    #                             to_rename = all_paths[new_neibor].popitem()
    #                             set_to_rename = set(to_rename[0])
    #                             logging.getLogger().debug("to_rename:"+str(set_to_rename))
    #                             if node in set_to_rename:
    #                                 set_to_rename.remove(node)
    #                                 set_to_rename.add(new_node)
    #                             renamed_paths[frozenset(set_to_rename)]=to_rename[1]
    #                             logging.getLogger().debug("renamed:"+str(set_to_rename))
    #                         all_paths[new_neibor].update(renamed_paths)
    #                     untangeled_neighbors_graph.add_node(new_node)
    #                     for neighbor in path_group:
    #                         logging.getLogger().debug(all_paths[node])
    #                         untangeled_neighbors_graph.add_edge(new_node,neighbor)
    #                         weight=0
    #                         for path, genes in all_paths[node].items():
    #                             if neighbor in path:
    #                                 for org, gene in genes:
    #                                     try:
    #                                         untangeled_neighbors_graph.node[new_node][org]=+" "+gene
    #                                     except:
    #                                         untangeled_neighbors_graph.node[new_node][org]=gene
    #                                     try:
    #                                         untangeled_neighbors_graph[new_node][neighbor][org]+=1
    #                                     except:
    #                                         untangeled_neighbors_graph[new_node][neighbor][org]=1
    #                                     inversed_dict.update({gene:new_node})
    #                                     weight+=1
    #                         untangeled_neighbors_graph[new_node][neighbor]["weight"]=weight
    #                 untangeled_neighbors_graph.remove_node(node)

    #         for organism, annot_contigs in self.annotations.items():
    #             for contig, contig_annot in annot_contigs.items():
    #                 if contig not in self.circular_contig:
    #                         continue
    #                 order = 0
    #                 for index, row in enumerate(contig_annot):
    #                     fam = row[FAMILLY_CODE]
    #                     if fam in high_degree_node.union(new_high_degree_node):
    #                         continue
    #                     try:
    #                         untangeled_neighbors_graph.node[fam]["order_"+contig.replace(".", "_")]=order
    #                     except:
    #                         try:
    #                             logging.getLogger().debug(inversed_dict[row[GENE]])
    #                             untangeled_neighbors_graph.node[inversed_dict[row[GENE]]]["order_"+contig.replace(".", "_")]=order
    #                         except:
    #                             logging.getLogger().warning(row)
    #                     order+=1
    #     return untangeled_neighbors_graph


#inspired from http://stackoverflow.com/a/9114443/7500030
    def merge_overlapping_path(data):
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


# class Node:
#     def __init__(self, name):
#         self.name = name
#         self.neibors_before = list()
#         self.neibors_after = list()

#     def split(self):
#         sets = set()
#         for index in range(0,len(self.neibors_before)):
#             sets.add(set([self.neibors_before,self.neibors_after]))
#         return merge_overlapping_path(sets)
#     def merge(self, other_node):

# class Syntenie(Node):
#     def __init__(self, name):
#         super(Syntenie, self).__init__(name):
#         self.subnode = list()
#     def split(self):
#         super(Syntenie, self).split()
#     def __str__:
#         pass

# class Familly(Node):
#     def __init__(self, name):
#         super(Familly, self).__init__(name,neibors_before,neibors_after,genes):
#         self.neibors_before = list()
#         self.neibors_after = list()
#         self.genes = list()

#     def __str__:
#         pass

# class Gene:
#     def __init__(self, name, org, contig, pos):
#         self.name=name
#         self.org=org
#         self.contig=contig
#         self.pos = pos
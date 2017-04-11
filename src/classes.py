#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

#import mysql.connector as mysql
import pandas as pd
import numpy as np
from scipy.sparse import coo_matrix
from collections import defaultdict
from collections import OrderedDict
from collections import Counter
from intspan import intspan
import networkx as nx
import subprocess
import os
import re
import math
import logging
import gffutils

NEM_LOCATION  = "../NEM/"
(GENE, TYPE, ORGANISM, CONTIG_ID, FAMILLY_CODE, START, END, STRAND) = range(0, 8)

class Pangenome:

    presences_absences = coo_matrix((0,0)).todense()
    annotations        = list()
    #remove_singleton   = None
    neighbors_graph    = None
    
    @staticmethod
    def clear_Pangenome():
        Pangenome.presences_absences = coo_matrix((0,0)).todense()
        Pangenome.annotations        = list()
        #Pangenome.remove_singleton   = None
        neighbors_graph              = None
    
    def __init__(self, init_from = "args", *args):
        self.nb_organisms         = 0
        self.organism_positions   = OrderedDict()
        self.familly_positions    = OrderedDict()
        self.annotation_positions = OrderedDict()
        self.circular_contig      = set()
        self.familly_ratio        = dict()
        self.core_list            = list()
        self.pan_size             = 0
        self.core_size            = 0
        self.classified_famillies = None
        self.class_ratio          = None
        self.k                    = None
        self.BIC                  = None # Bayesian Index Criterion
        self.weights              = None

        if init_from == "file":
            self.__initialize_from_files(*args)
        elif init_from == "args":
            (self.nb_organisms,
             self.organism_positions,
             self.familly_positions,
             self.annotation_positions) = args 
        else:
            raise ValueError("init_from parameter is required")
 
        self.pan_size = len(self.familly_positions)
        for fam in self.familly_positions.keys():
            conservation = (np.sum(np.count_nonzero(Pangenome.presences_absences[self.familly_positions[fam],list(self.organism_positions.values())]))/float(self.nb_organisms))
            if conservation == 1.00:
                self.core_size += 1
                self.core_list.append(fam)
            self.familly_ratio[fam] = conservation

    def __initialize_from_files(self, organisms_file, families_tsv_file):#, remove_singletons = False
        
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

        for line in organisms_file:
            elements = line.split()
            self.nb_organisms +=1
            Pangenome.annotations += self.__load_gff(elements[1], families, elements[0])
            if len(elements)>2:
                self.circular_contig += elements[2:len(elements)]

        self.circular_contig = set(self.circular_contig)
        #Pangenome.remove_singleton = remove_singletons

        Pangenome.presences_absences = coo_matrix((nb_families+1,self.nb_organisms), dtype=np.dtype('uint8')).todense()
        cpt_org = 0
        cpt_fam = 0
        cur_org = Pangenome.annotations[0][ORGANISM]
        self.organism_positions[cur_org]=cpt_org
        self.annotation_positions[cur_org]="0-"
        for i in range(0, len(Pangenome.annotations)):
            if Pangenome.annotations[i][ORGANISM] != cur_org:
                self.annotation_positions[cur_org] = intspan(self.annotation_positions[cur_org]+str(i-1))
                cur_org = Pangenome.annotations[i][ORGANISM]
                self.annotation_positions[cur_org] = str(i)+"-"
                cpt_org += 1
                self.organism_positions[cur_org]=cpt_org
            fam_id = Pangenome.annotations[i][FAMILLY_CODE]
            self.familly_positions[fam_id] = int(fam_id)
            Pangenome.presences_absences[int(fam_id),cpt_org] += 1   

        self.annotation_positions[cur_org] = intspan(self.annotation_positions[cur_org]+str(i))
        self.familly_positions             = OrderedDict(sorted(self.familly_positions.items(), key=lambda t: t[0]))

    def __load_gff(self, gff_file, groups, organism):

        db_gff = gffutils.create_db(gff_file, ':memory:')
        annot = []
        for row in db_gff.all_features(featuretype='CDS', order_by=('seqid','start')):
            protein = row.id
            annot.append([protein,"CDS",organism,row.seqid,groups[protein],row.start,row.end,row.strand])
        return(annot)

    def __str__(self):

        pan_str = ""
        pan_str += "----------- Statistics -----------\n"
        pan_str += "Number of organisms: "+str(self.nb_organisms)+"\n"
        pan_str += "Pan-genome size:"+str(self.pan_size)+"\n"
        pan_str += "Exact core-genome size:"+str(self.core_size)+"\n"
        pan_str += "Exact variable-genome size:"+str(self.pan_size-self.core_size)+"\n"
        if self.classified_famillies is not None:
            pan_str += "Classification in "+str(self.k)+" classes: \n"
            for class_i in self.classified_famillies.keys():
                pan_str += "\t> # of families in class "+class_i+": "+str(len(self.classified_famillies[class_i]))+"\n"
        else:
            pan_str += "No classification have been performed on this Pangenome instance\n"
        pan_str += "----------------------------------\n"

        return(pan_str)    

    def assign_weights(self, weights):
        """ weights must be a dictionary having organism names as key and weights > to 0.0 and <= to 1.0"""
        logging.getLogger().debug(len(self.organism_positions))
        logging.getLogger().debug(len(weights)) 
        if(len(weights) == len(self.organism_positions)):
            if(weights.keys() not in self.organism_positions.keys()):
                if all([True if value <=1.0 and value >0 else False for value in weights.values()]):
                    self.weights = weights
                else: 
                    raise ValueError("weights must be include in ]0.0;1.0]")
            else:
                raise ValueError("Organism keys in the weights arguments contained in this object")
        else:
            raise ValueError("weights argument must have the same length than the number of organisms in this object")

    def sub_pangenome(self, sub_organisms):
       
        if set(sub_organisms).issubset(set(self.organism_positions.keys())):
            if(len(set(sub_organisms)) != len(sub_organisms)):
                logging.getLogger().warning("sub_organisms contain duplicated organism names. Only unique organism names will be used")
                sub_organisms = set(sub_organisms)

            sub_organism_positions   = OrderedDict((org, self.organism_positions[org]) for org in sorted(sub_organisms))
            logging.getLogger().debug(sub_organism_positions) 
            sub_annotation_positions = OrderedDict((org, self.annotation_positions[org]) for org in sorted(sub_organisms))
            logging.getLogger().debug(sub_annotation_positions) 
            subset_familly_code = set()
            for org in sub_organisms:
                for pos in sub_annotation_positions[org]:
                    subset_familly_code.add(Pangenome.annotations[pos][FAMILLY_CODE])
            logging.getLogger().debug(subset_familly_code)
            subset_familly_code = subset_familly_code - set([None])
            logging.getLogger().debug(subset_familly_code)
            sub_familly_positions = OrderedDict((fam, self.familly_positions[fam]) for fam in sorted(subset_familly_code))
            sub_pangenome = Pangenome("args", len(sub_organisms), sub_organism_positions,sub_familly_positions,sub_annotation_positions)
            if self.weights is not None:
                sub_pangenome.weights = dict((org, self.weights[org]) for org in sorted(sub_organisms)) 

            return(sub_pangenome)
        else:
            raise ValueError("Organisms contained in the sub_organisms arguments are not a subset of the organisms contained in this object") 

    def classify(self, result_path, k = 3, use_neighborhood = True, write_graph = None, neighbor_jumps = 1):
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
            #NEM requires 5 files: file.index, file.str, file.dat, file.ck (optional) and file.nei
            os.makedirs(result_path)
            index = OrderedDict()
            index_file = open(result_path+"/file.index", "w")
            logging.getLogger().info("Writing file.index file")
            
            for i, fam in enumerate(self.familly_positions.keys(),1):
                index[fam]=i
                index_file.write(str(i)+"\t"+str(fam)+"\n")

            index_file.close()

            str_file = open(result_path+"/file.str", "w")
            logging.getLogger().info("Writing file.str file")
            str_file.write("S\t"+str(self.pan_size)+"\t"+str(self.nb_organisms)+"\n")
            str_file.close()

            nei_file = open(result_path+"/file.nei", "w")
            dat_file = open(result_path+"/file.dat", "w")
            ck_file  = open(result_path+"/file.ck", "w")
            logging.getLogger().info("Writing file.nei file.dat and file.ck files")
            if use_neighborhood:
                if Pangenome.neighbors_graph is None: 
                    logging.getLogger().info("Start neighborhood graph construction")
                    untangeled_neighbors_graph = self.__neighborhoodComputation(15)
                else:
                    logging.getLogger().info("Use already computed neighbors")
                nei_file.write("1\n")
            else:
                nei_file.write("0\n")

            ck_file.write(str(k)+"\n")
        else:
            raise ValueError("result_path already exist")

        threshold = round(float(1)/self.nb_organisms,4) if self.nb_organisms >1 else 0
        
        if write_graph is not None:
            if use_neighborhood == False:
                raise ValueError("use_neighborhood must be True to write graph")
        
        for fam, fam_id in self.familly_positions.items():
            #logging.getLogger().debug(self.organism_positions.values())
            if self.weights is None:
                dat_file.write("\t".join(["1" if Pangenome.presences_absences[fam_id,p_org]>0 else "0" for p_org in self.organism_positions.values()])+"\n")
            else:
                dat_file.write("\t".join([str(round(self.weights[org],3)) if Pangenome.presences_absences[fam_id,p_org]>0 else "0.0" for org, p_org in self.organism_positions.items()])+"\n")
            ck_value = 0
            if self.familly_ratio[fam] == float(1):
                ck_value = 1
            elif self.familly_ratio[fam] <= threshold:
                ck_value = k
            ck_file.write(str(ck_value)+"\n")
            if use_neighborhood:
                row_fam         = []
                row_dist_score  = []
                neighbor_number = 0
                try:
                    for neighbor in nx.all_neighbors(Pangenome.neighbors_graph, fam):
                        distance_score = Pangenome.neighbors_graph[fam][neighbor]["weight"]/self.nb_organisms
                        row_fam.append(str(index[neighbor]))
                        row_dist_score.append(str(round(distance_score,4)))
                        neighbor_number += 1
                    if neighbor_number>0:
                        nei_file.write("\t".join([str(item) for sublist in [[index[fam]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
                        logging.getLogger().debug("\t".join([str(item) for sublist in [[index[fam]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
                    else:
                        raise nx.exception.NetworkXError("no all_neighbors in selected organismss")
                except nx.exception.NetworkXError as nxe:
                    logging.getLogger().warning("The familly: "+fam+" is an isolated familly")
                    nei_file.write(str(index[fam])+"\t0\n")

        nei_file.close()
        dat_file.close()
        ck_file.close()
        logging.getLogger().info("Running NEM...")

        #bernouli -> no weight or normal -> weight
        model = "bern" if self.weights is None else "norm"
        print_log = " -l y" if logging.getLogger().getEffectiveLevel() < 20 else "" 
        command = NEM_LOCATION+"nem_exe "+result_path+"/file "+str(k)+" -a nem -i 2000 -m "+model+" pk sk_ -s r 10 -B fix -b "+("1" if use_neighborhood else "0")+" -T -O random"+print_log
        logging.getLogger().info(command)
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = proc.communicate()[1]
        # M = float(re.search("\(M = (.+?)\)",output).group(1))# M=markov ps-like
        # if model == "norm":
        #     self.BIC = 2 * M - (k * (self.nb_organisms + 1) + (1 if use_neighborhood else 0) + k- 1) * math.log(len(self.familly_positions)) 
        # elif model == "bern":
        #     self.BIC = 2 * M - (k * self.nb_organisms + (1 if use_neighborhood else 0) + k - 1) * math.log(len(self.familly_positions))
        # logging.getLogger().info("Based on "+str(k)+" classes, BIC = "+str(round(self.BIC,4)))
        if os.path.isfile(result_path+"/file.cf"):
            logging.getLogger().info("Reading NEM results")
        else:
            logging.getLogger().error("No NEM output found in file.cf")
        with open(result_path+"/file.cf","r") as classification_file:
            classification = classification_file.readline().split()
            print(index.keys())
            classification = {k: int(v) for k, v in zip(index.keys(), classification)}
            self.k = k
            if write_graph is not None:
                logging.getLogger().info("Writing graphML file")
                for node in list(untangeled_neighbors_graph.node):
                    try:
                        logging.getLogger().debug(node+": "+str(classification[elements[0]]))
                        untangeled_neighbors_graph.node[node]['nem_class'] = classification[node]
                    except:
                        elements = str(node).split("_")
                        logging.getLogger().debug(elements[0]+": "+str(classification[elements[0]]))
                        untangeled_neighbors_graph.node[node]['nem_class'] = classification[elements[0]]

                getattr(nx,'write_'+write_graph)(untangeled_neighbors_graph,result_path+"/file."+write_graph)
                # mst = nx.maximum_spanning_tree(graph)
                # for u,v,d in mst.edges(data=True):
                #     d["MST"]=True
                # getattr(nx,'write_'+write_graph)(nx.disjoint_union(graph,mst),result_path+"/file_mst."+write_graph)
        self.classified_famillies = dict(zip([str(i_k) for i_k in range(1,k+1)],[[]] * k))
        self.class_ratio          = dict(zip(self.classified_famillies.keys(),[float("0.0")] * k))
        for i_k in self.classified_famillies.keys():
                self.classified_famillies[i_k] = [fam for (fam, nem_class) in classification.items() if nem_class == i_k]
                self.class_ratio[i_k] = [self.familly_ratio[fam] for fam in self.classified_famillies[i_k]]
                self.class_ratio[i_k] = sum(self.class_ratio[str(i_k)])/len(self.class_ratio[str(i_k)]) if len(self.class_ratio[i_k]) > 0 else 0   
        #logging.getLogger().debug(sorted(self.class_ratio, key=self.class_ratio.__getitem__, reverse=True))        
        self.classified_famillies = OrderedDict((o, self.classified_famillies[o]) for o in sorted(self.class_ratio, key=self.class_ratio.__getitem__, reverse=True))
        #logging.getLogger().debug(self.classified_famillies)
        logging.getLogger().info("Writing stats summury")
        with open(result_path+"/nem.stat","w") as nem_stat_file:
            nem_stat_file.write(str(self.nb_organisms)+"\t"+"\t".join([str(len(fams)) for nem_class, fams in self.classified_famillies.items()])+"\t"+str(self.pan_size)+"\n")
        with open(result_path+"/exact.stat","w") as exact_stat_file:        
            exact_stat_file.write(str(self.nb_organisms)+"\t"+str(self.core_size)+"\t"+str(self.pan_size-self.core_size)+"\t"+str(self.pan_size)+"\n")    

        #return(graph)

    def __neighborhoodComputation(self, max_degree = float("Inf"), high_degree_node = set()):#initNeighborDistance, maxNeighborDistance,
        
        i = 0 # = compteur remis a zero chaque fois qu'on change de sequence
        neighbors_graph = nx.Graph()
        S_id_prev      = Pangenome.annotations[0][CONTIG_ID]
        O_id_nei       = Pangenome.annotations[0][ORGANISM]
        familly_id_nei = Pangenome.annotations[0][FAMILLY_CODE]
        gene_nei       = Pangenome.annotations[0][GENE]
        prec           = None
        circularize    = True if S_id_prev in self.circular_contig else False

        all_paths = defaultdict(lambda : defaultdict(list))  
        #path_organisms_gene = defaultdict(lambda : defaultdict(list)) 

        def add_neighbors(fam_id, fam_nei, prec,org,gene):
            neighbors_graph.add_node(fam_id)
            neighbors_graph.add_node(fam_nei)

            if prec is not None:
                # if gene == "1280.PRJNA239544.CP007447_03016":
                #     logging.getLogger().debug([prec,fam_nei,fam_id])
                #     exit()
                all_paths[fam_nei][frozenset([prec,fam_id])].append(tuple([org,gene]))
                #path_organisms_gene[([frozenset([prec,fam_id])])][fam_nei].append(O_id)

            try:
                neighbors_graph.node[fam_id][org]+=" "+gene
            except:
                neighbors_graph.node[fam_id][org]=gene

            if not neighbors_graph.has_edge(fam_id,fam_nei):
                neighbors_graph.add_edge(fam_id, fam_nei)
            try:
                neighbors_graph[fam_id][fam_nei]["weight"]+=1
            except KeyError:
                neighbors_graph[fam_id][fam_nei]["weight"]=1
            try:
                neighbors_graph[fam_id][fam_nei][org]+=1
            except:
                neighbors_graph[fam_id][fam_nei][org]=1

        projection = "1280.PRJNA239544.CP007447"
        i_first = None
        i_last = None
        logging.getLogger().debug(Pangenome.annotations[0])
        for index, row in enumerate(Pangenome.annotations[1:]):
            familly_id = row[FAMILLY_CODE]
            gene       = row[GENE]
            S_id       = row[CONTIG_ID]
            logging.getLogger().debug(row)
       

            if i_first is None and O_id_nei==projection:
                i_first = index
            if i_first is not None and i_last is None and O_id_nei!=projection:
                i_last = index
            if familly_id in high_degree_node or familly_id_nei in high_degree_node:
                logging.getLogger().debug(familly_id)
                i+=1
                continue
            if (S_id != S_id_prev):
                if circularize:
                    familly_id = Pangenome.annotations[index-i+1][FAMILLY_CODE]
                    add_neighbors(familly_id, familly_id_nei, prec, O_id_nei, gene_nei)
                circularize = True if S_id in self.circular_contig else False
                i=0
                prec = None
            else:
                i+=1
                add_neighbors(familly_id, familly_id_nei, prec, O_id_nei, gene_nei)
                prec = familly_id_nei
            
            familly_id_nei = familly_id
            gene_nei       = gene
            O_id_nei       = row[ORGANISM]
            S_id_prev      = S_id

        if circularize:
            add_neighbors(familly_id, familly_id_nei, prec, O_id_nei, gene_nei)
        if i_first is not None and i_last is None and O_id_nei!=projection:
            i_last = index

        new_high_degree_node = set([node for node, degree in neighbors_graph.degree_iter() if degree>max_degree])
        if len(new_high_degree_node)>max_degree:
            logging.getLogger().debug(len(new_high_degree_node))
            neighbors_graph = self.__neighborhoodComputation(max_degree, high_degree_node.union(new_high_degree_node))
            return neighbors_graph
        else:

            untangeled_neighbors_graph = neighbors_graph.copy()

            #inspired from http://stackoverflow.com/a/9114443/7500030
            def merge_overlapping_path(data):
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


            inversed_dict = dict()

            #untangling stage:
            for node in list(untangeled_neighbors_graph.node):
                path_groups = merge_overlapping_path(all_paths[node])
                logging.getLogger().debug(node)
                logging.getLogger().debug(all_paths[node])
                logging.getLogger().debug(path_groups)
                logging.getLogger().debug(len(path_groups))
                if len(path_groups)>1:
                    logging.getLogger().debug("split "+str(len(path_groups)))
                    for suffix, path_group in enumerate(path_groups):
                        new_node = node+"_"+str(suffix)
                        all_paths[new_node][frozenset(path_group)]= -1
                        logging.getLogger().debug("new_node:"+new_node)
                        #all_paths[new_node][frozenset(path_group)]+=1#change number
                        neibors_new_node = [neighbor for neighbor in nx.all_neighbors(untangeled_neighbors_graph, node) if neighbor in path_group]
                        
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
                        untangeled_neighbors_graph.add_node(new_node)
                        for neighbor in path_group:
                            print(all_paths[node])
                            
                            untangeled_neighbors_graph.add_edge(new_node,neighbor)
                            weight=0
                            for path, genes in all_paths[node].items():
                                if neighbor in path:
                                    for org, gene in genes:
                                        try:
                                            untangeled_neighbors_graph.node[new_node][org]+=" "+gene
                                        except:
                                            untangeled_neighbors_graph.node[new_node][org]=gene
                                        try:
                                            neighbors_graph[new_node][neighbor][org]+=1
                                        except:
                                            neighbors_graph[new_node][neighbor][org]=1
                                        inversed_dict.update({gene:new_node})
                                        weight+=1
                            untangeled_neighbors_graph[new_node][neighbor]["weight"]=weight
                    untangeled_neighbors_graph.remove_node(node)

            for index, row in enumerate(Pangenome.annotations[i_first:i_last-1]):
                fam = row[FAMILLY_CODE]
                if fam not in high_degree_node.union(new_high_degree_node):
                    try:
                        untangeled_neighbors_graph.node[fam]["order"]=int(index)
                    except:
                        try:
                            logging.getLogger().debug(inversed_dict[row[GENE]])
                            #untangeled_neighbors_graph.node[inversed_dict[row[GENE]]]["order"]=int(index)
                        except:
                            logging.getLogger().warning(row)

                        # prev_fam = Pangenome.annotations[i_first+index-1][FAMILLY_CODE]
                        # next_fam = Pangenome.annotations[i_first+index+1][FAMILLY_CODE]
                        # com = nx.neighbors(untangeled_neighbors_graph,prev_fam
                        # com2 = nx.neighbors(untangeled_neighbors_graph,next_fam)
                        # logging.getLogger().debug(com)
                        # logging.getLogger().debug(com2)
                        # exit()
                        #untangeled_neighbors_graph.node[inversed_dict[tuple([row[ORGANISM],row[GENE]])]]["order"]=int(index)
                    #if i_first+index-1 >= i_first and i_first+index+1 <=i_last-1:
                    #     prev_fam = Pangenome.annotations[i_first+index-1][FAMILLY_CODE]
                    #     next_fam = Pangenome.annotations[i_first+index+1][FAMILLY_CODE]
                    #     com = nx.common_neighbors(untangeled_neighbors_graph,prev_fam, next_fam)
                    #     untangeled_neighbors_graph.node[next(com)]["order"]=int(index)
                    #     if len(list(com))>2:
                    #         logging.getLogger().debug(list(com))
                    # else:
                    #     logging.getLogger().debug(index)
        Pangenome.neighbors_graph = neighbors_graph
        return untangeled_neighbors_graph

        #logging.getLogger().debug(list(neighbors_graph.degree_iter()))
        #logging.getLogger().debug(len(high_degree_node))
        #neighbors_graph.remove_nodes_from(high_degree_node)

        #logging.getLogger().debug(nx.number_of_nodes(neighbors_graph))
        #logging.getLogger().debug(nx.number_of_edges(neighbors_graph))
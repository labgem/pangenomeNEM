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
        self.class_ratio          = None
        self.k                    = None
        self.BIC                  = None # Bayesian Index Criterion

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

    def __load_gff(self, gff_file, groups, organism):
        logging.getLogger().info("Reading "+gff_file+" file ...")
        db_gff = gffutils.create_db(gff_file, ':memory:')
        annot = defaultdict(list)
        for row in db_gff.all_features(featuretype='CDS', order_by=('seqid','start')):
            logging.getLogger().debug(row)
            protein = row.id
            annot[row.seqid].append([protein,"CDS",organism,groups[protein],row.start,row.end,row.strand])
        return(annot)

    # def __str__(self):

    #     pan_str = ""
    #     pan_str += "----------- Statistics -----------\n"
    #     pan_str += "Number of organisms: "+str(self.nb_organisms)+"\n"
    #     pan_str += "Pan-genome size:"+str(self.pan_size)+"\n"
    #     pan_str += "Exact core-genome size:"+str(self.core_size)+"\n"
    #     pan_str += "Exact variable-genome size:"+str(self.pan_size-self.core_size)+"\n"
    #     if self.classified_famillies is not None:
    #         pan_str += "Classification in "+str(self.k)+" classes: \n"
    #         for class_i in self.classified_famillies.keys():
    #             pan_str += "\t> # of families in class "+class_i+": "+str(len(self.classified_famillies[class_i]))+"\n"
    #     else:
    #         pan_str += "No classification have been performed on this Pangenome instance\n"
    #     pan_str += "----------------------------------\n"

    #     return(pan_str)    

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

        else:
            raise ValueError("result_path already exist")

        logging.getLogger().info("Writing file.str file.index file.nei file.dat and file.ck files")
        index_file = open(result_path+"/file.index", "w")
        #ck_file  = open(result_path+"/file.ck", "w")
        nei_file = open(result_path+"/file.nei", "w")
        dat_file = open(result_path+"/file.dat", "w")

        if use_neighborhood:
            nei_file.write("1\n")
        else:
            nei_file.write("0\n")
        
        index = {node: index+1 for index, node in enumerate(self.neighbors_graph.nodes(data=False))}

        str_file = open(result_path+"/file.str", "w")
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
                        distance_score = self.neighbors_graph[node_name][neighbor]["weight"]/self.nb_organisms
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
                    nei_file.write(str(index[fam])+"\t0\n")

        index_file.close()
        #ck_file.close()
        nei_file.close()
        dat_file.close()             

        logging.getLogger().info("Running NEM...")

        #bernouli -> no weight or normal -> weight
        model = "bern" #if self.weights is None else "norm"
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
            for node, nem_class in zip(self.neighbors_graph.nodes(), classification):
                self.neighbors_graph.node[node]["nem_class"]=nem_class
            if write_graph is not None:
                logging.getLogger().info("Writing graphML file")
                getattr(nx,'write_'+write_graph)(self.neighbors_graph,result_path+"/file."+write_graph)
        # self.classified_famillies = dict(zip([str(i_k) for i_k in range(1,k+1)],[[]] * k))
        # self.class_ratio          = dict(zip(self.classified_famillies.keys(),[float("0.0")] * k))
        # for i_k in self.classified_famillies.keys():
        #         self.classified_famillies[i_k] = [fam for (fam, nem_class) in classification.items() if nem_class == i_k]
        #         self.class_ratio[i_k] = [self.familly_ratio[fam] for fam in self.classified_famillies[i_k]]
        #         self.class_ratio[i_k] = sum(self.class_ratio[str(i_k)])/len(self.class_ratio[str(i_k)]) if len(self.class_ratio[i_k]) > 0 else 0   
        # #logging.getLogger().debug(sorted(self.class_ratio, key=self.class_ratio.__getitem__, reverse=True))        
        # self.classified_famillies = OrderedDict((o, self.classified_famillies[o]) for o in sorted(self.class_ratio, key=self.class_ratio.__getitem__, reverse=True))
        # #logging.getLogger().debug(self.classified_famillies)
        # logging.getLogger().info("Writing stats summury")
        # with open(result_path+"/nem.stat","w") as nem_stat_file:
        #     nem_stat_file.write(str(self.nb_organisms)+"\t"+"\t".join([str(len(fams)) for nem_class, fams in self.classified_famillies.items()])+"\t"+str(self.pan_size)+"\n")
        # with open(result_path+"/exact.stat","w") as exact_stat_file:        
        #     exact_stat_file.write(str(self.nb_organisms)+"\t"+str(self.core_size)+"\t"+str(self.pan_size-self.core_size)+"\t"+str(self.pan_size)+"\n")    

        #return(graph)

    def __neighborhoodComputation(self, max_degree = float("Inf"), high_degree_node = set()):#initNeighborDistance, maxNeighborDistance,

        neighbors_graph = nx.Graph()
        all_paths = defaultdict(lambda : defaultdict(list))  

        def add_neighbors(fam_id, fam_nei, prec,org,gene):
            neighbors_graph.add_node(fam_id)
            neighbors_graph.add_node(fam_nei)

            if prec is None:
                logging.getLogger().debug([fam_id, fam_nei, prec,org,gene])

            if prec is not None:
                all_paths[fam_nei][frozenset([prec,fam_id])].append(tuple([org,gene]))
            else:
                all_paths[fam_nei][None].append(fam_id)
            try:
                neighbors_graph.node[fam_id][org]+=" "+gene
            except KeyError:
                neighbors_graph.node[fam_id][org]=gene
            if not neighbors_graph.has_edge(fam_id,fam_nei):
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
                prec  = None
                first = None
                shift = 0
                while True:
                    row = contig_annot[shift]
                    familly_id_nei = row[FAMILLY_CODE]
                    gene_nei = row[GENE]
                    if familly_id_nei not in high_degree_node:
                       first = shift
                       break
                    else:
                        shift+=1

                logging.getLogger().debug(row)
                
                for index, row in enumerate(contig_annot[first+1:]):
                    logging.getLogger().debug(row)        
                    gene       = row[GENE]
                    familly_id = row[FAMILLY_CODE]
                    skipped_node = set()
                    if familly_id not in high_degree_node:
                        add_neighbors(familly_id, familly_id_nei, prec, organism, gene_nei)
                        prec           = familly_id_nei
                        familly_id_nei = familly_id
                        gene_nei       = gene
                    else:
                        logging.getLogger().debug("skipped"+gene)
                        skipped_node.add(index)
                
                row = contig_annot[first]
                familly_id = row[FAMILLY_CODE]
                gene = row[GENE]
                
                if contig in self.circular_contig:
                    add_neighbors(familly_id, familly_id_nei, prec, organism, gene_nei)
                    logging.getLogger().debug("first2 "+familly_id)
                    logging.getLogger().debug(shift)
                    logging.getLogger().debug("prec "+str(all_paths[familly_id]))
                    prec = all_paths[familly_id].pop(None)[0]
                    all_paths[familly_id][frozenset([familly_id_nei,prec])].append(tuple([organism,gene]))
                else:
                    logging.getLogger().debug(contig)
                    all_paths[familly_id].pop(None)[0]

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
            logging.getLogger().debug(all_paths)
            #untangling stage:
            for node in list(untangeled_neighbors_graph.node):
                logging.getLogger().debug(node)
                logging.getLogger().debug(all_paths[node])
                path_groups = merge_overlapping_path(all_paths[node])

                logging.getLogger().debug(path_groups)
                logging.getLogger().debug(len(path_groups))
                if len(path_groups)>1:
                    logging.getLogger().debug("split "+str(len(path_groups)))
                    for suffix, path_group in enumerate(path_groups):
                        new_node = node+"_"+str(suffix)
                        all_paths[new_node][frozenset(path_group)]= -1
                        logging.getLogger().debug("new_node:"+new_node)
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
                                            untangeled_neighbors_graph.node[new_node][org]=+" "+gene
                                        except:
                                            untangeled_neighbors_graph.node[new_node][org]=gene
                                        try:
                                            untangeled_neighbors_graph[new_node][neighbor][org]+=1
                                        except:
                                            untangeled_neighbors_graph[new_node][neighbor][org]=1
                                        inversed_dict.update({gene:new_node})
                                        weight+=1
                            untangeled_neighbors_graph[new_node][neighbor]["weight"]=weight
                    untangeled_neighbors_graph.remove_node(node)

            for organism, annot_contigs in self.annotations.items():
                for contig, contig_annot in annot_contigs.items():
                    if contig not in self.circular_contig:
                            continue
                    order = 0
                    for index, row in enumerate(contig_annot):
                        fam = row[FAMILLY_CODE]
                        if fam in high_degree_node.union(new_high_degree_node):
                            continue
                        try:
                            untangeled_neighbors_graph.node[fam]["order_"+contig.replace(".", "_")]=order
                        except:
                            try:
                                logging.getLogger().debug(inversed_dict[row[GENE]])
                                untangeled_neighbors_graph.node[inversed_dict[row[GENE]]]["order_"+contig.replace(".", "_")]=order
                            except:
                                logging.getLogger().warning(row)
                        order+=1
        return untangeled_neighbors_graph
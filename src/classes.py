#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

#import mysql.connector as mysql
import pandas as pd
import numpy as np
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
import shlex

NEM_LOCATION  = "../NEM/"
(GENE, TYPE, ORGANISM, CONTIG_ID, FAMILLY_CODE, START, END, STRAND) = range(0, 8)

def __dbConnect():
    try:
        conn = mysql.connect(host="mysqlagcdb.genoscope.cns.fr",port=3306,user="agc",passwd="admagc21",db="pkgdb_dev")
    except:
        print("Connexion with pkgdb_dev failed.")
        exit(1)
    return conn

class Pangenome:

    presences_absences = OrderedDict()
    annotations        = tuple()
    remove_singleton   = False
    neighbors          = None
    neighbor_jumps     = None
    
    @staticmethod
    def clear_Pangenome():
        Pangenome.presences_absences = OrderedDict()
        Pangenome.annotations        = tuple()
        Pangenome.remove_singleton   = False
        neighbors                    = None
        neighbor_jumps               = None
    
    def __init__(self, init_from = "args", *args):
        self.nb_organisms         = 0
        self.organism_positions   = OrderedDict()
        self.familly_positions    = OrderedDict()
        self.annotation_positions = OrderedDict()
        self.familly_ratio        = dict()
        self.core_list            = list()
        self.pan_size             = 0
        self.core_size            = 0
        self.classified_famillies = None
        self.class_ratio          = None
        self.k                    = None
        self.BIC                  = None # Bayesian Index Criterion
        self.weights              = None

        if init_from == "progenome":
            logging.getLogger().info("Load progenome files...")
            Pangenome.clear_Pangenome()
            self.__initialize_from_progenome(*args)
        elif init_from == "microscope":
            Pangenome.clear_Pangenome()
            self.__initialize_from_microscope(*args)
        elif init_from == "prokka/roary":                
            Pangenome.clear_Pangenome()
            self.__initialize_from_prokka_roary(*args)
        elif init_from == "prokka/MMseqs2":
            Pangenome.clear_Pangenome()
            self.__initialize_from_prokka_mmseqs2(*args)
        elif init_from == "args":
            (self.nb_organisms,
             self.organism_positions,
             self.familly_positions,
             self.annotation_positions) = args 
        else:
            raise ValueError("init_from parameter is required")
 
        self.pan_size    = len(self.familly_positions)
        for fam in self.familly_positions.keys():
            conservation = (sum([1 if Pangenome.presences_absences[fam][p_org]>0 else 0 for p_org in self.organism_positions.values()])/float(self.nb_organisms))
            if conservation == 1.00:
                self.core_size += 1
                self.core_list.append(fam)
            self.familly_ratio[fam] = conservation

    def __initialize_from_microscope(self):
        
        Oids = [int(o) for o in self.options.organisms]
        Oids.sort()
        Oids_str = ','.join(str(Oid) for Oid in Oids)
        Oids_bracket_str = '),('.join(str(Oid) for Oid in Oids)
        Oids_size = len(Oids)

        conn = _dbConnect()
        cur = conn.cursor()

                
        if self.options.artefact:
            ArtefactOpt="0,1"
            print("> Artefacts are taken into account in this analysis.")
        else:
            ArtefactOpt="0"
            print("> Artefacts are excluded from the analysis.")

        try:
            
            cur.execute("""SELECT DISTINCT O_id FROM Organism INNER JOIN
                    Replicon USING(O_id) 
                    INNER JOIN Sequence USING (R_id) 
                    WHERE Replicon.O_id IN (%s) 
                    AND S_status = 'inProduction';""" % Oids_str)


        except mysql.Error as e:
            conn.close()
            print("Error %d: %s" % (e.args[0], e.args[1]))
            exit(1)
 
        rows = cur.fetchall() 
        
        Oids_bug = []
        Oids_ok = []
        if len(rows) < Oids_size:
            for row in rows:
                Oids_ok.append(str(row[0]))
            for Oid in Oids:
                if Oid not in Oids_ok:
                    Oids_bug.append(Oid)
            if len(Oids_bug)>1:
                Oids_bug_str = ','.join(str(Oid) for Oid in Oids_bug)
                print("/!\ Error: Following Oids: "+Oids_bug_str+" have S_status != inProduction")
            else:
                print("/!\ Error: Oid "+Oids_bug[0]+ " has S_status != inProduction")
            exit(1)

        
        self.organisms = set([row[0] for row in rows])
            
        try:    
            cur.execute("""CREATE TEMPORARY TABLE tmp_GO_all_NEMtest
            SELECT GO_ori_id as GO_id, GO_type, O_id, S_id, GO_begin, GO_end, SUBSTRING(GO_frame,1,1) as GO_strand, IF(GO_status='Artefact',1,0) AS isArtefact
            FROM Genomic_Object G
            INNER JOIN Sequence S USING (S_id) 
            INNER JOIN Replicon R USING (R_id)
            INNER JOIN Organism O USING (O_id)
            WHERE O_id IN (%s)
            AND S_status = 'inProduction'
            AND GO_type IN ('CDS', 'fCDS')
            AND GO_update = 'current';""" % Oids_str)
            cur.execute("""ALTER TABLE tmp_GO_all_NEMtest ADD PRIMARY KEY(GO_id);""")
            cur.execute("""ALTER TABLE tmp_GO_all_NEMtest ADD INDEX(O_id);""")

            cur.execute("""CREATE TEMPORARY TABLE tmp_GO_NEMtest
            SELECT * FROM tmp_GO_all_NEMtest
            WHERE GO_type IN ('CDS','fCDS');""")
            cur.execute("""ALTER TABLE tmp_GO_NEMtest ADD PRIMARY KEY(GO_id);""")
            cur.execute("""ALTER TABLE tmp_GO_NEMtest ADD INDEX(O_id);""")

            cur.execute("""CREATE TEMPORARY TABLE tmp_cluster_GO_NEMtest
            SELECT cluster_id, O_id, GO_id, isArtefact
            FROM tmp_GO_NEMtest
            INNER JOIN MICFAM_cluster USING(GO_id)
            WHERE MICFAM_param_id=%i;""" % self.options.micfamparam)
            cur.execute("""ALTER TABLE tmp_cluster_GO_NEMtest ADD PRIMARY KEY(GO_id);""")
            cur.execute("""ALTER TABLE tmp_cluster_GO_NEMtest ADD INDEX(cluster_id);""")

            cur.execute("""CREATE TEMPORARY TABLE tmp_clustO 
            SELECT cluster_id, O_id, if(sum(isArtefact)>0,1,0) as clustO_artefact
            FROM tmp_cluster_GO_NEMtest
            GROUP BY cluster_id, O_id;""")
            cur.execute("""ALTER TABLE tmp_clustO ADD PRIMARY KEY(cluster_id,O_id);""")

            cur.execute("""CREATE TEMPORARY TABLE tmp_clust 
            SELECT cluster_id, if(sum(clustO_artefact)>0,1,0) as cluster_artefact
            FROM tmp_clustO
            GROUP BY cluster_id;""")
            cur.execute("""ALTER TABLE tmp_clust ADD PRIMARY KEY(cluster_id);""")
                
            cur.execute("""CREATE TEMPORARY TABLE tmp_clust_final_NEMtest
            SELECT tmp_clust.cluster_id, tmp_clust.cluster_artefact
            FROM tmp_clust
            WHERE tmp_clust.cluster_artefact IN (%s);""" % ArtefactOpt)
            cur.execute("""ALTER TABLE tmp_clust_final_NEMtest ADD PRIMARY KEY(cluster_id);""")    

            cur.execute("""CREATE TEMPORARY TABLE tmp_Oid(O_id INT NOT NULL, PRIMARY KEY(O_id));""")

            cur.execute("""INSERT INTO tmp_Oid VALUES (%s)""" % Oids_bracket_str)
            if (self.options.remove_singleton):
                query = """SELECT tmp_Oid.O_id, tmp_clust_final_NEMtest.cluster_id, IF(tmp_clustO.O_id,1,0) AS Bool 
                FROM tmp_clust_final_NEMtest
                INNER JOIN tmp_Oid
                LEFT JOIN tmp_clustO USING(cluster_id,O_id)
                ORDER BY O_id, cluster_id;"""
            else:
                
                cur.execute("""CREATE TEMPORARY TABLE tmp_singletons 
                SELECT CONCAT("GO_",GO_id) as cluster_id, tmp_GO_NEMtest.O_id, GO_id
                FROM tmp_GO_NEMtest
                LEFT JOIN tmp_cluster_GO_NEMtest  USING(GO_id)
                WHERE tmp_cluster_GO_NEMtest .GO_id IS NULL
                AND tmp_GO_NEMtest.isArtefact IN (%s);""" % ArtefactOpt)

                cur.execute("""CREATE TEMPORARY TABLE tmp_clustO_singl
                SELECT cluster_id, O_id
                FROM tmp_clustO
                UNION
                SELECT cluster_id, O_id
                FROM tmp_singletons ;""")

                cur.execute("""ALTER TABLE tmp_clustO_singl ADD PRIMARY KEY(cluster_id,O_id);""")

                cur.execute("""CREATE TEMPORARY TABLE tmp_clust_singl_final
                SELECT cluster_id
                FROM tmp_clust_final_NEMtest
                UNION
                SELECT cluster_id
                FROM tmp_singletons ;""")

                query = """SELECT tmp_Oid.O_id, tmp_clust_singl_final.cluster_id, IF(tmp_clustO_singl.O_id,1,0) AS Bool 
                FROM tmp_clust_singl_final
                INNER JOIN tmp_Oid
                LEFT JOIN tmp_clustO_singl USING(cluster_id,O_id)
                ORDER BY O_id, cluster_id;"""

            rows_binary = pd.read_sql_query(query, conn)
        except mysql.Error as e:
            conn.close()
            print("Error %d: %s" % (e.args[0], e.args[1]))
            exit(1)

        
        if rows_binary.empty:
            print("Problem with data results")            
            exit(1)
        else:
            cluster_occurence_matrix = pd.pivot_table(rows_binary, index=['cluster_id'], columns=['O_id'],  aggfunc=np.sum, fill_value = 0).xs('Bool', axis=1, drop_level=True)
            
            self.families = set(cluster_occurence_matrix.index)

            #TODO passer les cases suprÃrieures a 1 a 0

            self.presence_absence_matrix = cluster_occurence_matrix
        try:

            cur.execute("""CREATE TEMPORARY TABLE tmp_GO
            SELECT * FROM tmp_GO_NEMtest
            WHERE O_id IN (%s);""" % Oids_str)
            cur.execute("""ALTER TABLE tmp_GO ADD PRIMARY KEY(GO_id);""")
        
            cur.execute("""CREATE TEMPORARY TABLE tmp_GO_all
            SELECT * FROM tmp_GO_all_NEMtest
            WHERE O_id IN (%s);""" % Oids_str)
            cur.execute("""ALTER TABLE tmp_GO_all ADD PRIMARY KEY(GO_id);""")

            cur.execute("""CREATE TEMPORARY TABLE tmp_cluster_GO 
            SELECT cluster_id, O_id, GO_id, isArtefact
            FROM tmp_cluster_GO_NEMtest
            WHERE O_id IN (%s);""" % Oids_str)
            cur.execute("""ALTER TABLE tmp_cluster_GO  ADD PRIMARY KEY(GO_id);""")
            cur.execute("""ALTER TABLE tmp_cluster_GO  ADD INDEX(cluster_id);""")

            cur.execute("""CREATE TEMPORARY TABLE tmp_clust_2 
            SELECT cluster_id, if(sum(isArtefact)>0,1,0) as cluster_artefact
            FROM tmp_cluster_GO
            GROUP BY cluster_id;""")
            cur.execute("""ALTER TABLE tmp_clust_2 ADD PRIMARY KEY(cluster_id);""")
        
            if (self.options.remove_singleton):
                cur.execute("""CREATE TEMPORARY TABLE tmp_clust_final_2
                SELECT GO_id, cluster_id
                FROM tmp_clust_2 
                INNER JOIN tmp_cluster_GO USING (cluster_id)
                INNER JOIN tmp_clust_final_NEMtest USING(cluster_id)
                WHERE tmp_clust_2.cluster_artefact IN (%s);""" % ArtefactOpt)
            else:
                cur.execute("""CREATE TEMPORARY TABLE tmp_singletons_2
                SELECT CONCAT("GO_",GO_id) as cluster_id, tmp_GO.O_id, GO_id
                FROM tmp_GO
                LEFT JOIN tmp_cluster_GO USING(GO_id)
                WHERE tmp_cluster_GO.GO_id IS NULL
                AND tmp_GO.isArtefact IN (%s);""" % ArtefactOpt)

                cur.execute("""CREATE TEMPORARY TABLE tmp_clust_final_2
                SELECT GO_id, cluster_id
                FROM tmp_clust_2 
                INNER JOIN tmp_cluster_GO USING (cluster_id)
                INNER JOIN tmp_clust_final_NEMtest USING(cluster_id)
                WHERE tmp_clust_2.cluster_artefact IN (%s)
                UNION
                SELECT GO_id, cluster_id
                FROM tmp_singletons_2;""" % ArtefactOpt)

            cur.execute("""ALTER TABLE tmp_clust_final_2 ADD PRIMARY KEY(GO_id);""")
            
            self.annotations = pd.read_sql_query("""SELECT GO_id, GO_type, O_id, S_id, cluster_id, GO_begin, GO_end, GO_strand 
            FROM tmp_GO_all 
            LEFT JOIN tmp_clust_final_2 USING(GO_id)
            WHERE isArtefact IN (%s)
            AND cluster_id IS NOT NULL
            ORDER BY O_id, S_id, GO_begin;""" % ArtefactOpt, conn) 

            self.annotations.columns = ['GENE','TYPE','ORGANISM','CONTIG_ID','FAMILLY_CODE','START','END','STRAND']
            self.annotations.set_index(keys = 'GENE', drop = False, inplace = True)
            
            
        except mysql.Error as e:
            print("Error %d: %s" % (e.args[0], e.args[1]))
            conn.close()
            exit(1)
    
    def __initialize_from_progenome(self,
                                    annotations_file,
                                    eggNOG_clusters_file,
                                    remove_singletons = False):
                
        useful_cols_annot       = ['GENE_ID', 'CONTIG_ID', 'TYPE', 'START', 'END', 'STRAND']
        type_cols_annot         = {'GENE_ID':"str",'CONTIG_ID':"str",'TYPE':"str",'START':"int",'END':"int",'STRAND':"str"}
        useful_cols_orthoGroups = ['PROTEIN_ID','EGGNOG_GROUP_CODE']
        type_cols_orthoGroups   = {'PROTEIN_ID':"str",'EGGNOG_GROUP_CODE':"str"}

        try:
            annotations    = pd.read_csv(annotations_file, sep="\t", usecols=useful_cols_annot, dtype=type_cols_annot)
            familly_groups = pd.read_csv(eggNOG_clusters_file, sep="\t", usecols=useful_cols_orthoGroups, dtype=type_cols_orthoGroups)
        except ValueError as verr:
            if str(verr) == "Usecols do not match names.":
                raise ValueError("Column names do not match expected progenome names. Did you swapped annatotation file with eggNOGclusters file ?")
            else:
                raise verr

        familly_groups.rename(columns={'PROTEIN_ID':'PROTEIN_ID','EGGNOG_GROUP_CODE':'FAMILLY_CODE'}, inplace=True)

        #split a column with 3 fields separated by a point in 3 other columns
        #expand = True make pandas to return a DataFrame with the 3 columns comming from the splitted GENE_ID column
        annotations = pd.concat([annotations,annotations.loc[:,('GENE_ID')].str.split('.', expand = True).rename(columns={0:'TAX_ID',1:'PROJET',2:'GENE'})],axis=1)    
        annotations = pd.concat([annotations,annotations.loc[:,('CONTIG_ID')].str.split('.', expand = True).loc[:,2].to_frame("CONTIG")],axis=1)
        annotations.loc[:,('ORGANISM')] = annotations.loc[:,('TAX_ID')] + "." + annotations.loc[:,('PROJET')]

        annotations    = annotations.loc[:,('GENE','TYPE','ORGANISM','CONTIG_ID','START','END','STRAND')]
        
        familly_groups = pd.concat([familly_groups,familly_groups.loc[:,('PROTEIN_ID')].str.split('.', expand = True).rename(columns={0:'TAX_ID',1:'PROJET',2:'GENE'})],axis=1)
        familly_groups = familly_groups.loc[:,('GENE','FAMILLY_CODE')]
        annotations    = annotations.set_index('GENE').join(familly_groups.loc[:,('GENE','FAMILLY_CODE')].set_index("GENE"))

        del familly_groups
        annotations.loc[:,('GENE')] = annotations.index

        self.organisms_set = set(annotations.loc[:,('ORGANISM')])

        annotations = annotations.loc[:,('GENE','TYPE','ORGANISM','CONTIG_ID','FAMILLY_CODE','START','END','STRAND')]
        annotations.sort_values(['ORGANISM', 'CONTIG_ID', 'START'], inplace = True)
        annotations.fillna(value="None", inplace=True)

        Pangenome.annotations = list(map(list, annotations.as_matrix()))

        familly_code = None
        if not remove_singletons:
            familly_code = [record[GENE] if record[TYPE] == "CDS" and record[FAMILLY_CODE]=="None" else record[FAMILLY_CODE] for record in Pangenome.annotations]

        Pangenome.remove_singleton=remove_singletons

        self.nb_organisms = len(set([record[ORGANISM] for record in Pangenome.annotations]))

        cpt_org = 0
        cpt_fam = 0
        cur_org = Pangenome.annotations[0][ORGANISM]
        self.organism_positions[cur_org]=cpt_org
        self.annotation_positions[cur_org]="0-"
        for i in range(0, len(Pangenome.annotations)):
            #logging.getLogger().debug(cur_org)
            if familly_code != None:
                Pangenome.annotations[i][FAMILLY_CODE] = familly_code[i] if familly_code[i]!="None" else None
            else:
                Pangenome.annotations[i][FAMILLY_CODE] = Pangenome.annotations[i][FAMILLY_CODE] if Pangenome.annotations[i][FAMILLY_CODE]!="None" else None
            Pangenome.annotations[i] = tuple(Pangenome.annotations[i])
            if Pangenome.annotations[i][ORGANISM] != cur_org:
                self.annotation_positions[cur_org] = intspan(self.annotation_positions[cur_org]+str(i-1))
                cur_org = Pangenome.annotations[i][ORGANISM]
                self.annotation_positions[cur_org] = str(i)+"-"
                cpt_org += 1
                self.organism_positions[cur_org]=cpt_org 
            #logging.getLogger().debug(str(Pangenome.annotations[i][FAMILLY_CODE] is None))
            if Pangenome.annotations[i][FAMILLY_CODE] is None:#case of singleton or non CDS genes
                continue
            #logging.getLogger().debug(len(Pangenome.presences_absences[Pangenome.annotations[i][FAMILLY_CODE]]))
            if Pangenome.annotations[i][FAMILLY_CODE] not in Pangenome.presences_absences:
                Pangenome.presences_absences[Pangenome.annotations[i][FAMILLY_CODE]] = [0] * self.nb_organisms
                self.familly_positions[Pangenome.annotations[i][FAMILLY_CODE]]=cpt_fam
                cpt_fam += 1
                logging.getLogger().debug(str(cpt_fam) +" " + str(self.familly_positions[Pangenome.annotations[i][FAMILLY_CODE]]))
            #logging.getLogger().debug(Pangenome.annotations[i][FAMILLY_CODE])
            #logging.getLogger().debug(self.organism_positions[cur_org])
            Pangenome.presences_absences[Pangenome.annotations[i][FAMILLY_CODE]][self.organism_positions[cur_org]]=1
        
        Pangenome.annotations = tuple(Pangenome.annotations)
        self.annotation_positions[cur_org] = intspan(self.annotation_positions[cur_org]+str(i))
        Pangenome.presences_absences = OrderedDict(sorted(Pangenome.presences_absences.items(), key=lambda t: t[0]))#optional line of code 
        self.familly_positions       = OrderedDict(sorted(self.familly_positions.items(), key=lambda t: t[0]))

        logging.getLogger().debug(Pangenome.annotations)
        logging.getLogger().debug(self.annotation_positions)
        logging.getLogger().debug(self.organism_positions)
        logging.getLogger().debug(self.familly_positions) 
        logging.getLogger().debug(Pangenome.presences_absences)
 
    def __initialize_from_prokka_roary(self):
        self.presence_absence_matrix = pd.read_csv(self.options.roary_csv_file[0], sep=",")

        self.presence_absence_matrix.drop(["Non-unique Gene name","Annotation","No. isolates","No. sequences","Avg sequences per isolate","Genome Fragment","Order within Fragment","Accessory Fragment","Accessory Order with Fragment","QC","Min group size nuc","Max group size nuc","Avg group size nuc"], axis=1, inplace=True)
        melted_presence_absence_matrix = pd.melt(self.presence_absence_matrix, id_vars = "Gene").dropna()
        melted_presence_absence_matrix.set_index("value",inplace=True)
        melted_presence_absence_matrix = melted_presence_absence_matrix.loc[:,"FAMILLY_CODE"]
        melted_presence_absence_matrix = melted_presence_absence_matrix.as_dict()
        
        self.presence_absence_matrix.set_index("Gene",inplace=True)
        self.presence_absence_matrix = self.presence_absence_matrix.notnull().astype('int')

        print(melted_presence_absence_matrix)
        print(self.presence_absence_matrix)

        self.annotations = pd.DataFrame(self.__load_prokka_gff(melted_presence_absence_matrix), columns = ['GENE','TYPE','ORGANISM','CONTIG_ID','FAMILLY_CODE','START','END','STRAND'])
        self.annotations.sort_values(['ORGANISM', 'CONTIG_ID', 'STRAND', 'START'], inplace = True)
        print(self.annotations)
        exit()    

    def __initialize_from_prokka_mmseqs2(self):
        
        groups = dict()

        self.options.tsv_file[0]

        nb_groups=1
        first_iter = True
        for line in self.options.tsv_file[0]:
            elements = line.split()
            if first_iter:
                prec = elements[0]
                first_iter = False
            if elements[0] == prec:
                groups[elements[1]]=str(nb_groups)
            else :
                prec = elements[0]
                nb_groups +=1
                groups[elements[1]]=str(nb_groups)
        
        self.annotations = pd.DataFrame.from_dict(self.__load_prokka_gff(groups), orient = "index")
        self.annotations.columns = ['GENE','TYPE','ORGANISM','CONTIG_ID','FAMILLY_CODE','START','END','STRAND']
        self.annotations.sort_values(['ORGANISM', 'CONTIG_ID', 'START'], inplace = True)

        print(self.annotations)

        self.organisms = list(self.annotations.loc[:,"ORGANISM"].unique()).sort()
        self.families = list(self.annotations.loc[:,"FAMILLY_CODE"].unique()).sort()
        cluster_occurence_matrix = pd.DataFrame(0,columns = self.organisms, index = self.families, dtype = int)
        
        for index, row in self.annotations.loc[:,('FAMILLY_CODE','ORGANISM')].iterrows():
            if cluster_occurence_matrix.loc[row[0],row[1]] == 0:
                cluster_occurence_matrix.loc[row[0],row[1]]=1#+= 1
        
        self.presence_absence_matrix = cluster_occurence_matrix
        print(self.presence_absence_matrix)

    def __load_prokka_gff(self, groups):

        annot_dict = dict()    
        for gff_file in self.options.gff:
            filename = os.path.splitext(os.path.basename(gff_file.name))[0]
            index_date = filename[::-1].find("_")
            if index_date != -1:
                organism = filename[0:len(filename)-index_date-1]
            else:
                organism = filename
            for line in gff_file:
                if line.startswith("##"):
                    if line.startswith("##FASTA"):
                        break
                    else:
                        continue
                elements = line.split("\t")
                index_protein_start = elements[8].find('ID=')
                if index_protein_start != -1:
                    index_protein_end = elements[8][index_protein_start:].find(";")
                else:
                    print("gene name (ID feild) not found in gff file: "+gff_file.name+" \nline:"+line)
                    exit(1)

                protein = elements[8][index_protein_start+3:index_protein_end]
                if elements[2] == "CDS":
                    try:
                        familly_code = groups[protein]
                    except:
                        if self.options.remove_singleton:
                            continue
                        else:
                            familly_code = protein
                    
                    annot_dict[protein] = [protein,"CDS",organism,elements[0],familly_code,elements[3],elements[4],elements[6]]
            print(organism)

        print(annot_dict)
        return(annot_dict)

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
        if(len(weights) != len(organism_positions)):
            if(weights.keys() in organism_positions.keys()):
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
            index = {} 
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
                if Pangenome.neighbors is None: 
                    Pangenome.neighbors      = self.__neighborhoodComputation(1,neighbor_jumps)
                    Pangenome.neighbor_jumps = neighbor_jumps
                    logging.getLogger().info("Start optimal neighborhood graph construction (neighbor distance iteration) with max distance = "+str(neighbor_jumps)+" jumps")
                else:
                    logging.getLogger().info("Use already computed neighbors with max distance ="+str(neighbor_jumps)+" jumps")
                nei_file.write("1\n")
            else:
                nei_file.write("0\n")

            ck_file.write(str(k)+"\n")
        else:
            raise ValueError("result_path already exist")

        threshold = round(float(1)/self.nb_organisms,4) if self.nb_organisms >1 else 0
        
        graph = None
        if write_graph is not None:
            if use_neighborhood == False:
                raise ValueError("use_neighborhood must be True to write graph")
            graph  = nx.Graph()
       
        for fam in self.familly_positions.keys():
            logging.getLogger().debug(self.organism_positions.values())
            if self.weights is None:
                dat_file.write("\t".join(["1" if Pangenome.presences_absences[fam][p_org]>0 else "0" for p_org in self.organism_positions.values()])+"\n")
            else:
                dat_file.write("\t".join([str(round(self.weights[p_org],3)) if Pangenome.presences_absences[fam][p_org]>0 else "0.0" for p_org in self.organism_positions.values()])+"\n")
            ck_value = 0
            if self.familly_ratio[fam] == float(1):
                ck_value = 1
            elif self.familly_ratio[fam] <= threshold:
                ck_value = k
            ck_file.write(str(ck_value)+"\n")
            if use_neighborhood:
                if write_graph:                
                    graph.add_node(index[fam], label = fam, conservation = round(self.familly_ratio[fam],2))
                row_fam         = []
                row_dist_score  = []
                neighbor_number = 0
                logging.getLogger().debug(Pangenome.neighbors[fam])
                fam_neighbors = {nei:{org:dis for org,dis in orgs_nei.iteritems() if org in self.organism_positions.keys()} for nei, orgs_nei in Pangenome.neighbors[fam].iteritems() if nei != None}
                fam_neighbors = {nei: orgs_nei for nei, orgs_nei in fam_neighbors.iteritems() if len(orgs_nei)>0}
                if len(fam_neighbors.keys()) == 0:
                    logging.getLogger().warning("The familly: "+fam+" is an isolated familly")
                    nei_file.write(str(index[fam])+"\t0\n")
                    continue
                for neighbor, orgs_nei in Pangenome.neighbors[fam].iteritems():
                    logging.getLogger().debug(neighbor)
                    logging.getLogger().debug(orgs_nei)
                    if neighbor in index.keys():
                        distance_penality = float(sum(orgs_nei.values()))/len(orgs_nei)
                        distance_score = float(len(orgs_nei)) / distance_penality / self.nb_organisms
                        if distance_score>threshold:
                            if write_graph is not None:
                                graph.add_node(index[neighbor], label = neighbor, conservation = round(self.familly_ratio[fam],2))
                                graph.add_edge(index[fam],index[neighbor], weight = distance_score)
                            row_fam.append(str(index[neighbor]))
                            row_dist_score.append(str(round(distance_score,4)))
                            neighbor_number += 1
                    else:
                        logging.getLogger().debug("familly: "+neighbor+" is not in index file")
                nei_file.write("\t".join([str(item) for sublist in [[index[fam]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
                logging.getLogger().debug("\t".join([str(item) for sublist in [[index[fam]],[neighbor_number],row_fam,row_dist_score] for item in sublist])+"\n")
        nei_file.close()
        dat_file.close()
        ck_file.close()
        logging.getLogger().info("Running NEM...")

        #bernouli -> no weight or normal -> weight
        model = "bern" if self.weights is None else "norm"
        print_log = " -l y" if logging.getLogger().getEffectiveLevel() < 20 else "" 
        command = NEM_LOCATION+"nem_exe "+result_path+"/file "+str(k)+" -a nem -i 2000 -m "+model+" pk sk_ -s r 10 -B fix -b "+("1" if use_neighborhood else "0")+" -T -O random"+print_log
        logging.getLogger().info(command)
        proc = subprocess.Popen(shlex.split(command), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = proc.communicate()[1]
        M = float(re.search("\(M = (.+?)\)",output).group(1))# M=markov ps-like
        if model == "norm":
            self.BIC = 2 * M - (k * (self.nb_organisms + 1) + 1 + k- 1) * math.log(len(self.familly_positions)) 
        elif model == "bern":
            self.BIC = 2 * M - (k * self.nb_organisms + 1 + k - 1) * math.log(len(self.familly_positions))
        logging.getLogger().info("Based on "+str(k)+" classes, BIC = "+str(round(self.BIC,4)))
        if os.path.isfile(result_path+"/file.cf"):
            logging.getLogger().info("Reading NEM results")
        else:
            logging.getLogger().error("No NEM output found in file.cf")
        with open(result_path+"/file.cf","r") as classification_file:
            classification = classification_file.readline().split()
            classification = {k: v for k, v in zip(index.keys(), classification)}
            self.k = k
            if write_graph is not None:
                logging.getLogger().info("Writing graphML file")
                for fam, nem_class in classification.items():
                    graph.node[index[fam]]['nem_class'] = nem_class
                getattr(nx,'write_'+write_graph)(graph,result_path+"/file."+write_graph)
                mst = nx.maximum_spanning_tree(graph)
                for u,v,d in mst.edges(data=True):
                    d["MST"]=True
                getattr(nx,'write_'+write_graph)(nx.disjoint_union(graph,mst),result_path+"/file_mst."+write_graph)
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


    def __neighborhoodComputation(self, initNeighborDistance, maxNeighborDistance):
        """Algo voisinage taille flexible
            Principe : Tant qu'il reste des noeuds (=MICFAM) isoles dans le graphe (mais potentiellement "voisinable"), la distance de voisinage considere augmente jusqu'au max
            En sortie : Une matrice d'adjacence scoree avec les Oids des paires de voisins (le nombre d'Oid pour 1 paire sera utilise par la suite) 
            + la distance de voisinage ou la paire de voisin a ete trouvee, pour l'Oid considere (la moyenne de ces distances pour 1 paire sera utilise par la suite)
            En finalite : Produit un graphe connexe du voisinange des MICFAM ou chaque arete a pour score egal a : (Nb_Oid_voisins / Moyenne_Distances_voisins / Nb_Oid_total) """
        neighbors_dict = defaultdict(lambda : defaultdict(dict))
        dist = initNeighborDistance # Distance (en nb de regions, cad de lignes) de recherche des voisins
        neighbor_rows = [] # Garde en temp les n lignes precedentes avec n = dist, remis a zero chaque fois qu'on change de brin et de sequence
        i = 0 # = compteur remis a zero chaque fois qu'on change de brin et de sequence
        list_isolated_famillies = [] # = clusters isoles mais possible de leur trouver un voisinage : tant qu'il en reste => recursion (jusqu'a maxNeighborDistance)
        list_allowed_isolated_famillies = [] # = clusters vraiment isoles, impossible de leur trouver un voisin (ils sont autorises a passer la recursion)
        nb_allowed_isolated_famillies = 0
        first_row = self.annotations[0]
        GO_id_prev, GO_type_prev, O_id_prev, S_id_prev, familly_id_prev, GO_begin_prev, GO_end_prev, GO_strand_prev = first_row
        valid_familly_id_prev = familly_id_prev    
        neighbor_rows.append(first_row)    
        last_index = len(self.annotations)-2
        
        for index, row in enumerate(self.annotations[1:]):
     
            GO_id, GO_type, O_id, S_id, familly_id, GO_begin, GO_end, GO_strand = row

            if (S_id != S_id_prev):
                if (i == 0):# Cas particulier ou 1 seul CDS dans 1 sequence + 1 strand : Forcement aucun voisins #
                       list_allowed_isolated_famillies.append(str(familly_id_prev))
                elif ((dist-j) >= i):   # Cas ou vraiment aucun voisins qqsoit la distance consideree #
                    list_allowed_isolated_famillies.append(str(valid_familly_id_prev))
                elif index == last_index: # Cas ou il ne peux pas y avoir de suivant car on a atteint la fin#
                    list_allowed_isolated_famillies.append(str(familly_id))
                i=0
                neighbor_rows = []
            else:
                i+=1
                tmp_list_isolated_famillies = []
                for j in range((i if i<dist else dist),0,-1):
                    GO_id_nei, GO_type_nei, O_id_nei, S_id_nei, familly_id_nei, GO_begin_nei, GO_end_nei, GO_strand_nei = neighbor_rows[j-1]
                    if (familly_id_nei is not None) and GO_type_nei in ("CDS","fCDS"):
                        if ((familly_id is not None) and GO_type in ("CDS","fCDS") and familly_id != familly_id_nei):
                            neighbors_dict[str(familly_id)][str(familly_id_nei)][O_id]=dist-j+1
                            neighbors_dict[str(familly_id_nei)][str(familly_id)][O_id]=dist-j+1
                            tmp_list_isolated_famillies = []
                            break
                        else:
                            tmp_list_isolated_famillies.append(str(familly_id_nei))
                    else:
                        if ((familly_id is not None) and GO_type in ("CDS","fCDS") and familly_id != familly_id_nei):
                            tmp_list_isolated_famillies.append(str(familly_id))
                list_isolated_famillies.extend(tmp_list_isolated_famillies)

            if i<dist:
                neighbor_rows.append(row)
            else:
                neighbor_rows.pop(0)
                neighbor_rows.append(row)
            GO_id_prev, GO_type_prev, O_id_prev, S_id_prev, familly_id_prev, GO_begin_prev, GO_end_prev, GO_strand_prev = row
            valid_familly_id_prev = familly_id if ((familly_id is not None) and GO_type in ("CDS","fCDS")) else valid_familly_id_prev

        # Cas des clusters reelement isoles dans le graphe sans aucun voisinage, meme a rayon tres eleve
        for clust in (set(list_allowed_isolated_famillies) - set(neighbors_dict.keys())):
            neighbors_dict[str(clust)][None] = 1
            nb_allowed_isolated_famillies+=1
        
        logging.getLogger().debug("Neighbor maximum distance: "+str(dist)+" -> Isolated cluster counting: "+str(len(set(list_isolated_famillies) - set(list(neighbors_dict.keys())))+nb_allowed_isolated_famillies))
        if len(set(list_isolated_famillies) - set(neighbors_dict.keys())) != 0:
            if dist < maxNeighborDistance:
                return self.neighborhoodComputation(dist+1, maxNeighborDistance)
            else:
                # Cas des clusters isoles avec potentiellement un voisinage a rayon superieur a maxNeighborDistance
                for clust in (set(list_isolated_famillies) - set(neighbors_dict.keys())):    
                    neighbors_dict[str(clust)][None] = 1
          
        return neighbors_dict
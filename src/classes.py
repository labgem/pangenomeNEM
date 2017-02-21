#!/usr/bin/python
# -*- coding: iso-8859-1 -*-

#import mysql.connector as mysql
import pandas as pd
import numpy as np
from collections import defaultdict
from collections import OrderedDict
from intspan import intspan
import networkx as nx
import subprocess
import os
#from sklearn import manifold
#import urllib2
import logging

NEM_LOCATION  = "../NEM/"
(GENE, TYPE, ORGANISM, CONTIG_ID, GROUP_CODE, START, END, STRAND) = range(0, 8)

class Pangenome:

    presences_absences = OrderedDict()
    annotations        = tuple()
    remove_singleton   = False
    neighbors_dict     = dict()

    @staticmethod
    def clear_Pangenome():
        Pangenome.presences_absences = OrderedDict()
        Pangenome.annotations        = tuple()
        Pangenome.remove_singleton   = False
        neighbors_dict               = dict()
    
    def __init__(self, init_from = "args", *args):
        self.nb_organisms           = 0
        self.organism_positions     = OrderedDict()
        self.familly_positions      = OrderedDict()
        self.annotation_positions   = OrderedDict()
        self.familly_ratio          = dict()
        self.core_list              = list()
        self.pan_size               = 0
        self.core_size              = 0
        self.cluster_classification = None
        self.classnumber            = None
        self.neighbor_jumps         = None

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

            self.annotations.columns = ['GENE','TYPE','ORGANISM','CONTIG_ID','GROUP_CODE','START','END','STRAND']
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

        annotations    = pd.read_csv(annotations_file, sep="\t", usecols=useful_cols_annot, dtype=type_cols_annot)
        ortholog_group = pd.read_csv(eggNOG_clusters_file, sep="\t", usecols=useful_cols_orthoGroups, dtype=type_cols_orthoGroups)

        ortholog_group.rename(columns={'PROTEIN_ID':'PROTEIN_ID','EGGNOG_GROUP_CODE':'GROUP_CODE'}, inplace=True)

        #split a column with 3 fields separated by a point in 3 other columns
        #expand = True make pandas to return a DataFrame with the 3 columns comming from the splitted GENE_ID column
        annotations = pd.concat([annotations,annotations.loc[:,('GENE_ID')].str.split('.', expand = True).rename(columns={0:'TAX_ID',1:'PROJET',2:'GENE'})],axis=1)    
        annotations = pd.concat([annotations,annotations.loc[:,('CONTIG_ID')].str.split('.', expand = True).loc[:,2].to_frame("CONTIG")],axis=1)
        annotations.loc[:,('ORGANISM')] = annotations.loc[:,('TAX_ID')] + "." + annotations.loc[:,('PROJET')]

        annotations    = annotations.loc[:,('GENE','TYPE','ORGANISM','CONTIG_ID','START','END','STRAND')]
        
        ortholog_group = pd.concat([ortholog_group,ortholog_group.loc[:,('PROTEIN_ID')].str.split('.', expand = True).rename(columns={0:'TAX_ID',1:'PROJET',2:'GENE'})],axis=1)
        ortholog_group = ortholog_group.loc[:,('GENE','GROUP_CODE')]
        annotations    = annotations.set_index('GENE').join(ortholog_group.loc[:,('GENE','GROUP_CODE')].set_index("GENE"))

        del ortholog_group
        annotations.loc[:,('GENE')] = annotations.index

        self.organisms_set = set(annotations.loc[:,('ORGANISM')])

        annotations = annotations.loc[:,('GENE','TYPE','ORGANISM','CONTIG_ID','GROUP_CODE','START','END','STRAND')]
        annotations.sort_values(['ORGANISM', 'CONTIG_ID', 'START'], inplace = True)
        annotations.fillna(value="None", inplace=True)

        Pangenome.annotations = list(map(list, annotations.as_matrix()))

        group_code = None
        if not remove_singletons:
            group_code = [record[GENE] if record[TYPE] == "CDS" and record[GROUP_CODE]=="None" else record[GROUP_CODE] for record in Pangenome.annotations]

        Pangenome.remove_singleton=remove_singletons

        self.nb_organisms = len(set([record[ORGANISM] for record in Pangenome.annotations]))

        cpt_org = 0
        cpt_fam = 0
        cur_org = Pangenome.annotations[0][ORGANISM]
        self.organism_positions[cur_org]=cpt_org
        self.annotation_positions[cur_org]="0-"
        for i in range(0, len(Pangenome.annotations)):
            #logging.getLogger().debug(cur_org)
            if group_code != None:
                Pangenome.annotations[i][GROUP_CODE] = group_code[i] if group_code[i]!="None" else None
            else:
                Pangenome.annotations[i][GROUP_CODE] = Pangenome.annotations[i][GROUP_CODE] if Pangenome.annotations[i][GROUP_CODE]!="None" else None
            Pangenome.annotations[i] = tuple(Pangenome.annotations[i])
            if Pangenome.annotations[i][ORGANISM] != cur_org:
                self.annotation_positions[cur_org] = intspan(self.annotation_positions[cur_org]+str(i-1))
                cur_org = Pangenome.annotations[i][ORGANISM]
                self.annotation_positions[cur_org] = str(i)+"-"
                cpt_org += 1
                self.organism_positions[cur_org]=cpt_org 
            #logging.getLogger().debug(str(Pangenome.annotations[i][GROUP_CODE] is None))
            if Pangenome.annotations[i][GROUP_CODE] is None:#case of singleton or non CDS genes
                continue
            #logging.getLogger().debug(len(Pangenome.presences_absences[Pangenome.annotations[i][GROUP_CODE]]))
            if Pangenome.annotations[i][GROUP_CODE] not in Pangenome.presences_absences:
                Pangenome.presences_absences[Pangenome.annotations[i][GROUP_CODE]] = [0] * self.nb_organisms
                self.familly_positions[Pangenome.annotations[i][GROUP_CODE]]=cpt_fam
                cpt_fam += 1
                logging.getLogger().debug(str(cpt_fam) +" " + str(self.familly_positions[Pangenome.annotations[i][GROUP_CODE]]))
            #logging.getLogger().debug(Pangenome.annotations[i][GROUP_CODE])
            #logging.getLogger().debug(self.organism_positions[cur_org])
            Pangenome.presences_absences[Pangenome.annotations[i][GROUP_CODE]][self.organism_positions[cur_org]]=1
        
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
        melted_presence_absence_matrix = melted_presence_absence_matrix.loc[:,"GROUP_CODE"]
        melted_presence_absence_matrix = melted_presence_absence_matrix.as_dict()
        
        self.presence_absence_matrix.set_index("Gene",inplace=True)
        self.presence_absence_matrix = self.presence_absence_matrix.notnull().astype('int')

        print(melted_presence_absence_matrix)
        print(self.presence_absence_matrix)

        self.annotations = pd.DataFrame(self.__load_prokka_gff(melted_presence_absence_matrix), columns = ['GENE','TYPE','ORGANISM','CONTIG_ID','GROUP_CODE','START','END','STRAND'])
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
        self.annotations.columns = ['GENE','TYPE','ORGANISM','CONTIG_ID','GROUP_CODE','START','END','STRAND']
        self.annotations.sort_values(['ORGANISM', 'CONTIG_ID', 'START'], inplace = True)

        print(self.annotations)

        self.organisms = list(self.annotations.loc[:,"ORGANISM"].unique()).sort()
        self.families = list(self.annotations.loc[:,"GROUP_CODE"].unique()).sort()
        cluster_occurence_matrix = pd.DataFrame(0,columns = self.organisms, index = self.families, dtype = int)
        
        for index, row in self.annotations.loc[:,('GROUP_CODE','ORGANISM')].iterrows():
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
                        group_code = groups[protein]
                    except:
                        if self.options.remove_singleton:
                            continue
                        else:
                            group_code = protein
                    
                    annot_dict[protein] = [protein,"CDS",organism,elements[0],group_code,elements[3],elements[4],elements[6]]
            print(organism)

        print(annot_dict)
        return(annot_dict)

    def __str__(self):

        pan_str    = ""
        pan_str += "----------- Statistics -----------\n"
        pan_str += "Number of organisms: "+str(len(self.nb_organisms))+"\n"
        pan_str += "Pan-genome size:"+str(self.pan_size)+"\n"
        pan_str += "Exact core-genome size:"+str(self.core_size)+"\n"
        pan_str += "Exact variable-genome size:"+str(self.pan_size-self.core_size)+"\n"
        pan_str += "----------------------------------\n"

        return(pan_str)    

    def sub_pangenome(self, sub_organisms):
       
        if set(sub_organisms).issubset(set(self.organism_positions.keys())):
            
            sub_organism_positions   = OrderedDict((org, self.organism_positions[org]) for org in sorted(sub_organisms))
            logging.getLogger().debug(sub_organism_positions) 
            sub_annotation_positions = OrderedDict((org, self.annotation_positions[org]) for org in sorted(sub_organisms))
            logging.getLogger().debug(sub_annotation_positions) 
            
            subset_group_code = set()
            for org in sub_organisms:
                for pos in sub_annotation_positions[org]:
                    subset_group_code.add(Pangenome.annotations[pos][GROUP_CODE])
            logging.getLogger().debug(subset_group_code)
            subset_group_code = subset_group_code - set([None])
            logging.getLogger().debug(subset_group_code)

            sub_familly_positions = OrderedDict((fam, self.familly_positions[fam]) for fam in sorted(subset_group_code))

            #logging.getLogger().debug(len(sub_organisms))

            sub_pangenome = Pangenome("args", len(sub_organisms), sub_organism_positions,sub_familly_positions,sub_annotation_positions)

            return(sub_pangenome)
        else:
            raise ValueError("Organisms contained in the sub_organisms arguments are not a subset of the organisms contained in this object") 

    def classify(self, inc, path_prefix, max_neighbordistance = 1, plot = False):
                
        outputdir = path_prefix+"nborg"+str(len(self.organisms))+"_k"+str(self.options.classnumber[0])+"_i"+str(inc)+"/"

        #find most contigous organism
        reference = self.annotations.loc[lambda annot: annot.ORGANISM == 'EMPTY',:].loc[:,("GROUP_CODE")]
    
        self.nem_location = outputdir

        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
            ## Indexation des MICFAM
            MICFAM_index = {}
            i=0
            index_file = open(outputdir+"/file.index", "w")
            for ortholog in self.presence_absence_matrix.index:
                    i+=1
                    MICFAM_index[ortholog]=i
                    index_file.write(str(i)+"\t"+str(ortholog)+"\n")

            index_file.close()
            if self.options.verbose:
                    print("> Creation of file.index file")

            neighbors_dict = {}
            if self.options.neighborcomputation:
                    ## Recuperation des MICFAM voisins :
                neighbors_dict = self.__neighborhoodComputation(1,max_neighbordistance)
            if self.options.verbose:
                print("---- Start optimal neighborhood graph construction (neighbor distance iteration) with max distance = %s regions ----" % str(max_neighbordistance))

            # Ecriture des fichiers .dat, .nei et .str
            str_file = open(outputdir+"/file.str", "w")
            str_file.write("S\t"+str(self.pan_genome_size)+"\t"+str(len(self.organisms))+"\n")
            str_file.close()
            nei_file = open(outputdir+"/file.nei", "w")
            if self.options.neighborcomputation:
                    nei_file.write("1\n")
            else:
                    nei_file.write("0\n")
            dat_file = open(outputdir+"/file.dat", "w")
        ck_file = open(outputdir+"/file.ck", "w")

        ck_file.write(str(self.options.classnumber[0])+"\n")

        threshold = round(float(1)/len(self.organisms),4) if len(self.organisms)>1 else 0
    
        if inc == 0:
            graph  = nx.Graph()
        
        presence_absence_matrix = self.presence_absence_matrix

        if self.weights is not None :
            presence_absence_matrix = self.presence_absence_matrix.multiply(self.weights, axis=1)

        for ortholog, row in presence_absence_matrix.iterrows():
    
            dat_file.write("\t".join([str(e) for e in row])+"\n")
    
            conservation = round(float(self.cluster_pan.loc[ortholog]),4)

            ck_value = 0
            if conservation == 1:
                ck_value = 1
            elif conservation <= threshold:
                ck_value = self.options.classnumber[0]
            ck_file.write(str(ck_value)+"\n")
            
            if inc == 0:                
                try:
                    order = pd.Index(reference).get_loc(ortholog)
                    if type(order) != int:
                        raise KeyError 
                    graph.add_node(MICFAM_index[ortholog],label = ortholog, conservation = conservation, order = order)
                except KeyError:
                    graph.add_node(MICFAM_index[ortholog],label = ortholog, conservation = conservation)
            if self.options.neighborcomputation:
                row = []        
                row_tmp = []
                row_tmp2 = []
                neighbor_number = 0
                row.append(str(MICFAM_index[ortholog]))
                if neighbors_dict[ortholog].keys() == [None]:
                    if self.options.verbose:
                            print("Warning: The MICFAM "+ortholog+" is an isolated cluster node")
                    continue
                for neighbor, Oids_nei in neighbors_dict[ortholog].iteritems():
                    if neighbor in MICFAM_index.keys():
                        distance_penality = float(sum(Oids_nei.values()))/len(Oids_nei.values())
                        distance_score = float(len(Oids_nei)) / distance_penality / len(self.organisms)
                        if distance_score>threshold:
                            conservation = round(float(self.cluster_pan.loc[neighbor]),2)
                            if inc == 0:
                                try:
                                    order = pd.Index(reference).get_loc(neighbor)
                                    if type(order) != int:
                                        raise KeyError
                                    graph.add_node(MICFAM_index[neighbor],label = neighbor, conservation = conservation, order = order)
                                except KeyError:
                                    graph.add_node(MICFAM_index[neighbor],label = neighbor, conservation = conservation)
                                graph.add_edge(MICFAM_index[ortholog],MICFAM_index[neighbor], weight = distance_score)

                            row_tmp.append(str(MICFAM_index[neighbor]))
                            row_tmp2.append(str(distance_score))
                            neighbor_number+=1
                    else:
                        print(neighbor+" <= Warning : Not in index file")
                row.append(str(neighbor_number))
                row.extend(row_tmp)
                row.extend(row_tmp2)
                nei_file.write("\t".join(row)+"\n")
        else:
            nei_file.write(str(MICFAM_index[ortholog])+"\t0\n")
        nei_file.close()
        dat_file.close()
        ck_file.close()
        if self.options.verbose:
            print("> Creation of file.nei, file.dat and file.str files")

        #bernouli -> no weight or normal -> weight
        model = "bern" if self.weights is None else "norm"
        print_log = " -l y" if self.options.verbose else ""    
        command = NEM_LOCATION+"nem_exe "+outputdir+"/file "+str(self.options.classnumber[0])+" -a nem -i 2000 -m "+model+" pk sk_ -s r 10 -B fix -b "+("1" if self.options.neighborcomputation else "0")+" -T -O random"+print_log
        print(command)
        
        proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.communicate()    
    
        with open(outputdir+"/file.cf","r") as classification_file:
            #TODO store classification in attribute
            classification = classification_file.readline().split()
            self.cluster_classification = pd.Series(index = self.presence_absence_matrix.index, data = classification, name = "classification")
            if inc == 0:
                for i, nem_class in enumerate(classification):
                    graph.node[int(i)+1]['nem_class'] = nem_class

                nx.write_graphml(graph,outputdir+"/file.graphml")
                mst = nx.maximum_spanning_tree(graph)
                for u,v,d in mst.edges(data=True):
                    d["MST"]=True
                nx.write_graphml(nx.disjoint_union(graph,mst),outputdir+"/file_mst.graphml")

        if plot:
            self.plot(outputdir+"/occurences.pdf")

        with open(outputdir+"/nem.stat","w") as nem_stat_file:
            conservation = dict()
            for k in range(1,self.options.classnumber[0]+1):
                clusters_classified_in_k = self.cluster_classification[self.cluster_classification == str(k)]
                conservation[str(k)]=self.cluster_pan.loc[clusters_classified_in_k.index].mean()
                conservation[str(k)] = 0 if np.isnan(conservation[str(k)]) else conservation[str(k)]                
            order = sorted(conservation, key=conservation.__getitem__, reverse=True)    
            stats = self.cluster_classification.value_counts().to_dict()
            for item in order:
                if item not in stats:
                    stats[item]=0
            nem_stat_file.write(str(len(self.organisms))+"\t"+"\t".join([str(stats[item]) for item in order])+"\t"+str(self.pan_genome_size)+"\n")
        with open(outputdir+"/exact.stat","w") as exact_stat_file:        
            exact_stat_file.write(str(len(self.organisms))+"\t"+str(self.core_genome_size)+"\t"+str(self.pan_genome_size-self.core_genome_size)+"\t"+str(self.pan_genome_size)+"\n")    
    def __neighborhoodComputation(self, initNeighborDistance, maxNeighborDistance):
        """Algo voisinage taille flexible
                Principe : Tant qu'il reste des noeuds (=MICFAM) isoles dans le graphe (mais potentiellement "voisinable"), la distance de voisinage considere augmente jusqu'au max
                En sortie : Une matrice d'adjacence scoree avec les Oids des paires de voisins (le nombre d'Oid pour 1 paire sera utilise par la suite) 
                + la distance de voisinage ou la paire de voisin a ete trouvee, pour l'Oid considere (la moyenne de ces distances pour 1 paire sera utilise par la suite)
                En finalite : Produit un graphe connexe du voisinange des MICFAM ou chaque arete a pour score egal a : (Nb_Oid_voisins / Moyenne_Distances_voisins / Nb_Oid_total) """
        print("here")
        neighbors_dict = defaultdict(lambda : defaultdict(dict))
        dist = initNeighborDistance # Distance (en nb de regions, cad de lignes) de recherche des voisins
        neighbor_rows = [] # Garde en temp les n lignes precedentes avec n = dist, remis a zero chaque fois qu'on change de brin et de sequence
        i = 0 # = compteur de lignes, remis a zero chaque fois qu'on change de brin et de sequence
        list_isolated_cluster = [] # = clusters isoles mais possible de leur trouver un voisinage : tant qu'il en reste => recursion (jusqu'a maxNeighborDistance)
        list_allowed_isolated_cluster = [] # = clusters vraiment isoles, impossible de leur trouver un voisin (ils sont autorises a passer la recursion)
        nb_allowed_isolated_cluster = 0
        first_row = tuple(self.annotations.iloc[0].values)
        
        GO_id_prev, GO_type_prev, O_id_prev, S_id_prev, cluster_id_prev, GO_begin_prev, GO_end_prev, GO_strand_prev = first_row
        valid_cluster_id_prev = cluster_id_prev    
        neighbor_rows.append(first_row)    

        last_index = self.annotations.tail(1).index
        

        for key,row in self.annotations.iloc[1:].iterrows():
     
            GO_id, GO_type, O_id, S_id, cluster_id, GO_begin, GO_end, GO_strand = row

            if (S_id != S_id_prev):
                if (i == 0):# Cas particulier ou 1 seul CDS dans 1 sequence + 1 strand : Forcement aucun voisins #
                       list_allowed_isolated_cluster.append(str(cluster_id_prev))
                elif ((dist-j) >= i):   # Cas ou vraiment aucun voisins qqsoit la distance consideree #
                    list_allowed_isolated_cluster.append(str(valid_cluster_id_prev))
                elif key == last_index: # Cas ou il ne peux pas y avoir de suivant car on a atteint la fin#
                    list_allowed_isolated_cluster.append(str(cluster_id))
                i=0
                neighbor_rows = []
            else:
                i+=1
                tmp_list_isolated_cluster = []
                for j in range((i if i<dist else dist),0,-1):
                    GO_id_nei, GO_type_nei, O_id_nei, S_id_nei, cluster_id_nei, GO_begin_nei, GO_end_nei, GO_strand_nei = neighbor_rows[j-1]
                    if ((not pd.isnull(cluster_id_nei)) and GO_type_nei in ("CDS","fCDS")):
                        if ((not pd.isnull(cluster_id)) and GO_type in ("CDS","fCDS") and cluster_id != cluster_id_nei):
                            neighbors_dict[str(cluster_id)][str(cluster_id_nei)][O_id]=dist-j+1
                            neighbors_dict[str(cluster_id_nei)][str(cluster_id)][O_id]=dist-j+1
                            tmp_list_isolated_cluster = []
                            break
                        else:
                            tmp_list_isolated_cluster.append(str(cluster_id_nei))
                    else:
                        if ((not pd.isnull(cluster_id)) and GO_type in ("CDS","fCDS") and cluster_id != cluster_id_nei):
                            tmp_list_isolated_cluster.append(str(cluster_id))
                list_isolated_cluster.extend(tmp_list_isolated_cluster)

            if i<dist:
                neighbor_rows.append(row)
            else:
                neighbor_rows.pop(0)
                neighbor_rows.append(row)
            GO_id_prev, GO_type_prev, O_id_prev, S_id_prev, cluster_id_prev, GO_begin_prev, GO_end_prev, GO_strand_prev = row
            valid_cluster_id_prev = cluster_id if ((not pd.isnull(cluster_id)) and GO_type in ("CDS","fCDS")) else valid_cluster_id_prev

        # Cas des clusters reelement isoles dans le graphe sans aucun voisinage, meme a rayon tres eleve
        for clust in (set(list_allowed_isolated_cluster) - set(neighbors_dict.keys())):
            print(clust)
            neighbors_dict[str(clust)][None] = 1
            nb_allowed_isolated_cluster+=1
        
        if self.options.verbose:
            print("Neighbor maximum distance: "+str(dist)+" -> Isolated cluster counting: "+str(len(set(list_isolated_cluster) - set(list(neighbors_dict.keys())))+nb_allowed_isolated_cluster))
        if len(set(list_isolated_cluster) - set(neighbors_dict.keys())) != 0:
            if dist < maxNeighborDistance:
                return self.neighborhoodComputation(dist+1, maxNeighborDistance)
            else:
                # Cas des clusters isoles avec potentiellement un voisinage a rayon superieur a maxNeighborDistance
                for clust in (set(list_isolated_cluster) - set(neighbors_dict.keys())):    
                    neighbors_dict[str(clust)][None] = 1
          
        return neighbors_dict

def __dbConnect():
    try:
        conn = mysql.connect(host="mysqlagcdb.genoscope.cns.fr",port=3306,user="agc",passwd="admagc21",db="pkgdb_dev")
    except:
        print("Connexion with pkgdb_dev failed.")
        exit(1)
    return conn

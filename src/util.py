#!/usr/bin/python
# -*- coding: iso-8859-1 -*-


from optparse import OptionParser # deprecie mais + compatible avec les install python des etnas
#import pandas as pd
import os
import sys
import errno
from optparse import OptionParser # deprecie mais + compatible avec les install python des etnas
from collections import defaultdict
import collections
import numpy as np
import warnings
import math
from multiprocessing import Pool
import itertools
import shutil
import random
import string
from decimal import Decimal
import scipy.misc
import networkx as nx
from math import cos, sin, pi, sqrt, atan2
from Bio import SeqIO
import subprocess
import pandas as pd
import logging
warnings.filterwarnings("ignore", "Unknown table.*")


####################################################################################
#                                                                       Fonctions                                                                          #
####################################################################################

# Calcul du nombre de combinaisons de k elements parmi n
def combinationNb(k,n):
        if (k == 0):
                return 1
        result = 1
        for i in range(0, k):
                result *= Decimal(n - i)/(i + 1);
        return int(round(result))

# Calcul du nombre total de combinaisons uniques de n elements
def combinationTotalNb(size):
        return pow(2,size)-1

# Generation d'une sous-liste al<C3><A9>atoire de taille n
def randomSublist(items,n):
        item_array = np.array(items)
        index_array = np.arange(item_array.size)
        np.random.shuffle(index_array)
        ordered_index_array = sorted(index_array[:n])
        return list(item_array[ordered_index_array])

# Generation de toutes les combinaisons uniques (sans notion d'ordre) d'elements donnes
def exactCombinations(items):
	
        len_item  = len(items);
        combinations = defaultdict(list)
        for i in range(1, 1<<len_item):
                c = []
                for j in range(0, len_item):
                        if(i & (1 << j)):
                                c.append(items[j])
                combinations[len(c)].append(c)
        return combinations

# Echantillonage proportionnel d'un nombre donne de combinaisons (sans remise)
def samplingCombinations(items, sample_thr, sample_min, sample_max=None):
        samplingCombinationList = defaultdict(list)
        item_size = len(items)
        combTotNb = combinationTotalNb(item_size)
        sample_coeff = (Decimal(combTotNb)/sample_thr)
        for k in range(1,item_size+1):
                tmp_comb = []
                combNb = Decimal(scipy.misc.comb(item_size,k))#combinationNb(k,item_size)
                combNb = sys.float_info.max if combNb>sys.float_info.max else combNb# to avoid to reach infinit values
                combNb_sample = math.ceil(Decimal(combNb)/sample_coeff)
                # Plus petit echantillonage possible pour un k donn<C3><A9> = sample_min
                if ((combNb_sample < sample_min) and k != item_size):
                        combNb_sample = sample_min
                # Plus grand echantillonage possible
                if (sample_max is not None and combNb_sample > sample_max):
                        combNb_sample = sample_max
                i = 0;
                while (i < combNb_sample):
                        comb = randomSublist(items,k)
                        # Echantillonnage sans remise
                        if (comb not in tmp_comb):
                                tmp_comb.append(comb)
                                samplingCombinationList[len(comb)].append(comb)
                                i+=1
        return samplingCombinationList

# Generation des combinaisons d'une liste d'Oids (toutes les combinaisons ou bien un echantillon selon les seuils fixes)
def oidsCombinations(Oids, nbOrgThr, sample_thr, sample_min,sample_max=None):
        if (len(Oids) <= nbOrgThr):
                        comb_list = exactCombinations(Oids)
        else:
                comb_list = samplingCombinations(Oids, sample_thr, sample_min, sample_max)
        return comb_list

def run(cpt, pan, k, organisms):
        subpan = pan.sub_pangenome(organisms)
        subpan.classify("/tmp/test_pangenome2"+"_k3"+"_nb"+str(len(organisms))+"_i"+str(cpt),k=3, use_neighborhood=True, write_graph = "gexf")
        print(subpan)

def circle_intersection(circle1, circle2):
        '''
        @summary: calculates intersection points of two circles
        @param circle1: tuple(x,y,radius)
        @param circle2: tuple(x,y,radius)
        @result: string if no intersection and mean point (which is (x,y) tuple) between the center of the 2 circles
        '''
        x1,y1,r1 = circle1
        x2,y2,r2 = circle2
        # http://stackoverflow.com/a/3349134/798588
        dx,dy = x2-x1,y2-y1
        d = sqrt(dx*dx+dy*dy)
        if d > r1+r2:
            return "separated"
        if d < abs(r1-r2):
            return "contained"
        if d == 0 and r1 == r2:
            return (x1,y1)

        a = (r1*r1-r2*r2+d*d)/(2*d)
        #h = sqrt(r1*r1-a*a)
        xm = x1 + a*dx/d
        ym = y1 + a*dy/d
        #xs1 = xm + h*dy/d
        #xs2 = xm - h*dy/d
        #ys1 = ym - h*dx/d
        #ys2 = ym + h*dx/d

        return (xm,ym)

def included_in_circle(point, circle1):
        res = circle_intersection(point, circle1)
        if res == "contained" or res == point:
                return True
        else:
                return False

def calc_mash_distance(fasta, OUTPUTDIR, num_thread):

        MASH_DIRECTORY = OUTPUTDIR+"/"+"mash/"
        if not os.path.exists(MASH_DIRECTORY):
                os.makedirs(MASH_DIRECTORY)
        
        mash_parameters = set()
        fasta_sequences = SeqIO.parse(fasta,'fasta')
        for fasta in fasta_sequences:
                elements = fasta.id.split(".")
                out_file_name = MASH_DIRECTORY+elements[0]+"."+elements[1]
                mash_parameters.add(out_file_name)

                with open(out_file_name+".fasta","a") as out_file:
                        SeqIO.write(fasta,out_file,"fasta")

        command = 'mash sketch -n -p ' + str(num_thread) + " -o "+MASH_DIRECTORY+"all_sketch.msh "+MASH_DIRECTORY+"*.fasta"
        print(command)
        proc = subprocess.Popen(command, shell=True)
        proc.communicate()

        command = 'mash dist -t '+(MASH_DIRECTORY+'all_sketch.msh ') * 2 + " > "+OUTPUTDIR+"/mash_distance.csv"
        print(command)
        proc = subprocess.Popen(command, shell=True,stdout=subprocess.PIPE)
        proc.communicate()
        


        distances = pd.read_csv(OUTPUTDIR+"/mash_distance.csv", sep="\t", index_col =0)
        organisms_names = {i : os.path.splitext(os.path.basename(i))[0] for i in distances.index}
        distances.rename(index=organisms_names,columns=organisms_names, inplace=True)

        logging.getLogger().debug(distances.values)
        logging.getLogger().debug(distances)
        return(distances)

def findCOG(pangenome):

        print("HERE")   
        headers = {
            'Origin': 'http://eggnogdb.embl.de',
            'Accept-Encoding': 'gzip, deflate',
            'Accept-Language': 'en-US,en;q=0.8,fr;q=0.6',
            'User-Agent': 'Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/53.0.2785.143 Safari/537.36',
            'Content-Type': 'application/json;charset=UTF-8',
            'Accept': 'application/json, text/plain, */*',
            'Referer': 'http://eggnogdb.embl.de/',
            'Connection': 'keep-alive',
            'DNT': '1',
        }
                
        ortho_2_COG_funcat = pd.DataFrame("S", index = pangenome.families, columns = ["funcat"])
        for ortho in sorted(pangenome.families):
                if len(ortho) == 5:
                        try:
                                data = '{"desc":"","seqid":"","target_species":"","level":"","nognames":"'+ortho+'","page":0}'
                                content = requests.post('http://eggnogapi.embl.de/meta_search', headers=headers, data=data)
                        except :
                                continue
                        if content.status_code == 200:
                                funcat = json.loads(content.text)["matches"][0]["funcat"]
                                ortho_2_COG_funcat.loc[ortho,"funcat"] = funcat
                        print(ortho+"  "+str(ortho_2_COG_funcat.loc[ortho,:]))
        return(ortho_2_COG_funcat)

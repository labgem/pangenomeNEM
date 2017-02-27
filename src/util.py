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
		
		combNb = sys.float_info.max if combNb>sys.float_info.max else combNb # to avoid to reach infinit values
			
                combNb_sample = math.ceil(Decimal(combNb)/sample_coeff)
                # Plus petit echantillonage possible pour un k donn<C3><A9> = sample_min
                if ((combNb_sample < sample_min) and k != item_size):
                        combNb_sample = sample_min
		# Plus grand echantillonage possible
		if (sample_max != None and (combNb_sample > sample_max)):
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
#!/usr/bin/env python3
# -*- coding: iso-8859-1 -*-

import sys
import gzip
from decimal import Decimal
from collections import defaultdict
import math
from random import sample
from io import TextIOWrapper

""" argument can be a file descriptor (compressed or not) or file path (compressed or not) and return a readable file descriptor"""
def read_compressed_or_not(file_or_file_path):
    file = file_or_file_path
    if type(file) == str:
        file = open(file,"rb")
    else:
        try:
            file = open(file.name,"rb")
        except:
            return(file)
    if file.read(2).startswith(b'\x1f\x8b'):
        file.seek(0)
        return(TextIOWrapper(gzip.open(filename=file, mode = "r")))
    else:
        file.close()
        file = open(file.name,"r")
        return(file)

""" The number of combinations of n things taken k at a time."""
def comb_k_n(k,n):
    if (k == 0):
            return 1
    result = 1
    for i in range(0, k):
            result *= float(n - i)/(i + 1);
    return sys.maxsize if result==float("Inf") else int(round(result)) 

# proportional sampling
def samplingCombinations(items, sample_ratio, sample_min, sample_max=100, step = 1):
    samplingCombinationList = defaultdict(list)
    item_size = len(items)
    combTotNb = pow(2,item_size)-1
    for k in range(1, item_size+1, step):
        tmp_comb = []
        combNb = Decimal(comb_k_n(item_size, k))
        combNb = sys.float_info.max if combNb>sys.float_info.max else combNb # to avoid to reach infinit values
        combNb_sample = math.ceil(Decimal(combNb)/Decimal(sample_ratio))
        # Plus petit echantillonage possible pour un k donn<C3><A9> = sample_min
        if ((combNb_sample < sample_min) and k != item_size):
            combNb_sample = sample_min
        # Plus grand echantillonage possible
        if (sample_max != None and (combNb_sample > sample_max)):
            combNb_sample = sample_max
        
        i = 0;
        while(i < combNb_sample):
            comb_sub = sample(items,k)
            # Echantillonnage sans remise
            if (comb_sub not in tmp_comb):
                tmp_comb.append(comb_sub)
                samplingCombinationList[len(comb_sub)].append(comb_sub)
            i+=1
    return samplingCombinationList

"""simple arithmetic mean"""
def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

"""simple median"""
def median(numbers):
    numbers = sorted(numbers)
    n = len(numbers)
    if n == 0:
        return(None)
    if n%2 == 1:
        return numbers[n//2]
    else:
        i = n//2
        return (numbers[i - 1] + numbers[i])/2
# !/usr/bin/python
# -*- coding: iso-8859-1 -*-

from util import *
import shlex

import threading
import os

import argparse
import time
from multiprocessing import Pool
from joblib import Parallel, delayed
import multiprocessing
import subprocess

TOOL = "/home/ggautrea/pangenome/pangenomeNEM/src/nem.py"

logging.basicConfig(level = logging.DEBUG, format = '\n%(asctime)s %(filename)s:l%(lineno)d %(levelname)s\t%(message)s', datefmt='%H:%M:%S')

def removeFinishedProcesses(processes):
    """ given a list of (commandString, process), 
        remove those that have completed and return the result 
    """
    newProcs = []
    for pollCmd, pollProc in processes:
        retCode = pollProc.poll()
        if retCode==None:
            # still running
            newProcs.append((pollCmd, pollProc))
        elif retCode!=0:
            # failed
            raise Exception("Command %s failed" % pollCmd)
        else:
            logging.info("Command %s completed successfully" % pollCmd)
    return newProcs

def runCommands(commands, maxCpu):
            processes = []
            for command in commands:
                logging.info("Starting process %s" % command)
                proc =  subprocess.Popen(shlex.split(command))
                procTuple = (command, proc)
                processes.append(procTuple)
                while len(processes) >= maxCpu:
                    time.sleep(.2)
                    processes = removeFinishedProcesses(processes)

            # wait for all processes
            while len(processes)>0:
                time.sleep(0.5)
                processes = removeFinishedProcesses(processes)
            logging.info("All processes completed")

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Resample and run multiple pangenome analyses to mesure evolution rate of partition metrics')
    parser.add_argument('-o', '--organisms', type=argparse.FileType('r'), nargs=1, help="The tsv file provided by progenome containing the gene annotations")
    parser.add_argument('-g', '--gene_families', type=argparse.FileType('r'), nargs=1, help="The tsv file provided by progenome containing eggNOG orthologous groups related to the annotated genomes")
    parser.add_argument('-d', '--outputdirectory', type=str, nargs=1, default="output.dir", help="The output directory", required=True)
    parser.add_argument("-t", "--num_thread", type=int, default=1, nargs=1, help="The number of thread to use, 0 for autodetection")
    parser.add_argument("-]", "--max_resampling", default = 30, nargs=1, help="Number max of subsamples in each combinaison of organisms")
    parser.add_argument("-[", "--min_resampling", default = 10, nargs=1, help="Number min of subsamples in each combinaison of organisms")

    options = parser.parse_args()

    families = options.gene_families[0].name

    OUTPUTDIR = options.outputdirectory[0]
    if not os.path.exists(OUTPUTDIR):
        os.makedirs(OUTPUTDIR)

    organisms = list()    
    arguments = list()     

    with options.organisms[0] as org_file:
        for row in org_file:
            organisms.append(row)

    total_combinations    = oidsCombinations(range(0,len(organisms)),options.min_resampling,options.max_resampling,options.min_resampling)
    nb_total_combinations = 0
    for comb_size in total_combinations:
        if comb_size>5:
            nb_total_combinations += len(total_combinations[comb_size])
            for combination in total_combinations[comb_size]:
                arguments.append([organisms[i] for i in combination])
    print("..... Preparing to classify "+str(nb_total_combinations)+" subsampled pangenomes .....")

    random.shuffle(arguments)
    arguments.insert(0, list(organisms))
    commands = []
    for i, arg in enumerate(arguments):
        sub_name = str(len(arg))+"_"+str(i)
        sub_organism_file = open(OUTPUTDIR+"/"+sub_name+".list","w")
        sub_organism_file.writelines(arg)
        commands.append("python3 "+TOOL+" -o "+sub_organism_file.name+" -g "+families+" -d "+OUTPUTDIR+"/"+sub_name)

    runCommands(commands,options.num_thread[0])
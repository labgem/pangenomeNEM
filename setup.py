#!/usr/bin/python3
# -*- coding: iso-8859-1 -*-

import os, sys
from setuptools import setup, find_packages
from distutils.extension import Extension
import logging
import subprocess
from distutils.command.install import install
from distutils.command.build import build
from distutils.command.clean import clean
from Cython.Build import cythonize

name = find_packages().pop()
NEM_dir_path = name+"/NEM/"

# print(NEM_dir_path)
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

# class pre_build(build):
#     def run(self):
#         print("pre_build")
#         proc = subprocess.Popen("make -C "+NEM_dir_path+" && mv "+NEM_dir_path+"nem_exe "+name, shell=True)
#         proc.communicate()
#         build.run(self)

# class post_clean(clean):
#     def run(self):
#         clean.run(self)
#         print("post_clean")
#         proc = subprocess.Popen("make -C "+NEM_dir_path+" remobj", shell=True)
#         proc.communicate()
#         os.remove(name+"/nem_exe")
#         os.remove(NEM_dir_path+"tcpu")
#         os.remove(NEM_dir_path+"randord")
#         os.remove(NEM_dir_path+"geo2nei")
#         os.remove(NEM_dir_path+"txt2hlp")

if __name__ == "__main__":

    # if (sys.argv[1]=="install" and (not os.path.exists("ppanggolin/nem_exe"))):
    #     print("run 'python setup.py build' before running 'python setup.py install' to compile libraries")
    #     exit(0)

    setup(
        name = name,
        version = read("VERSION"),
        author = "Guillaume GAUTREAU",
        author_email = "ggautrea@genoscope.cns.fr",
        description = "Depict microbial diversity via a partitioned pangenome graph",
        license = "CeCILL-2.1",
        keywords = "pangenome comparative-genomics bioinformatics microbiology",
        url = "https://github.com/ggautreau/PPanGGOLiN",
        packages=[name],
        #package_dir = {name: 'lib'},
        long_description=read('README.rst'),
        classifiers=[
            "Environment :: Console",
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Programming Language :: Python",
            "Programming Language :: C",
            "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
        ],
        entry_points={
            'console_scripts': [
            name+' = '+name+'.command_line:__main__'
          ]},
        extras_require= {'all' : ['futures', 'collections', 'ordered-set', 'networkx >= 2.0', 'numpy', 'scipy','community', 'tqdm', 'python-highcharts'],
                        ':python_version == "2.7"': ['futures']},
        #cmdclass={'build': pre_build, 'clean': post_clean},
        #ext_modules = cythonize(Extension("nem",sources = [NEM_dir_path+'nem.pyx'], libraries=[NEM_dir_path])))#,
        ext_modules = cythonize([Extension(name = "nem",sources =[NEM_dir_path+'nem.pyx',
                                                                  NEM_dir_path+'nem_exe.c',
                                                                  NEM_dir_path+'nem_alg.c',
                                                                  NEM_dir_path+'nem_nei.c',
                                                                  NEM_dir_path+'nem_mod.c',
                                                                  NEM_dir_path+'nem_rnd.c',
                                                                  NEM_dir_path+'lib_io.c',
                                                                  NEM_dir_path+'nem_hlp.c',
                                                                  NEM_dir_path+'genmemo.c'],
                                                        include_dirs=[NEM_dir_path])]))

# ],include_dirs=[NEM_dir_path])])#, 

# from ctypes import *
# import nem

# nem_functions = cdll.LoadLibrary(nem.__file__)

# def call_c(L):
#     arr = (c_char_p * len(L))()
#     arr[:] = L
#     nem_functions.mainfunc(c_int(len(L)), arr)

# a = [b"nem",b"tests/yersinia_amandine/test/NEM_results/nem_file" ,b"3", b"-a", b"ncem" ,b"-i" ,b"100" ,b"-m" ,b"bern" ,b"pk" ,b"skd" ,b"-s" ,b"m" ,b"tests/yersinia_amandine/test/NEM_results/nem_file.m", b"-b", b"0.5" ,b"-n" ,b"f", b"-c", b"clas", b"0.00001"]
# call_c(a)

    #NEM_dir_path+'geo2nei.c',
#NEM_dir_path+'err2.c',
#NEM_dir_path+'randord.c',
#NEM_dir_path+'exemain.c',
#NEM_dir_path+'txt2hlp.c',
#NEM_dir_path+'tcpu.c',
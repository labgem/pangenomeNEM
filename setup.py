import os
from setuptools import setup, Extension

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

NEM_dir_path = "NEM"

nem = Extension('nem',
                sources = [#NEM_dir_path+"/exemain.c",
                           NEM_dir_path+"/nem_hlp.c",
                           NEM_dir_path+"/nem_mod.c",
                           NEM_dir_path+"/nem_nei.c",
                           #NEM_dir_path+"/randord.c",
                           #NEM_dir_path+"/txt2hlp.c",
                           NEM_dir_path+"/exememo.c",
                           NEM_dir_path+"/lib_io.c",
                           NEM_dir_path+"/nem_alg.c",
                           NEM_dir_path+"/nem_arg.c",
                           NEM_dir_path+"/nem_exe.c",
                           NEM_dir_path+"/nem_rnd.c",
                           NEM_dir_path+"/nem_ver.c"],
                include_dirs = [NEM_dir_path+"headers/"],
                extra_compile_args = [''])

if __name__ == "__main__":

    setup(
        name = "ppanggolin",
        version = "0.1",
        author = "Guillaume GAUTREAU",
        author_email = "ggautrea@genoscope.cns.fr",
        description = (""),
        license = "CeCILL-2.1",
        keywords = "pangenome graph microbial diversity",
        #url = "http://packages.python.org/ppanggolin",
        packages=['ppanggolin'],
        long_description=read('README.md'),
        classifiers=[
            "Environment :: Console",
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Bio-Informatics",
            "Programming Language :: Python",
            "Programming Language :: C",
            "License :: OSI Approved :: CEA CNRS Inria Logiciel Libre License, version 2.1 (CeCILL-2.1)",
        ],
        ext_modules=[nem],
        extras_require= {'all' : [ 'collections', 'ordered_set', 'networkx>=1.0', 'gffutils', 'numpy', 'community']})
        
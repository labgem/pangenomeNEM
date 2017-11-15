import os
from setuptools import setup
import logging
import subprocess

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

NEM_dir_path = "ppanggolin/NEM"

proc = subprocess.Popen("cd "+NEM_dir_path+" && make && cp nem_exe .. && cd ..", shell=True)
proc.communicate()

if __name__ == "__main__":
    setup(
        name = "ppanggolin",
        version = "0.0.1",
        author = "Guillaume GAUTREAU",
        author_email = "ggautrea@genoscope.cns.fr",
        description = "Depict microbial diversity via a partitioned pangenome graph",
        license = "CeCILL-2.1",
        keywords = "pangenome comparative-genomics bioinformatics microbiology",
        url = "https://github.com/ggautreau/PPanGGOLiN",
        packages=['ppanggolin'],
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
        include_package_data=True,
        entry_points={
            'console_scripts': [
            'ppanggolin = ppanggolin.command_line:__main__'
          ]},
        extras_require= {'all' : [ 'collections', 'ordered-set', 'networkx >= 2.0', 'numpy', 'community', 'tqdm']})
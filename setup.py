import os
from setuptools import setup, Extension
import logging
import subprocess
from ppanggolin import command_line

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

NEM_dir_path = "ppanggolin/NEM"

proc = subprocess.Popen("cd "+NEM_dir_path+" && make && cp nem_exe .. && cd ..", shell=True)
proc.communicate()

# nem_source = Extension('nem_source',
#                 sources = [NEM_dir_path+"/nem_arg.c",
#                            NEM_dir_path+"/nem_hlp.c",
#                            NEM_dir_path+"/nem_mod.c",
#                            NEM_dir_path+"/nem_nei.c",
#                            # NEM_dir_path+"/txt2hlp.c",
#                            NEM_dir_path+"/lib_io.c",
#                            NEM_dir_path+"/exememo.c",
#                            NEM_dir_path+"/nem_alg.c",
#                            NEM_dir_path+"/nem_exe.c",
#                            NEM_dir_path+"/nem_rnd.c",
#                            NEM_dir_path+"/nem_ver.c"],#   NEM_dir_path+"/randord.c"
#                 include_dirs = [NEM_dir_path+"/headers/",NEM_dir_path+"/nem_arg.c",
#                            NEM_dir_path+"/nem_hlp.c",
#                            NEM_dir_path+"/nem_mod.c",
#                            NEM_dir_path+"/nem_nei.c",
#                            # NEM_dir_path+"/txt2hlp.c",
#                            NEM_dir_path+"/lib_io.c",
#                            NEM_dir_path+"/exememo.c",
#                            NEM_dir_path+"/nem_alg.c",
#                            NEM_dir_path+"/nem_exe.c",
#                            NEM_dir_path+"/nem_rnd.c",
#                            NEM_dir_path+"/nem_ver.c"],
#                 )


#nem = Extension(name = "nem", sources = )

# nem = Extension(name = "nem", sources=[NEM_dir_path+"/nem.pyx",
#                                        NEM_dir_path+"/nem_arg.c",
#                                        NEM_dir_path+"/nem_hlp.c",
#                                        NEM_dir_path+"/nem_mod.c",
#                                        NEM_dir_path+"/nem_nei.c",
#                                        # NEM_dir_path+"/txt2hlp.c",
#                                        NEM_dir_path+"/lib_io.c",
#                                        NEM_dir_path+"/exememo.c",
#                                        NEM_dir_path+"/nem_alg.c",
#                                        # NEM_dir_path+"/nem_exe.c",
#                                        NEM_dir_path+"/nem_rnd.c",
#                                        NEM_dir_path+"/nem_ver.c"],
#                                        extra_compile_args = ['-lm -Qunused-arguments -Wno-unused-function -Wno-unused-variable']) 

if __name__ == "__main__":
    setup(
        name = "ppanggolin",
        version = command_line.__version__,
        author = "Guillaume GAUTREAU",
        author_email = "ggautrea@genoscope.cns.fr",
        description = "Depict microbial diversity via a partitioned pangenome graph",
        license = "CeCILL-2.1",
        keywords = "pangenome graph microbial diversity",
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
        #data_files=[('nem', ['NEM/nem_exe'])],
        # ext_modules=cythonize(nem),
        entry_points={
            'console_scripts': [
            'ppanggolin = ppanggolin.command_line:__main__'
          ]},
        extras_require= {'all' : [ 'collections', 'ordered-set', 'networkx >= 2.0', 'numpy', 'community', 'tqdm']}
        )
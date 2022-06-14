from distutils.core import setup, Extension
import numpy as np
import os
import sys

### change the path to IPOPT
#path Ipopt
path_ipopt = '/Users/pfh/lib/Ipopt-3.12.4/'


### nothing to do

#path of current example, you don't need to change it
path_example = os.path.abspath(".")
#path AMMSoCK
path_ammsock = os.path.dirname(os.path.abspath(sys.argv[0]))


#set directions
IPOPT_ICLUDE_DIRS=[path_ammsock+'/mex',path_example + '/cpp',path_ipopt + '/build/include/coin',np.get_include()]
#IPOPT_LIB_DIRS=[path_ipopt + '/build/lib64', path_ipopt + '/build/lib']
IPOPT_LIB_DIRS=[path_ipopt + '/build/lib64', path_ipopt + '/build/lib','/usr/local/Cellar/gcc/5.3.0/lib/gcc/5/gcc/x86_64-apple-darwin15.0.0/5.3.0/../../../']
IPOPT_LIBS=['ipopt','coinhsl','lapack','blas','gfortran','stdc++']

spammodule = Extension(
    'performReduction', 
    sources=[path_ammsock+'/mex/performReduction.cpp',path_ammsock+'/mex/ammsockNLP.cpp'],
    include_dirs=IPOPT_ICLUDE_DIRS,
    libraries=IPOPT_LIBS,
    library_dirs=IPOPT_LIB_DIRS,
    language='c++',
    extra_compile_args=['-w']
)
setup (ext_modules=[spammodule])

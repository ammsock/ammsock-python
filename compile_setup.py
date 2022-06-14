import os

def compile_setup(path_python, path_ammsock):

    #create string with compilation arguments
    cmd = path_python + ' ' + path_ammsock + '/setup.py build_ext --inplace --force'
    print(cmd)
    os.system(cmd)	#call cmd in shell

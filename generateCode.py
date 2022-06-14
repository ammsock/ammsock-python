def generateCode(path_ammsock,cur_path,mechfile):
    import sys
    import os

    #command to compile codeGenerator.cpp
    cmd = 'g++ -w ' + str(path_ammsock) +'/mex/codegenerator.cpp -o generateCode -I' + str(os.getcwd())

    status = os.system(cmd)  #calls cmd in shell

    #check if compile was successful
    if status != 0:
        sys.exit('Error: Can not compile codeGenerator.')

    if not os.path.exists('cpp'):  # check if directory "cpp" exists
        os.makedirs('cpp')         # create "cpp" directory

    os.chdir('./cpp')               #change to directory "cpp" 

    cmd = '../generateCode ' + str(cur_path) + '/mech ' + str(mechfile) + ' c++'

    os.system(cmd)

    os.chdir('..')                  #change directory

    if not os.path.exists('python'):
        os.makedirs('python')

    os.chdir('./python')            #change to directory "python"

    cmd = '../generateCode ' + str(cur_path) + '/mech ' + str(mechfile) + ' python'
    os.system(cmd)		# calls cmd in shell

    os.chdir('..')
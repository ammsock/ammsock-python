def readData(path_ammsock,cur_path,mechfile):
    import os
    import sys
    os.chdir('./mech')    #change to directory "mech"

    cmd = 'g++ -w  -std=c++11 ' + str(path_ammsock) + '/mex/readData.cpp -o readData -I' + str(os.getcwd())

    status = os.system(cmd)  #calls cmd in shell

    if status != 0:
        sys.exit('Error: Can not compile readData.')

    if not os.path.exists('data'):   #check if directory "data" exists
        os.makedirs('data')          #create "data" directory

    os.chdir('..')                  #change directory

    if not os.path.exists('cpp'):
        os.makedirs('cpp')

    if not os.path.exists('python'):
        os.makedirs('python')

    cmd = './mech/readData ' + str(cur_path) + '/mech ' + str(mechfile)

    os.system(cmd)	#calls cmd in shell

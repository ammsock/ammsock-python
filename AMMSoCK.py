import numpy as np      #import mathematical constants and structures
import sys              #need to print errors

class AMMSoCK:
    #AMMSoCK as Pythonclass

    #properties

    """
    problemname:   	name of the current problem
    nspec:		number of species
    nreac:		number of ractions
    natom:		number of atoms
    ntb:		number of third bodys
    nrtb:		number of third body reactions
    nfr:		number of just forward reactions
    nstages:		number of stages
    nrpv:		number of reaction progess variables
    nint:		number of collocation intervals
    nnop:		number of elements in Inop
    nop:		number of elements in Iop
    indexRpv:		indexes of reaction progess variables
    Iop:		Iop indexes
    Inop:		Inop indexes
    Acoll:
    species:		list of all occuring species
    atoms:		list of all occuring atoms
    workingDir:		current working direction
    oxidizer:		used oxidizer
    fuel:		used fuel
    exactHessian:	boolean
    """


    #constructor
    def __init__(self,name):
        self.problemname = name

    #methods

    def setGenerationParameter(self):   #writes required variables in "generationStats.hpp"
        fid = open(self.workingDir+'/cpp/generationStats.hpp','w')
        fid.write('#define NSPEC ' + str(self.nspec) + '\n')
        fid.write('#define NREAC ' + str(self.nreac) + '\n')
        fid.write('#define NATOM ' + str(self.natom) + '\n')
        fid.write('#define NTB ' + str(self.ntb) + '\n')
        fid.write('#define NRTB ' + str(self.nrtb) + '\n')
        fid.write('#define NRPV ' + str(self.nrpv) + '\n')
        fid.write('#define INDRPV { ')
        for i in range(0,len(self.indexRpv)-1):
            fid.write(str(np.where(self.Iop == self.indexRpv[i])[0][0]) + ', ')
        fid.write(str(list(self.Iop).index(self.indexRpv[-1])) + ' }\n')

        fid.write('#define NSTAGES ' + str(self.nstages) + '\n')
        fid.write('#define NINT ' + str(self.nint) + '\n')
        fid.write('#define ACOLL {')
        for i in range(0,self.nstages):
            fid.write('{')
            for j in range(0,self.nstages-1):
                fid.write("%1.17f"% self.Acoll[i][j] + ', ')
            fid.write("%1.17f"% self.Acoll[i][self.nstages-1])
            if i+1<self.nstages:
                fid.write('}, ')
            else:
                fid.write('}')
        fid.write('}\n')
        fid.write('#define NNOP ' + str(self.nnop) + '\n')
        fid.write('#define INOP {')
        for i in range(0,len(self.Inop)-1):
            fid.write(str(self.Inop[i]) + ', ')
        if list(self.Inop):
            fid.write(str(self.Inop[-1]) + ' }\n')
        else:
            fid.write('}\n')

        fid.write('#define NOP ' + str(self.nop) + '\n')
        fid.write('#define IOP {')
        for i in range(0,len(self.Iop)-1):
            fid.write(str(self.Iop[i]) + ', ')
        fid.write(str(self.Iop[-1]) + ' }\n')

        fid.write('#define FUEL "' + str(self.fuel) + '"\n')
        fid.write('#define OXIDIZER "' + str(self.oxidizer) + '"\n')
        if self.exactHessian is 'true':
            fid.write('#define EXACTHESSIAN ' + str(1) + '\n')
        else:
            fid.write('#define EXACTHESSIAN ' + str(0) + '\n')
        fid.close()

    def getParserStatistics(self):
        #loads corresponding files of the chemical reaction equation and assign values to corresponding variables
        fid = open(self.workingDir+'/mech/data/readerStats.dat','r')

        self.nreac = int(fid.readline())
        self.nspec = int(fid.readline())
        self.species = []

        for i in range(0,self.nspec):
            self.species.append(fid.readline())
            self.species[i] = self.species[i][:-1]      #deletes blank at the end of the row

        self.natom = int(fid.readline())
        self.atoms = []

        for i in range(0,self.natom):
            self.atoms.append(fid.readline())
            self.atoms[i] = self.atoms[i][:-1]

        self.ntb = int(fid.readline())
        self.nrtb = int(fid.readline())
        self.nfr = int(fid.readline())

    def setCombustion(self,fuel,oxidizer):   #assign oxidation variable
        fuelContained = 0
        #check whether species contains selected fuel
        for i in range(0,self.nspec):
            if (fuel == self.species[i]):
                fuelContained = 1
        if not fuelContained:
            sys.exit('Error: Please choose a valid species as fuel.')
        self.fuel = fuel

        #check if selected oxidizer is allowed
        if not (oxidizer =='Air(N2,O2)' or oxidizer == 'Air(N2,O2,Ar)' or oxidizer == 'O2'):
            sys.exit('Error: Oxidizer must be Air(N2,O2), Air(N2,O2,Ar) or O2.')
        self.oxidizer = oxidizer

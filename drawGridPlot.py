import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

def drawGridPlot(y,ammsock,indrpv,style):
    nspec = ammsock.nspec
    species = ammsock.species
    indexSet = np.setxor1d(np.array(range(1,nspec+1)),indrpv)

    nplots = nspec-indrpv.size
    dim = np.floor(np.sqrt(nplots)+1)

    #check if figure still exists
    if not plt.fignum_exists(1):
        global fig
        fig = plt.figure(1)

        #2D Plot
        if indrpv.size == 1:
            #check if subplots still exists
            if not 'ax1' in globals():
                #create subplots for species
                for i in range(0,nplots):
                    globals()['ax'+str(i+1)] = fig.add_subplot(dim,dim,i+1)
                    exec('ax'+str(i+1)+'.set_xlabel("Y_{'+str(species[indrpv[0]-1])+'}")')
                    exec('ax'+str(i+1)+'.set_ylabel("Y_{'+str(species[indexSet[i]-1])+'}")')
                    exec('ax'+str(i+1)+'.grid(True)')
                    #create subplot for temperature
                    globals()['ax'+str(nplots+1)] = fig.add_subplot(dim,dim,nplots+1)

                    exec('ax'+str(nplots+1)+'.set_xlabel("Y_{'+str(species[indrpv[0]-1])+'}")')
                    exec('ax'+str(nplots+1)+'.set_ylabel("T")')
                    exec('ax'+str(nplots+1)+'.grid(True)')
                #plot values for species
                if isinstance(style,str):
                    for i in range(0,nplots):
                        exec('ax'+str(i+1)+'.plot(y[:,indrpv[0]-1],y[:,indexSet['+str(i)+']-1],style)')
                        #plot vlaues for temperature
                        exec('ax'+str(nplots+1)+'.plot(y[:,indrpv[0]-1],y[:,-1],style)')
                elif isinstance(style,dict):
                    for i in range(0,nplots):
                        exec('ax'+str(i+1)+'.plot(y[:,indrpv[0]-1],y[:,indexSet['+str(i)+']-1],color=style["color"],linestyle=style["linestyle"],marker=style["marker"])')
                        #plot vlaues for temperature
                        exec('ax'+str(nplots+1)+'.plot(y[:,indrpv[0]-1],y[:,-1],color=style["color"],linestyle=style["linestyle"],marker=style["marker"])')
                else:
                    sys.exit('style has to be a string or a dictionary')

            #3D Plot
        elif indrpv.size == 2:
            #check if subplots still exists
            if not 'ax1' in globals():
                #create subplots for species
                for i in range(0,nplots):
                    globals()['ax'+str(i+1)] = fig.add_subplot(dim,dim,i+1,projection='3d')
                    exec('ax'+str(i+1)+'.set_xlabel("Y_{'+str(species[indrpv[0]-1])+'}")')
                    exec('ax'+str(i+1)+'.set_ylabel("Y_{'+str(species[indrpv[1]-1])+'}")')
                    exec('ax'+str(i+1)+'.set_zlabel("Y_{'+str(species[indexSet[i]-1])+'}")')
                    exec('ax'+str(i+1)+'.grid(True)')
                    #create subplot for temperature
                    globals()['ax'+str(nplots+1)] = fig.add_subplot(dim,dim,nplots+1,projection='3d')

                    exec('ax'+str(nplots+1)+'.set_xlabel("Y_{'+str(species[indrpv[0]-1])+'}")')
                    exec('ax'+str(nplots+1)+'.set_ylabel("Y_{'+str(species[indrpv[1]-1])+'}")')
                    exec('ax'+str(nplots+1)+'.set_zlabel("T")')
                    exec('ax'+str(nplots+1)+'.grid(True)')
                #plot values for species
                if isinstance(style,str):
                    for i in range(0,nplots):
                        exec('ax'+str(i+1)+'.plot(y[:,indrpv[0]-1],y[:,indrpv[1]-1],y[:,indexSet['+str(i)+']-1],style)')
                        #plot values for temperature
                        exec('ax'+str(nplots+1)+'.plot(y[:,indrpv[0]-1],y[:,indrpv[1]-1],y[:,-1],style)')
                elif isinstance(style,dict):
                    for i in range(0,nplots):
                        exec('ax'+str(i+1)+'.plot(y[:,indrpv[0]-1],y[:,indrpv[1]-1],y[:,indexSet['+str(i)+']-1],color=style["color"],linestyle=style["linestyle"],marker=style["marker"])')
                        #plot values for temperature
                        exec('ax'+str(nplots+1)+'.plot(y[:,indrpv[0]-1],y[:,indrpv[1]-1],y[:,-1],color=style["color"],linestyle=style["linestyle"],marker=style["marker"])')
                else:
                    sys.exit('style has to be a string or a dictionary')
            #error	
        else:
            sys.exit('Please select a 2d or 3d projction')


        return

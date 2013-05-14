from msmbuilder import io, Trajectory
from mpl_toolkits.mplot3d import Axes3D
import optparse
import pylab
import numpy
import glob
import os
from numpy import linalg

colors=['black', 'red', 'magenta', 'orange', 'orange', 'green', 'blue', 'cyan', 'purple']
colors=colors*40

font = {'family' : 'sans-serif',
        'sans-serif':['Helvetica'],
                'size'   : 24}
pylab.rcParams['lines.linewidth'] = 2

pylab.rc('font', **font)
params = {'legend.fontsize': 16,
          'legend.linewidth': 2}
pylab.rcParams.update(params)

def forceUpdate(event):
    global scatCollection
    scatCollection.changed()


percents=[20,10]
for percent in percents:
    sys='bi'
    k=3000
    lag=15
    type='strict'
    dir='/home/mlawrenz/wontkill/new-results/%s-bindb-pgeomx/structural-data/' % sys
    pathways=io.loadh('%s/Paths.lh5' % dir)
    fluxes=pathways['fluxes']
    refs=['h36', 'npxxy']
    types=range(0,20)
    file='%s_paths.txt' % sys
    fhandle=open(file)
    for (i, path) in enumerate(fhandle.readlines()):
        path=numpy.array(path.split())
        fig=pylab.figure(figsize=(15,13))
        ax = fig.add_subplot(111, projection='3d')
        #ax=fig.gca(projection='3d')
        print i, colors[i]        
        print "path %s" % i
        stat=fluxes[i]/fluxes[0]
        if stat < 0.3:
            break
        print "STAT is %s" % stat
        data1=numpy.loadtxt('%s/path%s.h36.dat' % (dir, i), usecols=(0,))
        data2=numpy.loadtxt('%s/inactive-npxxy-path%s.dat' % (dir, i), usecols=(0,))
        data3=numpy.loadtxt('%s/inactive-conn-path%s.dat' % (dir, i), usecols=(0,))
        volumes=[]
        for state in path:
            data=numpy.loadtxt('volmap-data/%s_state%s_volume.dat' % (sys, state))
            volumes.append(data)
        cm=pylab.cm.get_cmap('RdYlBu')
        scatCollection=ax.scatter(data1, data2, data3, linewidth=2, alpha=0.6,
c=volumes,vmin=40, vmax=200, cmap=cm, s=400)
        fig.canvas.mpl_connect('draw_event',forceUpdate)
        cbar=pylab.colorbar(scatCollection)
        #cbar.set_ticks([0.3, 1, 2.8])
        #cbar.set_ticklabels(['Antagonists', 'Equal', 'Agonists'])
        #cbar = fig.colorbar(scatCollection, ticks=[min(volumes), max(volumes)])
        #cbar.ax.set_yticklabels(['Antagonists',  'Agonists'])# vertically oriented colorbar
        fig.canvas.mpl_connect('draw_event',forceUpdate)
        #ax.plot(data1, data2, data3, c=colors[i], alpha=0.6, linewidth=((stat*20)-1))
        ax.set_xlim3d(6.5,15.5)
        ax.set_ylim3d(0.2,1.2)
        ax.set_zlim3d(0.5,2.0)
        ax.set_xlabel('Helix 3-6 Distance ($\AA$)', linespacing=5.2)
        ax.set_ylabel('NPxxY RMSD ($\AA$)', linespacing=5.1)
        ax.set_zlabel('Connector RMSD ($\AA$)', linespacing=5.4)
        ax.dist = 10
        pylab.title('Path%s %s%%Flux Pocket Volume' % (i, int(stat*100)))
        #pylab.legend(loc=2)
        #pylab.savefig('%s/select-%s-path%s-3d-top%s%%vol.png' % (dir, sys, i, percent), dpi=300)
        pylab.show()
    

import pickle, numpy, glob
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
import pylab
import optparse

def make_plot(type, pathdata, paths, combo=None):
    print type
    cm=pylab.cm.get_cmap('RdYlBu')
    for (n, line) in enumerate(paths.readlines()): 
        pylab.figure()
        if combo!=None:
            fig=pylab.scatter(range(0, len(line.split())), pathdata['volumes'][n], cmap=cm, c=combo[n],vmin=0.1, vmax=2.0, s=200)
        else:
            fig=pylab.scatter(range(0, len(line.split())), pathdata['volumes'][n], cmap=cm, c=pathdata[type][n], s=200)
        cbar=pylab.colorbar(fig)
        cbar.set_ticks([0.1, 1, 2.0])
        cbar.set_ticklabels(['Antagonists', 'Equal', 'Agonists'])
        pylab.legend()
        pylab.xlabel('Path %s Progress' % n)
        pylab.ylabel('binding site volume')
        pylab.show()
    #slope, intercept, R, pval, std_err = stats.linregress(data1, data2)
    #print type, R, R**2, pval
    #(ar,br)=polyfit(data1, data2, 1)
    #xr=polyval([ar,br], data1)    
    #pylab.plot(chems,xr,'r-')


def main(sys, percent):
    percent=float(percent)*100
    types=['agonist', 'antagonist']
    pathdata=dict()
    combo_data=dict()
    paths=open('%s_paths.txt' % sys)
    volume_data = pathdata.setdefault('volumes', {})
    for (n, path) in enumerate(paths.readlines()):
        volume_path = volume_data.setdefault(n, [])
        for (m, state) in enumerate(path.split()):
            volume_path.append(numpy.loadtxt('volmap-data/%s_state%s_volume.dat' % (sys, state)))
    for type in types:
        type_data = pathdata.setdefault(type, {})
        chems=numpy.loadtxt('enrichment/top%s%s.txt' % (percent, type), usecols=(1,))
        states=numpy.loadtxt('enrichment/top%s%s.txt' % (percent, type), usecols=(0,), dtype=int)
        paths=open('%s_paths.txt' % sys)
        for (n, path) in enumerate(paths.readlines()):
            type_path= type_data.setdefault(n, [])
            for (m, state) in enumerate(path.split()):
                location=numpy.where(states==int(state))[0]
                type_path.append(chems[location])
        paths=open('%s_paths.txt' % sys)
        #make_plot(type, pathdata, paths)
    for (n, path) in enumerate(paths.readlines()):
        combo_path= combo_data.setdefault(n, [])
        for (m, state) in enumerate(path.split()):
            combo_path.append(pathdata['agonist'][n][m]/float(pathdata['antagonist'][n][m]))
    paths=open('%s_paths.txt' % sys)
    make_plot('combo', pathdata, paths, combo_data)

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--system', dest='sys',
                      help='system')
    parser.add_option('-p', '--percent', dest='percent',
            help='top %% of database to analyze')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(sys=options.sys, percent=options.percent)


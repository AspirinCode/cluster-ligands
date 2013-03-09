import numpy 
import optparse
import pylab
import os
import glob

font = {'family' : 'sans-serif',
        'sans-serif':['Helvetica'],
                        'size'   : 16}

pylab.rcParams['lines.linewidth'] = 2

pylab.rc('font', **font)
params = {'legend.fontsize': 14,
          'legend.linewidth': 2}
pylab.rcParams.update(params)


def main(state):
    dir='/home/mlawrenz/wontkill/new-results/apo-%s' % state
    files=glob.glob('%s/eon-only-molecule*rpt' % dir)
    data=dict()
    data['combo']=[]
    for file in files:
        name=os.path.dirname(file)+'/mod-'+os.path.basename(file)
        os.system('sed "1d" < %s > %s' % (file, name))
        tmp1=numpy.loadtxt(name, usecols=(3,))    
        tmp3=numpy.loadtxt(name, usecols=(6,))    
        for (i, j) in zip(tmp1, tmp3):
            data['combo'].append((i+j))

    pylab.figure()
    type='combo'
    results=pylab.hist(data[type], alpha=0.7, bins=20, normed=True) #, label='shape+PB combo')
    colors=['black', 'magenta', 'cyan']
    cutoffs=[0.05, 0.10, 0.2]
    for (color, cutoff) in zip(colors, cutoffs):
        num=len(data[type])*cutoff
        test=sorted(data[type])[::-1][:int(num)]
        numpy.savetxt('%s_cutoff%s_eon_only_tanimotos.txt' % (type, cutoff), test)
        print cutoff, test[-1]
        list=numpy.arange(0,2,0.1)
        pylab.plot([test[-1]]*len(list), list, linewidth=5, color=color, label='%s%%' % int(cutoff*100))
    pylab.legend()
    pylab.xlim(0,2)
    pylab.ylim(0,2)
    pylab.xlabel('Shape + PB Tanimoto Scores')
    pylab.ylabel('Normed Probability')
    pylab.title('State %s' % state)
    pylab.savefig('%s/eon_only_agonist_apo%s_combo_tanimotos.png' % (dir, state), dpi=300)
    pylab.show()

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-s', '--state', dest='state',
                      help='state')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(state=options.state)
                                    

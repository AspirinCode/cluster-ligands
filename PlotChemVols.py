import pickle, numpy, glob
from scipy import linspace, polyval, polyfit, sqrt, stats, randn
import pylab
import optparse

def make_plot(type, data1, data2, percent):
    slope, intercept, R, pval, std_err = stats.linregress(data1, data2)
    print type, R, R**2, pval
    (ar,br)=polyfit(data1, data2, 1)
    xr=polyval([ar,br], data1)    
    pylab.figure()
    pylab.scatter(data1, data2, label='%s' % type)
    #pylab.plot(chems,xr,'r-')
    pylab.legend()
    pylab.xlabel('%s recovery' % type)
    pylab.ylabel('binding site volume')
    pylab.title('R=%s, R$^2$=%s, pval=%s' % (round(R,2), round(R**2,2), round(pval,2)))
    pylab.savefig('./enrichment/%srecov_vol_top%s_summary.png' % (type, percent*100), dpi=300)
    pylab.show()


def main(sys, percent):
    percent=float(percent)*100
    types=['agonist', 'antagonist']
    combo=dict()
    for type in types:
        chems=numpy.loadtxt('enrichment/top%s%s.txt' % (percent, type), usecols=(1,))
        combo[type]=chems
        states=numpy.loadtxt('enrichment/top%s%s.txt' % (percent, type), usecols=(0,), dtype=int)
        volumes=[]
        for state in states:
            data=numpy.loadtxt('volmap-data/%s_state%s_volume.dat' % (sys, state))
            volumes.append(data)
        make_plot(type, chems, volumes, percent)
    combos=[]
    for (x, y) in zip(combo['agonist'], combo['antagonist']):
        combos.append(x/y)
    make_plot('combo', combos, volumes, percent)

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


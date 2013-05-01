import numpy 
import pylab
import glob

def main(state):
    files=glob.glob('mod-eon-only-molecule*rpt')
    data=dict()
data['pb']=[]
data['coul']=[]
data['shape']=[]
data['combo']=[]
for file in files:
    tmp1=numpy.loadtxt(file, usecols=(3,))    
    for i in tmp1:
        data['pb'].append(i)
    tmp2=numpy.loadtxt(file, usecols=(4,))    
    for i in tmp2:
        data['coul'].append(i)
    tmp3=numpy.loadtxt(file, usecols=(6,))    
    for i in tmp3:
        data['shape'].append(i)
    for (i, j) in zip(tmp1, tmp3):
        data['combo'].append((i+j))

types=['pb', 'coul', 'shape', 'combo']
types=['combo']
pylab.figure()
for type in types:
    if type=='combo':
        results=pylab.hist(data[type], alpha=0.7, bins=20) #, label='shape+PB combo')
    else:
        results=pylab.hist(data[type], alpha=0.7, bins=10, label=type)
    colors=['black', 'magenta', 'cyan']  
    cutoffs=[0.05, 0.10, 0.2]
    for (color, cutoff) in zip(colors, cutoffs):
        num=len(data[type])*cutoff
        test=sorted(data[type])[::-1][:int(num)]
        numpy.savetxt('%s_cutoff%s_eon_only_tanimotos.txt' % (type, cutoff), test)
        print cutoff, test[-1]
        pylab.plot([test[-1]]*max(results[0]), range(0, max(results[0])), linewidth=5, color=color, label='%s%%' % (cutoff*100))
pylab.legend()
pylab.xlabel('Shape + PB Tanimoto Scores')
pylab.ylabel('Population')
pylab.savefig('eon_only_agonist_apo2665_combo_tanimotos.png', dpi=300)
pylab.show()

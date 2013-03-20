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

def name_catch(line, index1, index2):
    if "CL" == line.split()[1] or "ICI"==line.split()[1]:
        index1+=1   
        index2+=1   
        catch=True
    elif "CL" == line.split()[0] or "ICI"==line.split()[0]:
        index1+=1   
        index2+=1   
        catch=True
        if "CL" == line.split()[2] or "ICI"==line.split()[2]:
            index1+=1   
            index2+=1   
            catch=True
        else:
            pass
    elif "CL" == line.split()[2] or "ICI" == line.split()[2]:
        index1+=1   
        index2+=1   
        catch=True
        if "CL" == line.split()[0] or "ICI"==line.split()[0]:
            index1+=1   
            index2+=1   
            catch=True
        else:
            pass
    else:
        catch=False
    return catch, index1, index2


def main(dir, prefix):
    state=dir.split('results/')[1].strip('/')
    files=glob.glob('%s/%s*rpt' % (dir, prefix))
    data=dict()
    data['combo']=[]
    for file in files:
        fhandle=open(file)
        for line in fhandle.readlines():
            catch=False
            index1=3
            index2=6
            if 'Rank' not in line:
                if catch is False:
                    catch, index1, index2=name_catch(line,index1, index2)  
                else:
                    pass
                sum=float(line.split()[index1])+float(line.split()[index2])
                if sum<0:
                    data['combo'].append(0)
                elif sum >0 and sum<= 2.0:
                    data['combo'].append((sum))
                elif sum>2.0:
                    import pdb
                    pdb.set_trace()
                    print sum
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
    pylab.title('%s' % state)
    pylab.savefig('%s/%s_%s_combo_tanimotos.png' % (dir, prefix, state), dpi=300)
    pylab.show()

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dir', dest='dir',
                      help='dir')
    parser.add_option('-p', '--prefix', dest='prefix',
                      help='prefix')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(dir=options.dir, prefix=options.prefix)
                                    

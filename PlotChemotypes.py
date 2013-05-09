import pickle, numpy, glob
import pylab
import optparse

def get_types(refdata):
    types=dict()
    for x in refdata.keys():
        if 'mlawrenz' not in refdata[x]:
            if x!='rac-2':
                if refdata[x] not in types.keys():
                    types[refdata[x]]=0
    return types

def get_counts(refdata):
    total_agonist=0
    total_antagonist=0
    counts=dict()
    for x in refdata.keys():
        if 'mlawrenz' not in refdata[x]:
            if x!='rac-2':
                if refdata[x] not in counts.keys():
                    counts[refdata[x]]=1
                    if 'antagonist' in refdata[x]:
                        total_antagonist+=1
                    else:
                        total_agonist+=1
                else:
                    counts[refdata[x]]+=1
                    if 'antagonist' in refdata[x]:
                        total_antagonist+=1
                    else:
                        total_agonist+=1
    return counts, total_agonist, total_antagonist

def main(sys, percent):
    percent=float(percent)
    file=open('/home/mlawrenz/wontkill/new-results/annotate_dbase.pickle', 'rb')
    data=pickle.load(file)
    # Clenbutarol
    data['rac-2']='beta2 agonist'
    exceptions=['CGP-20712A','rac-3','rac-2','pos-Zilpaterol','rac-Zilpaterol', 'omega-92', 'neg-Zilpaterol']
    refdata=dict()
    for x in data.keys():
        if x not in exceptions:
            key=x.split('-')[0]
            refdata[key]=data[x]
        else:
            refdata[x]=data[x]

    types=get_types(refdata)
    counts, total_agonist, total_antagonist=get_counts(refdata)
    print "total agonists %s" % total_agonist
    print "total antagonists %s" % total_antagonist
    #types=[x for x in set(types)]
    paths=open('%s_paths.txt' % sys)
    output=True
    for (n, path) in enumerate(paths.readlines()):
        print "--------path %s---------" % n
        antagonist=[]
        agonist=[]
        for (m, state) in enumerate(path.split()):
            types=get_types(refdata)
            ohandle=open('path%s_state%s_top%ssummary.txt' % (n, m, percent*100), 'w')
            print state
            file='ranked-state%s-docking-scores.txt' % state
            names=numpy.loadtxt(file, usecols=(0,), dtype=str)
            scores=numpy.loadtxt(file, usecols=(1,))    
            if output==True:
                print "top %s%% %s molecules" % (str(percent*100), int(len(names)*percent))
                output=False
            names=names[:int(len(names)*percent)]
            scores=scores[:int(len(names)*percent)]
            for x in names:
                if 'ICI' in x:
                    if '118' not in x:
                        x='ICI89406'
                if 'CL' in x:
                    x='CL316243'
                if '_' in x:
                    x=x.strip('_').replace('_', '-')
                if 'omega-92' in x:
                    continue
                if 'mlawrenz' not in refdata[x]:
                    types[refdata[x]]+=1
            ag=0
            ant=0
            for x in types:
                if 'antagonist' in x:
                    ant+=types[x]
                else:   
                    ag+=types[x]
                ohandle.write('%s\t%s\n' % (x, types[x]))
            #print types[x]/counts[x]
            agonist.append(float(ag)/float(total_agonist))
            antagonist.append(float(ant)/float(total_antagonist))
            #print "agonist %% recovered %0.2f" % (float(ag)/float(total_agonist))
            #print "antagonist %% recovered %0.2f" % (float(ant)/float(total_antagonist))    
        pylab.figure()
        pylab.scatter(antagonist, agonist, label='Path %s' % n)
        pylab.legend()
        pylab.xlabel('antagonist recovery')
        pylab.ylabel('agonist recovery')
        pylab.show()

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


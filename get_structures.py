import numpy, pickle
file='gen-dat/rocs-bw_car-595-22337144-4_combo_0.4_g3265.dat'
dbase=numpy.loadtxt('all-top-agonist-dbase.list', dtype=str)
sys=file.split('bw_')[1].split('-')[0]
state=file.split('bw_')[1].split('-')[1]
name=file.split('bw_')[1].split('%s-%s-' % (sys, state))[1].split('_combo')[0]
target='/nobackup/mlawrenz/cluster-agonist-dbase/%s-%s/stereo-%s_000.mol2' % (sys, state, name)
index=numpy.where(dbase==target)[0]
targindex=index[0]+1

fhandle=open(file)
hitfile=open('test-hits.mol2', 'w')
for line in fhandle.readlines():
    sys=line.split('-')[0]
    state=line.split('-')[1]
    name=line.split('%s-%s-' % (sys, state))[1].split()[0]
    target='/nobackup/mlawrenz/cluster-agonist-dbase/%s-%s/stereo-%s_000.mol2' % (sys, state, name)
    index=numpy.where(dbase==target)[0]
    index=index[0]+1
    target='%s-%s/stereo-%s_000.mol2' % (sys, state, name)
    output='%s-%s/rocs-bw-%s-%s-%s.out' % (sys, state, sys, state, targindex)
    ohandle=open(output)
    write=False
    for line in ohandle.readlines():
        if 'Database' in line:
            if target in line.split()[-1].split('agonist/')[-1]:
                write=True
                print "writing %s" % target
                continue
            else:
                write=False
        elif write==True:
            stereo=name.rstrip()[-1]
            if name.split('-%s' % stereo)[0] in line:
                hitfile.write('%s-%s-%s\n' % (sys, state, name))
            else:
                hitfile.write(line)
        elif write==False:
            pass


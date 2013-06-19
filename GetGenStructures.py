import numpy, pickle
import optparse
import os
import glob


def main(database, prefix): 
    dir=os.path.dirname(database)
    dbase=numpy.loadtxt(database, dtype=str) 
    files=glob.glob('%s/gen-files/%s_*_combo*_g*dat' % (dir, prefix))
    for file in files:
        sys=file.split('%s_' % prefix)[1].split('-')[0]
        state=file.split('%s_' % prefix)[1].split('-')[1]
        name=file.split('%s_' % prefix)[1].split('%s-%s-' % (sys, state))[1].split('_combo')[0]
        target='%s-%s/stereo-%s_000.mol2' % (sys, state, name)
        index=numpy.where(dbase==target)[0]
        targindex=index[0]+1
        fhandle=open(file)
        hitfile=open('%s/gen-files/%s-%s-%s-hits.mol2' % (dir, sys, state, name), 'w')
        for line in fhandle.readlines():
            sys=line.split('-')[0]
            state=line.split('-')[1]
            name=line.split('%s-%s-' % (sys, state))[1].split()[0]
            if sys=='apo':
                target='%s-recent-%s/stereo-%s_000.mol2' % (sys, state, name)
            else:
                target='%s-%s/stereo-%s_000.mol2' % (sys, state, name)
            #target='/nobackup/mlawrenz/cluster-agonist-dbase/%s-%s/stereo-%s_000.mol2' % (sys, state, name)
            #index=numpy.where(dbase==target)[0]
            #index=index[0]+1
            output='%s/%s-%s/rocs-bw-%s-%s-%s.out' % (dir, sys, state, sys, state, targindex)
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

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dbase', dest='dbase',
                      help='database list')
    parser.add_option('-p', '--prefix', dest='prefix',
                      help='prefix name for rpt files')
    parser.add_option('-w', action="store_true", dest="writepdb")
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(database=options.dbase, prefix=options.prefix)


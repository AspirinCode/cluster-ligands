#!/bin/python
import numpy
import glob
import optparse
import operator
import random
import os
import operator


def main(file):
    fhandle=open(file)
    stereo=dict()
    test=dict()
    for line in fhandle.readlines():
        name=line.split(':')[1].split('_000')[0].split('-')[1]
        full=line.split(':')[1].split('_000')[0]
        score=float(line.split()[1])
        if name not in test.keys():
            stereo[full]=score
            test[name]=score
        elif name in test.keys():
            if test[name]<score:
                test[name]=score
                stereo[full]=score
            else:
                pass
    print "final %s ligands" % (len(stereo.keys()))
    sorted_stereo=reversed(sorted(stereo.iteritems(), key=operator.itemgetter(1)))
    ofile=open('reduce-%s' % file, 'w')
    for x in sorted_stereo:
        os.system('sed "s/FILE/%s/g" < mol2pdb.tcl > tmp.tcl' % x[0])
        os.system('vmd -dispdev text -e tmp.tcl')
        print x[0], x[1]
        ofile.write('%s\t%s\n' % (x[0], x[1]))
    ofile.close()

        
def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-f', '--file', dest='file',
                      help='file')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    main(file=options.file)

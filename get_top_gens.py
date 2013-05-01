import glob
import os
import numpy

def main(gens, dbase, ass):
    #gen_inds=[249, 95, 231, 161, 195]
    dir=gens.split('gens')[0]
    bindb=numpy.loadtxt('bindb_combo_1.5_gen_names.dat', dtype=str)
    gdd=numpy.loadtxt('gdd_combo_1.5_gen_names.dat', dtype=str)
    gdd_ind=numpy.loadtxt('gdd_combo_1.5_gen_indices.dat', dtype=int)
    bindb_ind=numpy.loadtxt('bindb_combo_1.5_gen_indices.dat', dtype=int)
    assignments=numpy.loadtxt(ass)
    for x in assignments:
        x=int(x)
        if x!=-1:
            if x not in pops.keys():
                pops[x]=1
            else:
                pops[x]+=1
    gen_indices=sorted(pops.keys())
    pylab.figure()
    pylab.plot(sorted(pops.keys()), [pops[key] for key in sorted(pops.keys())])
    pylab.xlabel('Gens Index')
    pylab.ylabel('Counts')
    pylab.savefig('%s/histogram_gens.png' % dir, dpi=300)
    pylab.show()
    gens[gen_inds]
    names=gens[gen_inds]
    category=dict()
    for (n, x) in enumerate(names):
        x=os.path.basename(x)
        if '.mol2' in x:
            x=x.split('stereo-')[1].split('-')[0]
        location=numpy.where(bindb==x)[0]
        if not location.size:
            location=numpy.where(gdd==x)[0]
            if not location.size:
                print "problem for ", x
            else:
                category[n]='gdd%s' % x # + str(gdd_ind[location][0])
        else:
            x=x.split('molecule_')[1].split('.pdb')[0]
            category[n]='bindb%s' % x# + str(bindb_ind[location][0])
    
    for key in sorted(category.keys()):
        value=category[key]
        print key, value
        if 'bindb' in  value:
            os.system('grep %s ../../bindingdb/oe-output/*1.5*dat | grep -v name' % value.split('bindb')[1])
        else:
            os.system('grep %s ../gdd_combo_1.5_g*dat | grep -v name' % value.split('gdd')[1])

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-g', '--gens', dest='gens',
                      help='gens list')
    parser.add_option('-d', '--dbase', dest='dbase',
                      help='database list')
    parser.add_option('-a', '--assignments', dest='assignments',
                      help='assignments file')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    assignments=options.assignments
    gens=options.gens
    dbase=options.dbase


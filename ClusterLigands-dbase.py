from msmbuilder import metrics
import os
import sys
import optparse
import pickle
import numpy
from msmbuilder import Trajectory

def format_stereo(list):
    formatted=[]
    for x in list:
        if 'stereo' in x:
            base=os.path.basename(x)
            formatted.append(base.split('.mol2')[0].split('stereo-')[1])
        elif 'ZINC' in x:
            base=os.path.basename(x)
            formatted.append(base.split('.mol2')[0].split('ZINC')[1])
        elif 'molecule' in x:
            if 'omega' in x:
                base=os.path.basename(x)
                formatted.append(base.split('_92.pdb')[0].split('molecule_')[1])
            else:
                base=os.path.basename(x)
                formatted.append(base.split('.pdb')[0].split('molecule_')[1])
        else:
            print "parser does not recognize molecule file name"
    return formatted

def format_names(list):
    formatted=[]
    for x in list:
        if '_000' in x:
            if 'ZINC' in x:
                base=os.path.basename(x)
                formatted.append(base.split('_000.mol2')[0].split('new-')[1])
            elif 'stereo' in x:
                base=os.path.basename(x)
                formatted.append(base.split('_000.mol2')[0].split('stereo-')[1].split('-')[0])
        elif 'stereo' in x:
            base=os.path.basename(x)
            formatted.append(base.split('.mol2')[0].split('stereo-')[1].split('.')[0])
        elif 'ZINC' in x:
            base=os.path.basename(x)
            formatted.append(base.split('.mol2')[0].split('ZINC')[1])
        elif 'molecule' in x:
            if 'omega' in x:
                base=os.path.basename(x)
                formatted.append(base.split('_92.pdb')[0].split('molecule_')[1])
            else:
                base=os.path.basename(x)
                formatted.append(base.split('.pdb')[0].split('molecule_')[1])
    return numpy.array(formatted)

# make a metric, where you have a grid of all distances to all
# call evaluation by matrix order
def parse(dir, file, output, patterns):
    patterns=list(patterns)
    fhandle=open(file)
    write=False
    for line in fhandle.readlines():
        if 'CRYS' in line:
            pass
        elif 'query' in line:
            pass
        elif 'COMPND' in line:
            # below covers the first round
            pattern=line.split()[1]
            if pattern in patterns:
                write=True
                print "writing %s" % pattern
                new=open('%s/%s_%s.pdb' % (dir, output, pattern), 'w')
                new.write(line)
                index=patterns.index(pattern)
                patterns.pop(index)
            else:
                write=False
        elif 'CONECT' in line:
            pass
        elif write==True:
		    new.write(line)
        elif write==False:
            pass
    return


def get_matrix(dir, database, reference, column, max, add=False, prefix='eon-bindb'):
    list=range(0, len(database))
    matrix=numpy.zeros((len(database), len(reference)))
    rankmatrix=numpy.zeros((len(database), len(reference)), dtype=int)
    # loop over database
    for (n,element) in enumerate(list):
        m=element+1
        file='%s/%s-%s.rpt' % (dir, prefix, m)
        name=os.path.dirname(file)+'/mod-'+os.path.basename(file)
        os.system('sed "1d" < %s | sed "s/ICI 89406/ICI89406/g" > %s' % (file, name))
        file=open(name)
        score=-100*numpy.ones((len(reference)))
        indices=-100*numpy.ones((len(reference)))
        names=[]
        k=0
        index1=column[0]
        index2=column[1]
        for line in file.readlines():
            if 'Rank' not in line.split():
                # looping over files for each member of dbase, matched with gens
                name=str(line.split()[0])
                location=numpy.where(reference==name)[0]
                if location.size:
                    names.append(name)
                    indices[k]=location
                    if add==True:
                        score[k]=(max-(float(line.split()[index1])+float(line.split()[index2])))
                        if score[k]<0:
                            import pdb
                            pdb.set_trace()
                        k+=1
                    else:
                        score[k]=(max-float(line.split()[index1]))
                        k+=1
                else:
                    print "no location for %s", name
                    pass
        order=numpy.argsort(indices)
        names=numpy.array(names)
        for i in range(0, len(score)):
            matrix[n,i]=score[order][i]
            rankmatrix[n,i]=indices[i]
    return rankmatrix, matrix

def assign(matrix, database, cutoff):
    distances=-1*numpy.ones(len(database))
    assignments=-1*numpy.ones(len(database), dtype=int)
    for j in xrange(matrix.shape[0]):
        d=matrix[j,:]
        ind=numpy.argmin(d)
        if not d[ind] < cutoff:
            pass
        else:
            assignments[j] = int(numpy.argmin(d))
            distances[j] = d[assignments[j]]
    return assignments, distances
 
def cluster(distance_cutoff, matrix, database):
    distance_list = numpy.inf * numpy.ones(len(database), dtype=numpy.float32)
    assignments = -1 * numpy.ones(len(database), dtype=int)
    distance_cutoff=float(distance_cutoff)
    # set k to be the highest 32bit integer
    k = sys.maxint

    seed=0
    generator_indices = []    
        
    for i in xrange(k):
        new_ind = seed if i == 0 else numpy.argmax(distance_list)
        print "K-centers: Finding generator %i. Will finish when % .4f drops below % .4f" % (i, float(distance_list[new_ind]), float(distance_cutoff))
        if distance_list[new_ind] < distance_cutoff:
            break
        new_distance_list = matrix[new_ind, :]
        updated_indices = numpy.where(new_distance_list < distance_list)[0]
        distance_list[updated_indices] = new_distance_list[updated_indices]
        assignments[updated_indices] = new_ind
        generator_indices.append(new_ind)
    return numpy.array(generator_indices), numpy.array(assignments), numpy.array(distance_list)

def parse_commandline():
    parser = optparse.OptionParser()
    parser.add_option('-d', '--dbase', dest='dbase',
                      help='database list')
    parser.add_option('-p', '--prefix', dest='prefix',
                      help='prefix name for rpt files')
    parser.add_option('-c', '--cutoff', dest='cutoff',
                      help='tanimoto cutoff')
    parser.add_option('-w', action="store_true", dest="writepdb")
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    cutoff=float(options.cutoff)
    prefix=options.prefix
    dbase=options.dbase
    dir=os.path.dirname(dbase)
    max=2
    print "using combo PB and shape Tamimoto score"
    column=[3,6]
    add=True
    cutoff=max-cutoff
    database=numpy.loadtxt(dbase, dtype=str)
    format_dbase=format_names(database)
    # here pass in reference as the database
    rankmatrix, matrix=get_matrix(dir, format_dbase, format_dbase, column, max, add,  prefix=prefix)
    gens, assignments, distances=cluster(cutoff, matrix, format_dbase)
    frames=numpy.where(assignments!=-1)[0]
    distances[frames]=[(max-i) for i in distances[frames]]
    numpy.savetxt('%s-assignments.dat' % dbase.split('.list')[0], assignments)
    numpy.savetxt('%s-distances.dat' % dbase.split('.list')[0], distances)
    print "%s assigned to gens" % len(frames)
    for i in gens:
        frames=numpy.where(assignments==i)[0]
        ohandle=open('%s/%s_%s_%s_%s_g%s.dat' % (dir, prefix, format_dbase[i], tanimoto, (max-cutoff), i), 'w')
        for (name, val) in zip(format_dbase[frames], distances[frames]):
            ohandle.write('%s\t%s\n' % (name, (max-val)))
        if options.writepdb==True:
            j=i+1
            hitfile='%s/%s-%s_hits.pdb' % (dir, prefix, j)
            output='g%s_%s' % (i, max-cutoff)
            parse(dir, hitfile, output, format_dbase[frames])
    numpy.savetxt('%s/%s_%s_%s_gen_indices.dat' % (dir, prefix, tanimoto, (max-cutoff)), gens, fmt='%i')
    numpy.savetxt('%s/%s_%s_%s_gen_names.dat' % (dir, prefix, tanimoto, (max-cutoff)), format_dbase[gens], fmt='%i')


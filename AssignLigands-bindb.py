from msmbuilder import metrics
import os
import sys
import optparse
import pickle
import numpy
from msmbuilder import Trajectory

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


def get_matrix(dir, database, reference, column, max, add=False, frames=None, prefix='eon-bindb'):
    if frames is not None:
        list=frames
        reference=reference[frames]
        matrix=numpy.zeros((len(reference), len(database)))
        rankmatrix=numpy.zeros((len(reference), len(database)), dtype=int)
    else:
        list=range(0, len(reference))
        matrix=numpy.zeros((len(reference), len(database)))
        rankmatrix=numpy.zeros((len(reference), len(database)), dtype=int)
    for (n,element) in enumerate(list):
        m=element+1
        file=open('%s/%s%s.rpt' % (dir, prefix, m))
        score=-100*numpy.ones((len(database)))
        indices=-100*numpy.ones((len(database)))
        names=numpy.zeros((len(database)), dtype=str)
        k=0
        index1=column[0]
        index2=column[1]
        for line in file.readlines():
            if 'Rank' not in line.split():
                # look at GDD molecule name, score, index for BINDB molecule
                name=str(line.split()[0])
                location=numpy.where(database==name)[0]
                if location.size:
                    names[k]=name
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
                    pass
        order=numpy.argsort(indices)
        for i in range(0, len(score)):
            matrix[n,i]=score[order][i]
            rankmatrix[n,i]=indices[i]
    return rankmatrix, matrix

def assign(matrix, database, cutoff):
    distances=-1*numpy.ones(len(database))
    assignments=-1*numpy.ones(len(database), dtype=int)
    for j in xrange(len(distances)):
        d=matrix[:,j]
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
    parser.add_option('-g', '--genindices', dest='genindices',
                      help='gens list')
    parser.add_option('-d', '--dbase', dest='dbase',
                      help='database list')
    parser.add_option('-c', '--cutoff', dest='cutoff',
                      help='tanimoto cutoff')
    parser.add_option('-w', action="store_true", dest="writepdb")
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    cutoff=float(options.cutoff)
    dbase=options.dbase
    dir=os.path.dirname(dbase)
    genindices=options.genindices
    max=2
    prefix=os.path.basename(genindices).split('_combo')[0]
    print "using combo PB and shape Tamimoto score"
    column=[3,6]
    add=True
    cutoff=max-cutoff
    formatted=[]
    database=numpy.loadtxt(dbase, dtype=str)
    for x in database:
        if 'ZINC' in x:
            formatted.append(x.split('_000.mol2')[0])
        else:
            formatted.append(x.split('_000.mol2')[0].split('stereo-')[1].split('-')[0])
    database=numpy.array(formatted)
    genindices=numpy.loadtxt(genindices, dtype=int)
    reference=numpy.loadtxt('/home/mlawrenz/wontkill/new-results/bindingdb/dbase.list', dtype=str)
    rankmatrix, matrix=get_matrix(dir, database, reference, column, max, add, frames=genindices, prefix=prefix)
    assignments, distances=assign(matrix, database, cutoff)
    #save the above assignments, distances
    database_frames=numpy.where(assignments==-1)[0]
    rankmatrix, matrix=get_matrix(dir, database[database_frames], database, column, max, add, frames=database_frames, prefix='eon-only-molecule')
    newgenindices, newassignments, newdistances=cluster(cutoff, matrix, database_frames)
    newgenindices=database_frames[newgenindices]
    numpy.savetxt('%s/gdd_combo_%s_gen_indices.dat' % (dir, (max-cutoff)), newgenindices, fmt='%i')
    data_gens=database[newgenindices]
    numpy.savetxt('%s/gdd_combo_%s_gen_names.dat' % (dir, (max-cutoff)), data_gens, fmt='%s')
    numpy.savetxt('%s/bindb_combo_%s_gen_indices.dat' % (dir, (max-cutoff)), genindices, fmt='%i')
    bindb_gens=reference[genindices]
    numpy.savetxt('%s/bindb_combo_%s_gen_names.dat' % (dir, (max-cutoff)), bindb_gens, fmt='%s')
    prefix='bindb'
    n=1
    for i in genindices:
        frames=numpy.where(assignments==i)[0]
        if frames.size:
            ohandle=open('%s/%s_combo_%s_g%s.dat' % (dir, prefix, (max-cutoff), n), 'w')
            for (name, val) in zip(database[frames], distances[frames]):
                ohandle.write('%s\t%s\n' % (name, (max-val)))
            if options.writepdb==True:
                j=i+1
                file='%s/eon-%s%s_hits.pdb' % (dir, j)
                output='g%s_%s' % (n, max-cutoff)
                parse(dir, file, output, database[frames])
            n+=1
    prefix='gdd'
    for i in newgenindices:
        frames=numpy.where(newassignments==i)[0]
        if frames.size:
            ohandle=open('%s/%s_combo_%s_g%s.dat' % (dir, prefix,  (max-cutoff), n), 'w')
            for (name, val) in zip(database[frames], newdistances[frames]):
                ohandle.write('%s\t%s\n' % (name, (max-val)))
            if options.writepdb==True:
                k=i+1
                file='%s/eon-only-molecule%s_hits.pdb' % (dir, k)
                output='g%s_%s' % (n, max-cutoff)
                parse(dir, file, output, database[frames])
            n+=1
    print "%s gens total done" % n


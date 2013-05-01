from msmbuilder import metrics
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


def get_matrix(dir, prefix, database, column, max, add=False):
    matrix=numpy.zeros((len(database), len(database)))
    rankmatrix=numpy.zeros((len(database), len(database)), dtype=int)
    for n in range(0, len(database)):
        m=n+1
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
        order=numpy.argsort(indices)
        for i in range(0, len(score)):
            matrix[n,i]=score[order][i]
            rankmatrix[n,i]=indices[i]
    return rankmatrix, matrix

def assign(matrix):
    distances=numpy.zeros(len(database))
    assignments=numpy.zeros(len(database), dtype=int)
    for j in xrange(len(distances)):
        d=matrix[:,j]
        assignments[j] = int(numpy.argmin(d))
        distances[j] = d[assignments[j]]
    return assignments, distances
 
def cluster(distance_cutoff, matrix, database, gens=None):
    if gens!=None:
        assignments=numpy.arange(0, len(reference))
        distance_list = numpy.inf * numpy.ones(len(reference), dtype=numpy.float32)
    else:
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
    parser.add_option('-d', '--dir', dest='dir',
                      help='directory')
    parser.add_option('-p', '--prefix', dest='prefix',
                      help='prefix to rocs/eon output - note that need to read dbase.list file of all molecule pdb names')
    parser.add_option('-t', '--tanimoto', dest='tanimoto',
                      help='tanimoto types: pb, shape, combo (pb + shape)')
    parser.add_option('-c', '--cutoff', dest='cutoff',
                      help='tanimoto cutoff')
    parser.add_option('-w', action="store_true", dest="writepdb")
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    dir=options.dir
    tanimoto=options.tanimoto
    cutoff=float(options.cutoff)
    prefix=options.prefix
    if tanimoto=='pb':
        print "using PB Tanimoto score"
        column=3
        max=1
        add=False
    elif tanimoto=='shape':
        print "using shape Tamimoto score"
        column=6
        max=1
        add=False
    elif tanimoto=='combo':
        print "using combo PB and shape Tamimoto score"
        column=[3,6]
        max=2
        add=True
    cutoff=max-cutoff
    database=[x.split('molecule_')[1].split('.pdb')[0] for x in numpy.loadtxt('%s/dbase.list' % dir, dtype=str)]
    database=numpy.array(database)
    rankmatrix, matrix=get_matrix(dir, prefix, database, column, max, add)
    gens, assignments, distances=cluster(cutoff, matrix, database)
    for i in gens:
        frames=numpy.where(assignments==i)[0]
        ohandle=open('%s/%s_%s_%s_%s_g%s.dat' % (dir, prefix, database[i], tanimoto, (max-cutoff), i), 'w')
        for (name, val) in zip(database[frames], distances[frames]):
            ohandle.write('%s\t%s\n' % (name, (max-val)))
        if options.writepdb==True:
            j=i+1
            file='%s/mod-eon-bindb%s_hits.pdb' % (dir, j)
            output='g%s_%s' % (i, max-cutoff)
            parse(dir, file, output, database[frames])
    numpy.savetxt('%s/%s_%s_%s_gens.dat' % (dir, prefix, tanimoto, (max-cutoff)), gens, fmt='%i')
    #assignments, distances=assign(matrix)
    #for (n, i) in enumerate(database):
    #    name=i.split('.pdb')[0].split('_')[1]
    #    frames=numpy.where(assignments==n)[0]
    #    if len(frames!=0):
    #        ohandle=open('%s/%s_%s_%s_%s_g%s.dat' % (dir, prefix, tanimoto, (max-cutoff), name, n), 'w')
    #        for (name, val) in zip(database[frames], distances[frames]):
    #            ohandle.write('%s\t%s\n' % (name, (max-val)))
    #numpy.savetxt('%s/%s_%s_%s_scores.dat' % (dir, prefix,  tanimoto, (max-cutoff)), [1-d for d in distances])
    #numpy.savetxt('%s/%s_%s_%s_gens.dat' % (dir, prefix, tanimoto, (max-cutoff)), gens, fmt='%i')
    #numpy.savetxt('%s/%s_%s_%s_assignments.dat' % (dir, prefix, tanimoto, (max-cutoff)), assignments, fmt='%i')
    print "done"


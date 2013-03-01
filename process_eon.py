from msmbuilder import metrics
import sys
import optparse
import pickle
import numpy
from msmbuilder import Trajectory

# make a metric, where you have a grid of all distances to all
# call evaluation by matrix order
def parse(dir, file, output, pdbframes):
    fhandle=open(file)
    count=0
    pattern=None
    for line in fhandle.readlines():
        if 'CRYS' in line:
            pass
        elif 'query' in line:
            pass
        elif 'COMPND' in line:
            # below covers the first round
            if pattern==None and count in pdbframes:
                pattern=line.split()[1]
                print "writing dbase.list molecule %s" % count
                print pattern
                new=open('%s/%s_%s.pdb' % (dir, output, pattern), 'w')
                new.write(line)
            elif pattern==None and count not in pdbframes:
                pattern=line.split()[1]
                #don't increment bc this stays frame 0
                pass
            elif pattern!=None:
                count+=1
                if count in pdbframes:
                    pattern=line.split()[1]
                    print "writing dbase.list molecule %s" % count
                    print pattern
                    new=open('%s/%s_%s.pdb' % (dir, output, pattern), 'w')
                    new.write(line)
                else:
                    pass
        elif pattern!=None and count in pdbframes:
		    new.write(line)
        elif pattern!=None and count not in pdbframes:
            pass
        elif pattern==None:
            pass
    return
                

def get_matrix(dir, database):
    matrix=numpy.zeros((len(database), len(database)))
    rankmatrix=numpy.zeros((len(database), len(database)), dtype=int)
    for n in range(0, len(database)):
        m=n+1
        file=open('%s/eon-molecule%s.rpt' % (dir, m))
        pb_T=[]
        #comb_T=[]
        #shape_T=[]
        names=[]
        for line in file.readlines():
            if 'Rank' not in line.split():
                names.append(line.split()[0].split('-')[2])
                pb_T.append(1-float(line.split()[3]))
                #comb_T.append(float(line.split()[5]))
                #shape_T.append(float(line.split()[6]))
        names=numpy.array(names, dtype=int)
        pb_T=numpy.array(pb_T)
        order=numpy.argsort(names)
        for i in range(0, len(pb_T)):
            matrix[n,i]=pb_T[order][i]
            rankmatrix[n,i]=names[i]
    return rankmatrix, matrix

def cluster(distance_cutoff, matrix, database):
    distance_cutoff=float(distance_cutoff)
    assignments=numpy.arange(0, len(database))
    # set k to be the highest 32bit integer
    k = sys.maxint

    seed=0
    distance_list = numpy.inf * numpy.ones(len(database), dtype=numpy.float32)
    assignments = -1 * numpy.ones(len(database), dtype=numpy.int32)
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
    parser.add_option('-t', '--tanimotoc', dest='tanimotoc',
                      help='tanimoto cutoff (inverse)')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    dir=options.dir
    tanimotoc=options.tanimotoc
    database=numpy.loadtxt('%s/dbase.list' % dir, dtype=str)
    rankmatrix, matrix=get_matrix(dir, database)
    gens, assignments, distances=cluster(tanimotoc, matrix, database)
    #for i in gens:
    for i in range(0, 1):
        frames=numpy.where(assignments==i)[0]
        file=('%s/eon-molecule%s_hits.pdb' % (dir, (i+1))) 
        pdbframes=numpy.zeros((len(frames)))
        scores=numpy.zeros((len(frames)))
        for (z, j) in enumerate(frames):
            k=j+1
            location=numpy.where(rankmatrix[i]==k)[0]
            pdbframes[z]=location
            scores[z]=1-matrix[i, j]
        order=numpy.argsort(scores)[::-1]
        print "g%s dbase.list frames: " % i, pdbframes[order]
        print "g%s dbase.list molecules: " % i, [i+1 for i in frames[order]]
        print "g%s dbase.list scores: " % i, scores[order]
        parse(dir, file, output='g%s' % i, pdbframes=pdbframes)
    numpy.savetxt('%s/molecule_gens.dat' % dir, database[gens], fmt='%s')
    numpy.savetxt('%s/gens.dat' % dir, gens, fmt='%i')
    numpy.savetxt('%s/assignments.dat' % dir, assignments, fmt='%i')
    numpy.savetxt('%s/distances.dat' % dir, distances)
    print "done"


from msmbuilder import metrics
import sys
import optparse
import pickle
import numpy
from msmbuilder import Trajectory

# make a metric, where you have a grid of all distances to all
# call evaluation by matrix order
def parse(file, output, climit):
    fhandle=open(file)
    pattern=None
    count=0
    pattern=None
    for line in fhandle.readlines():
        if 'COMPND' in line:
            if count==climit:
                break
            elif pattern==None:
                pattern=line.split()[1]
                print pattern
                new=open('%s_%s.pdb' % (output, pattern), 'w')
                new.write(line)
            else:
                newpattern=line.split()[1]
                if newpattern==pattern:
                    new.write(line)
                else:
                    pattern=line.split()[1]
                    print pattern
                    count+=1
                    new=open('%s_%s.pdb' % (output, pattern), 'w')
                    new.write(line)
        else:
		    new.write(line)
                


def get_matrix(database):
    matrix=numpy.zeros((len(database), len(database)))
    for n in range(0, len(database)):
        m=n+1
        file=open('eon-molecule%s.rpt' % m)
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
    return matrix

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
    parser.add_option('-d', '--distance', dest='distance',
                      help='cluster distance')
    (options, args) = parser.parse_args()
    return (options, args)

if __name__ == "__main__":
    (options, args) = parse_commandline()
    database=numpy.loadtxt('dbase.list', dtype=str)
    matrix=get_matrix(database)
    gens, assignments, distances=cluster(options.distance, matrix, database)
    for i in gens:
        frames=numpy.where(assignments==i)[0]
        file=('eon-molecule%s_hits.pdb' % i) 
    numpy.savetxt('molecule_gens.dat', database[gens], fmt='%s')
    numpy.savetxt('gens.dat', gens, fmt='%i')
    numpy.savetxt('assignments.dat', assignments, fmt='%i')
    numpy.savetxt('distances.dat', distances)


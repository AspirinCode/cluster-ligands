from msmbuilder import metrics
import sys
import glob
import os
import sys
import optparse
import pickle
import numpy
from msmbuilder import Trajectory

def format_top(list):
    formatted=[]
    namesonly=[]
    for x in list:
        state=x.split('/')[0]
        name=x.split('stereo-')[1].split('-')[0]
        num=x.split('stereo-')[1].split('-')[1].split('_000')[0]
        if 'ICI-' in x.split()[0] or 'AH3474-' in x.split()[0]:
            name=x.split('stereo-')[1].split('-')[0]+x.split('stereo-')[1].split('-')[1]
            num=x.split('stereo-')[1].split('-')[2].split('_000.mol2')[0]
            formatted.append('%s-%s-%s' % (state, name, num))
            namesonly.append(name)
        else:
            formatted.append('%s-%s-%s' % (state, name, num))
            namesonly.append(name)
    return numpy.array(namesonly), numpy.array(formatted)

def format_stereo(list):
    formatted=[]
    for x in list:
        if 'stereo' in x:
            base=os.path.basename(x)
            tmp=base.split('.mol2')[0].split('stereo-')[1]
            if len(tmp.split('.'))>1:
                formatted.append('%s-%s' % (tmp.split('.')[0], tmp.split('.')[1]))
            else:
                formatted.append(tmp)
        elif 'ZINC' in x:
            base=os.path.basename(x)
            formatted.append(base.split('.mol2')[0].split('ZINC')[1])
        elif 'molecule' in x:
            if 'omega' in x:
                base=os.path.basename(x)
                formatted.append(base.split('_92.pdb')[0].split('molecule_')[1])
            else:
                base=os.path.basename(x)
                formatted.append(base.split('.mol2')[0].split('molecule_')[1].replace('_',
'-'))
        else:
            print "parser does not recognize molecule file name"
    return numpy.array(formatted)

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

def checkpoint(dir, prefix, n, matrix):
    if os.path.exists('%s/%s-matrix-chkpt.pickle' % (dir, prefix)):
        os.system('rm %s/%s-matrix-chkpt.pickle' % (dir, prefix))
    pfile=open('%s/%s-matrix-chkpt.pickle' % (dir, prefix), 'wb')
    pickle.dump(matrix, pfile)
    pfile.close()
    numpy.savetxt('%s/chkpt-dbaseindex.dat' % dir, [n,], fmt='%i')
    return

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


def get_matrix(dir, database, reference, column, max, add=False, start=0, restart=None, prefix='eon-bindb'):
    errfile=open('errors.txt', 'w')
    list=range(0, len(database))
    dirs=numpy.loadtxt('%s/directory-list.txt' % dir, dtype=str)
    if restart==None:
        matrix=numpy.zeros((len(database), len(reference)))
    else:
        matrix=restart
    # for each n, want to load each dir into n column of matrix 
    for n in range(start, len(database)):
        index=n+1
        item=database[n]
        for subdir in dirs:
            format_names, format_sub=format_top(numpy.loadtxt(subdir, dtype=str))
            subdir=subdir.rstrip('-top-dbase.list').split('/')[-1]
            file='%s/%s/%s-%s-%s_1.rpt' % (dir, subdir, prefix, subdir, index)
            print "on file %s" % file
            fhandle=open(file)
            score=-100*numpy.ones((len(format_names)))
            indices=-100*numpy.ones((len(format_names)), dtype=int)
            names=[]
            k=0
            if len(column) >1:
                index1=column[0]
                index2=column[1]
            else:
                index1=column[0]
            checkfile=False
            for line in fhandle.readlines():
                if 'Rank' in line.split():
                    continue
                # looping over files for each member of dbase, matched with gens
                if checkfile==False:
                    refname=str(line.split()[1])
                    location=numpy.where(reference==refname)[0]
                    checkfile=True
                if 'ICI-' in line.split()[0] or 'AH3474-' in line.split()[0]:
                    testname=line.split()[0].split('-')[0]+line.split()[0].split('-')[1]
                else:
                    testname=line.split()[0]
                location=numpy.where(format_names==testname)
                stereo=format_sub[location]
                location=numpy.where(reference==stereo)[0]
                if location.size:
                    check=numpy.where(indices==location)[0]
                    if check.size:
                        if add==True:
                            sum=(float(line.split()[index1])+float(line.split()[index2]))
                        else:
                            sum=float(line.split()[index1])
                        s=(max-sum)
                        try: 
                            s==score[check]
                            continue
                        except ValueError:
                            print "same molecule %s exists, different score in %s" % (name, file)
                    names.append(stereo)
                    indices[k]=int(location)
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
                    errfile.write("no location for %s\n" % stereo)
                    pass
            for i in range(0, len(reference)):
                location=numpy.where(indices==i)[0]
                if not location.size:
                    errfile.write("no report for %s in %s\n" % (reference[i], file))
            order=numpy.argsort(indices)
            names=numpy.array(names)
            if -100 in score:
                print file
                print "not all reports found"
                sys.exit(0)
            for (i,j) in zip(indices, score): #locations in the dbase file
                matrix[n,i]=j
            if index % 1000==0:
                checkpoint(dir, prefix, n, matrix)
    return matrix

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

def clean_matrix(matrix):
    for x in range(0, matrix.shape[0]):
        for y in range(0, matrix.shape[1]):
            if matrix[x,y]!=matrix[y,x]:
                maxval=numpy.max([matrix[x,y], matrix[y,x]])
                matrix[x,y]=maxval
                matrix[y,x]=maxval
            else:
                pass
    return matrix
        

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
    tanimoto='combo'
    if 'eon' in prefix:
        print "using combo PB and shape Tamimoto score"
        column=[3,6]
        add=True
    else:
        print "using combo ROCS color and shape Tamimoto score"
        column=[3,]
        add=False
    cutoff=max-cutoff
    database=numpy.loadtxt(dbase, dtype=str)
    format_names, format_dbase=format_top(database)
    if not os.path.exists('%s/directory-list.txt' % dir):
        dirs=glob.glob('%s/*-top-dbase.list' % dir)
        numpy.savetxt('%s/directory-list.txt' % dir, dirs, fmt='%s')
    numpy.savetxt('%s/formatted_dbase.list' % dir, format_dbase, fmt='%s')
    # here pass in reference as the database
    if os.path.exists('%s-assignments.dat' % dbase.split('.list')[0]):
        print "assignments exist in %s cutoff %s" % (dir, cutoff)
        assignments=numpy.loadtxt('%s-assignments.dat' % dbase.split('.list')[0])
        distances=numpy.loadtxt('%s-distances.dat' % dbase.split('.list')[0])
        gens=numpy.loadtxt('%s/%s_%s_%s_gen_indices.dat' % (dir, prefix, tanimoto, (max-cutoff)))
        if options.writepdb==True:
            print "writing structure files"
            for i in gens:
                frames=numpy.where(assignments==i)[0]
                name=format_dbase[i]
                sys=name.split('-')[0]
                state=name.split('-')[1]
                molecule=name.split('-')[2]
                stereo=name.split('-')[3]
                subdir='%s-%s' % (sys, state)
                if 'eon' in prefix:
                    hitfile='%s/%s-%s_hits.pdb' % (dir, prefix, name)
                else:
                    hitfile='%s/%s/%s-%s-%i.out' % (dir, subdir, prefix, subdir, i)
                output='g%s_%s' % (i, max-cutoff)
                parse(dir, hitfile, output, format_dbase[frames])
    else:
        if os.path.exists('%s/%s-matrix-chkpt.pickle' % (dir, prefix)):
            print "loading checkpoint score matrix from %s" % dir
            pfile=open('%s/%s-matrix-chkpt.pickle' % (dir, prefix), 'rb')
            matrix=pickle.load(pfile)
            pfile.close()
            index=numpy.loadtxt('%s/chkpt-dbaseindex.dat' % dir, dtype=int)
            matrix=get_matrix(dir, format_dbase, format_dbase, column, max, add, start=index, restart=matrix, prefix=prefix)
        elif os.path.exists('%s/%s-matrix.pickle' % (dir, prefix)):
            print "loading score matrix from %s" % dir
            pfile=open('%s/%s-matrix.pickle' % (dir, prefix), 'rb')
            matrix=pickle.load(pfile)
            pfile.close()

        elif not os.path.exists('%s/%s-matrix.pickle' % (dir, prefix)) and not os.path.exists('%s/%s-matrix-chkpt.pickle' % (dir, prefix)):
            print "getting score matrix"
            matrix=get_matrix(dir, format_dbase, format_dbase, column, max, add, start=0, restart=None,  prefix=prefix)
        matrix=clean_matrix(matrix)
        pfile=open('%s/%s-matrix.pickle' % (dir, prefix), 'wb')
        pickle.dump(matrix, pfile)
        pfile.close()
        print "clustering scores"
        gens, assignments, distances=cluster(cutoff, matrix, format_dbase)
        frames=numpy.where(assignments!=-1)[0]
        distances[frames]=[(max-i) for i in distances[frames]]
        numpy.savetxt('%s-assignments.dat' % dbase.split('.list')[0], assignments)
        numpy.savetxt('%s-distances.dat' % dbase.split('.list')[0], distances)
        print "%s assigned to %s gens" % (len(frames), len(gens))
        for i in gens:
            frames=numpy.where(assignments==i)[0]
            ohandle=open('%s/%s_%s_%s_%s_g%s.dat' % (dir, prefix, format_dbase[i], tanimoto, (max-cutoff), i), 'w')
            for (name, val) in zip(format_dbase[frames], distances[frames]):
                ohandle.write('%s\t%s\n' % (name, val))
            if options.writepdb==True:
                name=format_dbase[i]
                hitfile='%s/all-%s-%s_hits.pdb' % (dir, prefix, name)
                output='g%s_%s' % (i, max-cutoff)
                parse(dir, hitfile, output, format_dbase[frames])
        numpy.savetxt('%s/%s_%s_%s_gen_indices.dat' % (dir, prefix, tanimoto, (max-cutoff)), gens, fmt='%i')
        numpy.savetxt('%s/%s_%s_%s_gen_names.dat' % (dir, prefix, tanimoto,(max-cutoff)), format_dbase[gens], fmt='%s')


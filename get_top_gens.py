import glob
import os
import numpy
gen_inds=[249, 95, 231, 161, 195]
#gen_inds=[95, 161] #249, 95, 231, 161, 195]
print gen_inds
gens=numpy.loadtxt('./reports/gens.list', dtype=str)
bindb=numpy.loadtxt('../bindb_combo_1.5_gen_names.dat', dtype=str)
gdd=numpy.loadtxt('../gdd_combo_1.5_gen_names.dat', dtype=str)
gdd_ind=numpy.loadtxt('../gdd_combo_1.5_gen_indices.dat', dtype=int)
bindb_ind=numpy.loadtxt('../bindb_combo_1.5_gen_indices.dat', dtype=int)
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


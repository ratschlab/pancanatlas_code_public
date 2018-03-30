import sys
import os
import scipy as sp
import glob 
import pdb
import gzip
import h5py

def appendToHDF5(file, data, name):
    
    ### get current shape 
    tmp = file[name].shape
    ### resize
    if len(tmp) == 1:
        file[name].resize((tmp[0] + data.shape[0],))
        file[name][tmp[0]:] = data
    elif len(tmp) == 2:
        file[name].resize((tmp[0], tmp[1] + 1))
        if len(data.shape) < 2:
            file[name][:, tmp[1]:] = data[:, sp.newaxis]
        else:
            file[name][:, tmp[1]:] = data
    else:
        print >> sys.stderr, "cannot append data to HDF5 with more than 2 dimensions" 
        sys.exit(-1)


if len(sys.argv) < 3:
    print >> sys.stderr, 'Usage: <outfile>  <count_pattern>'
    sys.exit(1)
else:
    outfile = sys.argv[1]
    count_pattern = sys.argv[2]

verbose = ("-v" in sys.argv)

files = glob.glob(count_pattern)

OUT = h5py.File(outfile, 'w')

counts = []
header = ['feature']
labels = []
for f, file in enumerate(files):
    if not os.path.exists(file.rsplit('.', 1)[0] + '.done'):
        print >> sys.stderr, '(%i / %i) Skipping %s - incomplete' % (f + 1, len(files), file)
        continue
    if verbose:
        print >> sys.stderr, '(%i / %i) Loading %s' % (f + 1, len(files), file)
    data = sp.loadtxt(file, dtype='str', delimiter='\t')

    sid = file.split('/')[-2]

    if f == 0:
        OUT.create_dataset('sids', data=sp.array([sid]), dtype='|S128', chunks=True, compression='gzip', maxshape=(None,))
        OUT.create_dataset('gids', data=data[:, 0])
        OUT.create_dataset('counts', data=data[:, 1][:, sp.newaxis].astype('int'), chunks=True, compression='gzip', maxshape=(data.shape[0], None))
    else:
        assert(sp.all(OUT['gids'][:] == data[:, 0]))
        appendToHDF5(OUT, sp.array([sid], dtype='|S128'), 'sids')
        appendToHDF5(OUT, data[:, 1].astype('int'), 'counts')
    del data
OUT.close()


# This script read data from grav_ files and save to .npy
# Fortran unformatted files seem to be composed with several blocks of different datatype.
# Each block looks like this:
# |(int) Size of DATA || (dtype) DATA ||(int) Size of DATA |

import sys,glob,time
from struct import *
import numpy as np
import argparse
import pandas as pd

t1 = time.clock()

def get_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Sky-map creation settings."
    )
    parser.add_argument(
        "--snapnr",
        dest='snapnr',
        action="store",
        type=int,
        default=6,
        help="Snapshot number",
    )
    parser.add_argument(
        "--indir",
        dest='indir',
        action="store",
        type=str,
        default="/cosma7/data/dp004/dc-beck3/3_Proca/cvg_b3_000001_with_cbf/",
        help="Folder in which simulation output is saved.",
    )
    parser.add_argument(
        "--outdir",
        dest='outdir',
        action="store",
        type=str,
        default="/cosma7/data/dp004/dc-beck3/3_Proca/cvg_b3_000001_with_cbf/",
        help="Folder in which transformed data will be saved.",
    )
    parser.add_argument(
        "--levelmin",
        dest='levelmin',
        action="store",
        type=int,
        default=8,
        help="Minimum resolution level.",
    )
    parser.add_argument(
        "--levelmax",
        dest='levelmax',
        action="store",
        type=int,
        default=8,
        help="Maximum resolution level.",
    )
    parser.add_argument(
        "--ngrid",
        dest='n_grid',
        action="store",
        type=int,
        default=256,
        help="Number of particles used.",
    )
    args = parser.parse_args()
    return args

def run(args):
    """ """
    # simulation details
    levelmax = args.levelmax
    dimfac= 8  # 2^n, where n is # of dimensions

    # files to be read
    snap_files = glob.glob(
        args.indir+"output_%05d/grav_%05d.out?????" % (args.snapnr, args.snapnr)
    )
    indx = [int(fil.split('.')[-1][-5:]) for fil in snap_files]
    snap_files = [snap_files[ii] for ii in np.argsort(indx)]
    print("There are %d files to be transcribed" % len(snap_files))
    # stored quantities to convert data format of
    nr_quant = 11
    datlis = [[] for t in range(nr_quant)]

    # run through output files of snapshot
    count = 0
    Ngrid = 0
    t1 = time.clock()
    for snap_file in snap_files:
        with open(snap_file, "rb") as f:
            content = f.read()
        print("Reading "+snap_file.split('/')[-1]+" the %d-th file %f" % (count, time.clock()-t1))

        # header info
        pmin=0; pmax=48  # lines of header
        info = unpack("i"*3*4, content[pmin:pmax])
        [ncpu, ndim, _nlevelmax, nboundary] = [info[1], info[4], info[7], info[10]]
        
        # run through resolution levels
        for ilevel in range(args.levelmin, levelmax+1):
            
            # run through nboundary+ncpu
            for ibound in range(1, nboundary+ncpu+1):
                pmin0 = pmax
                pmax0 = pmin0 + 4*3*2
                info = unpack("i"*3*2, content[pmin0:pmax0])
                [currlevel, ncache] = [info[1], info[4]]
                Ngrid += ncache*dimfac

                if ncache == 0:
                    # if no data in that boundary/cpu
                    pmax = pmax0
                    continue
                
                # simulation dimensions
                for dim in range(1, dimfac+1):
                    j = 0
                    
                    # run through parameters
                    for N in range(1, nr_quant+1):
                        pmin = pmax0 + (8*N - 4) + (N - 1)*8*ncache
                        pmax = pmin + ncache*8
                        info = unpack("d"*ncache, content[pmin:pmax])
                        #datlis[j].append(float(info[0]))
                        for floatelem in info:
                            datlis[j].append(floatelem)
                        j += 1
                    
                    pmax0 = pmax + 4
                pmax = pmax0
        f.close()
        count += 1

    print("Transposing: " + str(time.clock() - t1) + " s")
    datlis = np.transpose(datlis)

    print("Mapping to tuple: " + str(time.clock() - t1) + " s")
    datlis = map(tuple, datlis)

    print("Finding set of values: " + str(time.clock() - t1) + " s")
    datlis = set(datlis)

    print("Mapping back to array: " + str(time.clock() - t1) + " s")
    datlis = np.asarray([list(dat) for dat in datlis])
    datlis = np.transpose(datlis)

    data_dic = {
        "x": datlis[0],
        "y": datlis[1],
        "z": datlis[2],
        "phi": datlis[3],
        "f1": datlis[4],
        "f2": datlis[5],
        "f3": datlis[6],
        "sf": datlis[7],
        "cbf1": datlis[8],
        "cbf2": datlis[9],
        "cbf3": datlis[10],
    }
    data_df = pd.DataFrame(data_dic)
    data_df.to_hdf(args.indir+"grav_%05d.h5" % args.snapnr, key='df', mode='w')
    print("Done : " + str(time.clock() - t1) + " s, N = %d" % (len(datlis)))
    sys.exit()

if __name__ == "__main__":
    args = get_args()

    run(args)

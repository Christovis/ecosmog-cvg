import os, sys, glob
import argparse
import numpy as np
import pandas as pd
import astropy
import astropy.units as u
from astropy.cosmology import LambdaCDM
from scipy import ndimage
from nbodykit.lab import *


def get_args(argv=None):
    parser = argparse.ArgumentParser(
        description="The sky is not the limit."
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
        "--boxsize",
        dest='boxsize',
        action="store",
        type=int,
        default=200,
        help="Box-size of simulation in Mpc/h.",
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


def run(args, snap_nrs, quantities):
    """
    """
    for snapnr in snap_nrs:
        print("Reading data of snapshot %d" % snapnr)
        # load snapshot
        infile = args.indir + "grav_%05d.h5" % snapnr
        fields = pd.read_hdf(infile, key="df")

        # Pk settings
        #Delta_k  = 1.0e-2   # size of k bins  (where k is the wave vector in Fourier Space)
        k_min    = 2*np.pi/args.boxsize  # smallest k value

        # value map
        x = (args.n_grid*fields["x"].values).astype(int)
        y = (args.n_grid*fields["y"].values).astype(int)
        z = (args.n_grid*fields["z"].values).astype(int)
        count = 0
        pk_dict = {}
        for quant in quantities:
            print("Power-spectrum of %s" % quant)
            value_map = np.zeros((args.n_grid, args.n_grid, args.n_grid))
            value_map[(x, y, z)] = fields[quant].values

            if quant is "sf":
                print("Take divergence of gradient of sf field")
                di_sf,dj_sf,dk_sf = np.gradient(
                    value_map,
                    1/args.n_grid,
                    1/args.n_grid,
                    1/args.n_grid,
                )
                didi_sf = np.gradient(di_sf, 1/args.n_grid, axis=0)
                djdj_sf = np.gradient(dj_sf, 1/args.n_grid, axis=1)
                dkdk_sf = np.gradient(dk_sf, 1/args.n_grid, axis=2)
                #djdj_sf = np.gradient(
                #    np.gradient(value_map,axis=1,varargs=1/512),
                #    axis=1,
                #    varargs=1/512,
                #)
                #dkdk_sf = np.gradient(
                #    np.gradient(value_map,axis=2,varargs=1/512),
                #    axis=2,
                #    varargs=1/512,
                #)
                if count == 0:
                    #value_map = ndimage.laplace(
                    #    np.gradient(
                    #        value_map,
                    #        axis=0,
                    #        varargs=1/512,
                    #    )
                    #)
                    value_map = (didi_sf + djdj_sf + dkdk_sf)*di_sf
                elif count == 1:
                    #value_map = ndimage.laplace(
                    #    value_map
                    #)
                    value_map = didi_sf + djdj_sf + dkdk_sf
                elif count == 2:
                    value_map = np.gradient(value_map, axis=0)
                elif count == 3:
                    pass
            else:
                pass

            # power-spectrum of density-fluctuations 
            mesh = ArrayMesh(
                value_map, Nmesh=args.n_grid, compensated=False, BoxSize=args.boxsize,
            )
            r = FFTPower(
                mesh,
                mode='1d',
                #dk=Delta_k,
                kmin=k_min
            )
            k = np.array(r.power['k'])                 # the k-bins
            Pk = np.array(r.power['power'].real)       # the power spectrum
            Pk_shotnoise = r.power.attrs['shotnoise']  # shot-noise
            Pk -= Pk_shotnoise
            pk_dict["Pk_%s" % quant+str(count)] = Pk
            #print("PkPkPk", np.mean(Pk))
            count += 1

        pk_dict["k"] = k
        pk_df = pd.DataFrame(pk_dict)
        pk_df.to_hdf(indir+"pk_%05d.h5" % snapnr, key='df', mode='w')

if __name__ == "__main__":
    args = get_args()
    snap_nrs = [6,8,10]  # snapshots
    # TODO: compare \partial^2\partial_i\chi to \partial^2B_i
    quantities = ["sf", "cbf1", "cbf2", "cbf3"]  # cvg fields
    #quantities = ["phi", "sf", "cbf1", "cbf2", "cbf3"]  # cvg fields
    
    run(args, snap_nrs, quantities)

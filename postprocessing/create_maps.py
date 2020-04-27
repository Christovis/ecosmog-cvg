import os, sys, glob
import argparse
import numpy as np
import pandas as pd
from scipy.interpolate import griddata
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import pyplot as plt, rc
from matplotlib.pyplot import figure
import figure_size
plt.style.use("../../publication.mplstyle")


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
        indir = "/cosma7/data/dp004/dc-beck3/3_Proca/cvg_b3_000001_with_cbf/"
        infile = indir + "grav_%05d.h5" % snapnr
        fields = pd.read_hdf(infile, key="df")
       
        frac_x = 1
        
        # take slice
        delta = 1/args.n_grid * (1 + 0.1)
        centre = 0.5
        fields = fields.loc[
            (fields["z"] > centre - delta/2) &\
            (fields["z"] < centre + delta/2)
        ]
        print(len(fields["x"].values))
        
        # Make 2d mesh and map of slice data
        bins = np.linspace(0, frac_x, int(args.n_grid*frac_x))
        bincenter = (bins[1:]+bins[:-1])/2
        xxi, yyi = np.mgrid[0.:1.:256j, 0.:1.:256j]  # depending on grid_size

        for quant in quantities:
            print("Create map of %s" % quant)
            
            if np.min(fields["x"].values) < 0.0:
                print("Map lifted by %.3f" % np.abs(fields["x"].values))
                values = fields[quant].values + np.abs(fields[quant].values)
            else:
                values = fields[quant].values

            value_map = griddata(
                (fields["x"].values, fields["y"].values),
                values,  # avoid <0 values
                (xxi, yyi),
                method='linear',
                fill_value=np.mean(values),
            )
            np.save(indir+"map_%s_%05d.h5" % (quant, snapnr), value_map.T)


if __name__ == "__main__":
    args = get_args()
    snap_nrs = [10, 8, 6]  # snapshots
    # TODO: compare \partial^2\partial_i\chi to \partial^2B_i
    quantities = ["phi", "sf", "cbf1", "cbf2", "cbf3"]  # cvg fields
    
    run(args, snap_nrs, quantities)

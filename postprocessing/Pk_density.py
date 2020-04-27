# Functions to handle power spectra
import os, sys, glob
import numpy as np
import pandas as pd

import astropy
import astropy.units as u
from astropy.cosmology import LambdaCDM

from scipy import ndimage

from nbodykit.lab import *


def read_pk_file(infile, boxsize, Np):
    ik  = np.loadtxt(
        infile,
        usecols=[0],
        skiprows=0,
        unpack=True,
        dtype=np.float32
    )
    P_z00 = np.loadtxt(
        infile,
        usecols=[3],
        skiprows=0,
        unpack=True,
        dtype=np.float32)
    W_z00 = np.loadtxt(
        infile,
        usecols=[4],
        skiprows=0,
        unpack=True,
        dtype=np.float32
    )
    k = ik*2*np.pi/boxsize
    Pk = (P_z00 - W_z00/Np**3)*boxsize**3
    return k, Pk

def read_pk(indir, boxsize, res, Np, snap_nrs, b3, version):
    """
    """
    pk_dic = {
        "k": [],
        "Pk_cvg": [],
        "Pk_qcdm": [],
        "delta": [],
    }
   
    snap_nrs = np.sort(np.array(snap_nrs))
    cc= 0
    for snap_nr in snap_nrs:

        infile = indir+"cvg_b3_%s_v%d/pk_%05d.dat" % (b3, version, snap_nr)
        print(infile)
        k, Pk_cvg = read_pk_file(infile, boxsize, Np)
        pk_dic["k"].append(k)
        pk_dic["Pk_cvg"].append(Pk_cvg)
            
        infile = indir+"ecosmog_qcdm_b3_%s/pk_%05d.dat" % (b3, snap_nr)
        k, Pk_qcdm = read_pk_file(infile, boxsize, Np)
        pk_dic["k"].append(k)
        pk_dic["Pk_qcdm"].append(Pk_qcdm)

        pk_dic["delta"].append((Pk_cvg - Pk_qcdm)/Pk_qcdm)
            
    return pk_dic


#def read_pk(temp_dir, boxsize, res, Np, snap_nrs, b3):
#    """
#    """
#    pk_dic = {
#        "k": [],
#        "Pk_cvg": [],
#        "Pk_cvg_ref": [],
#        "Pk_cvg_ref_lin": [],
#        "Pk_qcdm": [],
#        "Pk_qcdm_ref": [],
#        "delta": [],
#    }
#   
#    snap_nrs = np.sort(np.array(snap_nrs))
#    cc= 0
#    for snap_nr in snap_nrs:
#        indirs = glob.glob(temp_dir + "ecosmog_cvg_b3_%s_*/" % b3)
#        for indir in indirs:
#
#            if "refine" in indir:
#                cvg_file = indir+"%s_v4_b3_%s_refine_Pk_B%d_PM%d_out%05d.dat" % \
#                    ("cvg", b3, boxsize, res, snap_nr)
#                k, Pk_cvg = read_pk_file(cvg_file, boxsize, Np)
#                pk_dic["k"].append(k)
#                pk_dic["Pk_cvg_ref"].append(Pk_cvg)
#            
#                if "linear" in indir:
#                    cvg_file = indir+"%s_v4_b3_%s_refine_linear_Pk_B%d_PM%d_out%05d.dat" % \
#                        ("cvg", b3, boxsize, res, snap_nr)
#                    k, Pk_cvg = read_pk_file(cvg_file, boxsize, Np)
#                    pk_dic["k"].append(k)
#                    pk_dic["Pk_cvg_ref_lin"].append(Pk_cvg)
#            
#            elif "linear" in indir:
#                cvg_file = indir+"%s_v4_b3_%s_linear_Pk_B%d_PM%d_out%05d.dat" % \
#                    ("cvg", b3, boxsize, res, snap_nr)
#                k, Pk_cvg = read_pk_file(cvg_file, boxsize, Np)
#                pk_dic["k"].append(k)
#                pk_dic["Pk_cvg_lin"].append(Pk_cvg)
#            else:
#                cvg_file = indir+"%s_v4_b3_%s_Pk_B%d_PM%d_out%05d.dat" % \
#                    ("cvg", b3, boxsize, res, snap_nr)
#                k, Pk_cvg = read_pk_file(cvg_file, boxsize, Np)
#                pk_dic["k"].append(k)
#                pk_dic["Pk_cvg"].append(Pk_cvg)
#        
#        indirs = glob.glob(temp_dir + "ecosmog_qcdm_b3_%s_*/" % b3)
#        for indir in indirs:
#
#            if "refine" in indir:
#                cvg_file = indir+"%s_v4_b3_%s_refine_Pk_B%d_PM%d_out%05d.dat" % \
#                    ("qcdm", b3, boxsize, res, snap_nr)
#                k, Pk = read_pk_file(cvg_file, boxsize, Np)
#                pk_dic["k"].append(k)
#                pk_dic["Pk_qcdm_ref"].append(Pk)
#                
#                pk_dic["delta_ref"].append((pk_dic["Pk_cvg_ref"][cc] - pk_dic["Pk_qcdm_ref"][cc])/pk_dic["Pk_qcdm_ref"][cc])
#            
#            else:
#                cvg_file = indir+"%s_v4_b3_%s_Pk_B%d_PM%d_out%05d.dat" % \
#                    ("qcdm", b3, boxsize, res, snap_nr)
#                k, Pk = read_pk_file(cvg_file, boxsize, Np)
#                pk_dic["k"].append(k)
#                pk_dic["Pk_qcdm"].append(Pk)
#
#                pk_dic["delta"].append((pk_dic["Pk_cvg"][cc] - pk_dic["Pk_qcdm"][cc])/pk_dic["Pk_qcdm"][cc])
#        cc += 1
#    return pk_dic


def align_lin_nonlin(lin, nonlin, k):
    return lin[0] - np.mean(nonlin[(1e-2 < k) & (k<1e-1)])


def run(snap_nrs, params_b3, quantities):
    """
    """
    for param_b3 in params_b3:
        indir = "/cosma7/data/dp004/dc-beck3/3_Proca/ecosmog_cvg_b3_%s/" % param_b3

        for snapnr in snap_nrs:
            print("Reading data of snapshot %d" % snapnr)
            # load snapshot
            infile = indir + "grav_%05d.h5" % snapnr
            fields = pd.read_hdf(infile, key="df")

            # Pk settings
            Boxsize = 200  #[Mpc/h]
            grid_size  = 512  # default number of mesh cells per coordinate axis
            #Delta_k  = 1.0e-2   # size of k bins  (where k is the wave vector in Fourier Space)
            k_min    = 2*np.pi/Boxsize  # smallest k value

            # value map
            x = (grid_size*fields["x"].values).astype(int)
            y = (grid_size*fields["y"].values).astype(int)
            z = (grid_size*fields["z"].values).astype(int)

            pk_dict = {}
            for quant in quantities:
                value_map = np.zeros((grid_size,grid_size,grid_size))

                if quant is "lp_dx_sf":
                    value_map[(x, y, z)] = fields["sf"].values
                    di_sf = np.gradient(
                        value_map,1/512,axis=0,edge_order=2
                    )
                    di_di_sf,dj_di_sf,dk_di_sf = np.gradient(
                        di_sf,1/512,1/512,1/512
                    )
                    didi_di_sf = np.gradient(di_di_sf,1/512,axis=0,edge_order=2)
                    djdj_di_sf = np.gradient(dj_di_sf,1/512,axis=1,edge_order=2)
                    dkdk_di_sf = np.gradient(dk_di_sf,1/512,axis=2,edge_order=2)
                    value_map = (didi_di_sf + djdj_di_sf + dkdk_di_sf)
                    value_map[0:5, :] = value_map[5:10, :]
                    value_map[-6:-1, :] = value_map[-10:-5, :]
                    value_map[abs(value_map) > 5e5] = 0.
                    label = "lp_dx_sf"
                else:
                    value_map[(x, y, z)] = fields[quant].values
                    label = quant
                print("Power-spectrum of %s" % label)

                # power-spectrum of density-fluctuations 
                mesh = ArrayMesh(
                    value_map, Nmesh=grid_size, compensated=False, BoxSize=Boxsize,
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
                pk_dict["Pk_%s" % label] = Pk
                print("PkPkPk", np.mean(Pk))

            pk_dict["k"] = k
            pk_df = pd.DataFrame(pk_dict)
            pk_df.to_hdf(indir+"pk_%05d.h5" % snapnr, key='df', mode='w')

if __name__ == "__main__":
    snap_nrs = [13]  # snapshots [13, 20]
    b3 = ["000001"]  # param_b3 ["01", "0001", "000001"]
    quantities = ["phi", "sf", "lp_dx_sf", "cbf1", "cbf2", "cbf3"]  # cvg fields
    
    run(snap_nrs, b3, quantities)

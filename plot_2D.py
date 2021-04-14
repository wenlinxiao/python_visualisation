"""
Generating 2-D plots of energy spectra

Input: binary data file of En(k), En(n)

Usage: Change the input data file path f_name;
       Change the output figure file path in plt.savefig();
       Run the script in terminal
@author: Wenlin Xiao
"""

import netCDF4 as nc
import numpy as np
import logging
import os
import matplotlib.pylab as plt
import matplotlib.colors as colors
#from matplotlib.backends.backend_pdf import PdfPages

def plot_2D_contourf(E_knl,gtype):
    ''' plot contourf for all vertical levels, all gtype
        x axis: k
        y axis: n
    '''
    #fig_path=path_to_save+"/kin_en_kn_"+gtype
    #if not os.path.isdir(fig_path):
        #os.makedirs(fig_path)
    #for lev in range(len(datain["level"])):
    for ind in range(len(level_inds)):    #
        l=level_inds[ind]    # wanted level's index in the level list in data
        level_value=level[l]
        print("l",l,"level_value",level_value,np.shape(E_knl))
        with np.errstate(divide="ignore"):
            E_knl_log=np.log(E_knl[:,:,l])
        print(wave_short,np.amax(E_knl_log),np.amin(E_knl_log))

        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)    # set levels
        plt.contourf(range(nk),range(nn),np.transpose(E_knl_log),norm=norm,levels=bounds)
        cb=plt.colorbar()
#        cb.ax.set_yticklabels([""]+bounds[1:])
        # plt.title("Kinetic Energy Distribution at Level: \u03C3=%f" % datain['level'][level]+" (log(KE))")
        #plt.title(wave_name+": Kinetic energy distribution \napprox. at level: %d hPa (log(KE))" % level_value)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim(1,350)
        plt.ylim(1,350)
        plt.xlabel("Zonal wavenumber k")
        plt.ylabel("Total wavenumber n")
        plt.savefig("./results_202010_mean/20201123/E_2D_"+wave_short+"_"+gtype+"_l%03d_log.eps" %l)
        plt.close("all")


def plot_2D_contourf_largescale(E_knl,gtype):
    ''' plot contourf for all vertical levels, all gtype
        x axis: k
        y axis: n
    '''
    #fig_path=path_to_save+"/kin_en_kn_"+gtype
    #if not os.path.isdir(fig_path):
        #os.makedirs(fig_path)
    #for lev in range(len(datain["level"])):
    for ind in range(len(level_inds)):    #
        l=level_inds[ind]    # wanted level's index in the level list in data
        level_value=level[l]
        with np.errstate(divide="ignore"):
            E_knl_log=np.log(E_knl[:,:,l])

        norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)    # set levels
        plt.contourf(range(nk),range(nn),np.transpose(E_knl_log),norm=norm,levels=bounds)
        cb=plt.colorbar()
#        cb.ax.set_yticklabels([""]+bounds[1:])
        # plt.title("Kinetic Energy Distribution at Level: \u03C3=%f" % datain['level'][level]+" (log(KE))")
        #plt.title(wave_name+": Kinetic energy distribution \napprox. at level: %d hPa (log(KE))" % level_value)
        plt.xscale("log")
        plt.yscale("log")
        plt.xlim(1,50)
        plt.ylim(1,50)
        plt.xlabel("Zonal wavenumber k")
        plt.ylabel("Total wavenumber n")
        plt.savefig("./results_202010_mean/20201123/E_2D_"+wave_short+"_"+gtype+"_l%03d_log_xy50_2.eps" %l)    #png,dpi=300
        plt.close("all")


 
# main
# Read data
nk=351
nn=351
nl=137

# read the pressure levels / sigma levels
level=np.fromfile("levels_202010",sep="",dtype=np.float32)    # dtype must be included, don't know why

# levels=[1,2,3,5,7,10,20,30,50,70,100,150,200,250,300,400,500]
# level_inds=[13,17,19,23,25,28,35,40,47,53,59,67,73,78,82,89,95]
level_inds=[65]    #65,80
#level_inds=np.arange(0,137,5)

#wave_names=["IG modes","Rossby modes","Total modes"]    #,"Era5 u-v"
#wave_names=["Total modes","Era5 u-v"]    #
wave_names=["Total modes"]
i=-1
#for wave_short in ["IG","ROS","TOT"]:    #,"Era5"
for wave_short in ["TOT"]:    #,"Era5"
    i+=1
    wave_name=wave_names[i]
    irrot_knl_sum=np.zeros((nk,nl))
    nondiv_knl_sum=np.zeros((nk,nl))
    total_knl_sum=np.zeros((nk,nl))
    for d in range(15,22):
    # for d in [15]:    # only one day!!!
        for t in ["00","06","12","18"]:
            f_name="../../../../../../scratch/cen/mi-theo/u300924/sph/sph_pyspharm_oper/data/202010/results_data_202010_mean/"+wave_short+"_newnorm/e_spectra_irrot_knl_"+wave_short+"_202002"+str(d)+t
            inp_tmp=np.fromfile(f_name)    #IG
            irrot_knl_tmp=inp_tmp.reshape(351,351,137)    # zonal wavenumber 640, vertical level 137
            irrot_knl_sum=irrot_knl_sum+irrot_knl_tmp

            f_name="../../../../../../scratch/cen/mi-theo/u300924/sph/sph_pyspharm_oper/data/202010/results_data_202010_mean/"+wave_short+"_newnorm/e_spectra_nondiv_knl_"+wave_short+"_202002"+str(d)+t
            inp_tmp=np.fromfile(f_name)    #IG
            nondiv_knl_tmp=inp_tmp.reshape(351,351,137)    # zonal wavenumber 640, vertical level 137
            nondiv_knl_sum=nondiv_knl_sum+nondiv_knl_tmp

            f_name="../../../../../../scratch/cen/mi-theo/u300924/sph/sph_pyspharm_oper/data/202010/results_data_202010_mean/"+wave_short+"_newnorm/e_spectra_total_knl_"+wave_short+"_202002"+str(d)+t
            inp_tmp=np.fromfile(f_name)    #IG
            total_knl_tmp=inp_tmp.reshape(351,351,137)    # zonal wavenumber 640, vertical level 137
            total_knl_sum=total_knl_sum+total_knl_tmp
        
    irrot_knl=irrot_knl_sum/(7*4)
    nondiv_knl=nondiv_knl_sum/(7*4)
    total_knl=total_knl_sum/(7*4)
    print(nondiv_knl[3,3,3],total_knl[3,3,3])
    
       
    # same colorbar for irrot, nondiv, l65
    if wave_short=="TOT":
        bounds=[-42,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
    if wave_short=="Era5":
        bounds=[-16,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
    plot_2D_contourf(irrot_knl,"irrot")
    plot_2D_contourf_largescale(irrot_knl,"irrot")
    
    
    if wave_short=="TOT":
        bounds=[-42,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
    if wave_short=="Era5":
        bounds=[-16,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
    plot_2D_contourf(nondiv_knl,"nondiv")
    plot_2D_contourf_largescale(nondiv_knl,"nondiv")
#
#    
    if wave_short=="TOT":
        bounds=[-42,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
    if wave_short=="Era5":
        bounds=[-15,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5]
    plot_2D_contourf(total_knl,"total")
    plot_2D_contourf_largescale(total_knl,"total")

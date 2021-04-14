#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Averaging 1D energy spectra En(k), En(n); Plotting En(k), En(n), En(k)& En(n); Calculating irrot./total ratio

Input: binary data of En(k), En(n) (output of loop_plot.sh)
Output: 1-D mean energy spectra figures, irrot/total ratio figures

Usage: Comment/uncomment the functions at the end of the script as needed;
       Change the output figure path plt.savefig();
       Run the script

@author: Wenlin Xiao
"""

import os
import numpy as np
#import netCDF4 as nc
from matplotlib import pyplot as plt
import matplotlib.colors as colors
from matplotlib.backends.backend_pdf import PdfPages

def plot_spectra():
    if not os.path.isdir("./results_202010_mean/20210104/"+wave_short+"/"+kn+"_full/"):
        os.makedirs("./results_202010_mean/20210104/"+wave_short+"/"+kn+"_full/")
#    pp = PdfPages("./results_202010_mean/20201123/"+wave_short+"_newnorm_every5th/"+kn+"_clearer_slope_2/"+wave_short+"_"+kn+".pdf")
    pp = PdfPages("./results_202010_mean/20210104/"+wave_short+"/"+kn+"_full/"+wave_short+"_"+kn+"full.pdf")

    for ind in range(len(level_inds)):
        l=level_inds[ind]    # wanted level's index in the level list in data
    #     level_value=levels[ind]
        level_value=level[l]
        fig, ax = plt.subplots(figsize=(8, 8))
        # total energy
        ax.loglog(np.arange(0,350),total_kl[:350,l],"-",color="black",lw=5,label="total")    #5
        # irrot energy
    #if 'irrot' in data_kin_en_sp['sph'].keys():
        ax.loglog(np.arange(0,350),irrot_kl[:350,l],"-",color="blue",lw=3,label="irrotational")    #1.5
    #        ax.loglog(datain['k'],np.ones(len(datain['k']))*data_kin_en_sp['sph']['irrot']['div_vort']['total'][level],\
    #                "--",lw=2,color="red",\
    #                label="irrotational (total on level from DIV,VORT) = %.6f J/kg" % (data_kin_en_sp['sph']['irrot']['div_vort']['total'][level]))
        # nondiv energy
    #if 'nondiv' in data_kin_en_sp['sph'].keys():
        ax.loglog(np.arange(0,350),nondiv_kl[:350,l],"-",color="red",lw=3,label="nondivergent")    #1.5
    #        ax.loglog(datain['k'],np.ones(len(datain['k']))*data_kin_en_sp['sph']['nondiv']['div_vort']['total'][level],\
    #                "--",lw=3,color="blue",\
    #                label="nondivergent (total on level from DIV,VORT) = %.6f J/kg" % (data_kin_en_sp['sph']['nondiv']['div_vort']['total'][level]))
        ax.legend(loc='lower left',frameon=False,fontsize=15)
        # reference line: k=-5/3,-3
        x1=[10,100]
        x2=[10,100]    #[np.exp(10),np.exp(100)]
    #     y1=[10,10**(1/3)]    #[np.exp(10),np.exp(-140)]    # old normalisation
    #     y2=[10,10**-2]    #[np.exp(10),np.exp(-260)]
        y1=[100,10**-1]    # new normalisation, for k
        y2=[100,10**(1/3)]    #
        if wave_short=="ROS" or wave_short=="TOT":
            x1=[20,200]    # new normalisation, for ROS-n, Total
            x2=[20,200]    #[np.exp(10),np.exp(100)]
            y1=[100,10**-1]    # new normalisation, for ROS-n, Total
            y2=[100,10**(1/3)]    #
        ax.loglog(x1,y1,":",color="black",lw=1.5,label="reference slope: -5/3")
        ax.loglog(x2,y2,":",color="black",lw=1.5,label="reference slope: -3")
    #     ax.scatter(15,100,color="green",s=5)    # check the point position in loglog plot
        
        ax.set_xticks([1,10,100,1000],minor=False)
        ax.set_xticklabels([1,10,100,1000])    #,fontsize=15
        ax.tick_params(labelsize=15)
        ax.set_xlim(1,350)    #np.nanmax(np.arange(0,nk))
        # set limit in y axes
        #_max_en=data_kin_en_sp['sph']['total']['div_vort']['total'][level]
        #_min_en=1e-7
        #_min_en=min(np.amin(data_kin_en_sp['sph']['nondiv']['div_vort']['dep_k'][:,level]),np.amin(data_kin_en_sp['sph']['irrot']['div_vort']['dep_k'][:,level]))
#        if wave_short=="ROS":
#            ax.set_ylim(10e-11,10e3)    # for ROS
#        if wave_short=="IG":
#            ax.set_ylim(10e-6,10e1)    #for IG
##        if wave_short=="TOT":
##            ax.set_ylim(10e-4,10e3)    # for Total
#        if wave_short=="TOT":
#            ax.set_ylim(10e-4,10e2)    # for Total-k, in folder /k_clearer_slope_2
        if kn=="k":
            ax.set_xlabel("Zonal wavenumber k",fontsize=15)
        if kn=="n":
            ax.set_xlabel("Global wavenumber n",fontsize=15)
        ax.set_ylabel("Kinetic energy (J/kg)",fontsize=15)
        #plt.title("Kinetic energy spectrum log(J/kg) Spherepack vs. ECMWF\nvertical level: %d (might be Pa or model level id)" % datain['level'][level])
        #plt.title("Kinetic energy spectrum log(J/kg) of vertical level: \u03C3=%f" % datain['level'][level])
        #plt.title(wave_name+": Kinetic energy spectrum log(J/kg) of vertical level: approx. %d hPa" % level_value)
    #     plt.title(wave_name+": approx. %d hPa" % level_value)
#        plt.title(wave_name+": %.2f hPa" % level_value)
    #    if path_to_file==None:
#        plt.savefig("./results_202010_mean/20201123/"+wave_short+"_newnorm_every5th/"+kn+"_clearer_slope_2/"+wave_short+"_"+kn+"_l%03d.eps" % l,bbox_inches="tight")    #,dpi=300
#        plt.savefig("./results_202010_mean/20201123/"+wave_short+"_newnorm_every5th/"+kn+"_full/"+wave_short+"_"+kn+"_l%03d_full.eps" % l,bbox_inches="tight")    #,dpi=300
        plt.savefig("./results_202010_mean/20210104/"+wave_short+"/"+kn+"_full/"+wave_short+"_"+kn+"_l%03d_full.eps" % l,bbox_inches="tight")    #,dpi=300
#         plt.show()
    #     plt.close("all")
#         plt.close("all")
        pp.savefig(fig)
        plt.close()
    pp.close()
    
def plot_spectra_kn():
    #pp = PdfPages("./results_202010_mean/20201123/"+wave_short+"_newnorm_every5th/kn_clearer_slope/"+wave_short+"_kn.pdf")
#    pp = PdfPages("./results_202010_mean/20201123/"+wave_short+"_newnorm_every5th/kn_full/"+wave_short+"_kn_full.pdf")
    if not os.path.isdir("./results_202010_mean/20210104/"+wave_short+"/kn_full/"):
        os.makedirs("./results_202010_mean/20210104/"+wave_short+"/kn_full/")
    pp = PdfPages("./results_202010_mean/20210104/"+wave_short+"/kn_full/"+wave_short+"_kn_full.pdf")
    for ind in range(len(level_inds)):
        l=level_inds[ind]    # wanted level's index in the level list in data
    #     level_value=levels[ind]
        level_value=level[l]
        fig, ax = plt.subplots(figsize=(8, 8))
        # total energy
        ax.loglog(np.arange(1,350),Ekn["k"]["total"][1:350,l],"-",color="black",lw=5,label="total - k")    #5
        # irrot energy
    #if 'irrot' in data_kin_en_sp['sph'].keys():
        ax.loglog(np.arange(1,350),Ekn["k"]["irrot"][1:350,l],"-",color="blue",lw=3,label="irrotational - k")    #1.5
    #        ax.loglog(datain['k'],np.ones(len(datain['k']))*data_kin_en_sp['sph']['irrot']['div_vort']['total'][level],\
    #                "--",lw=2,color="red",\
    #                label="irrotational (total on level from DIV,VORT) = %.6f J/kg" % (data_kin_en_sp['sph']['irrot']['div_vort']['total'][level]))
        # nondiv energy
    #if 'nondiv' in data_kin_en_sp['sph'].keys():
        ax.loglog(np.arange(1,350),Ekn["k"]["nondiv"][1:350,l],"-",color="red",lw=3,label="nondivergent - k")    #1.5
        
        
        
        ax.loglog(np.arange(1,350),Ekn["n"]["total"][1:350,l],"--",color="black",lw=3,label="total - n")    #5
        ax.loglog(np.arange(1,350),Ekn["n"]["irrot"][1:350,l],"--",color="blue",lw=3,label="irrotational - n")    #1.5
        ax.loglog(np.arange(1,350),Ekn["n"]["nondiv"][1:350,l],"--",color="red",lw=3,label="nondivergent - n")    #1.5
        
        
    #        ax.loglog(datain['k'],np.ones(len(datain['k']))*data_kin_en_sp['sph']['nondiv']['div_vort']['total'][level],\
    #                "--",lw=3,color="blue",\
    #                label="nondivergent (total on level from DIV,VORT) = %.6f J/kg" % (data_kin_en_sp['sph']['nondiv']['div_vort']['total'][level]))
        ax.legend(loc='lower left',frameon=False,fontsize=15)
        # reference line: k=-5/3,-3
        x1=[10,100]
        x2=[10,100]    #[np.exp(10),np.exp(100)]
    #     y1=[10,10**(1/3)]    #[np.exp(10),np.exp(-140)]    # old normalisation
    #     y2=[10,10**-2]    #[np.exp(10),np.exp(-260)]
        y1=[100,10**-1]    # new normalisation, for k
        y2=[100,10**(1/3)]    #
        if wave_short=="ROS" or wave_short=="TOT" or wave_short=="Era5":
            x1=[20,200]    # new normalisation, for ROS-n, Total
            x2=[20,200]    #[np.exp(10),np.exp(100)]
            y1=[100,10**-1]    # new normalisation, for ROS-n, Total
            y2=[100,10**(1/3)]    #
        ax.loglog(x1,y1,":",color="black",lw=1.5,label="reference slope: -5/3")
        ax.loglog(x2,y2,":",color="black",lw=1.5,label="reference slope: -3")
    #     ax.scatter(15,100,color="green",s=5)    # check the point position in loglog plot
        
        ax.set_xticks([1,10,100,1000],minor=False)
        ax.set_xticklabels([1,10,100,1000])
        ax.tick_params(labelsize=15)
        ax.set_xlim(1,350)    #np.nanmax(np.arange(0,nk))
        # set limit in y axes
        #_max_en=data_kin_en_sp['sph']['total']['div_vort']['total'][level]
        #_min_en=1e-7
        #_min_en=min(np.amin(data_kin_en_sp['sph']['nondiv']['div_vort']['dep_k'][:,level]),np.amin(data_kin_en_sp['sph']['irrot']['div_vort']['dep_k'][:,level]))
#        if wave_short=="ROS":
#            ax.set_ylim(10e-11,10e3)    # for ROS
#        if wave_short=="IG":
#            ax.set_ylim(10e-6,10e1)    #for IG
#        if wave_short=="TOT":
#           ax.set_ylim(10e-4,10e3)    # for Total
        ##if wave_short=="TOT":
            ##ax.set_ylim(10e-4,10e2)    # for Total-k, in folder /k_clearer_slope_2
        #if kn=="k":
            #ax.set_xlabel("k (zonal wave number)")
        #if kn=="n":
        ax.set_xlabel("Zonal wavenumber k / global wavenumber n",fontsize=15)
        ax.set_ylabel("Kinetic energy (J/kg)",fontsize=15)
        #plt.title("Kinetic energy spectrum log(J/kg) Spherepack vs. ECMWF\nvertical level: %d (might be Pa or model level id)" % datain['level'][level])
        #plt.title("Kinetic energy spectrum log(J/kg) of vertical level: \u03C3=%f" % datain['level'][level])
        #plt.title(wave_name+": Kinetic energy spectrum log(J/kg) of vertical level: approx. %d hPa" % level_value)
    #     plt.title(wave_name+": approx. %d hPa" % level_value)
#        plt.title(wave_name+": %.2f hPa" % level_value)
    #    if path_to_file==None:
        #plt.savefig("./results_202010_mean/20201123/"+wave_short+"_newnorm_every5th/kn_clearer_slope/"+wave_short+"_kn_l%03d.eps" % l,bbox_inches="tight")    #,dpi=300
#        plt.savefig("./results_202010_mean/20201123/"+wave_short+"_newnorm_every5th/kn_full/"+wave_short+"_kn_l%03d_full.eps" % l,bbox_inches="tight")    #,dpi=300
        plt.savefig("./results_202010_mean/20210104/"+wave_short+"/kn_full/"+wave_short+"_kn_l%03d_full.eps" % l,bbox_inches="tight")    #,dpi=300
#         plt.show()
    #     plt.close("all")
#         plt.close("all")
        pp.savefig(fig)
        plt.close()
    pp.close()
        
        
def plot_ratio():
    fig,ax=plt.subplots()
#     bounds = np.arange(0,1.1,0.1)*100
    bounds = [0,5,10,15,20,30,40,50,60,70,80,90,100]
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=256)
    plt.contourf(range(1,351),level,np.transpose(irrot_ratio)*100,cmap="YlOrBr",norm=norm,levels=bounds)
    plt.gca().invert_yaxis()
    cb=plt.colorbar(spacing="proportional")    #label=
    # cb.ax.set_xticklabels(["0%","20%","40%","60%","80%","100%"])
#     cb.set_label("%", labelpad=-40, y=1.05, rotation=0)
    cb.ax.set_title("%")
    plt.xlim(1,200)
    plt.ylim(700,0)
    plt.xscale("log")
    # plt.yscale("log")
    # plt.title("Total modes: E(irrot) / [E(irrot)+E(nondiv)]")
#    plt.title(wave_name+": E(irrot) / [E(irrot)+E(nondiv)]")
    # plt.xlabel("Zonal wavenumber k")
    plt.xlabel("Global wavenumber "+kn)
    plt.ylabel("Pressure level (hPa)")
    plt.savefig("./results_202010_mean/20201123/irrot_to_total_"+wave_short+"_"+kn+"_xlog.eps")    #,dpi=300
#     plt.show()





# main
# Read data
nk=351
nl=137

# read the pressure levels / sigma levels
level=np.fromfile("levels_202010",sep="",dtype=np.float32)    # dtype must be included, don't know why

# levels=[1,2,3,5,7,10,20,30,50,70,100,150,200,250,300,400,500]
# level_inds=[13,17,19,23,25,28,35,40,47,53,59,67,73,78,82,89,95]
level_inds=np.arange(0,137,5)

wave_names=["IG modes","Rossby modes","Total modes","Era5 u-v"]
Ekn=dict()
i=-1
for wave_short in ["IG","ROS","TOT","Era5"]:
    i+=1
    wave_name=wave_names[i]
    for kn in ["k","n"]:    #
        irrot_kl_sum=np.zeros((nk,nl))
        nondiv_kl_sum=np.zeros((nk,nl))
        total_kl_sum=np.zeros((nk,nl))
        for d in range(15,22):
        # for d in [15]:    # only one day!!!
            for t in ["00","06","12","18"]:
                f_name="./results_data_202010_mean/"+wave_short+"_newnorm/e_spectra_irrot_"+kn+"l_"+wave_short+"_202002"+str(d)+t
                inp_tmp=np.fromfile(f_name)    #IG
                irrot_kl_tmp=inp_tmp.reshape(351,137)    # zonal wavenumber 640, vertical level 137
                irrot_kl_sum=irrot_kl_sum+irrot_kl_tmp

                f_name="./results_data_202010_mean/"+wave_short+"_newnorm/e_spectra_nondiv_"+kn+"l_"+wave_short+"_202002"+str(d)+t
                inp_tmp=np.fromfile(f_name)    #IG
                nondiv_kl_tmp=inp_tmp.reshape(351,137)    # zonal wavenumber 640, vertical level 137
                nondiv_kl_sum=nondiv_kl_sum+nondiv_kl_tmp

                f_name="./results_data_202010_mean/"+wave_short+"_newnorm/e_spectra_total_"+kn+"l_"+wave_short+"_202002"+str(d)+t
                inp_tmp=np.fromfile(f_name)    #IG
                total_kl_tmp=inp_tmp.reshape(351,137)    # zonal wavenumber 640, vertical level 137
                total_kl_sum=total_kl_sum+total_kl_tmp
        
        irrot_kl=irrot_kl_sum/(7*4)
        nondiv_kl=nondiv_kl_sum/(7*4)
        total_kl=total_kl_sum/(7*4)

        # En(k), En(n) spectra plotting
        plot_spectra()

        # irrot/total ratio calculation
        irrot_ratio=irrot_kl[1:,]/total_kl[1:,]
        plot_ratio()
        
        # k,n spectra
        Ekn[kn]=dict()
        Ekn[kn]["irrot"]=irrot_kl
        Ekn[kn]["nondiv"]=nondiv_kl
        Ekn[kn]["total"]=total_kl
    # k,n spectra on the same plot
    plot_spectra_kn()
        


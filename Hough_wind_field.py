#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plotting wind field of Hough harmonics
Usage: Change the output figure path in plt.savefig();
       Choose the projection: proj="rob" / "eqcyl";
       Run the script

Created on Thu Dec 17 15:00:13 2020

@author: Wenlin Xiao
"""

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import os
import netCDF4 as nc

if __name__=="__main__":
    proj="rob"    # "eqcyl"
    level_ind=9
    
    path='../../../../../../scratch/cen/mi-theo/Share/MODES_Outputs/Hough/'
    inp_fns=sorted(os.listdir(path))
    inp_fns=inp_fns[:]    # change the file index accordingly

    for f in inp_fns:    #[0:1]
        mode=f[9:12]
        if mode=="eig":
            mode_string="EIG"
        if mode=="wig":
            mode_string="WIG"
        if mode=="rot":
            mode_string="Rossby"
        kk=f[7]
        if f[9:14]=="rot10":
            nn=f[12:14]
            if f[17]=="_":
                mm=f[16]
            else:
                mm=f[16:18]
        else:
            nn=f[12]
            if f[16]=="_":
                mm=f[15]
            else:
                mm=f[15:17]
        ed=["10 km","6905 m","3040 m","1536 m","886 m","576 m","413 m","313 m","247 m","203 m","173 m","152 m","131 m","110 m","92 m","76 m","64 m","53 m","45 m","38 m","33 m","27 m","21 m","17 m","13 m","10 m","8 m","7 m","6 m","5 m","4 m","3.2 m","2.6 m","2.2 m","1.9 m","1.7 m","1.5 m","1.4 m","1.2 m","1.0 m","0.6 m","0.3 m","0.06 m"]
        ed_mm=ed[int(mm)-1]
        print(kk,nn,mm)
        inp=nc.Dataset(path+f,format="NETCDF4")
        u=inp.variables["u"][0,level_ind,:,:]    # [level_p,lat,lon] [l=10,128,256]
        v=inp.variables["v"][0,level_ind,:,:]    # [level_p,lat,lon] [l=10,128,256]
        Z=inp.variables["Z"][0,level_ind,:,:]
        lon=inp.variables["lon"][:]
        lat=inp.variables["lat"][:]
        wind_spd=np.sqrt(u**2+v**2)

        # Shifting the longitudes to centralize the wave
        lon_ct_ind=np.argmax(u[64,:])    # index of the center of the equatorial positive half, (0,255)
        lon_ct=lon[lon_ct_ind]    #???
        print(lon_ct_ind,lon_ct)
        
        if lon_ct_ind>=64:
            lon_pjct_ind=lon_ct_ind-64
        else:
            lon_pjct_ind=lon_ct_ind+192    #+192
        if lon_pjct_ind>=128:
            lon_pjst_ind=lon_pjct_ind-128
        else:
            lon_pjst_ind=lon_pjct_ind+128
        lon_pjst=lon[lon_pjst_ind]
        
        print(lon_pjst)
        Z_shft,lon_shft=shiftgrid(180.,Z,lon,start=False)    #lon_pjst
        lon_pjct=lon_pjst-180
        if lon_pjct>0:
            lon_pjst=lon_pjct
            lon_pjct=lon_pjct-180    # lon_0 is [-360,0] so E and W don't flip
            
        else:
            lon_pjst=lon_pjct+180
        Z_shft,lon_shft=shiftgrid(lon_pjst,Z,lon,start=False)
        
        
                
        
        # plot
        fig,ax=plt.subplots()
        fig.set_tight_layout(True)        
        
        #### Robinson
        if proj=="rob":
            fig.set_size_inches(8,4.5)
            m = Basemap(projection='robin',lon_0=lon_pjct,resolution='c')    #"lon_0 must be between -360.000000 and 720.000000 degrees"; BUT when lon_0>0, the map reverses horizantally    ,llcrnrlon=lon_pjct-180,urcrnrlon=lon_pjct+180,llcrnrlat=-90,urcrnrlat=90
        
        #### Equidistant Cylindrical
        if proj=="eqcyl":
            if mode=="wig":
                if int(mm)<10:
                    print("lat range 1")
                    m = Basemap(projection='cyl',llcrnrlon=lon_pjct-180,urcrnrlon=lon_pjct+180,llcrnrlat=-85,urcrnrlat=85,resolution='c')    #"lon_0 must be between -360.000000 and 720.000000 degrees"; BUT when lon_0>0, the map reverses horizantally    ,llcrnrlon=lon_pjct-180,urcrnrlon=lon_pjct+180,llcrnrlat=-90,urcrnrlat=90
                    fig.set_size_inches(8,4.5)
                if int(mm)>=10 and int(mm)<18:
                    print("lat range 2")
                    m = Basemap(projection='cyl',llcrnrlon=lon_pjct-180,urcrnrlon=lon_pjct+180,llcrnrlat=-60,urcrnrlat=60,resolution='c')    #"lon_0 must be between -360.000000 and 720.000000 degrees"; BUT when lon_0>0, the map reverses horizantally
                    fig.set_size_inches(8,3)
                if int(mm)>=18:
                    print("lat range 3")
                    m = Basemap(projection='cyl',llcrnrlon=lon_pjct-180,urcrnrlon=lon_pjct+180,llcrnrlat=-40,urcrnrlat=40,resolution='c')    #"lon_0 must be between -360.000000 and 720.000000 degrees"; BUT when lon_0>0, the map reverses horizantally
                    fig.set_size_inches(8,2.5)
            else:
                if int(mm)<10:
                    print("lat range 1")
                    m = Basemap(projection='cyl',llcrnrlon=lon_pjct-180,urcrnrlon=lon_pjct+180,llcrnrlat=-85,urcrnrlat=85,resolution='c')    #"lon_0 must be between -360.000000 and 720.000000 degrees"; BUT when lon_0>0, the map reverses horizantally    ,llcrnrlon=lon_pjct-180,urcrnrlon=lon_pjct+180,llcrnrlat=-90,urcrnrlat=90
                    fig.set_size_inches(8,4.5)
                if int(mm)>=10 and int(mm)<18:
                    print("lat range 2")
                    m = Basemap(projection='cyl',llcrnrlon=lon_pjct-180,urcrnrlon=lon_pjct+180,llcrnrlat=-45,urcrnrlat=45,resolution='c')    #"lon_0 must be between -360.000000 and 720.000000 degrees"; BUT when lon_0>0, the map reverses horizantally
                    fig.set_size_inches(8,2.5)
                if int(mm)>=18:
                    print("lat range 3")
                    m = Basemap(projection='cyl',llcrnrlon=lon_pjct-180,urcrnrlon=lon_pjct+180,llcrnrlat=-30,urcrnrlat=30,resolution='c')    #"lon_0 must be between -360.000000 and 720.000000 degrees"; BUT when lon_0>0, the map reverses horizantally
                    fig.set_size_inches(8,2)
        
        
        
        m.drawcoastlines(color='dimgrey',linewidth=0.2)
        m.drawmapboundary(linewidth=0.1)
        m.fillcontinents(color='bisque',lake_color='lightcyan',alpha=0.3)            X,Y=m(*np.meshgrid(lon_shft[1:],lat))
        Zmax=np.amax(Z_shft)
        Zmin=np.amax(Z_shft)*-1
        levels=np.linspace(Zmin,Zmax,18)
        np.delete(levels,9)
        ct=m.contourf(X,Y,Z_shft[:,1:],cmap="Blues_r",levels=levels)    #bwr
            
        # plot wind vectors   
        u_shft,lon_shft=shiftgrid(lon_pjst,u,lon,start=False)
        v_shft,lon_shft=shiftgrid(lon_pjst,v,lon,start=False)
        
        if int(mm)<10:
            lat_ind=np.concatenate((np.arange(4,64,8),np.flip(np.arange(len(lat)-5,63,-8))),axis=None)    # make the resolution coarser for wind vector plotting
            lon_ind=np.arange(0,len(lon),8)
        if int(mm)>=10:    # & int(mm)<=18
            lat_ind=np.concatenate((np.arange(2,64,4),np.flip(np.arange(len(lat)-3,63,-4))),axis=None)    # make the resolution coarser for wind vector plotting
            lon_ind=np.arange(0,len(lon),6)
        latc=lat[lat_ind]
        lonc_shft=lon_shft[lon_ind]
        uc_shft=u_shft[lat_ind,:][:,lon_ind]
        vc_shft=v_shft[lat_ind,:][:,lon_ind]
        
        # normalization
        uc_norm=uc_shft/np.max(wind_spd)    #uc_shft
        vc_norm=vc_shft/np.max(wind_spd)    #uc_shft
        Xc,Yc=m(*np.meshgrid(lonc_shft,latc))
               
        
        if mode=="eig":
            if int(mm)<10:
                q=m.quiver(Xc,Yc,uc_norm,vc_norm,color='k',scale=30,width=0.002,ax=ax)    #,scale=wind_scl[int(mm)],scale_units=,scale=8*wscl eig: 60 rot,n1: 30 rot,n2: 20 ###30
            if int(mm)>=10 and int(mm)<18:
                q=m.quiver(Xc,Yc,uc_norm,vc_norm,color='k',scale=35,width=0.0015,ax=ax)   # 90 60 40 ###45
    #            m.ax.set_ylim(-45,45)
            if int(mm)>=18:
                q=m.quiver(Xc,Yc,uc_norm,vc_norm,color='k',scale=40,width=0.0015,ax=ax)   # 110 90 60 ###55
                
                
        if mode=="wig":
            if int(nn)>=3 and int(mm)<=4:
                q=m.quiver(Xc,Yc,uc_norm,vc_norm,color='k',scale=15,width=0.002,ax=ax)    #,scale=wind_scl[int(mm)],scale_units=,scale=8*wscl eig: 60 rot,n1: 30 rot,n2: 20
            else:
                if int(mm)<10:
                    q=m.quiver(Xc,Yc,uc_norm,vc_norm,color='k',scale=40,width=0.002,ax=ax)    #,scale=wind_scl[int(mm)],scale_units=,scale=8*wscl eig: 60 rot,n1: 30 rot,n2: 20
                if int(mm)>=10 and int(mm)<18:
                    q=m.quiver(Xc,Yc,uc_norm,vc_norm,color='k',scale=35,width=0.0015,ax=ax)   # 90 60 40 ###45
                if int(mm)>=18:
                    q=m.quiver(Xc,Yc,uc_norm,vc_norm,color='k',scale=40,width=0.0015,ax=ax)   # 110 90 60 ###55
                
        if mode=="rot" and int(nn)==1:
            if int(mm)<10:
                q=m.quiver(Xc,Yc,uc_norm,vc_norm,color='k',scale=15,width=0.002,ax=ax)    #,scale=wind_scl[int(mm)],scale_units=,scale=8*wscl eig: 60 rot,n1: 30 rot,n2: 20
            if int(mm)>=10 and int(mm)<18:
                q=m.quiver(Xc,Yc,uc_norm,vc_norm,color='k',scale=30,width=0.0015,ax=ax)   # 90 60 40   ###30
            if int(mm)>=18:
                q=m.quiver(Xc,Yc,uc_norm,vc_norm,color='k',scale=40,width=0.0015,ax=ax)   # 110 90 60  ###45
                
        if mode=="rot" and int(nn)!=1:
            if int(mm)<10:
                q=m.quiver(Xc,Yc,uc_norm,vc_norm,color='k',scale=15,width=0.002,ax=ax)    #,scale=wind_scl[int(mm)],scale_units=,scale=8*wscl eig: 60 rot,n1: 30 rot,n2: 20
            if int(mm)>=10 and int(mm)<18:
                q=m.quiver(Xc,Yc,uc_norm,vc_norm,color='k',scale=40,width=0.0015,ax=ax)   # 90 60 40 ###30
            if int(mm)>=18:
                q=m.quiver(Xc,Yc,uc_norm,vc_norm,color='k',scale=40,width=0.0015,ax=ax)   # 110 90 60 ###30
                
                
                
        if int(mm)==1:
            if mode=="rot" and int(kk)==1 and int(nn)==1:
                ax.set_title("MRG, k=1, D="+ed_mm,pad=5,fontsize=15)
                plt.savefig("./draw_x3/Hough/png7_3/HoughHarmonic_MRG_k1_D10km.png",dpi=300)
            elif mode=="eig" and int(kk)==1 and int(nn)==1: 
                ax.set_title("Kelvin wave, k=1, D="+ed_mm,pad=5,fontsize=15)
                plt.savefig("./draw_x3/Hough/png7_3/HoughHarmonic_Kelvin_k1_D10km.png",dpi=300)
            else:
                ax.set_title(mode_string+", k=1, n="+str(int(nn)-1)+", D="+ed_mm,pad=5,fontsize=15)
                plt.savefig("./draw_x3/Hough/png7_3/HoughHarmonic_"+mode+"_k1_n"+str(int(nn)-1)+"_D10km.png",dpi=300)    #+f[:-3]    #The PostScript backend does not support transparency; partially transparent artists will be rendered opaque.
        
        else:
            if mode=="rot" and int(kk)==1 and int(nn)==1:
                ax.set_title("MRG, k=1, D="+ed_mm,pad=5,fontsize=15)
                plt.savefig("./draw_x3/Hough/png7_3/HoughHarmonic_MRG_k1_D"+ed_mm[:-2]+"m.png",dpi=300)
            elif mode=="eig" and int(kk)==1 and int(nn)==1: 
                ax.set_title("Kelvin wave, k=1, D="+ed_mm,pad=5,fontsize=15)
                plt.savefig("./draw_x3/Hough/png7_3/HoughHarmonic_Kelvin_k1_D"+ed_mm[:-2]+"m.png",dpi=300)
            else:
                ax.set_title(mode_string+", k=1, n="+str(int(nn)-1)+", D="+ed_mm,pad=5,fontsize=15)
                plt.savefig("./draw_x3/Hough/png7_3/HoughHarmonic_"+mode+"_k1_n"+str(int(nn)-1)+"_D"+ed_mm[:-2]+"m.png",dpi=300)    #+f[:-3]    #The PostScript backend does not support transparency; partially transparent artists will be rendered opaque.
        
        
        plt.show()

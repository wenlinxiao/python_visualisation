#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plotting waves (wind speed contourf) on the spherical projection, e.g. Rossby waves
input: 3D wind field (u,v) of waves, from MODES
output: multiple figures of wind
time step: 12h

Created on Tue Nov 10 02:52:41 2020
@author: u300924
"""

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid
import os
import netCDF4 as nc
from scipy.interpolate import interp2d


if __name__=="__main__":
    pressures=[200,150,100,50,30,20,10,5,1]   # list of levels
    date="20200729"
    level_ind=0
    mode="rot"    # can also be used for ig mode, Kelvin waves...
    mode_string="Rossby"    # to be used in the figure title
    levels=np.arange(10,100,5)    # contourf levels
    
    path='../../../../../../scratch/uni/u234/u236001/MODES/testdata/tmp/'
    inp_fns=sorted(os.listdir(path))    # get the file names in the aim directory 
    inp_fns=inp_fns[15:29]    # only keep the data files
    
    # plot
    i=-1    # for shifting the projection center over time lags
    for f in inp_fns:    #[0:1]
        i=i+1
        date=f[-15:-7]
        t=f[-6:-3]
        print("time:",date,t,"h")
        
        # read data from nc file
        inp=nc.Dataset(path+f,format="NETCDF4")
        u=inp.variables["u"][level_ind+1,:,:]    # +1 because skip the interpolated 850 hPa (level==1)
        v=inp.variables["v"][level_ind+1,:,:]    # [level_p,lat,lon] [10,256,512]
        lon=inp.variables["lon"][:]
        lat=inp.variables["lat"][:]
        wind_spd=np.sqrt(u**2+v**2)
        

        fig,ax=plt.subplots()
        # plot contourf of wind speed
        m = Basemap(projection='ortho',lon_0=0-i*5,lat_0=20,resolution='l',round=True)
        m.drawcoastlines(color='dimgrey',linewidth=0.2)
        m.drawmapboundary(fill_color='lightcyan')
        m.fillcontinents(color='bisque',lake_color='lightcyan',alpha=0.5)
        m.drawparallels(np.arange(-90.,91.,30.),color='grey',linewidth=0.5)    
        m.drawmeridians(np.arange(0.,360.,30.),color='grey',linewidth=0.5)
        X,Y=m(*np.meshgrid(lon,lat))
        ct=m.contourf(X,Y,wind_spd,cmap="hot_r",levels=levels)    #
        cb1=fig.colorbar(ct,pad=0.1)
        cb1.set_label('Wind speed \nm/s',size=10)
        cb1.ax.tick_params(labelsize=10)
        
        # plot wind vectors   
        lat_ind=np.arange(0,len(lat),8)    # make the resolution coarser for wind vector plotting
        lon_ind=np.arange(0,len(lon),8)
        latc=lat[lat_ind]
        lonc=lon[lon_ind]
        uc=u[lat_ind,:][:,lon_ind]
        vc=v[lat_ind,:][:,lon_ind]
        wind_spdc=wind_spd[lat_ind,:][:,lon_ind]
        Xc,Yc=m(*np.meshgrid(lonc,latc))
        q=m.quiver(Xc,Yc,uc,vc,color='k',scale=1000,ax=ax)
        
        
        ax.set_title(mode_string+": "+date+" "+t+" UTC",pad=25)
        plt.savefig("./gif/PolarJet_"+mode+"_"+str(pressures[level_ind])+"hPa_"+date+"_"+t+".png",dpi=300)
        plt.show()
            
            
            
            

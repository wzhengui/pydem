#!/usr/bin/env python3
'''
pre-process dem data, and save them. modify this code according to the dem data
'''
from pylib import *
import tifffile as tiff
close("all")


#input
sdir='JP'
fname='japan_dem_with_ocean.tif'
sname='jdem'
xll0=122; yll0=46; dxy=0.000277777784503963722;
subdomain_size=3e8
depth_limit=[1,1e4]
nodata=-99999

#read dem
if fname.endswith('.tif'):
    dem=tiff.imread('./{}/{}'.format(sdir,fname))
    #get deminfo
    ym,xm=dem.shape; nsize=ym*xm;  #xi=xll+dxy*arange(xm); yi=yll-dxy*arange(ym)

#divide domain
nsub=int(max(nsize/subdomain_size,1)); nx=int(max(sqrt(nsub),1)); ny=int(max(nsub/nx,1))
dy0=int(ym/ny); dx0=int(xm/nx)

#save data
for i in arange(ny):
    for k in arange(nx):
        #get index for subdomain
        iy=i*dy0; dy=dy0; yll=yll0-iy*dxy; ix=k*dx0; dx=dx0; xll=xll0+ix*dxy
        if i==(ny-1): dy=ym-(i*dy0)
        if k==(nx-1): dx=xm-(k*dx0)
        iy2=iy+dy; ix2=ix+dx

        #save data
        S=npz_data(); S.info=npz_data()
        S.dem=dem[iy:iy2,ix:ix2]
        fpn=(S.dem<depth_limit[0])|(S.dem>depth_limit[1]); S.dem[fpn]=nodata
        S.info.ds=[dy,dx]; S.info.yll=yll; S.info.xll=xll; S.info.dxy=dxy; S.info.nodata=nodata
        snamei='{}_{}'.format(sname,i*nx+k)
        print('save {}'.format(snamei))
        save_npz('./{}/{}'.format(sdir,snamei),S)


#!/usr/bin/env python
# coding: utf-8

# # SWE Gridded Products

# In this notebook, we illustrate each of the four commonly used SWE gridded products to demonstrate differences in gridded estimates differ across products.

# In[1]:


import math
import datetime
import xarray as xr
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt


# In[2]:


# Set date of interest
date = datetime.datetime(2003, 4, 1)


# In[3]:


era5_dir = '/global/cfs/cdirs/dasrepo/yum/swe/gridded_products/data/ERA5/'
era5_file = 'ERA5_snw_1979-2014.nc'


# In[4]:


livneh_dir = '/global/cfs/cdirs/dasrepo/yum/swe/gridded_products/data/L15/'
livneh_file = 'L15_Fluxes_1985-2005.nc'


# In[5]:


# snsr_dir = '/global/cfs/cdirs/dasrepo/yum/swe/gridded_products/data/SNSR/'    # original UCLA SR that was only for the Sierra Nevada
# snsr_file = 'SN_SWE_WY2003.h5'
ucla_dir = '/global/cfs/cdirs/dasrepo/yum/swe/gridded_products/data/Margulis/Margulis/'
ucla_file = 'WUS_UCLA_SR_v01_ALL_0_agg_16_WY2002_SD_SWE_SCA_POST.nc'


# In[6]:


ua_dir = '/global/cfs/cdirs/dasrepo/yum/swe/gridded_products/data/UofASWE/'
ua_file = '4km_SWE_Depth_WY2003_v01.nc'


# In[7]:


save_dir = '/global/cfs/cdirs/dasrepo/yum/swe/uaswe-error-decomposition/figures/'


# ## Load and Plot Products in Default Spatial Extent

# ### ERA5

# In[8]:


# Load ERA5 data
era5_ds = xr.open_dataset(era5_dir+era5_file)


# In[9]:


# extract SWE variable from xr dataset
era5_da = era5_ds.snw


# In[10]:


# subset era5_da to desired date
era5_date = era5_da.sel(time=slice(date, date + datetime.timedelta(days=1)))


# In[11]:


# convert longitude from 0-360 to -180-180 coordinate system 
# https://stackoverflow.com/questions/53121089/regridding-coordinates-with-python-xarray
era5_date = era5_date.assign_coords(longitude=(((era5_date.longitude + 180) % 360) - 180)).sortby('longitude')


# In[12]:


# convert swe from meters to mm
era5_date = era5_date * 1000


# In[13]:


# replace 0mm SWE with NAN
era5_date.values[era5_date.values == 0] = np.nan


# In[14]:


# plot
# era5_date.plot(cmap='cool', vmax=2000)
# plt.savefig(save_dir+'era5.png', dpi=300)


# ### Livneh

# In[15]:


# load livneh data
livneh_ds = xr.open_dataset(livneh_dir+livneh_file)


# In[ ]:


# extract SWE variable
livneh_da = livneh_ds.SWE


# In[ ]:


# subset livneh_da to desired date
livneh_date = livneh_da.sel(time=slice(date, date + datetime.timedelta(days=1)))


# In[ ]:


# plot
livneh_date[0,:,:].plot(cmap='cool', vmax=2000)
# plt.savefig(save_dir+'livneh.png', dpi=300)


# ### UCLA SR

# In[ ]:


# load UCLA_SR data
ucla_ds = xr.open_dataset(ucla_dir+ucla_file)


# In[ ]:


# extract SWE variable
ucla_da = ucla_ds.SWE_Post


# In[ ]:


# subset UCLA_SR to desired date
ucla_date = ucla_da.sel(time=slice(date, date + datetime.timedelta(days=1)))


# In[ ]:


# Convert from meters to milimeters
ucla_date = ucla_date * 1000


# In[ ]:


ucla_wus = ucla_date[0,2,:,:] * 1000
ucla_wus.plot(cmap='cool')


# ### UA SWE

# In[ ]:


# load UA data
ua_ds = xr.open_dataset(ua_dir+ua_file)


# In[ ]:


# extract SWE variable 
ua_da = ua_ds.SWE


# In[ ]:


# subset ua_da to desired date
ua_date = ua_da.sel(time=slice(date, date + datetime.timedelta(days=1)))


# In[ ]:


# plot
ua_date[0,:,:].plot(cmap='cool', vmax=2000)
# plt.savefig(save_dir+'ua.png', dpi=300)


# ## Plot Products for the UCRB Study Area

# Load spatial extent of UCRB

# In[ ]:


# Get shapefile for Upper Colorado Riber Basin (UCRB)
uc_shp = "../data/Upper_Colorado_River_Basin_Boundary/Upper_Colorado_River_Basin_Boundary.shp"

# Read UCRB shapefile
gm_poly_gdf = gpd.read_file(uc_shp, encoding="utf-8")

# Get bounds of UCRB
gm_poly_geom = gm_poly_gdf.iloc[0].geometry


# Get lat/lon boundaries of UCRB

# In[ ]:


yy, xx = gm_poly_geom.exterior.coords.xy
yy, xx = np.array(yy), np.array(xx)
# Get UCRB bounds
xmin, xmax = float(np.min(xx)), float(np.max(xx))
ymin, ymax = float(np.min(yy)), float(np.max(yy))


# ### Subset Products to UCRB

# In[ ]:


# subset ERA5 to UCRB long/lat
era5_sub = era5_date.sel(latitude=slice(xmax,xmin), longitude=slice(ymin,ymax))
# subset ERA5 to UCRB bounds (i.e. when UCRB is NAN, make ERA5 NAN)
era5_ucrb = era5_sub[0,:,:]


# In[ ]:


# subset Livneh to UCRB long/lat
livneh_sub = livneh_date.sel(lat=slice(xmin,xmax), lon=slice(ymin,ymax))
# subset Livneh to UCRB bounds (i.e. when UCRB is NAN, make Livneh NAN)
livneh_ucrb = livneh_sub[0,:,:]


# In[ ]:


# subset UCLA_SR to UCRB long/lat
ucla_sub = ucla_date.sel(Latitude=slice(xmax,xmin), Longitude=slice(ymin,ymax))
# subset UCLA_SR to UCRB bounds (i.e. when UCRB is NAN, make UA NAN)
ucla_ucrb = ucla_sub[0,2,:,:]   # Posterior median
# ucla_ucrb = ucla_sub[0,:,:]    # Posterior mean


# In[ ]:


# subset UA to UCRB long/lat
ua_sub = ua_date.sel(lat=slice(xmin,xmax), lon=slice(ymin,ymax))
# subset UA to UCRB bounds (i.e. when UCRB is NAN, make UA NAN)
ua_ucrb = ua_sub[0,:,:]


# ### Apply Log Scaling to SWE Gridded Values

# In[ ]:


log_era5_ucrb = np.log10(era5_ucrb)
log_era5_ucrb = log_era5_ucrb.where(log_era5_ucrb != -math.inf)


# In[ ]:


log_livneh_ucrb = np.log10(livneh_ucrb)
log_livneh_ucrb = log_livneh_ucrb.where(log_livneh_ucrb != -math.inf)


# In[ ]:


log_ucla_ucrb = np.log10(ucla_ucrb)
log_ucla_ucrb = log_ucla_ucrb.where(log_ucla_ucrb != -math.inf)


# In[ ]:


log_ua_ucrb = np.log10(ua_ucrb)
log_ua_ucrb = log_ua_ucrb.where(log_ua_ucrb != -math.inf)


# ### Plot Products for UCRB

# #### Indivially

# In[ ]:


# plot era5 clipped to ucrb region
f, ax = plt.subplots(figsize=(10,6))
log_era5_ucrb.plot(ax=ax, cmap='coolwarm', vmin=-4, vmax=4, cbar_kwargs={'label':'log SWE (mm)'})
# era5_ucrb.plot(ax=ax, cmap='Blues', vmax=1800) #vmax=3200)
gm_poly_gdf.plot(ax=ax, facecolor="none", edgecolor='black', lw=1.5)
plt.title('ERA5 (Clipped to UCRB)')
plt.xlim(-112.5,-105.5)
plt.ylim(35.5,43.5)
# plt.savefig(save_dir+'era5_ucrb.png', dpi=300)
# plt.savefig(save_dir+'log_era5_ucrb.png', dpi=300)


# In[ ]:


# plot livneh clipped to ucrb region
f, ax = plt.subplots(figsize=(10,6))
log_livneh_ucrb.plot(ax=ax, cmap='coolwarm', vmin=-4, vmax=4, cbar_kwargs={'label':'log SWE (mm)'})
# livneh_ucrb.plot(ax=ax, cmap='Blues', vmax=1800)#vmax=3200)
gm_poly_gdf.plot(ax=ax, facecolor="none", edgecolor='black', lw=1.5)
plt.title('Livneh (Clipped to UCRB)')
plt.xlim(-112.5,-105.5)
plt.ylim(35.5,43.5)
# plt.savefig(save_dir+'livneh_ucrb.png', dpi=300) 
# plt.savefig(save_dir+'log_livneh_ucrb.png', dpi=300)


# In[ ]:


# plot wus_ucla_sr clipped to ucrb region
f, ax = plt.subplots(figsize=(10,6))
log_ucla_ucrb.plot(ax=ax, cmap='coolwarm', vmin=-4, vmax=4, cbar_kwargs={'label':'log SWE (mm)'})#, vmax=1800)#vmax=3200)
# ucla_ucrb.plot(ax=ax, cmap='Blues', vmax=1800)#vmax=3200)
gm_poly_gdf.plot(ax=ax, facecolor="none", edgecolor='black', lw=1.5)
plt.title('UCLA SR (Clipped to UCRB)')
plt.xlim(-112.5,-105.5)
plt.ylim(35.5,43.5)
# plt.savefig(save_dir+'snsr_ucrb.png', dpi=300)
# plt.savefig(save_dir+'log_snsr_ucrb.png', dpi=300)


# In[ ]:


# plot ua clipped to ucrb region
f, ax = plt.subplots(figsize=(10,6))
log_ua_ucrb.plot(ax=ax, cmap='coolwarm', vmin=-4, vmax=4, cbar_kwargs={'label':'log SWE (mm)'})
# ua_ucrb.plot(ax=ax, cmap='Blues', vmax=1800)#vmax=3200)
gm_poly_gdf.plot(ax=ax, facecolor="none", edgecolor='black', lw=1.5)
plt.title('UA (Clipped to UCRB)')
plt.xlim(-112.5,-105.5)
plt.ylim(35.5,43.5)
# plt.savefig(save_dir+'ua_ucrb.png', dpi=300)
# plt.savefig(save_dir+'log_ua_ucrb.png', dpi=300)


# #### All Together

# In[ ]:


# plot era5 clipped to ucrb region
f, axes = plt.subplots(nrows=2, ncols=2, figsize=(10,12))
log_era5_ucrb.plot(ax=axes[0,0], cmap='coolwarm', vmin=-4, vmax=4, cbar_kwargs={'label':'log SWE (mm)'})
log_livneh_ucrb.plot(ax=axes[0,1], cmap='coolwarm', vmin=-4, vmax=4, cbar_kwargs={'label':'log SWE (mm)'})
log_ua_ucrb.plot(ax=axes[1,0], cmap='coolwarm', vmin=-4, vmax=4, cbar_kwargs={'label':'log SWE (mm)'})
log_ucla_ucrb.plot(ax=axes[1,1], cmap='coolwarm', vmin=-4, vmax=4, cbar_kwargs={'label':'log SWE (mm)'})
# era5_ucrb.plot(ax=axes[0,0], cmap='Blues', vmax=1800) #vmax=3200)
# livneh_ucrb.plot(ax=axes[0,1], cmap='Blues', vmax=1800) #vmax=3200)
# ua_ucrb.plot(ax=axes[1,0], cmap='Blues', vmax=1800) #vmax=3200)
# ucla_ucrb.plot(ax=axes[1,1], cmap='Blues', vmax=1800) #vmax=3200)
gm_poly_gdf.plot(ax=axes[0,0], facecolor="none", edgecolor='black', lw=1.5)
gm_poly_gdf.plot(ax=axes[0,1], facecolor="none", edgecolor='black', lw=1.5)
gm_poly_gdf.plot(ax=axes[1,0], facecolor="none", edgecolor='black', lw=1.5)
gm_poly_gdf.plot(ax=axes[1,1], facecolor="none", edgecolor='black', lw=1.5)
axes[0,0].set_title('ERA5 (Clipped to UCRB)')
axes[0,1].set_title('Livneh (Clipped to UCRB)')
axes[1,0].set_title('UA SWE (Clipped to UCRB)')
axes[1,1].set_title('UCLA SR (Clipped to UCRB)')
plt.xlim(-112.5,-105.5)
plt.ylim(35.5,43.5)
plt.tight_layout()
# plt.savefig(save_dir+'swe_products.png', dpi=300)
# plt.savefig(save_dir+'swe_products.png', dpi=300)


# In[ ]:





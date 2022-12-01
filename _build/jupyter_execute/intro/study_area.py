#!/usr/bin/env python
# coding: utf-8

# # Study Area

# ## Upper Colorado River Basin (UCRB)

# In[1]:


import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import contextily as ctx
import matplotlib.pyplot as plt


# In[ ]:


save_dir = '/global/cfs/cdirs/dasrepo/yum/swe/uaswe-error-decomposition/figures/'


# ### Load UCRB shapefile to get region bounds

# In[ ]:


# Get shapefile for Upper Colorado Riber Basin (UCRB)
uc_shp = "/global/cscratch1/sd/yum/swe/Upper_Colorado_River_Basin_Boundary/Upper_Colorado_River_Basin_Boundary.shp"

# Read UCRB shapefile
gm_poly_gdf = gpd.read_file(uc_shp, encoding="utf-8")

# Get bounds of UCRB
gm_poly_geom = gm_poly_gdf.iloc[0].geometry

# Determine sites in UCRB
# sites_idx = sites_gdf.intersects(gm_poly_geom)

# Subset df to sites in UCRB
# gm_snotel_sites = sites_gdf.loc[sites_idx]


# ### Plot UCRB boundary over map

# In[ ]:


f, ax = plt.subplots(figsize=(8,8))
gm_poly_gdf.plot(ax=ax, facecolor="none", edgecolor='black', lw=2.5)
# gm_snotel_sites.plot(ax=ax, column='elevation_m', markersize=100, edgecolor='k', cmap='inferno', legend=True, legend_kwds={'label':'Elevation (m)'})
ctx.add_basemap(ax=ax, crs=gm_poly_gdf.crs,source=ctx.providers.Stamen.Terrain)    # add/plot basemap
# ctx.add_basemap(ax=ax, crs=gm_snotel_sites.crs, source=ctx.providers.Stamen.Terrain)
ax.set_title('Study Area: \n Upper Colorado River Basin (UCRB)', fontsize=21)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()
plt.savefig(save_dir+'ucrb.png', dpi=300)


# In[ ]:





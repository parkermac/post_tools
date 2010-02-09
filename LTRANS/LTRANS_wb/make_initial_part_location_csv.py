#!/usr/bin/env python
import os
import netCDF4 as nc
from optparse import OptionParser
import numpy as np

def dim2vals(dims,ncf):
	dim_vals = []
	for i in range(len(dims)):
		dim_vals.append(len(ncf.dimensions[dims[i]]))
	return tuple(dim_vals)

parser = OptionParser()


(options, args) = parser.parse_args()

if len(args) == 2:
	ncfile_name = args[0]
	init_file_name = args[1]
elif len(args) == 1:
	ncfile_name = args[0]
	init_file_name = 'initial_part_location.csv'
else:
	ncfile_name = 'test.nc'
	init_file_name = 'initial_part_location.csv'

#print(ncfile_name)
#print(init_file_name)

ncfile = nc.Dataset(ncfile_name,'r')

lat = ncfile.variables['lat_rho'][:]
lon = ncfile.variables['lon_rho'][:]
#mask = ncfile.variables['wetdry_mask_rho'][:]
mask = ncfile.variables['mask_rho'][:]

#print(mask.shape)
for i in range(mask.shape[0]-20):
	for j in range(mask.shape[1]-20):
		if mask[i+10,j+10] == 1:
			print('%f, %f, 10.0' %(lon[i+10,j+10],lat[i+10,j+10]))
#for i in range(mask.size):
#	if mask.flat[i] == 1:
#		print('%f, %f, 10.0' %(lon.flat[i],lat.flat[i]))

# salt_var = ncfile.variables['salt']
# salt_dims = salt_var.dimensions
# 
# salt_dims_vals = dim2vals(salt_dims, ncfile)
# 
# nitro_var = ncfile.createVariable('nitrogen','f4', salt_dims)
# nitro_var.units = 'g/l'
# 
# nitro_var[:] = 0.01*np.ones(salt_dims_vals)
# 
# nitro_var.coordinates = salt_var.coordinates
# nitro_var.long_name = 'nitrogen concentration'
# nitro_var.time = salt_var.time
# 
# print(nitro_var[:])




ncfile.close()

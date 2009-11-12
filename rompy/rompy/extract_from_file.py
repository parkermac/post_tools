import os

import netCDF4 as nc
import numpy as np

import load_grid
import plot_surface
import utils

def extract_from_file(file='',varname='zeta',extraction_type='full',**kwargs):
	# Check the netCDF file for existence and if the variable is in it
	if not os.path.exists(file):
		raise IOError('File %s could not be located on the filesystem' %file)
	ncf = nc.Dataset(file,mode='r')
	
	if varname not in ncf.variables:
		raise IOError('File %s does not have a variable named %s' % (file, varname))
	
	# start getting data
	ncvar = ncf.variables[varname]
	dims = ncvar.dimensions
	ndims = len(dims)
	print('dims: %s' % str(dims))
	shape = ncvar.shape
	print('shape: %s' % str(shape))
	
	if not ndims == 3 and not ndims == 4:
		raise TypeError('ndims is neither 3 nor 4')
	if not dims[0] == 'ocean_time':
		raise TypeError('first dimension is not ocean_time')
	if not shape[0] == 1:
		raise TypeError('first dimension is not of length one')
	
	grid = load_grid.load_grid(file)
	coords = {}
	if ndims == 3:
		if dims[1] == 'eta_rho' and dims[2] == 'xi_rho':
#			print('G.lat %s' % str(grid['lat']))
#			print('G.lon %s' % str(grid['lon']))
			y2 = grid['lat']
			x2 = grid['lon']
			mask2 = grid['mask']
		elif dims[1] == 'eta_u' and dims[2] == 'xi_u':
#			print('G.latu %s' % str(grid['lat_u']))
#			print('G.lonv %s' % str(grid['lon_u']))
			y2 = grid['latu']
			x2 = grid['lonu']
			mask2 = grid['masku']
		elif dims[1] == 'eta_v' and dims[2] == 'xi_v':
#			print('G.latv %s' % str(grid['lat_v']))
#			print('G.lonv %s' % str(grid['lon_v']))
			y2 = grid['latv']
			x2 = grid['lonv']
			mask2 = grid['maskv']
		else:
			raise TypeError('Unable to determine which gird to use')
		
		if extraction_type == 'full' or extraction_type == 'surface':
			x2 = x2[:]
			y2 = y2[:]
			mask2 = mask2[:]
			
#			data = np.squeeze(ncvar[:])
#			data[mask2==0] = np.NaN
			data = np.ma.array(np.squeeze(ncvar[:]),mask=(mask2==0))
			
			coords['ym'] = y2
			coords['xm'] = x2
		
		if (extraction_type == 'profile' or extraction_type == 'profiles' or
		    extraction_type == 'point' or extraction_type == 'points') and \
		   (kwargs.has_key('y') and kwargs.has_key('x')):
			print('In profiles')
			xm = np.array(kwargs['x'])
			ym = np.array(kwargs['y'])
			if not (xm.ndim == ym.ndim) or not (xm.shape == ym.shape):
				if xm.size == 1:
					xm = np.repmat(xm,ym.shape)
				elif ym.size == 1:
					ym = np.repmat(ym,xm.shape)
				else:
					raise RuntimeError('The x and y chosen to extract a point or profile on this 2D variable are incompatible.')

			x2 = x2[:]
			y2 = y2[:]
			mask2 = mask2[:]
			data2 = np.ma.array(np.squeeze(ncvar[:]),mask=(mask2==0))
#			data = utils.interp_2d(lat=y2,lon=x2,data=data2,lati=ym,loni=xm)
			mask = utils.interp_2d_xy(y=y2,x=x2,data=mask2,yi=ym,xi=xm)
			data = utils.interp_2d_xy(y=y2,x=x2,data=data2,yi=ym,xi=xm)
			data = np.ma.array(data,mask=(mask < 1))
			coords['ym'] = ym
			coords['xm'] = xm

	if ndims == 4:
		Hu = utils.interp_2d(lat=y2,lon=x2,data=grid['H'][:],lati=grid['latu'][:],loni=grid['lonu'])
	return (data,coords)
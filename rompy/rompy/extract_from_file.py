import os

import netCDF4 as nc

import load_grid

def extract_from_file(file='',varname='zeta',extractionType='full',**kwargs):
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
	
	if ndims == 3:
		if dims[1] == 'eta_rho' and dims[2] == 'xi_rho':
			print('G.lat %s' % str(grid['lat']))
			print('G.lon %s' % str(grid['lon']))
			y2 = grid['lat']
			x2 = grid['lon']
			mask2 = grid['mask']
		elif dims[1] == 'eta_u' and dims[2] == 'xi_u':
			print('G.latu %s' % str(grid['lat_u']))
			print('G.lonv %s' % str(grid['lon_u']))
			y2 = grid['latu']
			x2 = grid['lonu']
			mask2 = grid['masku']
		elif dims[1] == 'eta_v' and dims[2] == 'xi_v':
			print('G.latv %s' % str(grid['lat_v']))
			print('G.lonv %s' % str(grid['lon_v']))
			y2 = grid['latv']
			x2 = grid['lonv']
			mask2 = grid['maskv']
		else:
			raise TypeError('Unable to determine which gird to use')
		
		
		
	return
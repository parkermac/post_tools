#!/usr/bin/env python
from rompy import rompy, plot_utils
import numpy as np

map1 = True
map2 = True
map3 = True
map4 = True

if map1:
	(data, coords) = rompy.extract('ocean_his_1000.nc',varname='zeta')
	plot_utils.plot_surface(coords['xm'],coords['ym'],data)
	plot_utils.plot_map(coords['xm'],coords['ym'],data,filename='/Users/lederer/tmp/rompy.map.png')
	
	del(data)
	del(coords)

if map2:
	# full domain
	#x = np.linspace(-127.,-122.,100)
	#y = np.linspace(44.,50.,100)
	
	#puget sound area
	x = np.linspace(-123.,-122.,500)
	y = np.linspace(47.,48.,500)
	
	(data, coords) = rompy.extract('ocean_his_1000.nc', x=x, y=y)
	plot_utils.plot_map(coords['xm'],coords['ym'],data,filename='/Users/lederer/tmp/rompy.map2.png',resolution='h')

if map3:
	(data, coords) = rompy.extract('ocean_his_1000.nc',varname='v',extraction_type='full')
	print(data.shape)
	for key in coords:
		print(key, coords[key].shape)
		
	plot_utils.plot_profile(data[:,20,20],coords['zm'][:,20,20],filename='/Users/lederer/tmp/rompy.profile.png')

if map4:
	(data, coords) = rompy.extract('ocean_his_1000.nc',varname='salt',extraction_type='surface')
	plot_utils.plot_map(coords['xm'],coords['ym'],data,filename='/Users/lederer/tmp/rompy.map4.png',resolution='h')
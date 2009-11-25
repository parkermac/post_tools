#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure


from rompy import rompy, plot_utils, utils

map1 = False
map2 = False
map3 = False
map4 = False
map5 = False
map6 = False
map7 = False
map8 = False
map9 = False
map10 = True

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
	#x = np.linspace(-123.,-122.,500)
	#y = np.linspace(47.,48.,500)
	
	# hood canal
	x = np.linspace(-123.25,-122.5,400)
	y = np.linspace(47.33,48.0,400)
	(data, coords) = rompy.extract('ocean_his_1000.nc', varname='zeta',extraction_type='points', x=x, y=y)
	plot_utils.plot_map(coords['xm'],coords['ym'],data,filename='/Users/lederer/tmp/rompy.map2.png',resolution='h')
#	plot_utils.plot_surface(coords['xm'],coords['ym'],data,filename='/Users/lederer/tmp/rompy.map2.png')

if map3:
	(data, coords) = rompy.extract('ocean_his_1000.nc',varname='v',extraction_type='full')
	print(data.shape)
	for key in coords:
		print(key, coords[key].shape)
		
	plot_utils.plot_profile(data[:,20,20],coords['zm'][:,20,20],filename='/Users/lederer/tmp/rompy.profile.png')

if map4:
	(data, coords) = rompy.extract('ocean_his_1000.nc',varname='salt',extraction_type='surface')
	plot_utils.plot_map(coords['xm'],coords['ym'],data,filename='/Users/lederer/tmp/rompy.map4.png',resolution='h')
	
if map5:
#	middle of pacific
#	x = np.linspace(-126.0,-125.0,1001)
#	y = np.linspace(45.0,46.0,1001)

#	hood canal PRISM Cruise February 2009
	x,y = utils.hood_canal_xy()
	
	#cs = np.linspace(-0.96103753,-0.00143376,10)
	(data, coords) = rompy.extract('ocean_his_1000.nc',varname='salt',extraction_type='profile',x=x,y=y)#,cs=cs)
	fig = Figure(facecolor='white')
	ax = fig.add_subplot(111)
	
#	my_plot = ax.pcolormesh(np.arange(data.shape[1]),coords['zm'],data,clim=(0,35),colorbar=True)
	my_plot = ax.contourf(np.tile(np.arange(data.shape[1]),(coords['zm'].shape[0],1)),coords['zm'],data,100)
	my_plot2 = ax.contour(np.tile(np.arange(data.shape[1]),(coords['zm'].shape[0],1)),coords['zm'],data,100,linewidths=1,linestyle=None)
	ax.fill_between(np.arange(data.shape[1]),coords['zm'][0,:],ax.get_ylim()[0],color='grey')
	fig.colorbar(my_plot,ax=ax)
	ax.set_title('Hood Canal Salinity from a ROMS run')
	ax.set_ylabel('depth in meters')
	ax.set_xticks(np.arange(data.shape[1]))
	ax.set_xticklabels(utils.hood_canal_station_list())
	ax.set_xlabel('station ID')
	
	FigureCanvas(fig).print_png('/Users/lederer/tmp/rompy.map5.png')

if map6:
#	middle of pacific
#	x = np.linspace(-126.0,-125.0,1001)
#	y = np.linspace(45.0,46.0,1001)

#	hood canal PRISM Cruise February 2009
	x,y = utils.main_basin_xy()
	
	#cs = np.linspace(-0.96103753,-0.00143376,10)
	(data, coords) = rompy.extract('ocean_his_1000.nc',varname='salt',extraction_type='profile',x=x,y=y)#,cs=cs)
	fig = Figure(facecolor='white')
	ax = fig.add_subplot(111)
	
#	my_plot = ax.pcolormesh(np.arange(data.shape[1]),coords['zm'],data,clim=(0,35),colorbar=True)
	my_plot = ax.contourf(np.tile(np.arange(data.shape[1]),(coords['zm'].shape[0],1)),coords['zm'],data,100)
	my_plot2 = ax.contour(np.tile(np.arange(data.shape[1]),(coords['zm'].shape[0],1)),coords['zm'],data,100,linewidths=1,linestyle=None)
	ax.fill_between(np.arange(data.shape[1]),coords['zm'][0,:],ax.get_ylim()[0],color='grey')
	fig.colorbar(my_plot,ax=ax)
	ax.set_title('Main Basin Salinity from a ROMS run')
	ax.set_ylabel('depth in meters')
	ax.set_xticks(np.arange(data.shape[1]))
	ax.set_xticklabels(utils.main_basin_station_list())
	ax.set_xlabel('station ID')
	
	FigureCanvas(fig).print_png('/Users/lederer/tmp/rompy.map6.png')

if map7: # Main Basin
	n = 4
	x,y = utils.high_res_main_basin_xy(n=n)
	
	# Salinity
	(data, coords) = rompy.extract('ocean_his_1000.nc', varname='salt', extraction_type='profile', x=x, y=y)
	
	plot_utils.plot_mickett(coords=coords, data=data, varname='Salinity', region='Main Basin', filename='/Users/lederer/tmp/rompy.mickett_main_salt.png', n=n, x_axis_offset=utils.offset_region(coords), clim=[0,20,30,35],cmap='banas_hsv_cm')
	
	# Temperature
	(data, coords) = rompy.extract('ocean_his_1000.nc',varname='temp',extraction_type='profile',x=x,y=y)
	
	plot_utils.plot_mickett(coords=coords, data=data, varname='Temperature', region='Main Basin', filename='/Users/lederer/tmp/rompy.mickett_main_temp.png', n=n,  x_axis_offset=utils.offset_region(coords), clim=[6,8,11,18], cmap='banas_hsv_cm')

if map8: # Hood Canal
	n=3
	x,y = utils.high_res_hood_canal_xy(n=n)
	# Salinity
	(data, coords) = rompy.extract('ocean_his_1000.nc', varname='salt', extraction_type='profile', x=x, y=y)
	
	plot_utils.plot_mickett(coords=coords, data=data, varname='Salinity', region='Hood Canal', filename='/Users/lederer/tmp/rompy.mickett_hood_salt.png', n=n,  x_axis_offset=utils.offset_region(coords), clim=[10,25,30,35], cmap='banas_hsv_cm')
	
	# Temperature
	(data, coords) = rompy.extract('ocean_his_1000.nc', varname='temp', extraction_type='profile', x=x, y=y)
	
	plot_utils.plot_mickett(coords=coords, data=data, varname='Temperature', region='Hood Canal', filename='/Users/lederer/tmp/rompy.mickett_hood_temp.png', n=n, x_axis_offset=utils.offset_region(coords), clim=[6,18], cmap='banas_hsv_cm')

if map9: # velocity in Hood Canal
	n=1
	x,y = utils.high_res_hood_canal_xy(n=n)
	(u, coords) = rompy.extract('ocean_his_1000.nc',varname='u',extraction_type='profile',x=x,y=y)
	(v, coords) = rompy.extract('ocean_his_1000.nc',varname='v',extraction_type='profile',x=x,y=y)
	data = np.zeros(u.shape)

	for i in range(u.shape[1]):
		if i == u.shape[1]-1:
			x_vec = np.array([x[i] - x[i-1], y[i] - y[i-1]])
		else:
			x_vec = np.array([x[i+1] - x[i], y[i+1] - y[i]])
		for j in range(u.shape[0]):
			u_vec = np.array([u[j,i], v[j,i]])
			data[j,i] = np.dot(x_vec,u_vec)/(np.sqrt(np.dot(x_vec,x_vec)))
	
	plot_utils.plot_mickett(coords=coords,data=data,varname='U', region='Hood Canal', filename='/Users/lederer/tmp/rompy.mickett_hood_U.png', n=n, clim=[-2,2], x_axis_offset=utils.offset_region(coords),cmap='red_blue')

if map10: # velocity in Main Basin
	n=3
	x,y = utils.high_res_main_basin_xy(n=n)
	(u, coords) = rompy.extract('ocean_his_1000.nc',varname='u',extraction_type='profile',x=x,y=y)
	(v, coords) = rompy.extract('ocean_his_1000.nc',varname='v',extraction_type='profile',x=x,y=y)
	data = np.zeros(u.shape)

	for i in range(u.shape[1]):
		if i == u.shape[1]-1:
			x_vec = np.array([x[i] - x[i-1], y[i] - y[i-1]])
		else:
			x_vec = np.array([x[i+1] - x[i], y[i+1] - y[i]])
		for j in range(u.shape[0]):
			u_vec = np.array([u[j,i], v[j,i]])
			data[j,i] = np.dot(x_vec,u_vec)/(np.sqrt(np.dot(x_vec,x_vec)))
	
	plot_utils.plot_mickett(coords=coords,data=data,varname='U', region=' Main Basin', filename='/Users/lederer/tmp/rompy.mickett_main_U.png', n=n, clim=[-2,2], x_axis_offset=utils.offset_region(coords),cmap='red_blue')
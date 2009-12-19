#!/usr/bin/env python
import os
import glob
from optparse import OptionParser

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure


from rompy import rompy, plot_utils, utils, extract_utils

def surface_map(file,img_file=None,varname='salt',clim=None):
	(data, coords) = rompy.extract(file,varname=varname,extraction_type='surface')
#	plot_utils.plot_surface(coords['xm'],coords['ym'],data)
	
	
	title = '%s %s %s %s' % ( extract_utils.run_id(file), file, var_title_map[var], extract_utils.file_time(file).strftime(title_time_fmt) )
	plot_utils.plot_map(coords['xm'],coords['ym'],data,filename=img_file, clim=clim, title=title, resolution=whole_domain_coastline_res)

def main_basin_curtain(file,img_file,varname,n=4,clim=None): # Main Basin
	if var == 'U':
		main_basin_U_curtain(file,img_file,n,clim)
	else:
		x,y = utils.high_res_main_basin_xy(n=n)
		
		(data, coords) = rompy.extract(file, varname=varname, 	extraction_type='profile', x=x, y=y)
		
		title = 'Main Basin %s %s' % (var_title_map[var], extract_utils.file_time(file).strftime(title_time_fmt))
		
		plot_utils.plot_parker(coords=coords, data=data, varname=varname, 	region='Main Basin', filename=img_file, n=n, x_axis_offset=utils.offset_region(coords), clim=clim,cmap='banas_hsv_cm',labeled_contour_gap=2, title=title, resolution=inset_coastline_resolution)
	

def hood_canal_curtain(file,img_file,varname,n=1,clim=None): # Hood Canal
	if var == 'U':
		hood_canal_U_curtain(file,img_file,n,clim)
	else:
		x,y = utils.high_res_hood_canal_xy(n=n)
		
		(data, coords) = rompy.extract(file, varname=varname, extraction_type='profile', x=x, y=y)
		
		title = 'Hood Canal %s %s' % (var_title_map[var], extract_utils.file_time(file).isoformat())
		
		plot_utils.plot_parker(coords=coords, data=data, varname=varname, region='Hood Canal', filename=img_file, n=n,  x_axis_offset=utils.offset_region(coords), clim=clim, cmap='banas_hsv_cm',labeled_contour_gap=2, title=title, resolution=inset_coastline_resolution)

def hood_canal_U_curtain(file,img_file,n=1,clim=None): # velocity in Hood Canal
	x,y = utils.high_res_hood_canal_xy(n=n)
	(u, coords) = rompy.extract(file,varname='u',extraction_type='profile',x=x,y=y)
	(v, coords) = rompy.extract(file,varname='v',extraction_type='profile',x=x,y=y)
	data = np.zeros(u.shape)

	for i in range(u.shape[1]):
		if i == u.shape[1]-1:
			x_vec = np.array([x[i] - x[i-1], y[i] - y[i-1]])
		else:
			x_vec = np.array([x[i+1] - x[i], y[i+1] - y[i]])
		for j in range(u.shape[0]):
			u_vec = np.array([u[j,i], v[j,i]])
			data[j,i] = np.dot(x_vec,u_vec)/(np.sqrt(np.dot(x_vec,x_vec)))
	
	data = np.ma.array(data, mask=np.abs(data) > 100)
	
	title = 'Hood Canal %s %s' % (var_title_map[var], extract_utils.file_time(file).isoformat())
	
	hood_U_clim = (np.array(clim)/2.0).tolist()
	
	plot_utils.plot_parker(coords=coords,data=data,varname='U', region='Hood Canal', filename=img_file, n=n, clim=clim, x_axis_offset=utils.offset_region(coords), cmap='red_blue', title=title, resolution=inset_coastline_resolution)

def main_basin_U_curtain(file,img_file,n=1,clim=None): # velocity in Main Basin
	x,y = utils.high_res_main_basin_xy(n=n)
	(u, coords) = rompy.extract(file,varname='u',extraction_type='profile',x=x,y=y)
	(v, coords) = rompy.extract(file,varname='v',extraction_type='profile',x=x,y=y)
	data = np.zeros(u.shape)

	for i in range(u.shape[1]):
		if i == u.shape[1]-1:
			x_vec = np.array([x[i] - x[i-1], y[i] - y[i-1]])
		else:
			x_vec = np.array([x[i+1] - x[i], y[i+1] - y[i]])
		for j in range(u.shape[0]):
			u_vec = np.array([u[j,i], v[j,i]])
			data[j,i] = np.dot(x_vec,u_vec)/(np.sqrt(np.dot(x_vec,x_vec)))
	
	data = np.ma.array(data, mask=np.abs(data) > 100)
	
	title = 'Main Basin %s %s' % (var_title_map[var], extract_utils.file_time(file).strftime(title_time_fmt))
	
	plot_utils.plot_parker(coords=coords,data=data,varname='U', region=' Main Basin', filename=img_file, n=n, clim=clim, x_axis_offset=utils.offset_region(coords),cmap='red_blue', title=title, resolution=inset_coastline_resolution)

# begin actual code that runs.

parser = OptionParser()
(options, args) = parser.parse_args()

if args == []:
	fl = glob.glob('ocean_his*.nc')
	print(fl)
	file_list = [fl[0]]
else:
	file_list = args

img_dir = '/Users/lederer/Repositories/PSVS/rompy/image_sequence'

var_list = ['salt','temp','U']

var_title_map = {'salt':'Salinity','temp':'Temperature','U':'Velocity'}
title_time_fmt = '%Y-%m-%d %H:%M UTC'

clims = {'salt':[0, 21,33, 33], 'temp': [8, 20], 'U':[-2,2]}

inset_coastline_resolution = 'f'
whole_domain_coastline_res = 'f'

for file in file_list:
	ncf_index = file[:-3]
	print('ocean_his_%s' %ncf_index)
	if not os.path.exists(img_dir):	
		os.makedirs(img_dir)
	
	for var in var_list:
		hood_img_file = '%s/%s_hood_%s.png' %(img_dir, ncf_index,var)
		main_img_file = '%s/%s_main_%s.png' %(img_dir, ncf_index,var)
		surface_img_file = '%s/%s_surface_%s.png' % (img_dir, ncf_index, var)
		
		print('making hood canal %s' % var)
#		hood_canal_curtain(file, hood_img_file, var, n=8, clim=clims[var])
		print('making main basin %s' % var)
		main_basin_curtain(file, main_img_file, var, n=8, clim=clims[var])
		if not var == 'U':
			print('making surface %s' % var)
#			surface_map(file,surface_img_file,var,clim=clims[var])




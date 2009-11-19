from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

import utils

def plot_surface(x,y,data,filename='/Users/lederer/tmp/rompy.tmp.png'):
	print('Making plot')
	fig = Figure(facecolor='white',figsize=(12.0,12.0))
	ax = fig.add_subplot(111)
	
	ax.pcolormesh(x,y,data)
#	ax.contour(x,y,data,20)
	ax.axis('tight')
	ax.set_aspect('equal')
	ax.grid()
	FigureCanvas(fig).print_png(filename)

def plot_map(lon,lat,data,filename='/Users/lederer/tmp/rompy.map.png',resolution='i'):
	fig = Figure(facecolor='white',figsize=(24.0,24.0))
#	ax = fig.add_subplot(111)
	longest_side_size = 24.0
	ax = fig.add_axes((0.,0.,1.,1.),axisbg='grey')
	
	lllat = np.min(lat)
	urlat = np.max(lat)
	lllon = np.min(lon)
	urlon = np.max(lon)
#	print(lllat,urlat,lllon,urlon)
	m = Basemap(projection='merc',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon,resolution=resolution,ax=ax)
	x,y = m(*(lon,lat))

# Code to make the map fit snuggly with the png
#print(np.max(x), np.min(x), np.max(y),np.min(y))
#	width = np.max(x) - np.min(x)
#	height = np.max(y) - np.min(y)
	
#	if width >= height:
#		fig.set_size_inches(longest_side_size, (height/width)*longest_side_size)
#	else:
#		fig.set_size_inches((width/height)*longest_side_size, longest_side_size)
#	ax.set_position([0.,0.,1.,1.])




# 	bbox = ax.get_position()
# 	print(bbox.xmin, bbox.xmax, bbox.ymin, bbox.ymax)
# 	
# 	fig.set_size_inches((bbox.xmax - bbox.xmin)*longest_side_size, (bbox.ymax - bbox.ymin)*longest_side_size)
# 	ax.set_position([0.,0.,1.,1.])
# 	bbox = ax.get_position()
# 	print(bbox.xmin, bbox.xmax, bbox.ymin, bbox.ymax)
# 	
# 	
	m.pcolormesh(x,y,data)
	m.drawcoastlines(linewidth=0.5)

	FigureCanvas(fig).print_png(filename)

def plot_profile(data,depth,filename='/Users/lederer/tmp/rompy.profile.png'):
	fig = Figure()
	ax = fig.add_subplot(111)
	
	ax.plot(data,depth)
	
	ax.grid()
	
	FigureCanvas(fig).print_png(filename)

def plot_mickett(coords,data,varname='',region='',filename='/Users/lederer/tmp/rompy.mickett.png',n=1):
	fig = Figure(facecolor='white')
	
	ax1 = fig.add_axes([0.1, 0.5, 0.75, 0.4])
	ax2 = fig.add_axes([0.1, 0.1, 0.75, 0.4])
	cax = fig.add_axes([0.9,0.1,0.02,0.8],frameon=False)
	
	my_plot11 = ax1.contourf(np.tile(np.arange(data.shape[1]),(coords['zm'].shape[0],1)),coords['zm'],data,100)
	my_plot12 = ax1.contour(np.tile(np.arange(data.shape[1]),(coords['zm'].shape[0],1)),coords['zm'],data,100,linewidths=1,linestyle=None)
	ax1.fill_between(np.arange(data.shape[1]),coords['zm'][0,:],ax1.get_ylim()[0],color='grey')
	ax1.set_ylim((-20,0))
	ax1.set_xlim((0,data.shape[1]-1))
	
	my_plot21 = ax2.contourf(np.tile(np.arange(data.shape[1]),(coords['zm'].shape[0],1)),coords['zm'],data,100)
	my_plot22 = ax2.contour(np.tile(np.arange(data.shape[1]),(coords['zm'].shape[0],1)),coords['zm'],data,100,linewidths=1,linestyle=None)
	ax2.fill_between(np.arange(data.shape[1]),coords['zm'][0,:],ax2.get_ylim()[0],color='grey')
	ax2.set_ylim(ax2.get_ylim()[0],-20)
	ax2.set_xlim((0,data.shape[1]-1))
	
	fig.colorbar(my_plot11,cax=cax)
	ax1.set_title('%s %s from a ROMS run' % (region,varname))
	ax1.set_ylabel('depth in meters',position=(0.05,0))
	ax1.set_xticks(np.arange(data.shape[1]))
	ax1.set_xticklabels('')
	
	ax2.set_xticks(np.arange(data.shape[1]))
	
	print(ax2.get_xticks())
	if region == 'Hood Canal':
		ax2.set_xticks(n*np.arange(len(utils.hood_canal_station_list())))
		ax2.set_xticklabels(utils.hood_canal_station_list())
	elif region == 'Main Basin':
		ax2.set_xticks(n*np.arange(len(utils.main_basin_station_list())))
		ax2.set_xticklabels(utils.main_basin_station_list())
	else:
		ax2.set_xticklabels('')
	ax2.set_xlabel('station ID')
	
	FigureCanvas(fig).print_png(filename)
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

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
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.colors import Normalize, ListedColormap, LinearSegmentedColormap, hsv_to_rgb
from matplotlib.cm import ScalarMappable

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt

import utils

def red_blue_cm():
	cdict = {'red':		[(0.0, 0.0, 0.0),
						(0.5,1.0,1.0),
						(1.0, 1.0, 1.0)],
			
			'green':	[(0.0, 0.0, 0.0),
						(0.5, 1.0, 1.0),
						(1.0, 0.0, 0.0)],
			
			'blue':		[(0.0, 1.0, 1.0),
						(0.5, 1.0, 1.0),
						(1.0, 0.0, 0.0)]
			}
	return LinearSegmentedColormap('red_blue_cm',cdict,N=256)
#	return ListedColormap(['b','w','r'],name='red_blue',N=None)

def banas_cm(a,b,c,d):
	norm = Normalize(vmin=a,vmax=d,clip=False)
	cdict = {'red':[],'green':[],'blue':[]}
	
	if not a==b:
		# add dark blue
		cdict['red'].append((0., 0., 0.))
		cdict['green'].append((0., 0., 0.))
		cdict['blue'].append((0., 0., 0.25))
		
		# add blue
		cdict['red'].append((norm(b), 0., 0.))
		cdict['green'].append((norm(b), 0., 0.))
		cdict['blue'].append((norm(b), 1.0, 1.0))
	else:
		cdict['red'].append((0., 0., 0.))
		cdict['green'].append((0., 0., 0.))
		cdict['blue'].append((0., 0., 1.0))

	# add green between blue and yellow
	cdict['red'].append((norm(b + (c-b)/4.0), 0., 0.))
	cdict['green'].append((norm(b + (c-b)/4.0), 1.0, 1.0))
	cdict['blue'].append((norm(b + (c-b)/4.0), 0., 0.))

	# add yellow in the middle
	cdict['red'].append((norm((b+c)/2.0), 1.0, 1.0))
	cdict['green'].append((norm((b+c)/2.0), 1.0, 1.0))
	cdict['blue'].append((norm((b+c)/2.0), 0., 0.))

	if not c==d:
		# add red
		cdict['red'].append((norm(c), 1.0, 1.0))
		cdict['green'].append((norm(c), 0., 0.))
		cdict['blue'].append((norm(c), 0., 0.))
		
		# add dark red
		cdict['red'].append((1.0, 0.25, 0.25))
		cdict['green'].append((1.0, 0., 0.))
		cdict['blue'].append((1.0, 0., 0.))
	else:
		cdict['red'].append((1.0, 1.0, 1.))
		cdict['green'].append((1.0, 0., 0.))
		cdict['blue'].append((1.0, 0., 0.))
	
	return LinearSegmentedColormap('banas_cm',cdict,N=100)

def banas_hsv_cm(a,b,c,d,n=100):
	norm = Normalize(vmin=a,vmax=d,clip=False)
	cdict = {'red':[],'green':[],'blue':[]}
	
	aa = norm(a) # 0.0
	bb = norm(b)
	cc = norm(c)
	yy = 0.5*(bb+cc) # yellow is half way between blue and red
	dd = norm(d) # 1.0
	
	center_value = 0.87
	end_value = 0.65
	tail_end_value = 0.3
	
	blue_hue = 0.55
	yellow_hue = 1./6.
	red_hue = 0.04
	green_hue = 1./3.
	
	gg = ((green_hue - blue_hue)/(yellow_hue - blue_hue))*(yy-bb) + bb
	green_desaturation_width = 0.67
	green_desaturation_amount = 0.5
	
	ii = np.linspace(0.,1.,n)
	hue = np.zeros(ii.shape)
	sat = np.ones(ii.shape)
	val = np.zeros(ii.shape)
	hsv = np.zeros((1,n,3))
	
	for i in range(len(ii)):
		if ii[i] < bb: # if true then aa is less than bb
			#hue[i] = blue_hue
			hsv[0,i,0] = blue_hue
			#val[i] = tail_end_value*(1 - (ii[i]-aa)/(bb-aa) ) + end_value*( (ii[i]-aa)/(bb-aa) )
			hsv[0,i,2] = tail_end_value*(1 - (ii[i]-aa)/(bb-aa) ) + end_value*( (ii[i]-aa)/(bb-aa) )
		elif ii[i] <= yy:
			hsv[0,i,0] = blue_hue*(1 - (ii[i]-bb)/(yy-bb) ) + yellow_hue*( (ii[i]-bb)/(yy-bb) )
			hsv[0,i,2] = end_value*(1 - (ii[i]-bb)/(yy-bb) ) + center_value*( (ii[i]-bb)/(yy-bb) )
		elif ii[i] <= cc:
			hsv[0,i,0] = yellow_hue*(1 - (ii[i]-yy)/(cc-yy) ) + red_hue*( (ii[i]-yy)/(cc-yy) )
			hsv[0,i,2] = center_value*(1 - (ii[i]-yy)/(cc-yy) ) + end_value*( (ii[i]-yy)/(cc-yy) )
		elif ii[i] <= dd:
			hsv[0,i,0] = red_hue
			hsv[0,i,2] = end_value*(1 - (ii[i]-cc)/(dd-cc) ) + tail_end_value*( (ii[i]-cc)/(dd-cc) )
		hsv[0,i,1] = 1.0 - green_desaturation_amount * np.exp(-np.power(ii[i]/(gg*green_desaturation_width),2.0))
	
	rgb = hsv_to_rgb(hsv)
	cdict['red'].append((0.,0.,rgb[0,0,0]))
	cdict['green'].append((0.,0.,rgb[0,0,1]))
	cdict['blue'].append((0.,0.,rgb[0,0,2]))
	
	for j in range(len(ii)-2):
		i = j+1
		cdict['red'].append((ii[i],rgb[0,i,0],rgb[0,i+1,0]))
		cdict['green'].append((ii[i],rgb[0,i,1],rgb[0,i+1,1]))
		cdict['blue'].append((ii[i],rgb[0,i,2],rgb[0,i+1,2]))

	cdict['red'].append((1.0,rgb[0,-1,0],rgb[0,-1,0]))
	cdict['green'].append((1.0,rgb[0,-1,1],rgb[0,-1,1]))
	cdict['blue'].append((1.0,rgb[0,-1,2],rgb[0,-1,2]))
	
	return LinearSegmentedColormap('banas_cm',cdict,N=n)

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

def plot_mickett(coords,data,varname='',region='',filename='/Users/lederer/tmp/rompy.mickett.png',n=1,x_axis_style='kilometers',x_axis_offset=0,clim=None,cmap=None):
	fig = Figure(facecolor='white')
	fontsize = 8
	
	if cmap == 'red_blue':
		cmap = red_blue_cm()
	if cmap == 'banas_cm':
		if len(clim) == 2:
			cmap = banas_cm(clim[0],clim[0],clim[1],clim[1])
		elif len(clim) == 4:
			cmap = banas_cm(clim[0],clim[1],clim[2],clim[3])
	elif cmap == 'banas_hsv_cm':
		if len(clim) == 2:
			cmap = banas_hsv_cm(clim[0],clim[0],clim[1],clim[1])
		elif len(clim) == 4:
			cmap = banas_hsv_cm(clim[0],clim[1],clim[2],clim[3])		
	
	ax1 = fig.add_axes([0.1, 0.55, 0.75, 0.4])
	ax2 = fig.add_axes([0.1, 0.1, 0.75, 0.4])
	cax = fig.add_axes([0.9, 0.1, 0.02, 0.8],frameon=False)
	
	x_axis_as_km = utils.coords_to_km(coords)

	if not clim == None:
		norm = Normalize(vmin=clim[0],vmax=clim[-1],clip=False)
		sm = ScalarMappable(norm=norm,cmap=cmap)
		sm.set_clim(vmin=clim[0],vmax=clim[-1])
		sm.set_array(np.array([0]))
	else:
		norm = None	

	my_plot11 = ax1.contourf(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,100,norm=norm,cmap=cmap)
	my_plot12 = ax1.contour(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,100,linewidths=1,linestyle=None,norm=norm,cmap=cmap)
	
	ax1.fill_between(x_axis_as_km,coords['zm'][0,:],ax1.get_ylim()[0],color='grey')
#	ax1.set_ylim((-20,ax1.get_ylim()[1]))
	ax1.set_ylim((-20,2))
	ax1.set_xlim((0,x_axis_as_km[-1]))
	for yticklabel in ax1.get_yticklabels():
		yticklabel.set_fontsize(fontsize)
	
	my_plot21 = ax2.contourf(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,100,norm=norm,cmap=cmap)
	my_plot22 = ax2.contour(np.tile(x_axis_as_km,(coords['zm'].shape[0],1)),coords['zm'],data,100,linewidths=1,linestyle=None,norm=norm,cmap=cmap)
	ax2.fill_between(x_axis_as_km,coords['zm'][0,:],ax2.get_ylim()[0],color='grey')
#	ax2.set_ylim(ax2.get_ylim()[0],-20)
	ax2.set_xlim((0,x_axis_as_km[-1]))
	for yticklabel in ax2.get_yticklabels():
		yticklabel.set_fontsize(fontsize)
	if clim == None:
		sm = my_plot11
	
	fig.colorbar(sm,cax=cax)
	ax1.set_title('%s %s from a ROMS run' % (region,varname))
	ax1.set_ylabel('depth in meters',position=(0.05,0))
	ax1.set_xticks(10*np.arange(x_axis_as_km[-1]/10))
	ax1.set_xticklabels('')
	
	if x_axis_style == 'kilometers' or x_axis_style == 'kilometer':
			#tick_list = x_axis_as_km[::n]
			#ax2.set_xticks(tick_list)
			#ax2.set_xticklabels([int(tick) for tick in tick_list],size=fontsize)
			td = 10 #tick_distance
			
			ax2.set_xticks(td*np.arange(x_axis_as_km[-1]/td) + (x_axis_offset % td))
			ax2.set_xticklabels([int(num) for num in np.arange(-int(x_axis_offset - x_axis_offset % td),x_axis_as_km[-1],td)])
			for xticklabel in ax2.get_xticklabels():
				xticklabel.set_fontsize(fontsize)
			ax2.set_xlabel('Kilometers')

	elif x_axis_style == 'stations' or x_axis_style == 'station':
		if region == 'Hood Canal':
			tick_list = x_axis_as_km[::n]
			ax2.set_xticks(tick_list)
			ax2.set_xticklabels(utils.hood_canal_station_list(),size=fontsize)
			ax2.set_xlabel('Station ID')
			
		elif region == 'Main Basin':
			tick_list = x_axis_as_km[::n]
			ax2.set_xticks(tick_list)
			ax2.set_xticklabels(utils.main_basin_station_list(),size=fontsize)
	 		ax2.set_xlabel('Station ID') 			
 	else:
 		ax2.set_xticks(x_axis_as_km)
 		ax2.set_xticklabels('')
		ax2.set_xlabel('Kilometers')
	
	FigureCanvas(fig).print_png(filename)
from matplotlib.mlab import griddata

def interp_2d(lat,lon,data,lati,loni):
	return griddata(lat.reshape(lat.size),lon.reshape(lon.size),data.reshape(data.size),lati,loni)
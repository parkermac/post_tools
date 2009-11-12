from matplotlib.mlab import griddata

def interp_2d(lat,lon,data,lati,loni):
	return griddata(lat.reshape(lat.size),lon.reshape(lon.size),data.reshape(data.size),lati,loni)
	
def interp_2d_xy(x,y,data,xi,yi):
	return griddata(x.reshape(x.size),y.reshape(y.size),data.reshape(data.size),xi,yi)
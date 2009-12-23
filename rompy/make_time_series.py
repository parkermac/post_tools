#!/usr/bin/env python
import os
import glob
from optparse import OptionParser
import datetime as dt

from rompy import plot_utils
from rompy.extract_from_series import extract_from_series,extract_from_two_datetimes

def plot_time_series(varname='salt', x=-123.1126, y=47.4218, dir='.', seconds=3600, filelist=None, imgfile='~/tmp/rompy.time_series.png',title=None):
	
	if filelist == None or filelist == []:
		filelist = glob.glob(os.path.join(dir,'ocean_his*.nc'))
	
	# lat lon points for Hoodsport ORCA buoy
	x = [x]
	y = [y]
	
	clim_map = {
		'temp':[0,8,20,20],
		'salt':[0,21,33,33],
		'U':[-2,-2,2,2]
		}
	clabel_map = {
		'temp': u'\u00B0 C',
		'salt': 'psu',
		'U': 'm/s'
		}
	#(data, coords) = rompy.extract('ocean_his_1000.nc', varname='salt', extraction_type='profile', x=x, y=y)
#	(data, time,z) = extract_from_series(filelist,varname=varname,extraction_type='profile',x=x,y=y,freq=dt.timedelta(seconds=seconds))
	
	dt0 = dt.datetime(2006,06,01,0,0,0)
	dt1 = dt.datetime(2006,07,01,0,0,0)
	
	(data, time,z) = extract_from_two_datetimes(x=x,y=y,dt0=dt0,dt1=dt1,varname=varname,interval=seconds)
	
	
	plot_utils.plot_time_series_profile(time,z,data,filename=imgfile,clim=clim_map[varname],varname=varname,title=title, caxis_label=clabel_map[varname])
	
parser = OptionParser()

parser.add_option('-f', '--file', dest='imgfile', action='store', default='time_series.png', help='output name for image')

parser.add_option('-v', '--var', dest='varname', action='store', default='salt', help='select which variable do you want to plot as a time series, e.g. salt or temp')

parser.add_option('-y', '--lat', dest='y', action='store', default=47.4218, help='select latitude for the point')

parser.add_option('-x', '--lon', dest='x', action='store', default=-123.1126, help='select longitude for the point')

parser.add_option('-d','--dir', dest='dir', action='store', default='.', help='select the source directory of ocean_his_####.nc files. Defaults to the current directory.')

parser.add_option('-s', '--seconds', dest='seconds', action='store', default=3600, help='select frequency for data extraction')

parser.add_option('-t', '--title', dest='title', action='store', default='Time Series Plot', help='select title for image')

(options, args) = parser.parse_args()


plot_time_series(varname=options.varname, x=options.x, y=options.y, dir=options.dir, seconds=int(options.seconds), filelist=args, imgfile=options.imgfile, title=options.title)
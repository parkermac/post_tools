import os
import datetime as dt
import time

import pytz
import netCDF4 as nc
import numpy as np

import extract_from_file

def file_time(f):
	UTC = pytz.timezone('UTC')
	ncf = nc.Dataset(f,mode='r')
	ot = ncf.variables['ocean_time']
	base_time = dt.datetime.strptime(ot.units,'seconds since %Y-%m-%d %H:%M:%S').replace(tzinfo=UTC)
	offset = dt.timedelta(seconds=ot[0][0])
	ncf.close()
	return base_time + offset

def calc_num_time_slices(file1,file2,td):
	d1 = file_time(file1)
	d2 = file_time(file2)
	
	gap = d2-d1
	gap_seconds = gap.days*60*60*24 + gap.seconds
	td_seconds = td.days*60*60*24 + td.seconds
	
	t1 = time.mktime(d1.timetuple())
	t2 = time.mktime(d2.timetuple())
	tgap = t2 - t1
	
	print(gap_seconds, tgap)
	
	return int(gap_seconds/td_seconds + 1)
	

def extract_from_series(file_list,extraction_type='point',varname='zeta',freq=dt.timedelta(seconds=600),**kwargs):
	UTC = pytz.timezone('UTC')
	pacific = pytz.timezone('US/Pacific')
	
	print('Hello from extractFromSeries')
	
	if extraction_type == 'point':
		# assume var is a 2d var for now
		file_list.sort()
		ntimes = calc_num_time_slices(file_list[0],file_list[-1],freq)
		
		x = kwargs['x']
		y = kwargs['y']
		
		print(x,y,ntimes)
		data = np.zeros((len(file_list),len(x)))
		for i in range(len(x)):	
			data[0,i],junk = extract_from_file.extract_from_file(file_list[0],varname=varname,extraction_type='point',x=x[i],y=y[i])
		
		f1 = file_list[0]
		f2 = file_list[1]
		time_list = [file_time(file_list[0])]
		
		for i in range(1,len(file_list)):
			for j in range(len(x)):
				data[i,j],junk = extract_from_file.extract_from_file(file_list[i],varname=varname,extraction_type='point',x=x[j],y=y[j])
			time_list.append(file_time(file_list[i]))
		return (data,time_list)
		
#		for i in range(1,ntimes):
#			time_list.append(time_list[-1]+freq)
		
	
	if extraction_type == 'profile':
		
		file_list.sort()
		ntimes = calc_num_time_slices(file_list[0],file_list[-1],freq)
		
		ocean_time = None
		
		data = None
		z = None
		
		x = kwargs['x']
		y = kwargs['y']
		
		if len(x) > 1:
			x = np.array(x[0])
		else:
			x = np.array(x)
		if len(y) > 1:
			y = np.array(y[0])
		else:
			y = np.array(y)

		for i in range(len(file_list)):
			file = file_list[i]
			
			if not os.path.exists(file):
				raise IOError('File %s could not be located on the filesystem' %file)
			ncf = nc.Dataset(file,mode='r')
			
			if varname not in ncf.variables:
				raise IOError('File %s does not have a variable named %s' % (file, varname))
						
			#dstart = ncf.variable['dstart']
			#print(ocean_time)
			#print(dstart)
			
			# start getting data
			d,junk = extract_from_file.extract_from_file(file_list[i],varname=varname,extraction_type='profile',x=x,y=y)
			if data == None:
				data = np.zeros((d.shape[0],len(file_list)))
			data[:,i] = d.T
			if z == None:
				z = np.zeros(data.shape)
			z[:,i] = junk['zm'].T
			if ocean_time == None:
				ocean_time = []
#				ocean_time = np.zeros(data.shape)
			tmp_time = file_time(file)
			for j in range(data.shape[1]):
				ocean_time[j][i] = tmp_time
		
		return (data, ocean_time,z)
	return
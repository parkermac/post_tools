#!/usr/bin/env python
import os
import glob
import datetime as dt
import platform
from optparse import OptionParser

import threading
import Queue
from subprocess import Popen,PIPE

class Worker(threading.Thread):
	
	def __init__(self,queue):
		self.__queue = queue
		threading.Thread.__init__(self)
		
	def run(self):
		while True:
			try:
				file = self.__queue.get(True,4)
				msi_cmd = cmd + file
#				a = os.popen(msi_cmd)
#				a.close()
				p = Popen(msi_cmd, shell=True, stdout=PIPE)
				sts = os.waitpid(p.pid, 0)[1]
				self.__queue.task_done()
				print ("msi_cmd: %s has completed" % msi_cmd)
			except Queue.Empty:
				break

parser = OptionParser()
parser.add_option('-C', '--crude',
					dest='crude_coast',
					action='store_true',
					default=False,
					help='this option will set the coastline resolution to crude')

parser.add_option('-L', '--low',
					dest='low_coast',
					action='store_true',
					default=False,
					help='this option will set the coastline resolution to low')

parser.add_option('-I', '--intermediate',
					dest='int_coast',
					action='store_true',
					default=False,
					help='this option will set the coastline resolution to intermediate')

parser.add_option('-H', '--high',
					dest='high_coast',
					action='store_true',
					default=False,
					help='this option will set the coastline resolution to high')

parser.add_option('-F', '--full',
					dest='full_coast',
					action='store_true',
					default=False,
					help='this option will set the coastline resolution to full')

WORKERS = os.sysconf('SC_NPROCESSORS_ONLN')
(options, args) = parser.parse_args()

cmd_base = './make_standard_images.py '

if options.crude_coast:
	cmd = cmd_base + '-C '
elif options.low_coast:
	cmd = cmd_base + '-L '
elif options.int_coast:
	cmd = cmd_base + '-I '
elif options.high_coast:
	cmd = cmd_base + '-H '
elif options.full_coast:
	cmd = cmd_base + '-F '
else:
	cmd = cmd_base + '-F '

today   = dt.datetime.today()

queue = Queue.Queue(WORKERS + 1)

for i in range(WORKERS):
	Worker(queue).start()

if args == []:
	file_list = glob.glob('ocean_his*.nc')
else:
	file_list = args
	
try:
	for file in file_list:
		print ("adding %s to the queue" % file)
		queue.put(file)
		
except (KeyboardInterrupt, SystemExit):
	sys.exit()
	
else:
	queue.join()
	
	if platform.system() == 'Darwin':
		print('Making the movies and putting them on metoc1')
		p = Popen('./make_movie.py /Volumes/lederer/Sites/rompy/movies image_sequence 1 6',shell=True,stdout=PIPE)
		sts = os.waitpid(p.pid, 0)[1]
	
		print('making Hoodsport time series for salinity')
		p = Popen('./make_time_series.py -s 10800 -t "Salinity at Hoodsport ORCA Buoy" -v salt -f /Volumes/lederer/Sites/rompy/movies/hoodsport_salt.png',shell=True,stdout=PIPE)
		sts = os.waitpid(p.pid, 0)[1]
		
		print('making Hoodsport time series for temperature')
		p = Popen('./make_time_series.py -s 10800 -t "Temperature at Hoodsport ORCA Buoy" -v temp -f /Volumes/lederer/Sites/rompy/movies/hoodsport_temp.png',shell=True,stdout=PIPE)
		sts = os.waitpid(p.pid, 0)[1]
		
	elif platform.system() == 'Linux':
		print('making Hoodsport time series for salinity')
		p = Popen('./make_time_series.py -s 10800 -t "Salinity at Hoodsport ORCA Buoy" -v salt -f image_sequence/hoodsport_salt_200606.png -a 200606010000 -b 200607010000',shell=True,stdout=PIPE)
		sts = os.waitpid(p.pid, 0)[1]
		
		print('making Hoodsport time series for temperature')
		p = Popen('./make_time_series.py -s 10800 -t "Temperature at Hoodsport ORCA Buoy" -v temp -f image_sequence/hoodsport_temp_200606.png -a 200606010000 -b 200607010000',shell=True,stdout=PIPE)
		sts = os.waitpid(p.pid, 0)[1]
		
		print('making Hoodsport time series for salinity')
		p = Popen('./make_time_series.py -s 10800 -t "Salinity at Hoodsport ORCA Buoy" -v salt -f image_sequence/hoodsport_salt_200607.png -a 200607010000 -b 200608010000',shell=True,stdout=PIPE)
		sts = os.waitpid(p.pid, 0)[1]
		
		print('making Hoodsport time series for temperature')
		p = Popen('./make_time_series.py -s 10800 -t "Temperature at Hoodsport ORCA Buoy" -v temp -f image_sequence/hoodsport_temp_200607.png -a 200607010000 -b 200608010000',shell=True,stdout=PIPE)
		sts = os.waitpid(p.pid, 0)[1]
		
		print('making Hoodsport time series for salinity')
		p = Popen('./make_time_series.py -s 10800 -t "Salinity at Hoodsport ORCA Buoy" -v salt -f image_sequence/hoodsport_salt_200608.png -a 200608010000 -b 200609010000',shell=True,stdout=PIPE)
		sts = os.waitpid(p.pid, 0)[1]
		
		print('making Hoodsport time series for temperature')
		p = Popen('./make_time_series.py -s 10800 -t "Temperature at Hoodsport ORCA Buoy" -v temp -f image_sequence/hoodsport_temp_200608.png -a 200608010000 -b 200609010000',shell=True,stdout=PIPE)
		sts = os.waitpid(p.pid, 0)[1]
	
	print ('threaded_make_standard task completed. Total time elapsed: %d seconds' %(dt.datetime.today() - today).seconds)

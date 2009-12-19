#!/usr/bin/env python
import os
import glob
import datetime as dt
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
WORKERS = os.sysconf('SC_NPROCESSORS_ONLN')
(options, args) = parser.parse_args()

cmd = './make_standard_images.py '

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
	
	print('Making the movies and putting them on metoc1')
	p = Popen('./make_movie.py /Volumes/lederer/Sites/rompy/movies image_sequence 1 6',shell=True,stdout=PIPE)
	sts = os.waitpid(p.pid, 0)[1]
	
	print ('threaded_make_standard task completed. Total time elapsed: %d seconds' %(dt.datetime.today() - today).seconds)



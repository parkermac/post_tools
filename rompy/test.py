#!/usr/bin/env python
from rompy import rompy, plot_surface
import numpy as np

(data, coords) = rompy.extract('ocean_his_1000.nc')
plot_surface.plot_surface(coords['xm'],coords['ym'],data)
plot_surface.plot_map(coords['xm'],coords['ym'],data,filename='/Users/lederer/tmp/rompy.map.png')

del(data)
del(coords)
# full domain
#x = np.linspace(-127.,-122.,100)
#y = np.linspace(44.,50.,100)

#puget sound area
x = np.linspace(-123.,-122.,500)
y = np.linspace(47.,48.,500)

(data, coords) = rompy.extract('ocean_his_1000.nc', x=x, y=y)
plot_surface.plot_map(coords['xm'],coords['ym'],data,filename='/Users/lederer/tmp/rompy.map2.png',resolution='h')
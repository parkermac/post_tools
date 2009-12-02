#!/usr/bin/env python
import glob
import datetime as dt

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

from rompy.extract_from_series import extract_from_series
from rompy import utils

filelist = glob.glob('ocean_his*.nc')
x,y = utils.station_to_lat_lon([24])
#(data, coords) = rompy.extract('ocean_his_1000.nc', varname='salt', extraction_type='profile', x=x, y=y)
(data, time,z) = extract_from_series(filelist,varname='salt',extraction_type='profile',x=x,y=y,freq=dt.timedelta(seconds=1800))

print(time)
print(data)
print(z)

plt.contour(data)
plt.show()
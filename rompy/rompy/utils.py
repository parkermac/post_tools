import numpy as np
from matplotlib.mlab import griddata

def interp_2d(lat,lon,data,lati,loni):
	return griddata(lat.reshape(lat.size),lon.reshape(lon.size),data.reshape(data.size),lati,loni)
	
def interp_2d_xy(x,y,data,xi,yi):
	return griddata(x.reshape(x.size),y.reshape(y.size),data.reshape(data.size),xi,yi)

def interp_3d_point(x,y,z,d,xi,yi,zi):
	if not x.ndim == 1 or not y.ndim == 1 or not z.ndim == 1:
		raise(TypeError,'interp_3d_from_point needs the x, y, and z to be vectors')
	if not xi.size == 1 or not yi.size == 1 or not zi.size ==1:
		print(xi,yi,zi)
		raise(TypeError, 'interp_3d_from_point needs xi, yi, and zi to be a single value')
	xl = np.nonzero(x < xi)[0][-1]
	xh = np.nonzero(xi <= x)[0][0]
	yl = np.nonzero(y < yi)[0][-1]
	yh = np.nonzero(yi <= y)[0][0]
	zl = np.nonzero(z < zi)[0][-1]
	zh = np.nonzero(zi <= z)[0][0]
#	print((xl,xi, xh),(yl, yi, yh),( zl, zi, zh))
	
	xd = (xi-x[xh])/(x[xl]-x[xh])
	yd = (yi-y[yh])/(y[yl]-y[yh])
	zd = (zi-z[zh])/(z[zl]-z[zh])
	
	i0 = d[zl,yl,xl]*(1-zd) + d[zh,yl,xl]*zd
	i1 = d[zl,yh,xl]*(1-zd) + d[zh,yh,xl]*zd
	
	j0 = d[zl,yl,xh]*(1-zd) + d[zh,yl,xh]*zd
	j1 = d[zl,yh,xh]*(1-zd) + d[zh,yh,xh]*zd
	 
	w0 = i0*(1-yd) + i1*yd
	w1 = j0*(1-yd) + j1*yd
	
	return w0*(1-xd) + w1*xd

def interp_3d_from_list_of_points(x,y,z,d,p_list):
	di = np.zeros(len(p_list),1)
	for i in range(len(p_list)):
		di[i] = interp_3d_point(x,y,z,d,p_list[i][0],p_list[i][1],p_list[i][2])
	return di

	
def interp_3d(x,y,z,data,xi,yi,zi):
	# we make a lot of assumptions about the incoming data. this is not a universal interpn
	if x.shape == y.shape and x.shape == z.shape and x.shape == data.shape:
		print('Do this the hard way')
		# assume x, y, and z are the same everywhere in their respective dimension
		if x.ndim == 3:
			x_vec = x[0,0,:]
		else:
			x_vec = x
		if y.ndim == 3:
			y_vec = y[0,:,0]
		else:
			y_vec = y
		if z.ndim == 3:
			z_vec = z[:,0,0]
		else:
			z_vec = z
		
		# assume xi, yi, and zi are vectors
		if xi.ndim == 1 and yi.ndim == 1 and zi.ndim == 1:
			di = np.zeros((len(zi), len(yi),len(xi)))
			for i in range(len(xi)):
				for j in range(len(yi)):
					for k in range(len(zi)):
						di[k,j,i] = interp_3d_point(x_vec,y_vec,z_vec,data,xi[i],yi[j],zi[k])
		elif xi.shape == yi.shape and xi.shape == zi.shape:
			di = np.zeros(xi.shape)
			for i in range(xi.size):
				di.flat[i] = interp_3d_point(x_vec,y_vec,z_vec,data,xi.flat[i],yi.flat[i],zi.flat[i])
		return di
		
	elif (len(x),len(y),len(z)) == data.shape:
		print('Do this the other way')


def meshgrid(*xi,**kwargs):
    """
    Return coordinate matrices from one or more coordinate vectors.

    Make N-D coordinate arrays for vectorized evaluations of
    N-D scalar/vector fields over N-D grids, given
    one-dimensional coordinate arrays x1, x2,..., xn.

    Parameters
    ----------
    x1, x2,..., xn : array_like
        1-D arrays representing the coordinates of a grid.
    indexing : 'xy' or 'ij' (optional)
        cartesian ('xy', default) or matrix ('ij') indexing of output
    sparse : True or False (default) (optional)
         If True a sparse grid is returned in order to conserve memory.
    copy : True (default) or False (optional)
        If False a view into the original arrays are returned in order to
        conserve memory

    Returns
    -------
    X1, X2,..., XN : ndarray
        For vectors `x1`, `x2`,..., 'xn' with lengths ``Ni=len(xi)`` ,
        return ``(N1, N2, N3,...Nn)`` shaped arrays if indexing='ij'
        or ``(N2, N1, N3,...Nn)`` shaped arrays if indexing='xy'
        with the elements of `xi` repeated to fill the matrix along
        the first dimension for `x1`, the second for `x2` and so on.

    See Also
    --------
    index_tricks.mgrid : Construct a multi-dimensional "meshgrid"
                     using indexing notation.
    index_tricks.ogrid : Construct an open multi-dimensional "meshgrid"
                     using indexing notation.

    Examples
    --------
    >>> x = np.linspace(0,1,3)   # coordinates along x axis
    >>> y = np.linspace(0,1,2)   # coordinates along y axis
    >>> xv, yv = meshgrid(x,y)   # extend x and y for a 2D xy grid
    >>> xv
    array([[ 0. ,  0.5,  1. ],
           [ 0. ,  0.5,  1. ]])
    >>> yv
    array([[ 0.,  0.,  0.],
           [ 1.,  1.,  1.]])
    >>> xv, yv = meshgrid(x,y, sparse=True)  # make sparse output arrays
    >>> xv
    array([[ 0. ,  0.5,  1. ]])
    >>> yv
    array([[ 0.],
           [ 1.]])

    >>> meshgrid(x,y,sparse=True,indexing='ij')  # change to matrix indexing
    [array([[ 0. ],
           [ 0.5],
           [ 1. ]]), array([[ 0.,  1.]])]
    >>> meshgrid(x,y,indexing='ij')
    [array([[ 0. ,  0. ],
           [ 0.5,  0.5],
           [ 1. ,  1. ]]),
     array([[ 0.,  1.],
           [ 0.,  1.],
           [ 0.,  1.]])]

    >>> meshgrid(0,1,5)  # just a 3D point
    [array([[[0]]]), array([[[1]]]), array([[[5]]])]
    >>> map(np.squeeze,meshgrid(0,1,5))  # just a 3D point
    [array(0), array(1), array(5)]
    >>> meshgrid(3)
    array([3])
    >>> meshgrid(y)      # 1D grid; y is just returned
    array([ 0.,  1.])

    `meshgrid` is very useful to evaluate functions on a grid.

    >>> x = np.arange(-5, 5, 0.1)
    >>> y = np.arange(-5, 5, 0.1)
    >>> xx, yy = meshgrid(x, y, sparse=True)
    >>> z = np.sin(xx**2+yy**2)/(xx**2+yy**2)
    """
    copy = kwargs.get('copy',True)
    args = np.atleast_1d(*xi)
    if not isinstance(args, list):
        if args.size>0:
            return args.copy() if copy else args
        else:
            raise TypeError('meshgrid() take 1 or more arguments (0 given)')

    sparse = kwargs.get('sparse',False)
    indexing = kwargs.get('indexing','xy') # 'ij'


    ndim = len(args)
    s0 = (1,)*ndim
    output = [x.reshape(s0[:i]+(-1,)+s0[i+1::]) for i, x in enumerate(args)]

    shape = [x.size for x in output]

    if indexing == 'xy':
        # switch first and second axis
        output[0].shape = (1,-1) + (1,)*(ndim-2)
        output[1].shape = (-1, 1) + (1,)*(ndim-2)
        shape[0],shape[1] = shape[1],shape[0]

    if sparse:
        if copy:
            return [x.copy() for x in output]
        else:
            return output
    else:
        # Return the full N-D matrix (not only the 1-D vector)
        if copy:
            mult_fact = np.ones(shape,dtype=int)
            return [x*mult_fact for x in output]
        else:
            return np.broadcast_arrays(*output)


def ndgrid(*args,**kwargs):
    """
    Same as calling meshgrid with indexing='ij' (see meshgrid for
    documentation).
    """
    kwargs['indexing'] = 'ij'
    return meshgrid(*args,**kwargs)

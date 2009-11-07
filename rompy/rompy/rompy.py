import extract_from_file
import extract_from_series

def extract(files,**kwargs):
	if files.__class__ == str:
		if kwargs.has_key('x') and kwargs.has_key('y'):
			(data, coords) = extract_from_file.extract_from_file(file=files, varname='zeta', extraction_type='points', **kwargs)
		else:
			(data, coords) = extract_from_file.extract_from_file(files)
	elif files.__class__ == list:
		(data, coords) = extract_from_series.extract_from_series(files)
	else:
		raise TypeError('wrong type sent to the extractor')
		return False
	return (data,coords)

def test():
	print('rompy.test() works')
	extract('a.nc')
	#extract(['a.nc','b.nc'])
	try:
		extract({'first':'a.nc'})
	except TypeError, e:
		print(e)
	return

if __name__ == '__main__':
	test()
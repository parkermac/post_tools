import extract_from_file
import extract_from_series

def extract(files,**kwargs):
	if files.__class__ == str:
		extract_from_file.extract_from_file(files)
	elif files.__class__ == list:
		extract_from_series.extract_from_series(files)
	else:
		raise TypeError('wrong type sent to the extractor')
		return False
	return True

def test():
	print('rompy.test() works')
	extract('a.nc')
	extract(['a.nc','b.nc'])
	try:
		extract({'first':'a.nc'})
	except TypeError, e:
		print(e)
	return

if __name__ == '__main__':
	test()
function [data,varargout] = roms_extractFromSeries(series, varname, t, theType, varargin);

% [data,zm,ym,xm] = ...
%     roms_extractFromSeries(series, varname, t, 'full');
%                                           ..., 'surface');
%                                           ..., 'zslice', z);
%                                           ..., 'profile', y, x);
%                                           ..., 'point', z, y, x);
% [data,tm,zm,ym,xm] = ...
%     roms_extractFromSeries(series, varname, tvector, ...
%
% general routine for extracting data from a ROMS netcdf file series
% (e.g., ocean_his_*.nc) in data units.
%
% _series_ is a seriesDef: see roms_createSeriesDef.m.
% see roms_extract.m for more info.
%
% neil banas feb 2009


% note: for code simplicity, in this version all coordinate outputs are
% calculated regardless of nargout.

if ~isstruct(series) | ~isfield(series,'ncn') | ~isfield(series,'nctime')
	error('need to supply a seriesDef with ncn and nctime defined.');
elseif isempty(series.ncn) | isempty(series.nctime)
	error('file series timebase is empty.');
end

if length(t) > 1 % if t is a vector, recurse over each time in the vector

	for n = 1:length(t)
		% extract one time slice
		[data1,zm1,ym1,xm1] = roms_extract(series, varname, t(n), theType, varargin);
		if n==1 % intialize output variables
			data = repmat(nan,[length(t) size(data1)]);
			tm = data;
			zm = data;
			ym = data;
			xm = data;
		end
		% place the time slice in the output variables
		data(n,:) = data1;
		tm(n,:) = repmat(t(n),[1 size(data1)]);
		zm(n,:) = zm1;
		ym(n,:) = ym1;
		xm(n,:) = xm1;
	end
	
	switch nargout
		case 2, varargout = {tm};
		case 3, varargout = {tm,zm};
		case 4, varargout = {tm,zm,ym};
		case 5, varargout = {tm,zm,ym,xm};
	end
		
else % t is scalar

	% file numbers bracketing t
	n = interp1(series.nctime, series.ncn, t);
	n0 = floor(n);
	file0 = roms_filename([series.dirname series.basename],n0);
	n1 = ceil(n);
	file1 = roms_filename([series.dirname series.basename],n1);
	
	% make sure the variable exists
	info = nc_info(file0);
	if ~strmatch(varname, {info.Dataset.Name})
		error(['variable ''' varname ''' not found in ' file0 '.']);
	end
	
	% do the extraction
	[data,zm,ym,xm] = roms_extractFromFile(file0, varname, varargin);
	if n0 ~= n1 % interpolate between frames if necessary
		[data1,zm1,ym1,xm1] = roms_extractFromFile(file1, varname, varargin);
		fr = (n-n0)/(n1-n0);
		data = data + (data1-data).*fr;
		zm = zm + (zm1-zm).*fr;
	end
	
	switch nargout
		case 2, varargout = {zm};
		case 3, varargout = {zm,ym};
		case 4, varargout = {zm,ym,xm};
	end
		
end


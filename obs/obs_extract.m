function data = obs_extract(files, vars, timeRange, varargin)
%------------------------------------------------------------
% data = obs_extract(filename, vars, timeRange);
%                                         ..., LatLonVector);
%                                         ..., z, LatLonVector);
% data = obs_extract(directory, ...
% data = obs_extract(cell array of files and directories, ...
%
% looks inside a netcdf file for one or more observational variables and
% returns them in a structure _data_, along with coordinate variables.
%
% if a directory is given instead of a filename, recurses over 
% files inside it and concatenates the output.
%
% vars can be a string containing one variable name, or a cell array of strings,
% or 'all'.
%
% timeRange (if given) is a 2 element vector of matlab datenums or 'all'.
%
% LatLonVector can be:
%    1) geographic range specified either by a 2-column lat-long polygon or
%    2) a lat/lon pair and a radius (in km) as [lat lon radius]
%
% if z is given, the data is returned just at that depth
% 
% in the case of 1-D variables, data contains each variable along with
% t (time), y (lat), x (lon), depth, and if appropriate, cast.
% in the case of 2-D (depth-time) variables, data contains each variable along
% with plaid matrices t, depth and also lon, lat matching the data variables
% in size. (If a directory of files is requested, t and depth may be non-monotonic,
% but still plaid.)
%
% written by D. Sutherland and N. Banas, UW, Spring 2009
%
%
% rewritten by NSB, jun 09, to accept a cell array of names as input; this
% this makes obs_extractDir.m unnecessary. Now recurses over subdirectories too.
%------------------------------------------------------------

verbose = 1;
data = [];

% input = cell array of files and directories
if iscell(files)
	for i=1:length(files(:))
		data2 = obs_extract(files{i}, vars, timeRange, varargin{:});
		data = catstruct(data,data2);
	end
else
	type = exist(files, 'file');	
	if type == 7 % input = directory
		if files(end) ~= '/', files = [files '/']; end %make sure there's a backslash at end
		if verbose, disp(['obs_extract: looking inside ' files]); end
		filelist = dir(files);
		for i=1:length(filelist)
			rec = filelist(i);
			if rec.isdir
				if rec.name(1) ~= '.'
					% directory within the directory
					data2 = obs_extract(rec.name, vars, timeRange, varargin{:});
					data = catstruct(data,data2);
				end
			elseif length(rec.name)>3 & strcmp(rec.name(end-2:end),'.nc')
				% .nc file within the directory
				fullname = [files rec.name];
				if verbose, disp(['obs_extract: reading from ' fullname]); end
				data2 = obs_extractFile(fullname, vars, timeRange, varargin{:});
				data = catstruct(data,data2);
			end
		end
	elseif type == 2 % input = file
		if length(files)>3 & strcmp(files(end-2:end),'.nc')
			if verbose, disp(['obs_extract: reading from ' files]); end
			data = obs_extractFile(files, vars, timeRange, varargin{:});
		end
	else
		error(['Cannot find the file(s): ' files ' as requested']);
	end
end


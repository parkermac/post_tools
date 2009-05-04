function data = obs_extract(files, vars, timeRange, varargin)
%------------------------------------------------------------
% data = obs_extract(filename, vars);
%                              ..., timeRange);
%                                         ..., LatLonVector);
%                                         ..., z, LatLonVector);
% data = obs_extract(directory, ...
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
%------------------------------------------------------------

%check to see if file or directory of files exist
type = exist(files, 'file');
if type == 7 %then a directory exists
    if files(end) ~= '/', files = [files '/']; end %make sure there's a backslash at end
        filelist = dir([files, '*.nc']); %pick out netcdf files in that directory
elseif type == 2 %then files is just one file
        filelist = files;
else
    error(['Cannot find the file(s): ' files ' as requested']);
end

%do data extraction from other mfiles...
switch type 
    case 7
        data = obs_extractDir(filelist, vars, timeRange, varargin{:});
    case 2
        data = obs_extractFile(filelist, vars, timeRange, varargin{:});
end

function [data,varargout] = roms_extract(varargin);

% [data,zm,ym,xm] = roms_extract(series, varname, t, 'full');
%                                               ..., 'surface');
%                                               ..., 'zslice', z);
%                                               ..., 'profile', y, x);
%                                               ..., 'point', z, y, x);
%
% [data,tm,zm,ym,xm] = roms_extract(series, varname, tvector, ...
%
% [data,zm,ym,xm] = roms_extract(filename, varname, 'full');
%                                              ..., 'surface');
%                                              ..., 'zslice', z);
%                                              ..., 'profile', y, x);
%                                              ..., 'point', z, y, x);
%
% general routine for extracting data from a ROMS netcdf file series or a
% single file (e.g., ocean_his_*.nc) in data units.
%
% If a seriesDef _series_ is given (see roms_createSeriesDef.m), t can be
% either a scalar or a vector. If a filename is given, t isn't given at all.
% if the extraction type is 'full' and the coordinate variables zm,ym,xm aren't
% requested, then the variable _varname_ can be anything. Otherwise, it's assumed
% to be a 4D field of size [1 K J I], where K,J,I match the rho, u, v, or w grid.
%
% z, y, x can be scalars or vectors.
% zm,ym,xm are plaid matrices the same size as _data_.

%
% neil banas feb 2009

if ischar(varargin{1})
	[data,varargout] = roms_extractFromFilename(varargin{:});
elseif isstruct(varargin{1})
	[data,varargout] = roms_extractFromSeries(varargin{:});
end

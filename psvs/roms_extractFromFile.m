function [data,varargout] = roms_extractFromFile(filename, varname, theType, varargin);

% [data,zm,ym,xm] = ...
%     roms_extractFromFile(filename, varname, 'full');
%                                        ..., 'surface');
%                                        ..., 'zslice', z);
%                                        ..., 'profile', y, x);
%                                        ..., 'point', z, y, x);
% [data,tm,zm,ym,xm] = ...
%     roms_extractFromFile(series, varname, tvector, ...
%
% general routine for extracting data from a single ROMS netcdf file
% (e.g., ocean_his_*.nc) in data units. See roms_extract.m for more info.
%
% neil banas feb 2009


% make sure the file and the variable both exist
if ~exist(filename)
	error([filename ' not found.']);
end
info = nc_info(filename);
if ~strmatch(varname, {info.Dataset.Name})
	error(['variable ''' varname ''' not found in ' filename '.']);
end

% if the extraction type is 'full' and no coordinate variables are requested,
% then just extract it without checking dimensions or any such
if strcmp(theType,'full') & nargout==1
	data = nc_varget(filename, varname);
	return;
end

% from here on, we assume a 4D variable of size [1 K J I], where [K J I]
% match the rho, u, v, or w grids

% check variable size
sz = nc_varsize(filename, varname);
ndims = length(sz);
if ndims ~= 4
	error([varname ' is size [' num2str(sz) '] and this the wrong number of dimensions.']);
elseif sz(1) ~= 1
	error([varname ' is size [' num2str(sz) '] and the first dimension should be length 1.']);
end
K = sz(2);
J = sz(3);
I = sz(4);

% load grid and make full coordinate variables
% Some of this may not need to be done depending on which type of extraction is
% specified and whether coordinate variables were requested as outputs. But for
% code simplicity we just do the whole thing every time.
G = roms_loadGrid(filename);
% now K,J,I are the dimensions of _varname_, while G.K, G.J, G.I are the
% dimensions of the rho grid
% pick the right x,y,cs axes from the model grid
if I==G.I & J==G.J
	x2 = G.lon;
	y2 = G.lat;
	mask2 = G.mask;
elseif I==G.I-1 & J==G.J
	x2 = G.lonu;
	y2 = G.latu;
	mask2 = G.masku;
elseif I==G.I & J==G.J-1
	x2 = G.lonv;
	y2 = G.latv;
	mask2 = G.maskv;
else
	error([varname ' is size [' num2str(sz) '] and J,I are a bad size.']);
end
zeta2 = interp2(G.lon,G.lat,squeeze(nc_varget(filename,'zeta')),x2,y2);
H2 = interp2(G.lon,G.lat,G.H,x2,y2);
if K==G.K
	cs = G.cs;
elseif K==G.K+1
	cs = G.csw;
elseif K==1
	cs = 0;
else
	error([varname ' is size [' num2str(sz) '] and K is a bad length.']);
end
% make full 3d coordinate variables
x3 = repmat(reshape(x2, [1 J I]), [K 1 1]);
y3 = repmat(reshape(y2, [1 J I]), [K 1 1]);
mask3 = repmat(reshape(mask2, [1 J I]), [K 1 1]);
zeta3 = repmat(reshape(zeta2,[1 J I]),[K 1 1]);
H3 = repmat(reshape(H2,[1 J I]),[K 1 1]);
cs3 = repmat(cs(:),[1 J I]);
z3 = zeta3 + cs3.*(zeta3+H3);
% z3,y3,x3,etc. match the full variable in size
% zm,ym,xm,etc. (not yet created) match the extraction in size

% do the extraction. In this version, the extractions that require interpolation
% ('zslice', 'profile', and 'point') read the full variable in and do the
% interpolation in memory. If this is too slow, they could be optimized to
% read only the bits of the variable from the file
switch theType
	case 'full' % extract the whole variable
		data = nc_varget(filename,varname);
		zm = z3;
		ym = y3;
		xm = x3;
		data(mask3==0) = nan;
		
	case 'surface' % extract the surface field
		data = nc_varget(filename, varname, [0 K-1 0 0], [1 1 -1 -1]);
		zm = zeta2;
		ym = y2;
		xm = x2;
		data(mask2==0) = nan;
		
	case 'zslice' % extract a slice at one or more z values
		z = varargin{1};
		zm = repmat(z(:),[1 J I]);
		ym = repmat(reshape(y2,[1 J I]),[length(z(:)) 1 1]);
		xm = repmat(reshape(x2,[1 J I]),[length(z(:)) 1 1]);
		zetam = repmat(reshape(zeta2,[1 J I]),[length(z(:)) 1 1]);
		Hm = repmat(reshape(H2,[1 J I]),[length(z(:)) 1 1]);
		csm = (zm - zetam) ./ (zetam + Hm);
		bad = isnan(csm) | csm > 0 | csm < -1;
		csm(isnan(csm)) = 0;
		data3 = squeeze(nc_varget(filename,varname));
		data = interpn(cs3,y3,x3,data3,csm,ym,xm);
		data(bad) = nan;
		
	case 'profile' % extract vertical profiles at one or more (y,x) pairs
		y = varargin{1}(:);
		x = varargin{2}(:);
		if length(x)==1 & length(y)>1
			x = repmat(x,size(y));
		elseif length(x)>1 & length(y)==1
			y = repmat(y,size(x));
		end
		L = length(x);
		cs = unique([-1; sort(cs); 0]); % extend to top & bottom of water column
		K = length(cs);
		ym = repmat(reshape(y,[1 L]),[K 1]);
		xm = repmat(reshape(x,[1 L]),[K 1]);
		zetam = interp2(x2,y2,zeta2,xm,ym);
		Hm = interp2(x2,y2,H2,xm,ym);
		csm = repmat(cs,[1 L]);
		zm = zetam + csm .* (zetam + Hm);
		data3 = squeeze(nc_varget(filename,varname));
		data = interpn(cs3,y3,x3,data3,csm,ym,xm);
	
	case 'point' % extract data at one or more (z,y,x) triplets
		z = varargin{1}(:);
		y = varargin{2}(:);
		x = varargin{3}(:);
		L = max([length(x) length(y) length(z)]);
		if L > 1
			if length(x)==1, x = repmat(x,[L 1]); end
			if length(y)==1, y = repmat(y,[L 1]); end
			if length(z)==1, z = repmat(z,[L 1]); end
		end
		zm = z;
		ym = y;
		xm = x;
		zetam = interp2(x2,y2,zeta2,xm,ym);
		Hm = interp2(x2,y2,H2,xm,ym);
		csm = (zm - zetam) ./ (zetam + Hm);
		bad = isnan(csm) | csm > 0 | csm < -1;
		csm(isnan(csm)) = 0;
		data3 = squeeze(nc_varget(filename,varname));
		data = interpn(cs3,y3,x3,data3,csm,ym,xm);
		data(bad) = nan;
		
	otherwise
		error(['''' theType ''' is not a valid extraction type.']);			
end


% clean up outputs
data = squeeze(data);
switch nargout
	case 2, varargout = {squeeze(zm)};
	case 3, varargout = {squeeze(zm), squeeze(ym)};
	case 4, varargout = {squeeze(zm), squeeze(ym), squeeze(xm)};
end
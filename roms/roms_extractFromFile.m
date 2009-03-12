function [data,coords] = roms_extractFromFile(filename, varname, theType, varargin);

% [data,coords] = ...
%     roms_extractFromFile(filename, 3Dvarname, 'full');
%                                          ..., 'point', y, x);
%     roms_extractFromFile(filename, 4Dvarname, 'full');
%                                          ..., 'surface');
%                                          ..., 'zslice', z);
%                                          ..., 'profile', y, x);
%                                          ..., 'point', z, y, x);
%
% general routine for extracting data from a single ROMS netcdf file
% (e.g., ocean_his_*.nc) in data units. There should be no reason for the
% user to call this directly: use roms_extract.m.
%
% neil banas feb 2009


% make sure the file and the variable both exist
if ~exist(filename)
	error([filename ' not found.']);
end
if ~nc_isvar(filename,varname)
	error(['variable ''' varname ''' not found in ' filename '.']);
end

% if theType is 'full' and no coordinate variables were requested, just
% return the whole variable without any checking of dimensions
if strcmp(theType,'full') & nargout==1
	data = nc_varget(filename,varname);
	return;
end

% get variable size & check dimensions
info = nc_getvarinfo(filename, varname);
dims = info.Dimension;
sz = info.Size;
ndims = length(sz);
if ~(ndims==3 | ndims==4), error([varname ' has the wrong number of dimensions.']); end
if ~strcmp(info.Dimension{1},'ocean_time'), error(['the first dimension of ' varname ' isn''t ocean_time.']); end
if sz(1) ~= 1, error(['the first dimension of ' varname ' isn''t length 1.']); end


if ndims==3
	% 3D variables ---------------------------------------------------------------------------------
	% grid and coordinates -------------------------------
	G = roms_loadGrid(filename);
	if strcmp(dims{2},'eta_rho') & strcmp(dims{3},'xi_rho')
		y2 = G.lat;
		x2 = G.lon;
		mask2 = G.mask;
	elseif strcmp(dims{2},'eta_v') & strcmp(dims{3},'xi_v')
		y2 = G.latv;
		x2 = G.lonv;
		mask2 = G.maskv;
	elseif strcmp(dims{2},'eta_u') & strcmp(dims{3},'xi_u')
		y2 = G.latu;
		x2 = G.lonu;
		mask2 = G.masku;
	else
		error(['the 2nd & 3rd dimensions of ' varname ' should be (eta_something, xi_something).']);
	end
	
	% do the extraction ----------------------------------
	switch theType
		case {'full','surface'} % ----- 3D, full
			data = nc_varget(filename,varname);
			data(mask2==0) = nan;
			ym = y2;
			xm = x2;
			
		case {'profile','profiles','point','points'} % ----- 3D, points
			ym = varargin{1}(:);
			xm = varargin{2}(:);
			L = max([length(ym) length(xm)]);
			if length(xm)==1, xm = repmat(xm,[L 1]); end
			if length(ym)==1, ym = repmat(ym,[L 1]); end
			data2 = nc_varget(filename,varname);
			data2(mask2==0) = nan;
			data = interp2(x2,y2,data2,xm,ym);
			
		otherwise
			error(['don''t know how to extract ''' theType ''' for 3D variables.']);
	end
	
	% clean up outputs -----------------------------------
	data = squeeze(data);
	coords.ym = squeeze(ym);
	coords.xm = squeeze(xm);

else
	% 4D variables --------------------------------------------------------------------------------
	% grid and coordinates -------------------------------
	K = sz(2);
	J = sz(3);
	I = sz(4);
	G = roms_loadGrid(filename);
	% make (x2,y2,mask2,zeta2,H2) matching the variable in size in the y and x directions
	if strcmp(dims{3},'eta_rho') & strcmp(dims{4},'xi_rho')
		y2 = G.lat;
		x2 = G.lon;
		mask2 = G.mask;
	elseif strcmp(dims{3},'eta_v') & strcmp(dims{4},'xi_v')
		y2 = G.latv;
		x2 = G.lonv;
		mask2 = G.maskv;
	elseif strcmp(dims{3},'eta_u') & strcmp(dims{4},'xi_u')
		y2 = G.latu;
		x2 = G.lonu;
		mask2 = G.masku;
	else
		error(['the 3rd & 4th dimensions of ' varname ' should be (eta_something, xi_something).']);
	end
	zeta2 = interp2(G.lon,G.lat,squeeze(nc_varget(filename,'zeta')),x2,y2);
	H2 = interp2(G.lon,G.lat,G.H,x2,y2);
	% make (cs) matching the variable in size in the z direction
	if strcmp(dims{2},'s_rho')
		cs = G.cs;
	elseif strcmp(dims{2},'s_w')
		cs = G.csw;
	elseif K==1 % i.e, some constructed variable of size [1 1 J I]
		cs = 0;
	else
		error(['the 2nd dimension of ' varname ' should be s_rho, s_w, or length 1.']);
	end
	% make (x3,y3,z3,cs3,mask3) matching the variable in all dimensions
	x3 = repmat(reshape(x2, [1 J I]), [K 1 1]);
	y3 = repmat(reshape(y2, [1 J I]), [K 1 1]);
	mask3 = repmat(reshape(mask2, [1 J I]), [K 1 1]);
	zeta3 = repmat(reshape(zeta2,[1 J I]),[K 1 1]);
	H3 = repmat(reshape(H2,[1 J I]),[K 1 1]);
	cs3 = repmat(cs(:),[1 J I]);
	z3 = zeta3 + cs3.*(zeta3+H3);
	% below, we'll make zm,ym,xm,... to match the _extraction_ in size
	
	% do the extraction ----------------------------------
	switch theType
		case 'full' % ----- 4D, whole variable
			data = nc_varget(filename,varname);
			zm = z3;
			ym = y3;
			xm = x3;
			data(mask3==0) = nan;
			
		case 'surface' % ----- 4D, surface field
			data = nc_varget(filename, varname, [0 K-1 0 0], [1 1 -1 -1]);
			zm = zeta2;
			ym = y2;
			xm = x2;
			data(mask2==0) = nan;
			
		case {'zslice','zslices'} % ----- 4D, slice at one or more z values
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
			
		case {'profile','profiles'} % ----- 3D, vertical profiles at one or more (y,x) pairs
			y = varargin{1}(:);
			x = varargin{2}(:);
			if length(x)==1 & length(y)>1
				x = repmat(x,size(y));
			elseif length(x)>1 & length(y)==1
				y = repmat(y,size(x));
			end
			L = length(x);
			cs = unique([-1; sort(cs); 0]); % extend to top & bottom of water column
			Km = length(cs);
			ym = repmat(reshape(y,[1 L]),[Km 1]);
			xm = repmat(reshape(x,[1 L]),[Km 1]);
			zetam = interp2(x2,y2,zeta2,xm,ym);
			Hm = interp2(x2,y2,H2,xm,ym);
			csm = repmat(cs,[1 L]);
			zm = zetam + csm .* (zetam + Hm);
			data3 = squeeze(nc_varget(filename,varname));
			data = interpn(cs3,y3,x3,data3,csm,ym,xm);
		
		case {'point','points'} % ----- 4D, data at one or more (z,y,x) triplets
			zm = varargin{1}(:);
			ym = varargin{2}(:);
			xm = varargin{3}(:);
			L = max([length(xm) length(ym) length(zm)]);
			if length(xm)==1, xm = repmat(xm,[L 1]); end
			if length(ym)==1, ym = repmat(ym,[L 1]); end
			if length(zm)==1, zm = repmat(zm,[L 1]); end
			zetam = interp2(x2,y2,zeta2,xm,ym);
			Hm = interp2(x2,y2,H2,xm,ym);
			csm = (zm - zetam) ./ (zetam + Hm);
			bad = isnan(csm) | csm > 0 | csm < -1;
			csm(isnan(csm)) = 0;
			data3 = squeeze(nc_varget(filename,varname));
			data = interpn(cs3,y3,x3,data3,csm,ym,xm);
			data(bad) = nan;
			
		otherwise
			error(['don''t know how to extract ''' theType ''' for 4D variables.']);
	end
	
	% clean up outputs ------------------------------------
	data = squeeze(data);
	coords.zm = squeeze(zm);
	coords.ym = squeeze(ym);
	coords.xm = squeeze(xm);
	
end
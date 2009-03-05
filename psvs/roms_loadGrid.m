function G = roms_loadGrid(A);

% G = roms_loadGrid(filename);
% G = roms_loadGrid(seriesDef);

if ischar(A)
	fname = A;
elseif isfield(A,'dirname') & isfield(A,'basename') & isfield(A,'ncn')
	fname = roms_filename([A.dirname A.basename],A.ncn(1));
else
	error('roms_loadGrid requires either a filename or a seriesDef.');
end

G.lon = nc_varget(fname, 'lon_rho');
G.lat = nc_varget(fname, 'lat_rho');
G.lonu = nc_varget(fname, 'lon_u');
G.latu = nc_varget(fname, 'lat_u');
G.lonv = nc_varget(fname, 'lon_v');
G.latv = nc_varget(fname, 'lat_v');
G.cs = nc_varget(fname, 'Cs_r');
G.csw = nc_varget(fname, 'Cs_w');

G.mask = nc_varget(fname,'mask_rho');
G.masku = nc_varget(fname,'mask_u');
G.maskv = nc_varget(fname,'mask_v');
G.H = nc_varget(fname, 'h');
G.H(G.mask==0) = nan;

G.K = length(G.cs);
[G.J, G.I] = size(G.lon);
function F = roms_extractRiseForcing(SD);

% forcing = roms_extractRiseForcing(seriesdef);
%
% returns a structure containing forcing time series relevant to RISE 
% for a series of ocean_his_#.nc files
%
% neil banas oct 2006



% where to extract the tide timeseries
tide_lon = -124.1;
tide_lat = 46.25;


for ni=1:length(SD.ncn) % for each nc file in the list
	n = SD.ncn(ni);
	nc = netcdf(roms_filename([SD.dirname SD.basename],n));
	
	% grid info; round tide_lon,lat -> tide_i,j
	if ni==1
		lon = nc{'lon_rho'}(:);
		lat = nc{'lat_rho'}(:);
		mask = nc{'mask_rho'}(:);
		[ii,jj] = meshgrid(1:size(lon,2), 1:size(lon,1));
		tide_i = round(interp2(lon,lat,ii,tide_lon,tide_lat));
		tide_j = round(interp2(lon,lat,jj,tide_lon,tide_lat));
	end
	
	% tidal height
	eta1 = nc{'zeta'}(1,tide_j,tide_i);
	F.eta(ni) = eta1;
	
	% u wind stress
	sustr1 = nc{'sustr'}(:);
	fwet = find( (mask(:,1:end-1)+mask(:,2:end))./2 > 0 );
	F.sustr(ni) = mean(sustr1(fwet));
	F.sustr_rms(ni) = std(sustr1(fwet));

	% v wind stress
	svstr1 = nc{'svstr'}(:);
	fwet = find( (mask(1:end-1,:)+mask(2:end,:))./2 > 0 );
	F.svstr(ni) = mean(svstr1(fwet));
	F.svstr_rms(ni) = std(svstr1(fwet));
	
	% shortwave radiation
	swrad1 = nc{'swrad'}(:);
	fwet = find(mask == 1);
	F.swrad(ni) = mean(swrad1(fwet));
	F.swrad_rms(ni) = std(swrad1(fwet));
		
	close(nc);
end
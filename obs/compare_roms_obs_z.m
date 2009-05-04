function [Obs, Roms, Forcing] = compare_roms_obs_z(seriesDef, obsFile, pos, z)
%------------------------------------------------
%
% script to compare roms output and ctd data
%------------------------------------------------

out_dir = seriesDef.dirname; %OUT directory with _his files
r_dir = out_dir(1:end-4); % run directory with grid, forcing files

if nargin<3
    % pick out locations near ORCA buoy at Duckabush on shelf
    pos = [-123 47.5533];
    z = 20; 
elseif nargin < 4
    z = 20;
end

% create time span of interest, variables of interest 
tspan = [min(seriesDef.nctime) max(seriesDef.nctime)];
v = [{'salinity'},{'temperature'}];
radius = 3; %km radius to search around positions

% extract obs from locations above 
n = length(pos(:,1)); % number of stations
for i = 1:n
    obs_1 = obs_extract(obsFile, v, tspan, z, [pos(i,2) pos(i,1) radius]);
    good = 1:length(obs_1.time);
    Obs(i).salt = obs_1.salinity; Obs(i).temp = obs_1.temperature; 
    Obs(i).td = obs_1.time; Obs(i).pos = [obs_1.longitude(1) obs_1.latitude(1)];
    % use first cast only if more than 1 comes back
end
Obs.tdepth = z;

% extract t/s series from model output
for i = 1:n
    [salt_1,coords_1] = roms_extract(seriesDef,'salt',obs_1.time,'point',z,pos(i,2),pos(i,1));
    [temp_1,coords_1] = roms_extract(seriesDef,'temp',obs_1.time,'point',z,pos(i,2),pos(i,1));
    Roms(i).salt = salt_1; Roms(i).temp = temp_1; Roms(i).coords = coords_1; 
    Roms(i).pos = [Roms(i).coords.xm(1) Roms(i).coords.ym(1)]; Roms(i).td = coords_1.tm;
    Roms(i).dirname = seriesDef.dirname;  Roms(i).basename = seriesDef.basename;
end

% extract forcing from model output
%%%rivers
 riverfile = [r_dir,'rivers.nc'];
 rt = nc_varget(riverfile, 'river_time');
    tvec = datevec(Obs(1).td); tyear = tvec(1);
 rt = datenum(tyear,1,1)-1+rt/86400;
 rgood = find(rt>tspan(1) & rt<tspan(2));
 Forcing.rivertime = rt(rgood);
 Q = nc_varget(riverfile, 'river_transport');
 %%% determine which river to use
  if mean(pos(:,1))<-123.9; rindex = 4; rname = 'Columbia';
  elseif mean(pos(:,2))>48.24; rindex = 12; rname = 'Fraser';
  elseif mean(pos(:,2))<47.6 && mean(pos(:,1))<-122.9; rindex = 10; rname = 'Skokomish';
  else rindex = 1; rname = 'Skagit';
  end
 Forcing.riverQ = abs(Q(rgood,rindex));    
 Forcing.rivername = rname;
%%% winds
 windfile = [r_dir,'Atm/Vwind.nc'];
 vtime = nc_varget(windfile, 'wind_time');
 vtime = datenum(tyear,1,1)+vtime/86400;
 vpos = mean(pos,1);
 latr = nc_varget([r_dir 'grid.nc'],'lat_rho');
 lonr = nc_varget([r_dir 'grid.nc'],'lon_rho');
 Vwind = nc_varget(windfile,'Vwind');
 for i=1:length(vtime)
     vwind(i) = interp2(lonr,latr,squeeze(Vwind(i,:,:)),vpos(1),vpos(2));
 end
 %now interpolate to time period of interest (tspan)- 6-hourly for now
 vwind = interp1(vtime, vwind, tspan(1):1/4:tspan(2));
 vstress = stresslp(vwind, 10); %10-m wind height-> stress using Large and Pond
 vstress(vwind<0) = -vstress(vwind<0);
 Forcing.vstress = vstress;
 Forcing.time = tspan(1):1/4:tspan(2);
 Forcing.pos = vpos; %lon/lat of wind stress measurements
 
if(0) %plot option
    compare_roms_obs_z_plot(Obs, Roms, Forcing);
end


    
    
            
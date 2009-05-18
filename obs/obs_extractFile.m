function data = obs_extractFile(filename, vars, timeRange, varargin)
%------------------------------------------------------------
% data = obs_extractFile(filename, vars);
%                                  ..., timeRange);
%                                           ..., LatLonVector);
%                                           ..., z, LatLonVector);
%
% routing called by obs_extract to get out observational data. see that
% mfile for full syntax...
%
% written by D. Sutherland and N. Banas, UW, Spring 2009
%------------------------------------------------------------

do_all = 0; do_depth =0;
vinfo=nc_info(filename); ndim = length(vinfo.Dimension); %2-D if ndim=2 (e.g., ORCA)
if ischar(vars); %then either 1 variable or all
    if strcmp(vars,'all')
        do_all = 1;
        varlist =[{'temperature'};{'salinity'};{'density'};{'oxygen'};{'oxygen_saturation'};...
            {'fluorescence'};{'nitrate'};{'par'}];
        for i=1:length(vinfo.Dataset) %get all variables out of netcdf file
            varlist2(i)={vinfo.Dataset(i).Name};
        end
        for i=1:length(varlist2);
            tf=strcmp(varlist,varlist2(i)); goodvar(i) = sum(tf);
        end
        varlist = varlist2(find(goodvar==1));
    else
        varlist = {vars}; 
    end;
else
    varlist = vars;
end

% read in station time and location first
    data.time = nc_varget(filename, 'time');
    data.latitude = nc_varget(filename, 'latitude');
    data.longitude = nc_varget(filename, 'longitude');
    data.pressure = nc_varget(filename, 'pressure');
    [Mm,Nn]=size(data.latitude); %Nn is 1 if 1-D, >1 if 2-D, etc. 
    data.time = repmat(data.time,1,Nn); 
    if(Nn>1);data.pressure = repmat(data.pressure(:)',Mm,1); end
    good = 1:Mm; % set this in case no timerange or lat/lon's given
    
% now truncate according to time and/or distance
if nargin > 2
    if ~ischar(timeRange) %i.e., not 'all', then find data only within time range specified
        if(length(timeRange) ~= 2) 
            error('time range needs to be a 2-element vector in MATLAB Datenum format');
        end
        intime = data.time(:) < timeRange(2) & data.time(:) > timeRange(1);
        good = find(intime == 1);
        if(isempty(good)); error(['NO Data in this file for time range given!']); end
    else  
        intime = ones(data.time(:));
        good = (1:(Mm*Nn))'; %otherwise, get all timesteps
    end
    
    % now see if position data exists
    if nargin > 3
        nvarargin = length(varargin); 
        if nvarargin == 2; type = 3;  %do both depth and position
             npos = length(varargin{2}); 
         elseif nvarargin == 1 & length(varargin{1}(1,:))==1; type = 2; %do depth only
         else type = 1; %do pos only
             npos = length(varargin{1}(1,:));
        end  
       switch type
            case 1 %given position only
                poly = varargin{1};
                if npos == 2 %if 2 then a polygon is given
                  inpoly = inpolygon(data.longitude(:), data.latitude(:), poly(:,2), poly(:,1));
                elseif npos == 3 %if 3 then a lat/lon/radius is given
                  ring = make_range_ring(poly(2), poly(1), poly(3));
                  inpoly = inpolygon(data.longitude(:), data.latitude(:), ring(:,1), ring(:,2));
                end
                good = find(inpoly == 1 & intime == 1);
            case 2 %this means depth given only
                do_depth = 1;
                tdepth = varargin{1};
                pp = data.pressure; 
                dd = sw_dpth(pp,mean(data.latitude(:))); 
            case 3 %given position and depth
                do_depth = 1;
                tdepth = varargin{1};
                pp = data.pressure; 
                dd = sw_dpth(pp,mean(data.latitude(:)));
                poly = varargin{2};
                if npos == 2 %if 2 then a polygon is given
                  inpoly = inpolygon(data.longitude(:), data.latitude(:), poly(:,2), poly(:,1));
                elseif npos == 3 %if 3 then a lat/lon/radius is given
                  ring = make_range_ring(poly(2), poly(1), poly(3));
                  inpoly = inpolygon(data.longitude(:), data.latitude(:), ring(:,1), ring(:,2));
                end
                good = find(inpoly == 1 & intime == 1);    
            otherwise
                error(['need to specify one or both of depth, position'])
       end %end switch type
    end %end nargin > 3
end %end nargin > 2
% add info variables to list if not doing all
do_castid = nc_isvar(filename,'castid');
if ~do_all
    varlist = ([varlist(:);{'latitude'};{'longitude'};{'pressure'}]); 
    if(do_castid); %get cast id's if there, i.e. in ctd files only
        varlist = ([varlist(:);{'castid'}]);
    end
end
nv = length(varlist);

% now loop over each variable
for j = 1:nv
    varname = char(varlist(j));
    if ~nc_isvar(filename,varname) % check to see if variable in the netcdf file
        error(['variable ''' varname ''' not found in ' filename '.']);
    end
    % now put each variable into data structure, all for now, then truncate later
    disp(['    extracting ' varname ' from ' filename])
    datavec = nc_varget(filename, varname);
    if(strcmp(varname,'pressure')); %if variable is pressure, use from before
        datavec = data.pressure;
    end
    %%% interpolate to depth
    if(do_depth)
      data_int = []; tint = [];
      datavec = datavec(good); tgood = data.time(good); dint = dd(good);
      if(Nn>1); ntime = length(tgood)/Nn;
         tgood = reshape(tgood,ntime,Nn); 
         datavec = reshape(datavec,ntime,Nn);
         dint = reshape(dint,ntime,Nn);
         for jj = 1:ntime
            ind = find(~isnan(datavec(jj,:))); if(isempty(ind));ind=1:Nn;end
            data_int(jj)=interp1(dint(jj,ind),datavec(jj,ind),tdepth);              
            tint(jj) = tgood(jj,1);
         end
      else
        [junk, first, junk2] = unique(tgood, 'first');
        [junk, last, junk2] = unique(tgood, 'last');
        data.time_index = sort([first(:) last(:)]);
        ntime = length(data.time_index(:,1));
        for jj = 1:ntime
            ind = data.time_index(jj,1):data.time_index(jj,2);
            data_int(jj)=interp1(dint(ind),datavec(ind),tdepth);              
            tint(jj) = tgood(ind(1));
        end
      end
      datavec = data_int; 
      data.(varname) = datavec;
    else
      datavec = datavec(good);
      data.(varname) = datavec;
    end
end
if do_depth; data.time = tint; 
 else data.time =data.time(good); 
end
disp(['...DONE'])

% now make indices for easier looking at data
if(do_castid)
    [junk, data.cast_index, junk2] = unique(data.castid, 'first');
end
    [junk, data.time_index, junk2] = unique(data.time, 'first');
    position = [data.longitude data.latitude];
    [junk, data.pos_index, junk2] = unique(position, 'rows', 'first');
    data.pos_index = sort(data.pos_index);




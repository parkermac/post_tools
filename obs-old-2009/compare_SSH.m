function compareSSH = compare_SSH(seriesDef, startTime, numDays)
%****************************************************
%
% compareSSH = compare_SSH(seriesDef, startTime, numDays)
%
% uses: to compare ROMS model output (zeta) to tidal heights calculated from PS_tides
%  model at 5 segments for some arbitrary time range (must be covered by model time range)
%
% needs: PS_tides toolbox, roms_extract toolbox (and coastline file for plotting)
%
% INPUT 
%    seriesDef: series definition created by roms_createSeriesDef.m
%    startTime: starting date (datenum format)- will start on previous day to be able
%               to interpolate to model time
%    numDays: number of days to compare (14 is default)
%
% OUTPUT
%    compareSSHis structure containing:
%       zeta: model SSH at 5 locations
%       pos: lon/lat of stations interpolated to
%       tide_ps: SSH from PS_tides (on same time base as model, assume PST)
%       tt: time vector (datenum format)    
%       WS: Wilmott score based on time period, 5 x 1 vector
%       region: name of region WS refers to (fixed for now)
%
% written by: D. Sutherland 4/17/2009
%****************************************************

if(nargin < 3)
    numDays = 14; %set default to 14 days
end

% make time vector, assume hourly outputs
tt = floor(startTime-1):1/24:(floor(startTime)+numDays);
tt(end) = []; % get rid of last time because PS_tides doesn't use it

% get PS_tides predictions and locations of segments
segs = [416; % Admiralty Inlet
        360; % Whidbey Basin
        512; % Hood Canal
        403; % Main Basin
        118];% South 
lats = [48.0824; 48.2367; 47.7551; 47.6739; 47.1704];
lons = [-122.6477; -122.5872; -122.7385; -122.4632; -122.7904];
region=[{'Admiralty Inlet'};{'Whidbey Basin'};{'Hood Canal'};{'Main Basin'};{'South Sound'}];
for i=1:5
    Tide(i,:) = PStide_at_seg(segs(i), datestr(tt(1)), numDays+1, 1);
    Tide(i,:) = Tide(i,:)-mean(Tide(i,:)); %demean for comparison to model SSH
end

% extract zeta from model for time specified
t_model = startTime:1/24:(startTime+numDays); 
t_model(end) = [];
[zeta_full, coords] = roms_extract(seriesDef, 'zeta', t_model, 'full');
% now interpolate to segments
zeta = NaN*ones(5, length(t_model)); %pre-allocate zeta
lonr = squeeze(coords.xm(1,:,:)); latr = squeeze(coords.ym(1,:,:)); 
for i=1:5
    for j = 1:length(t_model)
        zeta(i,j) = interp2(lonr, latr, squeeze(zeta_full(j,:,:)), lons(i), lats(i));
    end
end
    
% interpolate onto same time base now
for i=1:5
    tide_ps(i,:)=interp1(tt, Tide(i,:), t_model-7/24);
end

% calculate Wilmott scores
for i=1:5;
    MSE(i)=mean((zeta(i,:)-tide_ps(i,:)).^2);
    denom(i)=mean((abs(zeta(i,:)-mean(tide_ps(i,:)))+abs(tide_ps(i,:)-mean(tide_ps(i,:)))).^2);
end
WS = 1 - MSE./denom;

% make structure and plot some things
compareSSH.zeta = zeta; compareSSH.tide_ps = tide_ps;
compareSSH.region = region; compareSSH.pos = [lons lats];
compareSSH.tt = t_model; compareSSH.WS = WS;

if(1) %to plot things
    % plot comparison over time of interest
    figure; set(gcf,'position',[163 278 1114 678]);
    for i=1:5
        subplot(5,1,i)
        plot(compareSSH.tt, compareSSH.tide_ps(i,:))
        hold on; grid on;
        plot(compareSSH.tt, compareSSH.zeta(i,:),'r')
        axis([compareSSH.tt(1) ceil(compareSSH.tt(end)) -3 3])
        set(gca,'xtick',compareSSH.tt(1):ceil(compareSSH.tt(end)));
        datetick('x',6,'keeplimits','keepticks')
        title(['PStides (blue) vs. Model SSH (red), ' char(compareSSH.region(i)) ', WS = '...
            num2str(compareSSH.WS(i),'%2.2f')], 'fontsize', 14);
    end
    % plot map to show stations
    figure; set(gcf,'position',[167 249 643 705]);
        plot_WAcoast('detailed','linewidth', 2);
        hold on; grid on; 
        axis([-123.2 -122.1 47 48.4])
        scatter(compareSSH.pos(:,1), compareSSH.pos(:,2), 100, compareSSH.WS, 'filled');
        ecolorbar([0.5:.05:1]);
        set(gca,'fontsize',14,'box','on','tickdir','out')
end



function plot_WAcoast(type, varargin)
%-----------------------------
%
% plot_WAcoast(type)
%
% this adds a coastline to the current plot of type:
% 'detailed': loads coastline_detailed.mat which is good for Puget Sound
% 'regional': loads coastline_regional.mat which is good for WA coast and S of Georgia
%
% varargin are line specs to pass to plot 
%
% ex) 'linewidth',2,'color','k' (default is black)
%-----------------------------

if nargin < 1 %default to regional
    type = 'regional';
end

% define where these coastline files exist
dirname = '/Users/daves/Desktop/PugetSound/R3_tools/In_data/Topo/coastlines/';
switch type
    case 'detailed'
        filename = [dirname, 'coastline_detailed.mat'];
        load(filename);
    case 'regional'
        filename = [dirname, 'coastline_regional.mat'];
        load(filename);
    otherwise
        error(['specify either detailed or regional or nothing'])
end

hold on;
plot(lon_coast, lat_coast, 'color', 'k', varargin{:})



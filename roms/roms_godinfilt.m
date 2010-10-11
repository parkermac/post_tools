function Sout = roms_godinfilt(Sin, outputTimebase, suffix)

% seriesDefOut = roms_godinfilt(seriesDefIn, outputTimebase, suffix);
% seriesDefOut = roms_godinfilt(seriesDefIn, suffix);
%
% filters all the variables in a file series using a 24-24-25 hr window and saves
% a new series of the results (Godin filter).
%
% NOTE: Sin must be hourly saves and default is the 24-24-25 hr window, with one output file
% per day, centered on noon.
%
% neil banas feb 2009, modified by DAS, Oct. 2010
%-------------------------------------------------------------------------

if nargin<3
	suffix = outputTimebase;
	outputTimebase = [];	
end

% Godin filter is 71 hours long
filterWindow = 71/24; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build weights for Godin filter- start with 24 hr. filter
 filter24 = ones(24,1); 
 filter24 = filter24 ./ sum(filter24);
 
% then build 25 hr filter
 filter25 = ones(25,1); 
 filter25 = filter25 ./ sum(filter25);
 
% covolve filters together, works because conv is associative
 temp_filter = conv(filter24, filter24);
 filter = conv(temp_filter, filter25);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if length(Sin.ncn)<2, error('no file series specified.'); end
t0 = Sin.nctime(1) + filterWindow/2;
t1 = Sin.nctime(end) - filterWindow/2;
if isempty(outputTimebase)
	% default: one file per day at noon
	outputTimebase = ceil(t0)+0.5 : t1;
else
	% make sure the output times requested fall within the file series,
	% with a margin for half the filter window on each end
	outputTimebase = outputTimebase(outputTimebase >= t0 & outputTimebase <= t1);
end

iAvg = 1;
for i = 1:length(outputTimebase)
	ti = outputTimebase(i);
	t0 = ti - filterWindow/2;
	t1 = ti + filterWindow/2;
	f = find(Sin.nctime >= ti-filterWindow/2 & Sin.nctime <= ti+filterWindow/2);
	weights = interp1(linspace(t0,t1,length(filter)), filter, Sin.nctime(f));
	nn = Sin.ncn(f);
	outname = roms_filename([Sin.dirname Sin.basename suffix '_'], iAvg);
	roms_averageFileSet(Sin, nn, outname, weights);
	iAvg = iAvg+1;
end

Sout = roms_createSeriesDef(Sin.dirname, [Sin.basename suffix '_']);	
function roms_averageFileSet(seriesIn, ncn, filename, weights);

% roms_averageFileSet(seriesIn, frameNums, outputFilename);
%                                                    ...,weights);
%
% averages a set of netcdf files into a single netcdf file with the same format.
% looks for files with numbers _ncn_, which should be a subset of seriesIn.ncn.
% All variables are averaged, except for those of type nc_char.
% 
% a vector of weights the same length as frameNums can be passed, in order to
% use this routine for filtering.
%
% neil banas feb 2009

if nargin < 4
	weights = ones(size(frameNums)) ./ length(frameNums);
end

% assemble filenames
inNames = cell(length(ncn),1);
for ni = 1:length(ncn)
	n = ncn(ni);
	inNames{ni} = roms_filename([seriesIn.dirname seriesIn.basename], n);
	if ~exist(inNames{ni}), error(['can''t find ' inNames{ni}]); end
end
disp(['averaging ' seriesIn.dirname seriesIn.basename ' #' num2str(ncn(1)) ' to ' num2str(ncn(end))]);

% duplicate the first input file to create a template for the output
tmp = [seriesIn.dirname 'temp.nc'];
system(['cp ' inNames{1} ' ' tmp]);

% loop through variables inside
info = nc_info(tmp);
vars = info.Dataset;
for vi=1:length(vars)
	if length(vars(vi).Size) >= 2, disp(['    ' vars(vi).Name]); end
	if vars(vi).Nctype ~= nc_char
		% loop through files and make a weighted average of that variable
		A = weights(1) .* nc_varget(inNames{1}, vars(vi).Name);
		for ni = 2:length(ncn)
			A = A + weights(ni) .* nc_varget(inNames{ni}, vars(vi).Name);
		end
		% write the average into the output file
		nc_varput(tmp, vars(vi).Name, A);
	end
end

% if the output filename is a relative path, prepend the series' dirname
if filename(1) ~= '/' & filename(1) ~= '~'
	filename = [seriesIn.dirname filename];
end
% name the output file correctly to indicate that we're done
system(['mv ' tmp ' ' filename]);
disp(['    done: created ' filename]);
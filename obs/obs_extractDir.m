function data = obs_extractDir(files, vars, timeRange, varargin)
%------------------------------------------------------------
% data = obs_extractDir(files, vars);
%                             ..., timeRange);
%                                        ..., LatLonPoly);
%                                        ..., lat, lon, radius);
%
% routing called by obs_extract to get out observational data. see that
% mfile for full syntax...
%
% this file just runs through the files listed and use obs_extractFile
%
% written by D. Sutherland and N. Banas, UW, Spring 2009
%------------------------------------------------------------

nfiles = length(files);

for j =1:nfiles
    filename = files(j).name;
    if j == 1  %initialize array
        data = obs_extractFile(filename, vars, timeRange, varargin{:});
    else
        data2 = obs_extractFile(filename, vars, timeRange, varargin{:});
        data = catstruct(data, data2);
    end
end

   





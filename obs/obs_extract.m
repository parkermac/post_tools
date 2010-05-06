

function data = obs_extract(files, vars, timeRange, varargin)
%------------------------------------------------------------
% data = obs_extract(filename, vars, timeRange
%                                         ..., 'section',x,y,range);
%                                         ..., 'polygon',x,y);
%                                         ..., 'all');
% data = obs_extract(directory, ...
% data = obs_extract(cell array of files and directories, ...
%
% looks inside a netcdf file for one or more observational variables and
% returns them in a structure _data_, along with coordinate variables.
%
% if a directory is given instead of a filename, recurses over
% files inside it and concatenates the output.
%
%Only includes data from files with single dimensions, skips mooring data
%
% vars can be a string containing one variable name, or a cell array of strings,
% or 'all'.
%
% timeRange is a 2 element vector of matlab datenums or 'all'.
%

%
% written by C. Bassin ,D. Sutherland and N. Banas, UW, Jan 2010
%
%------------------------------------------------------------


[netcdfpaths,files]=dirwalker(files);

x=1;
for i=1:length(netcdfpaths);

    Info_file=nc_info(netcdfpaths{i});
   
    if size(Info_file.Dimension,1)==1
        new_netcdfpaths{x}=netcdfpaths{i};
        x=x+1;
    else
        disp(['skipping ' netcdfpaths{i} ' due to incorrect dimension size'])
    end

end

if exist('new_netcdfpaths','var')==0
    data=[];
    return
end

N=length(new_netcdfpaths);
dataC=cell(N,1);


for j=1:N;
  disp(['reading data from ' new_netcdfpaths{j}])

  dataC{j,1}= obs_extractFromFile(new_netcdfpaths{j}, vars, timeRange, varargin{:});
end 

if N>1

    data=obs_merge(dataC,'any');
else
    data= dataC{1,1};
end

% get rid of nan data in coordinate variables
data=obs_omit(data,'x == nan');
data=obs_omit(data,'y == nan');
data=obs_omit(data,'z == nan');
data=obs_omit(data,'t == nan');


data.cast=obs_identifyCasts(data);



end


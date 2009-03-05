function RD = roms_createRunDef(shortname, dirname, hisbasename, diabasename);

% RD = roms_createRunDef(shortname, dirname);
% RD = roms_createRunDef(shortname, dirname, hisbasename, diabasename);
%
% makes a structure containing the run shortname, dirname, seriesDefs
% for the _his and _dia files (if any _his or _dia files were found),
% and model grid
%
% neil banas feb 2009

if nargin==2
	hisbasename = 'ocean_his_';
	diabasename = 'ocean_dia_';
end

RD.shortname = shortname;
if dirname(end) ~= '/', dirname = [dirname '/']; end
RD.dirname = dirname;
his = roms_createSeriesDef(dirname, hisbasename);
if ~isempty(his.ncn), RD.his = his; end
dia = roms_createSeriesDef(dirname, diabasename);
if ~isempty(dia.ncn), RD.dia = dia; end
RD.grid = roms_loadGrid(RD.his);
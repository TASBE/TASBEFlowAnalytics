% function VERSION = tasbe_version()
%   returns a cell array of {major minor patch}
%   major and minor are numbers, patch is a string
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function version = tasbe_version()

try 
    file = fopen([fileparts(mfilename('fullpath')) '/../version.txt'],'r');
    raw_version = textscan(file,'%s');
    if(numel(raw_version) ~= 1), TASBESession.error('TASBE:Version','NotSingleVersion','Could not read single version number from version file.'); return; end;
    raw_version = char(raw_version{1});
    fclose(file);
catch e
    TASBESession.error('TASBE:Version','CannotReadVersionFile','Could not read version information from TASBE distribution version.txt');
    throw e;
end

strs = strsplit(raw_version,'.');
if(numel(strs) ~= 3), TASBESession.error('TASBE:Version','VersionNotSemVer','Version does not follow SemVer convention: %s',raw_version); end;

version = {str2double(strs{1}) str2double(strs{2}) strs{3}};

if(numel(version{1})+numel(version{2}) ~= 2), TASBESession.error('TASBE:Version','VersionNotSemVer','Version does not follow SemVer convention: %s',raw_version); end;

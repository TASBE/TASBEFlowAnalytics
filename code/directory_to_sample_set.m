% DIRECTORY_TO_SAMPLE_SET returns a set of sample FCS files from inputted
% directory 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function set = directory_to_sample_set(directory,exclusions,prefix)

if nargin<3,
    prefix = '';
end

fileset = dir(directory);
set = cell(0);
for i=1:numel(fileset),
    name = fileset(i).name;
    % filter out unwanted
    skip = false;
    if(isempty(strmatch(fliplr('.fcs'),fliplr(name)))), skip = true; end;
    for j=1:numel(exclusions),
        if(~isempty(strfind(name,exclusions{j}))), skip = true; end;
    end
    if skip, continue; end;
    % add to collection
    set{end+1,1} = [prefix fileset(i).name(1:(end-4))];
    set{end,2} = { [directory '/' fileset(i).name] };
end
if numel(set)==0,
    TASBESession.warn('TASBE:Utilities','NoFilesFound','Warning: no FCS files found in directory: %s',directory);
end

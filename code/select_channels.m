% SELECT_CHANNELS makes sure that the inputted raw fcs data and header
% directory matches with one of the inputted channel names.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [selected, selected_par] = select_channels(channels,rawfcs,fcshdr)

selected = zeros(size(rawfcs,1),numel(channels));
selected_par = cell(numel(channels),1);
for i=1:numel(channels)
    found = false;
    for j=1:numel(fcshdr.par)
        name=getName(channels{i});
        if(strcmp(name,fcshdr.par(j).name) || strcmp(name,fcshdr.par(j).rawname))
            selected(:,i)=rawfcs(:,j); 
            selected_par{i}=fcshdr.par(j);
            found=true; continue;
        end
    end
    if(~found), TASBESession.error('TASBE:Analysis','MissingChannel','Could not find channel %s',getPrintName(channels{i})); end
end

end

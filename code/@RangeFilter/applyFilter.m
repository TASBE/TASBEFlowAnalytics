% APPLYFILTER overrides applyFilter method in Filter class. Corresponds to
% a RangeFilter object. 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function data = applyFilter(RF,fcshdr,rawfcs)
% Discard data that doesn't pass the specified filters

if(numel(rawfcs)==0), data = rawfcs; return; end

switch RF.mode
    case 'And'
        retained = ones(size(rawfcs,1),1);
    case 'Or'
        retained = zeros(size(rawfcs,1),1);
    otherwise
        TASBESession.error('TASBE:RangeFilter','BadRangeFilterMode','Unrecognized range filter mode: %s',value); 
end

for i=1:numel(RF.channels)
    found = false;
    % switch based on fcshdr vs. Channel objects:
    if isstruct(fcshdr) % it's a header, not Channel objects:
        for j=1:numel(fcshdr.par)
            if(strcmp(RF.channels{i},fcshdr.par(j).name))
                channeldata=rawfcs(:,j);
                found=true; continue;
            end
        end
    else % assume it's a cell array of channel objects:
        for j=1:numel(fcshdr)
            if(strcmp(RF.channels{i},getName(fcshdr{j})))
                channeldata=rawfcs(:,j);
                found=true; continue;
            end
        end
    end
    if(~found), 
        TASBESession.error('TASBE:RangeFilter','MissingChannel','Could not find range filter channel: %s',RF.channels{i}); 
    else
        pass_values = channeldata >= RF.ranges(i,1) & channeldata <= RF.ranges(i,2);
        switch RF.mode
            case 'And'
                retained = retained & pass_values; % keep only those that pass all tests
            case 'Or'
                retained = retained | pass_values; % keep anything that passes any test
            otherwise
                TASBESession.error('TASBE:RangeFilter','BadRangeFilterMode','Unrecognized range filter mode: %s',value); 
        end;
    end
    
end

data = rawfcs(retained,:);

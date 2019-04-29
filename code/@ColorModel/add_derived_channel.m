% ADD_DERIVED_CHANNEL adds inputted filter to postfilters property for
% ColorModel object
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function CM=add_derived_channel(CM,channel_name,print_name,units)
% make sure there's not a duplicated name
if ~isempty(find(cellfun(@(c)(strcmp(channel_name,getName(c))),CM.Channels),1)),
    TASBESession.error('TASBE:ColorModel','DuplicateChannel','Drived channel given duplicate name %s',channel_name);
end

% add the channel as an unprocessed channel
C = Channel(channel_name,0,0,0,false,true);
C = setUnits(C,units);
C = setPrintName(C,print_name);
% Add the new channels
CM.Channels{end+1} = C; 
% add empty placeholders for all n-channel structures
CM.autofluorescence_model{end+1} = [];
CM.ColorFiles{end+1} = [];
CM.noise_model.noisemin(end+1) = 0;
CM.noise_model.noisemean(end+1) = 0;
CM.noise_model.noisestd(end+1) = 0;
CM.noise_model.detailcounts{end+1} = {};
CM.noise_model.detailmeans{end+1} = {};
CM.noise_model.detailstds{end+1} = {};

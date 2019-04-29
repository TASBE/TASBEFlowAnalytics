% REMAPCHANNELS creates a new ColorModel object with an updated Channels
% property with new print names. Also updates the ERF_channel. 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function CM = remapChannels(CM,new_names)
for i=1:numel(new_names)
    reset_ERF = (CM.Channels{i}==CM.ERF_channel);
    CM.Channels{i} = remapChannelName(CM.Channels{i},new_names{i});
    if reset_ERF, CM.ERF_channel = CM.Channels{i}; end
end

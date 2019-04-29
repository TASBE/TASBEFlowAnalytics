% getChannel returns channel with given an inputted AnalysisParameter and
% channel name
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function channel = getChannel(AP, name)
    channel_id = find(strcmp(AP.ChannelLabels(:,1),name));
    if(numel(channel_id)>1)
        TASBESession.error('TASBE:ColorModel','DuplicateChannels','Analysis specifies multiple "%s" channels',name);
    end
    if(numel(channel_id)==0)
        TASBESession.error('TASBE:ColorModel','MissingChannel','Analysis does not specify a "%s" channel',name);
    end
    channel = AP.ChannelLabels{channel_id,2};

% CHANNEL_NAMED returns the channel of a ColorModel object with a certain
% print name.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function channel = channel_named(CM,name)
    found=false;
    for i=1:numel(CM.Channels)
        % FIXME not sure if we want to count on these strings of
        % what?
        if(strcmp(getPrintName(CM.Channels{i}),name))
            channel = CM.Channels{i}; found = true; break;
        end;
    end;
    if(~found), TASBESession.error('TASBE:ColorModel','MissingChannel','Unable to find channel %s',name); end;
       

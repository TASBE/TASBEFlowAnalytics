% set_ERF_channel_name is a setter function for the ColorModel class.
% 
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function CM=set_um_channel_name(CM, v)
    CM.um_channel_name=v;
    found=false;
    for i=1:numel(CM.Channels), 
        if(strcmp(CM.um_channel_name,getName(CM.Channels{i}))), 
            CM.um_channel = CM.Channels{i}; 
            CM.Channels{i} = setUnits(CM.Channels{i},'Eum');
            found=true; 
            break; 
        end;
    end;
    if(~found), TASBESession.error('TASBE:ColorModel','MissingumChannel','Unable to find um channel %s',CM.um_channel_name); end;


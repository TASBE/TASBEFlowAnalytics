% FIND returns index of inputted channel within the Channels property of a
% ColorModel object.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function index = find(CM, channel)

foundset=[]; nameeqset=[];
for i=1:numel(CM.Channels)
    if(CM.Channels{i} == channel)
        foundset(end+1) = i; 
        if(strcmp(getName(CM.Channels{i}),getName(channel)))
            nameeqset(end+1) = i;
        end
    end;
end;

if(isempty(foundset)), 
    TASBESession.error('TASBE:ColorModel','MissingChannel','Unable to find channel %s',getName(channel)); 
elseif numel(foundset)==1
    index = foundset(1);
elseif numel(nameeqset)==1 % more than one match, but only one matches name precisely
    % Turning this warning off, since it gets called too much and isn't very useful
    %TASBESession.warn('TASBE:ColorModel','DisambiguateChannel','Multiple channels match %s, discriminating by name',getName(channel)); 
    index = nameeqset(1);
else
    TASBESession.error('TASBE:ColorModel','MultipleChannels','Multiple channels match %s, and cannot discriminate by name',getName(channel));
end;
       
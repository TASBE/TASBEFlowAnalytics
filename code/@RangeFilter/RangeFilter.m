% Constructor of RangeFilter class
% Optional discarding of data outside of a certain range
% Specified as a sequence of 'ChannelA',[minA maxA],'ChannelB',[minB maxB]
% Also, argument 'Mode' can be 'And' (a value is excluded if any channel is outside), or 'Or' (a value is excluded if all are outside
% The default mode is 'And'
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function RF = RangeFilter(varargin)
RF.mode = 'And';
RF.channels = {};
RF.ranges = [];

for i=1:2:numel(varargin),
    arg = varargin{i};
    value = varargin{i+1};
    if ~ischar(arg), TASBESession.error('TASBE:RangeFilter','BadRangeFilterArgument','Range filter argument %i was not a string',i); end;
    
    % Either 'Mode' or a channel name
    if strcmp(arg,'Mode')
        if ~(strcmp(value,'And') || strcmp(value,'Or')), TASBESession.error('TASBE:RangeFilter','BadRangeFilterMode','Unrecognized range filter mode: %s',value); end;
        RF.mode = value;
    else
        if numel(value)~=2, TASBESession.error('TASBE:RangeFilter','BadRange','Range is not a 2-element vector: %s',num2str(value)); end;
        RF.channels{end+1} = arg;
        RF.ranges(end+1,1:2) = value;
    end
end

RF = class(RF,'RangeFilter',Filter());

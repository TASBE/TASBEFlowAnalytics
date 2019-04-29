function [names, hdr] = fcs_channel_names(filename)
% [names, hdr] = fcs_channel_names(filename): 
%   Returns a cell-array of channel names in an FCS file
%   The second return is the set of headers that these names have been extracted from 

% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

[~, hdr] = fca_read(filename);
names = {hdr.par(:).name};


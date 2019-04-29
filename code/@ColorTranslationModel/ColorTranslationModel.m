% Constructor for the ColorTranslationModel class with properties: 
%  Channels        % array of Channel objects
%  scales          % Channel x Channel scalars: Color_j = scales(i,j)*Color_i
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function CTM = ColorTranslationModel(channels, scales)
    if nargin == 0
        channels{1} = Channel();
        CTM.Channels = channels;
        CTM.scales = zeros(3,3);
    elseif nargin == 2
        CTM.Channels = channels;
        CTM.scales = scales;
    end
    CTM=class(CTM,'ColorTranslationModel');
    


% EQ returns true if inputted channels are equal to each other.
%
% The nice thing about this equals method is that it does not look at
% the name thus we can find matches for channels even if the channel
% has a different name (as long as we have the laser and filter
% information).
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function y = eq(C, b)
    y = (strcmp('Channel', class(b)) && (C.FilterCenter == b.FilterCenter) && ...
        (C.FilterWidth == b.FilterWidth) && (C.Laser==b.Laser));

% GEOSTD returns the geometric standard deviation of inputted data with
% optional weights.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function s = geostd(data,weights)
    if (isempty(data))
        s = NaN;
    else
        if nargin < 2
            s = 10.^std(log10(data));
        else
            valid = weights>0;
            s = 10.^sqrt(var(log10(data(valid)),weights(valid)));
        end
    end
end

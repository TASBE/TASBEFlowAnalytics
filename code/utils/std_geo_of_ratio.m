% Geometric std.dev. of ratio of independent random variables
% (i.e. assuming covariance = 0)
% No place for covariance here because I'm not sure how to translate it correctly
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function sr = std_geo_of_ratio(num_mean,denom_mean,num_std,denom_std)
sr = 10.^(std_of_ratio(log10(num_mean),log10(denom_mean),log10(num_std),log10(denom_std)));


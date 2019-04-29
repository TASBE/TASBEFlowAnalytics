% std.dev. of ratio of independent random variables
% (i.e. assuming covariance = 0)
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function sr = std_of_ratio(num_mean,denom_mean,num_std,denom_std)
% Taylor expansion for moments of random variables, taken from:
% http://en.wikipedia.org/wiki/Taylor_expansions_for_the_moments_of_functions_of_random_variables
var_of_ratio = @(num_mean,denom_mean,num_var,denom_var,cov) ...
    (num_var ./ (denom_mean.^2)) - ...
    ((2*num_mean ./ (denom_mean.^3)) .* cov) + ...
    ((num_mean.^2 ./ (denom_mean.^4)) .* denom_var);

sr = sqrt(var_of_ratio(num_mean,denom_mean,num_std.^2,denom_std.^2,0));

% P_GAUSSIAN(low,high,mean,std): compute the probability of falling into the interval
%  [low,high] on a gaussian distribution with the specifiec mean
%  and standard deviation (default 0,1)
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function prob = p_gaussian(low,high,mean,std) 
if nargin <= 2, mean=0; end;
if nargin <= 3, std=1; end;

cum_low = 0.5*(1+erf((low-mean)/(std*sqrt(2))));
cum_high = 0.5*(1+erf((high-mean)/(std*sqrt(2))));
prob = cum_high-cum_low;

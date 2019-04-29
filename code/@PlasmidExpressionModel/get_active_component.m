% GET_ACTIVE_COMPONENT returns information (mean, variance and weight) of
% the active component of a PlasmidExpressionModel object.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function y = get_active_component(PEM)
   y.mu = PEM.fp_dist.mu(PEM.active_component);
   y.Sigma = PEM.fp_dist.Sigma(PEM.active_component);
   y.weight = PEM.fp_dist.weight(PEM.active_component);



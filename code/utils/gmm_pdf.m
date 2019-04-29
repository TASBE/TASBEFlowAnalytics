% GMM_PDF returns the probability density function for a fp_dist obtained
% from a PlasmidExpressionModel object.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function p = gmm_pdf(gmm,x)
% Note: "Sigma" in the GMM indicates the *variance* and not the std.dev, as would usually be the case.
p = zeros(size(x));
for i=1:numel(gmm.mu),
    p = p + gmm.weight(i) * 1/(sqrt(gmm.Sigma(1,1,i))*sqrt(2*pi)) * exp( -0.5 * (((x-gmm.mu(i))/sqrt(gmm.Sigma(1,1,i))).^2));
end

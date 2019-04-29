%%%% WARNING: ESTIMATES EXPERIMENTAL AND UNTRUSTWORTHY
%%%% DISABLE AND DEPRECATE
% ESTIMATE_PLASMIDS is currently kept to run Gaussian mixture bimodal fit
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function plasmids = estimate_plasmids(PEM,ERFs)
    plasmids = zeros(size(ERFs));
    %%%% WARNING: ESTIMATES EXPERIMENTAL AND UNTRUSTWORTHY
    %%%% DISABLE AND DEPRECATE
%     for i=1:numel(ERFs)
%         if isnan(ERFs(i)), plasmids(i) = NaN; continue; end;
%         binCenters = get_bin_centers(PEM.bins);
%         idx = find(binCenters(1:end-1)<=ERFs(i) & binCenters(2:end)>=ERFs(i),1);
%         ratio = log10(ERFs(i)/binCenters(idx))/log10(get_bin_widths(PEM.bins));
%         plasmids(i) = (1-ratio)*PEM.estimated_plasmids(idx) + ratio*PEM.estimated_plasmids(idx+1);
%     end

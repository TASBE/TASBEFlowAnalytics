% ESTIMATE_FRACTION_ACTIVE computes the fraction of active cells in a data
% sample for a PlasmidExpressionModel object.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function active = estimate_fraction_active(PEM,ERFs)
    active = zeros(size(ERFs));
    if numel(PEM.fraction_active)<2
        TASBESession.warn('TASBE:PEM','NotBimodalGaussian','Cannot compute fraction active: distribution did not fit bimodal gaussian');
        return
    end
    for i=1:numel(ERFs)
        if isnan(ERFs(i)), active(i) = NaN; continue; end;
        binCenters = get_bin_centers(PEM.bins);

	 if binCenters(1) > ERFs(i)
		idx = 1;
	elseif binCenters(end) < ERFs(i)
		idx = length(binCenters);
	else
		idx = find(binCenters(1:end-1)<=ERFs(i) & binCenters(2:end)>=ERFs(i),1);
	end
	
	ratio = log10(ERFs(i)/binCenters(idx))/log10(get_bin_widths(PEM.bins));
        active(i) = (1-ratio)*PEM.fraction_active(idx) + ratio*PEM.fraction_active(idx+1);
    end
    

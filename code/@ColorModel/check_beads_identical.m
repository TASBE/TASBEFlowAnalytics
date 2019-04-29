%  [ok ratios] = check_beads_identical(CM, file, tolerance)
%  Loading 'file' as FCS data, tests whether the bead peaks found are
%  within tolerance of the bead peaks for the color model.
%  Tolerance is an optional argument, defaulting to 0.5
%
%  Returns a boolean 'ok' indicating whether all chanels passed
%  and a set 'ratios' of the relative shift between peak sets
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.
function [ok ratios] = check_beads_identical(CM, file, tolerance)
if nargin < 4, tolerance = 0.05; end; % 5 percent tolerance by default
alt_units = beads_to_ERF_model(CM, file);
ok = true;
ratios = zeros(size(CM.Channels));
for i=1:numel(CM.Channels),
    peaks = get_peak(CM.unit_translation, i);
    alt_peaks = get_peak(alt_units, i);
    if(numel(peaks) ~= numel(alt_peaks)),
        TASBESession.warn('TASBE:CompareBeads','PeakNumberMismatch','Number of peaks does not match on channel %s: %d vs. %d',getPrintName(CM.Channels{i}),numel(peaks),numel(alt_peaks));
        ok=false;
        ratios(i) = NaN;
    else
        ratioset = zeros(size(peaks));
        for j=1:numel(peaks),
            ratio = alt_peaks(j)/peaks(j);
            ratioset(j) = ratio;
            if (10^abs(log10(ratio)) - 1) > tolerance
                TASBESession.warn('TASBE:CompareBeads','PeakValueMismatch','Peak %d does not match on channel %s: ratio %f',j,getPrintName(CM.Channels{i}),ratio);
                ok=false;
            end
        end
        ratios(i) = mean(ratioset);
    end
end

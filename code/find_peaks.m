% FIND_PEAKS returns peak statistics (maxima, means, and counts) obtained
% from bead histogram.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [peak_maxima, peak_means, peak_counts] = find_peaks(bin_counts,bin_centers,peak_min_threshold,peak_ratio_threshold)

if nargin<4, peak_ratio_threshold = 2; end;
if nargin<3, peak_min_threshold = 100; end;

% initial set of maxima and minima by walking 
plateau_start = [];
maxima = []; minima = [];
for i=2:(numel(bin_counts)-1),
    if numel(plateau_start),
        if bin_counts(i) ~= bin_counts(i+1),
            plateau_middle = round(mean([plateau_start i]));
            if bin_counts(plateau_start)>bin_counts(plateau_start-1) && bin_counts(i)>bin_counts(i+1),
                maxima(end+1) = plateau_middle; 
            end;
            if bin_counts(plateau_start)<bin_counts(plateau_start-1) && bin_counts(i)<bin_counts(i+1),
                minima(end+1) = plateau_middle; 
            end;
            plateau_start = [];
        end
    else
        if (bin_counts(i)>bin_counts(i-1)) && (bin_counts(i)>bin_counts(i+1)), maxima(end+1) = i; end;
        if (bin_counts(i)<bin_counts(i-1)) && (bin_counts(i)<bin_counts(i+1)), minima(end+1) = i; end;
        if bin_counts(i)==bin_counts(i+1), plateau_start = i; end;
    end
end
% maxima = find(bin_counts(2:(end-1)) > bin_counts(1:(end-2)) & bin_counts(2:(end-1)) >= bin_counts(3:end));
% minima = find(bin_counts(2:(end-1)) <= bin_counts(1:(end-2)) & bin_counts(2:(end-1)) < bin_counts(3:end));

if(numel(maxima) == 0), 
    TASBESession.warn('TASBE:FindPeaks','NoPeaksFound','Problematic distribution: no peaks found in histogram'); 
    maxima(1) = 1; bin_counts(1) = NaN;% put in fake data
end;
if(abs(numel(maxima)-numel(minima))>1), TASBESession.error('TASBE:FindPeaks','ImpossibleDistribution','Internal error: impossible distribution of %i maxima and %i minima',numel(maxima),numel(minima)); end;

% enhance with first/last minima to get minima surrounding each maximum
if(numel(minima)==0 || minima(1)>maxima(1)), minima = [1 minima]; end;
if(minima(end)<maxima(end)), minima = [minima numel(bin_counts)]; end;

% filter out peaks with too small a count:
too_low_peaks = find(bin_counts(maxima)<peak_min_threshold);
for i=1:numel(too_low_peaks)
    % walk in reverse order to preserve indices of too_low_peaks 
    bad_peak = too_low_peaks(end-(i-1));
    % remove maxima
    maxima = [maxima(1:(bad_peak-1)) maxima((bad_peak+1):end)];
    % consolidate minima
    if bin_counts(minima(bad_peak))<bin_counts(minima(bad_peak+1)), bestmin = minima(bad_peak);
    else bestmin = minima(bad_peak+1); end;
    minima = [minima(1:(bad_peak-1)) bestmin minima((bad_peak+2):end)];
end

for i=1:numel(maxima),
    upper_ratios = bin_counts(maxima)./bin_counts(minima(1:(end-1)));
    lower_ratios = bin_counts(maxima)./bin_counts(minima(2:end));
    worst_ratio = min([upper_ratios lower_ratios]);
    % If the worst ratio is still good, terminate the search
    if(worst_ratio >= peak_ratio_threshold) break; end;
    % Otherwise, find and remove worst peak ...
    bad_peak = find(upper_ratios == worst_ratio | lower_ratios == worst_ratio,1);
    maxima = [maxima(1:(bad_peak-1)) maxima((bad_peak+1):end)];
    % ... and consolidate minima
    if bin_counts(minima(bad_peak))<bin_counts(minima(bad_peak+1)), bestmin = minima(bad_peak);
    else bestmin = minima(bad_peak+1); end;
    minima = [minima(1:(bad_peak-1)) bestmin minima((bad_peak+2):end)];
end

% finally, compute statistics for each peak region
peak_maxima = bin_centers(maxima);
peak_means = zeros(size(maxima)); peak_counts = peak_means;
for i=1:numel(maxima)
    % identify peak region
    lo = minima(i); hi = minima(i+1);
    count_set = bin_counts(lo:hi);
    % split edge bins to avoid double-counting
    if(lo ~= 1), count_set(1) = count_set(1)/2; end;
    if(hi ~= numel(bin_counts)), count_set(end) = count_set(end)/2; end;
    % get stats
    peak_means(i) = geomean(bin_centers(lo:hi),count_set);
    peak_counts(i) = sum(count_set);
end

% testing:
% centers = 10.^(3:0.1:8.2);
% counts = [1472          32          15          20          29 ...
%           46          50          80          75         108         137 ...
%          136         143         220         212         201         109 ...
%           72          28          12          12          11          11 ...
%           14          10          16          15          18          24 ...
%           38          48          63          78         105         147 ...
%          210         273         474         621        1027        1535 ...
%         1861        2168        1622         909         443         164 ...
%           53           8           0           0           0           0];
% [maxima pmeans pcounts] = find_peaks(counts,centers,100,2);
% figure; loglog(centers,counts,'b-'); hold on;
% plot(pmeans,pcounts,'r+');
% for i=1:numel(maxima), plot([maxima(i) maxima(i)],[1 max(counts)],'g-'); end;

function UT = beads_to_mefl_model(CM, settings, beadfile, makePlots, path)
% BEADS_TO_MEFL_MODEL: Computes a linear function for transforming FACS 
% measurements on the FITC channel into MEFLs, using a calibration run of
% RCP-30-5A.
% 
% Takes the name of the FACS file of bead measurements, plus optionally the
% name of the channel to be used (if not FITC-A) and a flag for whether to
% record the calibration plot.
%
% Returns:
% * k_MEFL:  MEFL = k_MEFL * FITC
% * first_peak: what is the first peak visible?
% * fit_error: residual from the linear fit

% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

FITC_channel = CM.FITC_channel;
if (nargin < 4)
    makePlots = CM.bead_plot;
end
if (nargin < 5)
    path = getSetting(settings, 'path', './');
end
force_peak = getSetting(settings,'force_first_bead_peak',[]);
if ~isempty(force_peak)
    warning('TASBE:Beads','Forcing interpretation of first detected peak as peak number %i',force_peak);
end


peak_threshold = CM.bead_peak_threshold;
bin_min = CM.bead_min;
bin_max = CM.bead_max;

nameFC=getName(FITC_channel);
i_FITC = find(CM,FITC_channel);

% SpheroTech RCP-30-5A beads (8 peaks) - only 7 are used here, since the
% first is not given a MEFL value in the tech notes
% fprintf('Matching to FITC values for %s beads\n', CM.bead_model);
% fprintf('Assuming Lot AA01, AA02, AA03, AA04, AB01, AB02, AC01, or GAA01-R\n');
% PeakMEFLs = [692 2192 6028 17493 35674 126907 290983];
%PeakRelative = [77.13 108.17 135.42 164.11 183.31 217.49 239.84];
%warning('Substituting a different RCP set');
%PeakMEFLs = [791	2083	6562	16531	47575	136680	271771];

[PeakMEFLs,units,actualBatch] = get_bead_peaks(CM.bead_model,CM.bead_channel,CM.bead_batch);

totalNumPeaks = numel(PeakMEFLs);
numQuantifiedPeaks = sum(~isnan(PeakMEFLs));
quantifiedPeakMEFLs = PeakMEFLs((end-numQuantifiedPeaks+1):end);

% identify peaks
bin_increment = 0.02;
bin_edges = 10.^(bin_min:bin_increment:bin_max);
n = (size(bin_edges,2)-1);
bin_centers = bin_edges(1:n)*10.^(bin_increment/2);

% option of segmenting FITC on a separate secondary channel
segment_secondary = hasSetting(settings,'SecondaryBeadChannel');
if segment_secondary
    segmentName = getSetting(settings,'SecondaryBeadChannel');
else
    segmentName = nameFC;
end

[fcsraw fcshdr fcsdat] = fca_readfcs(beadfile);
bead_data = get_fcs_color(fcsdat,fcshdr,nameFC);
segment_data = get_fcs_color(fcsdat,fcshdr,segmentName);

% The full range of bins (for plotting purposes) covers everything from 1 to the max value (rounded up)
range_max = max(bin_max,ceil(log10(max(bead_data(:)))));
range_bin_edges = 10.^(0:bin_increment:range_max);
range_n = (size(range_bin_edges,2)-1);
range_bin_centers = range_bin_edges(1:range_n)*10.^(bin_increment/2);

% All of this should be converted to use bin sequence classes
bin_counts = zeros(size(bin_centers));
for i=1:n
    which = segment_data(:)>bin_edges(i) & segment_data(:)<=bin_edges(i+1);
    bin_counts(i) = sum(which);
end
range_bin_counts = zeros(size(bin_centers));
for i=1:range_n
    which = segment_data(:)>range_bin_edges(i) & segment_data(:)<=range_bin_edges(i+1);
    range_bin_counts(i) = sum(which);
end

n_peaks = 0;
segment_peak_means = []; % segmentation channel (normally FITC)
peak_means = [];
peak_counts = [];
if (isempty(peak_threshold)),
    peak_threshold = 0.2 * max(bin_counts); % number of points to start a peak
end
if numel(peak_threshold)==1
    peak_threshold = peak_threshold*ones(numel(CM.Channels),1);
else if numel(peak_threshold)~=numel(CM.Channels)
        error('Bead calibration requires 0,1, or n_channels thresholds');
    end
end

in_peak = 0;
for i=1:n
    if in_peak==0 % outside a peak: look for start
        if(bin_counts(i) >= peak_threshold(i_FITC))
            peak_min = bin_edges(i);
            in_peak=1;
        end
    else % inside a peak: look for end
        if(bin_counts(i) < peak_threshold(i_FITC))
            peak_max = bin_edges(i);
            in_peak=0;
            % compute peak statistics
            n_peaks = n_peaks+1;
            which = segment_data(:)>peak_min & segment_data(:)<=peak_max;
            segment_peak_means(n_peaks) = mean(segment_data(which)); % arithmetic q. beads are primarily measurement noise
            peak_means(n_peaks) = mean(bead_data(which));
            peak_counts(n_peaks) = sum(which);
        end
    end
end

% Gather the peaks from all the channels
peak_sets = cell(numel(CM.Channels),1);
for i=1:numel(CM.Channels),
    alt_bead_data = get_fcs_color(fcsdat,fcshdr,getName(CM.Channels{i}));
    alt_bin_counts = zeros(size(bin_centers));
    for j=1:n
        which = alt_bead_data(:)>bin_edges(j) & alt_bead_data(:)<=bin_edges(j+1);
        alt_bin_counts(j) = sum(which);
    end
    alt_range_bin_counts = zeros(size(range_bin_centers));
    for j=1:range_n
        which = alt_bead_data(:)>range_bin_edges(j) & alt_bead_data(:)<=range_bin_edges(j+1);
        alt_range_bin_counts(j) = sum(which);
    end

    % identify peaks
    alt_n_peaks = 0;
    alt_peak_means = [];
    alt_peak_counts = [];
    in_peak = 0;
    for j=1:n
        if in_peak==0 % outside a peak: look for start
            if(alt_bin_counts(j) >= peak_threshold(i))
                alt_peak_min = bin_edges(j);
                in_peak=1;
            end
        else % inside a peak: look for end
            if(alt_bin_counts(j) < peak_threshold(i))
                alt_peak_max = bin_edges(j);
                in_peak=0;
                % compute peak statistics
                alt_n_peaks = alt_n_peaks+1;
                which = alt_bead_data(:)>alt_peak_min & alt_bead_data(:)<=alt_peak_max;
                alt_peak_means(alt_n_peaks) = mean(alt_bead_data(which)); % arithmetic q. beads only have measurement noise
                alt_peak_counts(alt_n_peaks) = sum(which);
            end
        end
    end
    peak_sets{i} = alt_peak_means;

    % Make plots for all peaks, not just FITC
    if makePlots >= 2
        graph_max = max(alt_range_bin_counts);
        h = figure('PaperPosition',[1 1 5 3.66]);
        set(h,'visible','off');
        semilogx(range_bin_centers,alt_range_bin_counts,'b-'); hold on;
        for j=1:alt_n_peaks
            semilogx([alt_peak_means(j) alt_peak_means(j)],[0 graph_max],'r-');
        end
        % show range where peaks were searched for
        plot(10.^[bin_min bin_min],[0 graph_max],'k:');
        text(10.^(bin_min),graph_max/2,'peak search min value','Rotation',90,'FontSize',7,'VerticalAlignment','top','FontAngle','italic');
        plot(10.^[bin_max bin_max],[0 graph_max],'k:');
        text(10.^(bin_max),graph_max/2,'peak search max value','Rotation',90,'FontSize',7,'VerticalAlignment','bottom','FontAngle','italic');
        xlabel(sprintf('a.u. for %s channel',getPrintName(CM.Channels{i}))); ylabel('Beads');
        title(sprintf('Peak identification for %s for %s beads',getPrintName(CM.Channels{i}), CM.bead_model));
        outputfig(h, sprintf('bead-calibration-%s',getPrintName(CM.Channels{i})),path);
    end
end

% look for the best linear fit of log10(peak_means) vs. log10(PeakMEFLs)
if(n_peaks>numQuantifiedPeaks)
    warning('Bead calibration found unexpectedly many bead peaks: truncating to use top peaks only');
    n_peaks = numQuantifiedPeaks;
    segment_peak_means = segment_peak_means((end-numQuantifiedPeaks+1):end);
    peak_means = peak_means((end-numQuantifiedPeaks+1):end);
    peak_counts = peak_counts((end-numQuantifiedPeaks+1):end);
end

% Use log scale for fitting to avoid distortions from highest point
if(n_peaks>=2)
    fit_error = Inf;
    first_peak = 0;
    if(n_peaks>2)
        best_i = -1;
        for i=0:(numQuantifiedPeaks-n_peaks),
          [poly,S] = polyfit(log10(peak_means),log10(quantifiedPeakMEFLs((1:n_peaks)+i)),1);
          if S.normr <= fit_error, fit_error = S.normr; model = poly; first_peak=i+2; best_i = i; end;
        end
        % Warn if setting to anything less than the top peak, since top peak should usually be visible
        fprintf('Bead peaks identified as %i to %i of %i\n',first_peak,first_peak+n_peaks-1,totalNumPeaks);
        if best_i < (numQuantifiedPeaks-n_peaks) && n_peaks < 5,
            warning('TASBE:Beads','Few bead peaks and fit does not include highest: error likely');
        end
    else % 2 peaks
        warning('TASBE:Beads','Only two bead peaks found, assuming brightest two');
        [poly,S] = polyfit(log10(peak_means),log10(quantifiedPeakMEFLs(end-1:end)),1);
        fit_error = S.normr; model = poly; first_peak = numQuantifiedPeaks;
    end
    if ~isempty(force_peak), first_peak = force_peak; end
    constrained_fit = mean(log10(quantifiedPeakMEFLs((1:n_peaks)+first_peak-2)) - log10(peak_means));
    cf_error = mean(10.^abs(log10((quantifiedPeakMEFLs((1:n_peaks)+first_peak-2)./peak_means) / 10.^constrained_fit)));
    % Final fit_error should be close to zero / 1-fold
    if(cf_error>1.05), warning('TASBE:Beads','Bead calibration may be incorrect: fit more than 5 percent off: error = %.2d',cf_error); end;
    %if(abs(model(1)-1)>0.05), warning('TASBE:Beads','Bead calibration probably incorrect: fit more than 5 percent off: slope = %.2d',model(1)); end;
    k_MEFL = 10^constrained_fit;
elseif(n_peaks==1) % 1 peak
    warning('TASBE:Beads','Only one bead peak found, assuming brightest');
    fit_error = 0; first_peak = totalNumPeaks;
    if ~isempty(force_peak), first_peak = force_peak; end
    k_MEFL = PeakMEFLs(first_peak)/peak_means;
else % n_peaks = 0
    warning('Bead calibration failed: found no bead peaks; using single dummy peak');
    k_MEFL = 1;
    fit_error = Inf;
    first_peak = NaN;
end;

% Plot fitted channel
if makePlots
    graph_max = max(range_bin_counts);
    h = figure('PaperPosition',[1 1 5 3.66]);
    set(h,'visible','off');
    semilogx(range_bin_centers,range_bin_counts,'b-'); hold on;
    % Show identified peaks
    for i=1:n_peaks
        semilogx([segment_peak_means(i) segment_peak_means(i)],[0 graph_max],'r-');
        text(peak_means(i),graph_max,sprintf('%i',i+first_peak-1),'VerticalAlignment','top');
    end
    % show range where peaks were searched for
    plot(10.^[bin_min bin_min],[0 graph_max],'k:');
    text(10.^(bin_min),graph_max/2,'peak search min value','Rotation',90,'FontSize',7,'VerticalAlignment','top','FontAngle','italic');
    plot(10.^[bin_max bin_max],[0 graph_max],'k:');
    text(10.^(bin_max),graph_max/2,'peak search max value','Rotation',90,'FontSize',7,'VerticalAlignment','bottom','FontAngle','italic');
    plot(10.^[0 range_max],[peak_threshold(i_FITC) peak_threshold(i_FITC)],'k:');
    text(1,peak_threshold(i_FITC),'clutter threshold','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','bottom','FontAngle','italic');
    title(sprintf('Peak identification for %s beads', CM.bead_model));
    xlim(10.^[0 range_max]);
    ylabel('Beads');
    if segment_secondary
        xlabel([segmentName ' a.u.']); 
        outputfig(h,'bead-calibration-secondary',path);
    else
        xlabel([CM.bead_channel ' a.u.']); 
        outputfig(h,'bead-calibration',path);
    end
end


% Plot bead fit curve
if makePlots>1
    h = figure('PaperPosition',[1 1 5 3.66]);
    set(h,'visible','off');
    loglog(peak_means,quantifiedPeakMEFLs((1:n_peaks)+first_peak-2),'b*-'); hold on;
    %loglog([1 peak_means],[1 peak_means]*(10.^model(2)),'r+--');
    loglog([1 peak_means],[1 peak_means]*k_MEFL,'go--');
    for i=1:n_peaks
        text(peak_means(i),quantifiedPeakMEFLs(i+first_peak-2)*1.3,sprintf('%i',i+first_peak-1));
    end
    xlabel([CM.bead_channel ' a.u.']); ylabel('Beads MEFLs');
    title(sprintf('Peak identification for %s beads', CM.bead_model));
    %legend('Location','NorthWest','Observed','Linear Fit','Constrained Fit');
    legend('Observed','Constrained Fit','Location','NorthWest');
    outputfig(h,'bead-fit-curve',path);
end

% Plog 2D fit
if makePlots
    % plot FITC linearly, since we wouldn't be using a secondary if the values weren't very low
    % there is probably much negative data
    if segment_secondary
        h = figure('PaperPosition',[1 1 5 3.66]);
        set(h,'visible','off');
        pos = segment_data>0;
        smin = log10(percentile(segment_data(pos),0.1)); smax = log10(percentile(segment_data(pos),99.9));
        bmin = percentile(bead_data(pos),0.1); bmax = percentile(bead_data(pos),99.9);
        % range to 99th percentile
        smoothhist2D([bead_data(pos) log10(segment_data(pos))],10,[200, 200],[],'image',[bmin smin; bmax smax]); hold on;
        for i=1:n_peaks
            semilogy([min(bead_data) max(bead_data)],log10([segment_peak_means(i) segment_peak_means(i)]),'r-');
            semilogy(peak_means(i),log10(segment_peak_means(i)),'k+');
            text(peak_means(i),log10(segment_peak_means(i))+0.1,sprintf('%i',i+first_peak-1));
        end
        xlabel([CM.bead_channel ' a.u.']); ylabel([segmentName ' a.u.']);
        title(sprintf('Peak identification for %s beads', CM.bead_model));
        outputfig(h,'bead-calibration',path);
    end
end

UT = UnitTranslation([CM.bead_model ':' CM.bead_channel ':' CM.bead_batch],k_MEFL, first_peak, fit_error, peak_sets);

end


% BEADS_TO_ERF_MODEL: Computes a linear function for transforming FACS 
% measurements on the ERF channel into ERFs, using a calibration run of
% the bead model.
% 
% Takes the name of the FACS file of bead measurements, plus optionally the
% name of the channel to be used (if not FITC-A) and a flag for whether to
% record the calibration plot.
%
% Returns:
% * k_ERF:  ERF = k_ERF * ERF_channel_AU
% * first_peak: what is the first peak visible?
% * fit_error: residual from the linear fit
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [UT CM] = beads_to_ERF_model(CM, beadfile)
ERF_channel = CM.ERF_channel;
if(isUnprocessed(ERF_channel)) 
    TASBESession.error('TASBE:Beads','ERFUnprocessed','ERF channel %s cannot be unprocessed',getName(ERF_Channel));
end

makePlots = TASBEConfig.get('beads.plot');
visiblePlots = TASBEConfig.get('beads.visiblePlots');
plotPath = TASBEConfig.get('beads.plotPath');
plotSize = TASBEConfig.get('beads.plotSize');
beadModel = TASBEConfig.get('beads.beadModel');
beadChannel = TASBEConfig.get('beads.beadChannel');
beadBatch = TASBEConfig.getexact('beads.beadBatch',[]);

force_peak = TASBEConfig.getexact('beads.forceFirstPeak',[]);
if ~isempty(force_peak)
    TASBESession.warn('TASBE:Beads','ForcedPeak','Forcing interpretation of first detected peak as peak number %i',force_peak);
end


peak_threshold = TASBEConfig.getexact('beads.peakThreshold',[]);
bin_min = TASBEConfig.get('beads.rangeMin');
bin_max = TASBEConfig.get('beads.rangeMax');
bin_increment = TASBEConfig.get('beads.binIncrement');


erfChannelName=getName(ERF_channel);
i_ERF = find(CM,ERF_channel);

[PeakERFs,units,actualBatch] = get_bead_peaks(beadModel,beadChannel,beadBatch);
CM.standardUnits = units;
if ~strcmp(units,'MEFL')
    TASBESession.warn('TASBE:Beads','NonMEFL','MEFL units are recommended, rather than %s',units);
end
ERF_channel_idx = indexof(CM.Channels,CM.ERF_channel);
CM.Channels{ERF_channel_idx} = setUnits(CM.Channels{ERF_channel_idx},units);


% NOTE: Calculations are done against the QuantifiedPeaks not PeakERFs.
% The value of first_peak is the first valid peak in QuantifiedPeaks not
% PeakERFs.  At the end of function, the peakOffset is added to first_peak
% and is used in calculation of UnitTranslation.

% NOTE: Thus, if reported messages or labels on plots are supposed to be
% based on PeakERFS, then add peakOffset to first_peak and
% numQuantifiedPeaks in messages and labels. Do not add it to num_peaks.
totalNumPeaks = numel(PeakERFs);
numQuantifiedPeaks = sum(~isnan(PeakERFs));
quantifiedPeakERFs = PeakERFs((end-numQuantifiedPeaks+1):end);
peakOffset = totalNumPeaks - numQuantifiedPeaks;

TASBESession.succeed('TASBE:Beads','ObtainBeadPeaks','Found specified bead model and lot');


% identify peaks
bin_edges = 10.^(bin_min:bin_increment:bin_max);
n = (size(bin_edges,2)-1);
bin_centers = bin_edges(1:n)*10.^(bin_increment/2);

% option of segmenting ERF on a separate secondary channel
segment_secondary = TASBEConfig.isSet('beads.secondaryBeadChannel');
if segment_secondary
    segmentName = TASBEConfig.get('beads.secondaryBeadChannel');
else
    segmentName = erfChannelName;
end

[fcsraw fcshdr fcsdat] = fca_read(beadfile);
if(isempty(fcshdr)), TASBESession.error('TASBE:Beads','MissingBeadFile','Cannot calibrate without beads'); end;
bead_data = get_fcs_color(fcsdat,fcshdr,erfChannelName);
segment_data = get_fcs_color(fcsdat,fcshdr,segmentName);

TASBESession.succeed('TASBE:Beads','ObtainBeadData','Successfully read bead data');

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
segment_peak_means = []; % segmentation channel (normally ERF channel)
peak_means = [];
peak_counts = [];
if (isempty(peak_threshold)),
    peak_threshold = 0.2 * max(bin_counts); % number of points to start a peak
end
if numel(peak_threshold)==1
    peak_threshold = peak_threshold*ones(numel(CM.Channels),1);
else if numel(peak_threshold)~=numel(CM.Channels)
        TASBESession.error('TASBE:Beads','ThresholdCountMismatch','Bead calibration requires 0,1, or n_channels thresholds');
    end
end

in_peak = 0;
for i=1:n
    if in_peak==0 % outside a peak: look for start
        if(bin_counts(i) >= peak_threshold(i_ERF))
            peak_min = bin_edges(i);
            in_peak=1;
        end
    else % inside a peak: look for end
        if(bin_counts(i) < peak_threshold(i_ERF))
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
for i=1:numel(CM.Channels)
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
    alt_peak_means = []; % x-value of peaks
    alt_peak_counts = [];
    alt_peak_starts = []; % starting points of peaks
    alt_peak_ends = []; % ending points of peaks
    alt_peak_maximas = []; % y-value of peaks
    in_peak = 0;
    peak_max = 0; % used to find y-value of peak
    for j=1:n
        if in_peak==0 % outside a peak: look for start
            if(alt_bin_counts(j) >= peak_threshold(i)) % found a start of a peak
                alt_peak_min = bin_edges(j);
                in_peak=1;
                peak_max = alt_bin_counts(j);
            end
        else % inside a peak: look for end
            if peak_max < alt_bin_counts(j)
                peak_max = alt_bin_counts(j);
            end
            if(alt_bin_counts(j) < peak_threshold(i)) % found the end of the peak
                alt_peak_max = bin_edges(j);
                in_peak=0;
                % compute peak statistics
                alt_n_peaks = alt_n_peaks+1;
                which = alt_bead_data(:)>alt_peak_min & alt_bead_data(:)<=alt_peak_max;
                alt_peak_means(alt_n_peaks) = mean(alt_bead_data(which)); % arithmetic q. beads only have measurement noise
                alt_peak_counts(alt_n_peaks) = sum(which);
                alt_peak_starts(alt_n_peaks) = alt_peak_min;
                alt_peak_ends(alt_n_peaks) = alt_peak_max;
                alt_peak_maximas(alt_n_peaks) = peak_max;
                peak_max = 0;
            end
        end
    end
    peak_sets{i} = alt_peak_means;
    % replace the ERF channel peak-set if we're doing a secondary segmentation
    if segment_secondary && strcmp(getName(CM.Channels{i}),erfChannelName);
        peak_sets{i} = peak_means;
    end
    
    % Running some initial tests on the peak statistics to generate some
    % warnings, either run only on channel used for ERF conversion or on all channels (determined by validateAllChannels):
    if alt_n_peaks > 0 && ((TASBEConfig.get('beads.validateAllChannels') && i == i_ERF) || ~TASBEConfig.get('beads.validateAllChannels'))
        % Check to see if peaks in ascending order and whether an extra peak of
        % combined beads is identified
        if ~issorted(alt_peak_maximas)
            TASBESession.warn('TASBE:Beads','PeaksNotAscending','Peaks are not in ascending order. May need to adjust rangeMin or rangeMax for %s.',clean_for_latex(getPrintName(CM.Channels{i})));
            if alt_peak_maximas(alt_n_peaks) < max(alt_peak_maximas) %if the last peak is lower than the maximum peak
                TASBESession.warn('TASBE:Beads','PotentialBeadClump','Last peak may consist of beads stuck together. May need to adjust rangeMax or peakThresholdfor %s.',clean_for_latex(getPrintName(CM.Channels{i})));
            end
        end
        % Check to make sure that the true peak was properly identified
        if max(alt_range_bin_counts) > max(alt_peak_maximas)
            TASBESession.warn('TASBE:Beads','HighestPeakUndetected','Did not detect highest peak for %s.',clean_for_latex(getPrintName(CM.Channels{i})));
        end
        % Check to see if a deceptive peak very close to rangeMin was
        % identified
        if abs(bin_min - log10(alt_peak_means(1))) < 0.05
            TASBESession.warn('TASBE:Beads','QuestionableFirstPeak','First peak very close to rangeMin. May need to increase rangeMin or peakThreshold for %s.',clean_for_latex(getPrintName(CM.Channels{i})));
        end
    end

    % Make plots for all peaks, not just ERF
    if makePlots
        graph_max = max(alt_range_bin_counts);
        h = figure('PaperPosition',[1 1 plotSize]);
        if(~visiblePlots), set(h,'visible','off'); end;
        semilogx(range_bin_centers,alt_range_bin_counts,'b-'); hold on;
        for j=1:alt_n_peaks
            semilogx([alt_peak_means(j) alt_peak_means(j)],[0 graph_max],'r-');
            % input the two dashed start and end lines for each peak
            semilogx([alt_peak_starts(j) alt_peak_starts(j)],[0 graph_max],'r:');
            semilogx([alt_peak_ends(j) alt_peak_ends(j)],[0 graph_max],'r:');
        end
        % show range where peaks were searched for
        plot(10.^[bin_min bin_min],[0 graph_max],'k:');
        text(10.^(bin_min),graph_max/2,'peak search min value','Rotation',90,'FontSize',7,'VerticalAlignment','top','FontAngle','italic');
        plot(10.^[bin_max bin_max],[0 graph_max],'k:');
        text(10.^(bin_max),graph_max/2,'peak search max value','Rotation',90,'FontSize',7,'VerticalAlignment','bottom','FontAngle','italic');
        xlabel(sprintf('a.u. for %s channel',clean_for_latex(getPrintName(CM.Channels{i})))); ylabel('Beads');
        title(sprintf('Peak identification for %s for %s beads',clean_for_latex(getPrintName(CM.Channels{i})), beadModel));
        outputfig(h, sprintf('bead-peak-identification-%s',clean_for_latex(getPrintName(CM.Channels{i}))),plotPath);
    end
end

% look for the best linear fit of log10(peak_means) vs. log10(PeakERFs)
if(n_peaks>numQuantifiedPeaks)
    TASBESession.warn('TASBE:Beads','PeakDetection','Bead calibration found unexpectedly many bead peaks: truncating to use top peaks only');
    n_peaks = numQuantifiedPeaks;
    segment_peak_means = segment_peak_means((end-numQuantifiedPeaks+1):end);
    peak_means = peak_means((end-numQuantifiedPeaks+1):end);
    peak_counts = peak_counts((end-numQuantifiedPeaks+1):end);
else
    TASBESession.succeed('TASBE:Beads','PeakDetection','Bead calibration found %i bead peaks',n_peaks);
end

% Use log scale for fitting to avoid distortions from highest point
if(n_peaks>=2)
    fit_error = Inf;
    first_peak = 0;
    if(n_peaks>2)
        best_i = -1;
        for i=0:(numQuantifiedPeaks-n_peaks),
          [poly,S] = polyfit(log10(peak_means),log10(quantifiedPeakERFs((1:n_peaks)+i)),1);
          if S.normr <= fit_error, fit_error = S.normr; model = poly; first_peak=i+1; best_i = i; end;
        end
        % Warn if setting to anything less than the top peak, since top peak should usually be visible
        fprintf('Bead peaks identified as %i to %i of %i\n',first_peak+peakOffset,first_peak+n_peaks-1+peakOffset,numQuantifiedPeaks+peakOffset);
        if best_i < (numQuantifiedPeaks-n_peaks) && n_peaks < 5,
            TASBESession.warn('TASBE:Beads','PeakIdentification','Few bead peaks and fit does not include highest: error likely');
        else
            TASBESession.succeed('TASBE:Beads','PeakIdentification','Matched multiple peaks in reasonable range');
        end
    else % 2 peaks
        TASBESession.warn('TASBE:Beads','PeakIdentification','Only two bead peaks found, assuming brightest two');
        [poly,S] = polyfit(log10(peak_means),log10(quantifiedPeakERFs(end-1:end)),1);
        fit_error = S.normr; model = poly; first_peak = numQuantifiedPeaks-1;
    end
    if ~isempty(force_peak), first_peak = force_peak-peakOffset; end
    constrained_fit = mean(log10(quantifiedPeakERFs((1:n_peaks)+first_peak-1)) - log10(peak_means));
    cf_error = mean(10.^abs(log10((quantifiedPeakERFs((1:n_peaks)+first_peak-1)./peak_means) / 10.^constrained_fit)));
    % Final fit_error should be close to zero / 1-fold
    if(cf_error>1.05), 
        TASBESession.warn('TASBE:Beads','PeakFitQuality','Bead calibration may be incorrect: fit more than 5 percent off: error = %.2d',cf_error); 
    else
        TASBESession.succeed('TASBE:Beads','PeakFitQuality','Bead fit quality acceptable: error = %.2d',cf_error);
    end;
    %if(abs(model(1)-1)>0.05), TASBESession.warn('TASBE:Beads','Bead calibration probably incorrect: fit more than 5 percent off: slope = %.2d',model(1)); end;
    k_ERF = 10^constrained_fit;
elseif(n_peaks==1) % 1 peak
    TASBESession.warn('TASBE:Beads','PeakIdentification','Only one bead peak found, assuming brightest');
    TASBESession.skip('TASBE:Beads','PeakFitQuality','Fit quality irrelevant for single peak');
    fit_error = 0; first_peak = numQuantifiedPeaks;
    if ~isempty(force_peak), first_peak = force_peak-peakOffset; end
    k_ERF = quantifiedPeakERFs(first_peak)/peak_means;
else % n_peaks = 0
    TASBESession.warn('TASBE:Beads','PeakIdentification','Bead calibration failed: found no bead peaks; using single dummy peak');
    TASBESession.skip('TASBE:Beads','PeakFitQuality','Fit quality irrelevant for single peak');
    k_ERF = 1;
    fit_error = Inf;
    first_peak = NaN;
    CM.standardUnits = 'arbitrary units';
end;

% Plot fitted channel
if makePlots
    graph_max = max(range_bin_counts);
    h = figure('PaperPosition',[1 1 plotSize]);
    if(~visiblePlots), set(h,'visible','off'); end;
    semilogx(range_bin_centers,range_bin_counts,'b-'); hold on;
    % Show identified peaks
    for i=1:n_peaks
        semilogx([segment_peak_means(i) segment_peak_means(i)],[0 graph_max],'r-');
        text(segment_peak_means(i),graph_max,sprintf('%i',i+first_peak-1+peakOffset),'VerticalAlignment','top');
    end
    % show range where peaks were searched for
    plot(10.^[bin_min bin_min],[0 graph_max],'k:');
    text(10.^(bin_min),graph_max/2,'peak search min value','Rotation',90,'FontSize',7,'VerticalAlignment','top','FontAngle','italic');
    plot(10.^[bin_max bin_max],[0 graph_max],'k:');
    text(10.^(bin_max),graph_max/2,'peak search max value','Rotation',90,'FontSize',7,'VerticalAlignment','bottom','FontAngle','italic');
    plot(10.^[0 range_max],[peak_threshold(i_ERF) peak_threshold(i_ERF)],'k:');
    text(1,peak_threshold(i_ERF),'clutter threshold','FontSize',7,'HorizontalAlignment','left','VerticalAlignment','bottom','FontAngle','italic');
    title(sprintf('Peak identification for %s beads', clean_for_latex(beadModel)));
    xlim(10.^[0 range_max]);
    ylabel('Beads');
    if segment_secondary
        xlabel([clean_for_latex(segmentName) ' a.u.']); 
        outputfig(h,'bead-peak-identification-secondary',plotPath);
    else
        xlabel([clean_for_latex(beadChannel) ' a.u.']); 
        outputfig(h,'bead-peak-identification',plotPath);
    end
end


% Plot bead fit curve
if makePlots
    h = figure('PaperPosition',[1 1 plotSize]);
    if(~visiblePlots), set(h,'visible','off'); end;
    loglog(peak_means,quantifiedPeakERFs((1:n_peaks)+first_peak-1),'b*-'); hold on;
    %loglog([1 peak_means],[1 peak_means]*(10.^model(2)),'r+--');
    loglog([1 peak_means],[1 peak_means]*k_ERF,'go--');
    for i=1:n_peaks
        text(peak_means(i),quantifiedPeakERFs(i+first_peak-1)*1.3,sprintf('%i',i+first_peak-1+peakOffset));
    end
    xlabel(clean_for_latex([beadChannel ' a.u.'])); ylabel('Beads ERFs');
    title(sprintf('Peak identification for %s beads', beadModel));
    %legend('Location','NorthWest','Observed','Linear Fit','Constrained Fit');
    legend('Observed','Constrained Fit','Location','NorthWest');
    outputfig(h,'bead-fit-curve',plotPath);
end

% Plog 2D fit
if makePlots
    % plot ERF linearly, since we wouldn't be using a secondary if the values weren't very low
    % there is probably much negative data
    if segment_secondary
        h = figure('PaperPosition',[1 1 plotSize]);
        if(~visiblePlots), set(h,'visible','off'); end;
        pos = segment_data>0;
        smin = log10(percentile(segment_data(pos),0.1)); smax = log10(percentile(segment_data(pos),99.9));
        bmin = percentile(bead_data(pos),0.1); bmax = percentile(bead_data(pos),99.9);
        % range to 99th percentile
        smoothhist2D([bead_data(pos) log10(segment_data(pos))],10,[200, 200],[],'image',[bmin smin; bmax smax]); hold on;
        for i=1:n_peaks
            semilogy([min(bead_data) max(bead_data)],log10([segment_peak_means(i) segment_peak_means(i)]),'r-');
            semilogy(peak_means(i),log10(segment_peak_means(i)),'k+');
            text(peak_means(i),log10(segment_peak_means(i))+0.1,sprintf('%i',i+first_peak-1+peakOffset));
        end
        xlabel(clean_for_latex([beadChannel ' a.u.'])); ylabel(clean_for_latex([segmentName ' a.u.']));
        title(sprintf('Peak identification for %s beads', beadModel));
        outputfig(h,'bead-peak-identification',plotPath);
    end
end

UT = UnitTranslation([beadModel ':' beadChannel ':' actualBatch],k_ERF, first_peak+peakOffset, fit_error, peak_sets);

end


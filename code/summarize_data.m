% SUMMARIZE_DATA Given the data from each sample, summarize it into a single
% collection of statistics per condition
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function experimentresults = summarize_data( colorModel, experiment, analysisParams, sampleresults)
filenames = getInducerLevelsToFiles(experiment); % array of file names
n_conditions = numel(filenames);

n_channels = numel(getChannels(colorModel));
n_bins = get_n_bins(getBins(analysisParams));
n_components = getNumGaussianComponents(analysisParams);

% Combine all material together across replicates
bincounts = zeros(n_bins,n_conditions);
repcounts = zeros(n_bins,n_conditions);
tmp_expt_means = zeros(n_bins,n_conditions,n_channels);
tmp_expt_stds = zeros(n_bins,n_conditions,n_channels);
tmp_expt_popmeans = zeros(n_conditions,n_channels);
tmp_expt_popstds = zeros(n_conditions,n_channels);
tmp_expt_histograms = zeros(n_conditions,n_bins,n_channels);
tmp_expt_excluded = zeros(n_conditions, n_channels);
expt_nonexpressing = zeros(n_conditions);
tmp_expt_popcmeans = zeros(n_components,n_conditions,n_channels);
tmp_expt_popcstds = zeros(n_components,n_conditions,n_channels);
tmp_expt_popcweights = zeros(n_components,n_conditions,n_channels);
expt_plasmids = zeros(n_bins,n_conditions);
tmp_expt_stdofmeans = zeros(n_bins,n_conditions,n_channels);
tmp_expt_stdofstds = zeros(n_bins,n_conditions,n_channels);
tmp_expt_stdofpopmeans = zeros(n_conditions,n_channels);
tmp_expt_stdofpopstds = zeros(n_conditions,n_channels);
tmp_expt_stdofhistograms = zeros(n_conditions,n_bins,n_channels);
tmp_expt_stdofexcluded = zeros(n_conditions,n_channels);
expt_stdofnonexpressing = zeros(n_conditions);
tmp_expt_stdofpopcmeans = zeros(n_components,n_conditions,n_channels);
tmp_expt_stdofpopcstds = zeros(n_components,n_conditions,n_channels);
tmp_expt_stdofpopcweights = zeros(n_components,n_conditions,n_channels);
expt_stdofplasmids = zeros(n_bins,n_conditions);

max_events = 0; max_sample = [];
min_events = Inf; min_sample = [];

for i=1:n_conditions,
    sampleset = [sampleresults{i}{:}];
    % Only use bins with at least the minimum number of data elements
    valid = [sampleset.BinCounts]>getMinValidCount(analysisParams);
    % Add all valid bins to counts
    repcounts(:,i) = sum(valid,2);
    % pull sample data
    sample_bc = [sampleset.BinCounts];
    sample_means = reshape([sampleset.Means],n_bins,n_channels,numel(sampleset));
    sample_stds = reshape([sampleset.StandardDevs],n_bins,n_channels,numel(sampleset));
    
    sample_popmeans = reshape([sampleset.PopMeans],n_channels,numel(sampleset));
    sample_popstds = reshape([sampleset.PopStandardDevs],n_channels,numel(sampleset));
    sample_histograms = reshape([sampleset.Histograms],n_bins,n_channels,numel(sampleset));
    sample_excluded = reshape([sampleset.Excluded],n_channels,numel(sampleset));
    sample_nonexpressing = [sampleset.NonExpressing];
    
    sample_popcmeans = reshape([sampleset.PopComponentMeans],n_components,n_channels,numel(sampleset));
    sample_popcstds = reshape([sampleset.PopComponentStandardDevs],n_components,n_channels,numel(sampleset));
    sample_popcweights = reshape([sampleset.PopComponentWeights],n_components,n_channels,numel(sampleset));
    sample_plas = [sampleset.PlasmidEstimates];
    sample_active = [sampleset.FractionActive];
    
    % track min/max total number of events in a sample
    for j=1:numel(sampleset)
        event_count = sampleresults{i}{j}.Excluded(1) + sum(sampleresults{i}{j}.Histograms(:,1));
        if event_count>max_events, max_events = event_count; max_sample = [i j]; end;
        if event_count<min_events, min_events = event_count; min_sample = [i j]; end;
    end

    % first the bin statistics
    for j=1:n_bins
        bincounts(j,i) = sum(sample_bc(j,valid(j,:)));
        % TODO: allow arithmetic statistics as well
        % Note: Matlab handled all of these statistics using the general
        % case of valud(j,:)>1; Octave has more safeties on statistics
        % over higher dimensions, and in this case they backfire
        if sum(valid(j,:))==0,
            tmp_expt_means(j,i,:) = NaN;
            tmp_expt_stdofmeans(j,i,:) = NaN;
            tmp_expt_stds(j,i,:) = NaN;
            tmp_expt_stdofstds(j,i,:) = NaN;
            expt_plasmids(j,i) = NaN;
            expt_stdofplasmids(j,i) = NaN;
            expt_activity(j,i) = NaN;
            expt_stdofactivity(j,i) = NaN;
        else if sum(valid(j,:))==1,
                tmp_expt_means(j,i,:) = sample_means(j,:,valid(j,:));
                tmp_expt_stdofmeans(j,i,:) = 0;
                tmp_expt_stds(j,i,:) = sample_stds(j,:,valid(j,:));
                tmp_expt_stdofstds(j,i,:) = 0;
                expt_plasmids(j,i) = sample_plas(j,valid(j,:));
                expt_stdofplasmids(j,i) = 0;
                expt_activity(j,i) = sample_active(j,valid(j,:));
                expt_stdofactivity(j,i) = 0;
            else
                for k=1:n_channels,
                    tmp_expt_means(j,i,k) = geomean(sample_means(j,k,valid(j,:)));
                    tmp_expt_stdofmeans(j,i,k) = geostd(sample_means(j,k,valid(j,:)));
                    tmp_expt_stds(j,i,k) = geomean(sample_stds(j,k,valid(j,:)));
                    tmp_expt_stdofstds(j,i,k) = geostd(sample_stds(j,k,valid(j,:)));
                end
                expt_plasmids(j,i) = geomean(sample_plas(j,valid(j,:)));
                expt_stdofplasmids(j,i) = geostd(sample_plas(j,valid(j,:)));
                expt_activity(j,i) = geomean(sample_active(j,valid(j,:)));
                expt_stdofactivity(j,i) = geostd(sample_active(j,valid(j,:)));
            end
        end
    end
    
    % and now the bulk statistics
    for j=1:n_channels,
        tmp_expt_popmeans(i,j) = geomean(sample_popmeans(j,:));
        tmp_expt_stdofpopmeans(i,j) = geostd(sample_popmeans(j,:));
        tmp_expt_popstds(i,j) = geomean(sample_popstds(j,:));
        tmp_expt_stdofpopstds(i,j) = geostd(sample_popstds(j,:));
        for k=1:n_bins,
            tmp_expt_histograms(i,k,j) = geomean(sample_histograms(k,j,:));
            tmp_expt_stdofhistograms(i,k,j) = geostd(sample_histograms(k,j,:));
        end
        tmp_expt_excluded(i,j) = geomean(sample_excluded(j,:));
        tmp_expt_stdofexcluded(i,j) = geostd(sample_excluded(j,:));
    end
    expt_nonexpressing = geomean(sample_nonexpressing);
    expt_stdofnonexpressing = geostd(sample_nonexpressing);
    
    % and now the k-component gaussian population statistics
    for j=1:n_components
        for k=1:n_channels,
            tmp_expt_popcmeans(j,i,k) = geomean(sample_popcmeans(j,k,:));
            tmp_expt_stdofpopcmeans(j,i,k) = geostd(sample_popcmeans(j,k,:));
            tmp_expt_popcstds(j,i,k) = geomean(sample_popcstds(j,k,:));
            tmp_expt_stdofpopcstds(j,i,k) = geostd(sample_popcstds(j,k,:));
            tmp_expt_popcweights(j,i,k) = mean(sample_popcweights(j,k,:));
            tmp_expt_stdofpopcweights(j,i,k) = std(sample_popcweights(j,k,:));
		end
    end	    
end

% Test for sanity on variation in sample/sample event counts - should be less than 10x
if max_events/min_events > TASBEConfig.get('flow.replicateEventRatioWarning')
    TASBESession.warn('TASBE:SummarizeData','HighSampleSizeVariation','High variation in events per sample:\n  Condition %i, sample %i = %i\n  Condition %i, sample %i = %i',...
        max_sample(1), max_sample(2), max_events, min_sample(1), min_sample(2), min_events);
end

% Fill in per-channel outputs
expt_means = cell(n_channels,2); expt_stds = cell(n_channels,2);
expt_stdofmeans = cell(n_channels,2); expt_stdofstds = cell(n_channels,2);
expt_popmeans = cell(n_channels,2); expt_popstds = cell(n_channels,2);
expt_histograms = cell(n_channels,2); expt_excluded = cell(n_channels,2);
expt_stdofpopmeans = cell(n_channels,2); expt_stdofpopstds = cell(n_channels,2);
expt_stdofhistograms = cell(n_channels,2); expt_stdofexcluded = cell(n_channels,2);
expt_popcmeans = cell(n_channels,2); expt_popcstds = cell(n_channels,2); expt_popcweights = cell(n_channels,2);
expt_stdofpopcmeans = cell(n_channels,2); expt_stdofpopcstds = cell(n_channels,2); expt_stdofpopcweights = cell(n_channels,2);

% First entry of each output table is the Channel, second is the statistic array
for i=1:n_channels,
    expt_means{i,1} = getChannel(colorModel,i);
    expt_stds{i,1} = getChannel(colorModel,i);
    expt_stdofmeans{i,1} = getChannel(colorModel,i);
    expt_stdofstds{i,1} = getChannel(colorModel,i);
    expt_means{i,2} = tmp_expt_means(:,:,i);
    expt_stds{i,2} = tmp_expt_stds(:,:,i);
    expt_stdofmeans{i,2} = tmp_expt_stdofmeans(:,:,i);
    expt_stdofstds{i,2} = tmp_expt_stdofstds(:,:,i);

    expt_popmeans{i,1} = getChannel(colorModel,i);
    expt_popstds{i,1} = getChannel(colorModel,i);
    expt_histograms{i,1} = getChannel(colorModel,i);
    expt_excluded{i,1} = getChannel(colorModel,i);
    expt_stdofpopmeans{i,1} = getChannel(colorModel,i);
    expt_stdofpopstds{i,1} = getChannel(colorModel,i);
    expt_stdofhistograms{i,1} = getChannel(colorModel,i);
    expt_stdofexcluded{i,1} = getChannel(colorModel,i);
    expt_popmeans{i,2} = tmp_expt_popmeans(:,i);
    expt_popstds{i,2} = tmp_expt_popstds(:,i);
    expt_histograms{i,2} = tmp_expt_histograms(:,:,i);
    expt_excluded{i,2} = tmp_expt_excluded(:,i);
    expt_stdofpopmeans{i,2} = tmp_expt_stdofpopmeans(:,i);
    expt_stdofpopstds{i,2} = tmp_expt_stdofpopstds(:,i);
    expt_stdofhistograms{i,2} = tmp_expt_stdofhistograms(:,:,i);
    expt_stdofexcluded{i,2} = tmp_expt_stdofexcluded(:,i);

    expt_popcmeans{i,1} = getChannel(colorModel,i);
    expt_popcstds{i,1} = getChannel(colorModel,i);
    expt_popcweights{i,1} = getChannel(colorModel,i);
    expt_stdofpopcmeans{i,1} = getChannel(colorModel,i);
    expt_stdofpopcstds{i,1} = getChannel(colorModel,i);
    expt_stdofpopcweights{i,1} = getChannel(colorModel,i);
    expt_popcmeans{i,2} = tmp_expt_popcmeans(:,:,i);
    expt_popcstds{i,2} = tmp_expt_popcstds(:,:,i);
    expt_popcweights{i,2} = tmp_expt_popcweights(:,:,i);
    expt_stdofpopcmeans{i,2} = tmp_expt_stdofpopcmeans(:,:,i);
    expt_stdofpopcstds{i,2} = tmp_expt_stdofpopcstds(:,:,i);
    expt_stdofpopcweights{i,2} = tmp_expt_stdofpopcweights(:,:,i);
end

experimentresults = ExperimentResults(experiment, analysisParams, repcounts, bincounts, ...
    expt_means, expt_stds, expt_plasmids, expt_activity, ...
    expt_popmeans, expt_popstds, expt_histograms, expt_excluded, expt_nonexpressing, ...
    expt_popcmeans, expt_popcstds, expt_popcweights, ...
    expt_stdofmeans, expt_stdofstds, expt_stdofplasmids, expt_stdofactivity, ...
    expt_stdofpopmeans, expt_stdofpopstds, expt_stdofhistograms, expt_stdofexcluded, expt_stdofnonexpressing, ...
    expt_stdofpopcmeans, expt_stdofpopcstds, expt_stdofpopcweights);

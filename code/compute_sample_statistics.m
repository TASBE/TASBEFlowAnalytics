% COMPUTE_SAMPLE_STATISTICS computes distribution structure of a data sample
% This can be applied to either real data (from FCS files) or simulated
% data.
% Data is an array of ERF levels (events x channels), in colormodel order of channels
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function sampleresults = compute_sample_statistics(colorModel,experiment,samplename,analysisParams,data)
n_components = getNumGaussianComponents(analysisParams);
n_channels = numel(getChannels(colorModel));

% Get the constitutive
constitutive_channel = getChannel(analysisParams, 'constitutive');
if isempty(constitutive_channel), TASBESession.error('TASBE:Analysis','NoConstitutiveChannel','A constitutive channel is required!'); end;
c_index = find(colorModel,constitutive_channel); % Determine data column from ColorModel

% Lookup distribution models
CFP_af = get_autofluorescence_model(colorModel,c_index);
CFP_noise_model = get_noise_model(colorModel);
CFP_noise = CFP_noise_model.noisemin(c_index);

bins = getBins(analysisParams);
bin_centers = get_bin_centers(bins);
bin_edges = get_bin_edges(bins);
n_bins = numel(bin_centers);

% compute population bulk statistics
histograms = zeros(n_bins,n_channels);
popmeans = zeros(n_channels,1);
popstds = zeros(n_channels,1);
poppeaks = cell(n_channels,1);
excluded = zeros(n_channels,1);
nonexpressing_set = ones(size(data,1),1);
popcmeans = zeros(n_components,n_channels);
popcstds = zeros(n_components,n_channels);
popcweights = zeros(n_components,n_channels);
drop_thresholds = zeros(n_channels,1);
for k=1:n_channels
    channel = getChannel(colorModel,k);
    
    % TODO: generalize this for other types of derived channel
    if(strcmp(getUnits(channel),'Boolean')) % linear stats for Boolean channels
        % nothing gets excluded from a derived channel
        drop_thresholds(k) = -1;
        excluded(k) = 0;
        
        % compute bulk statistics
        tmp_hist = histc(data(:,k),log10(bin_edges)); % note: histcounts is better, but missing in octave & older Matlab
        histograms(:,k) = tmp_hist(1:(end-1));

        popmeans(k) = mean(data(:,k));
        popstds(k) = std(data(:,k));
        
    else
        % compute channel-specific drop threshold:
        drop_threshold = au_to_ERF(colorModel,channel,get_pem_drop_threshold(analysisParams));
        drop_thresholds(k) = drop_threshold;

        % divide into included/excluded set
        pos = data(:,k)>drop_threshold;
        excluded(k) = sum(~pos);
        nonexpressing_set = nonexpressing_set & ~pos; % remove non-excluded set from non-expressing

        % compute bulk statistics
        tmp_hist = histc(log10(data(pos,k)),log10(bin_edges)); % note: histcounts is better, but missing in octave & older Matlab
        histograms(:,k) = tmp_hist(1:(end-1));

        popmeans(k) = geomean(data(pos,k));
        popstds(k) = geostd(data(pos,k));
    end
    poppeaks{k} = find_peaks(histograms(:,k)',bin_centers); % need to include peak detection params in analysisParams

    % For now, disable k-component gaussian statisics computation:
    % Needs to be rebuilt as a multi-dimensional gaussian mixture model
    if k~=c_index, popcmeans(:,k) = NaN; popcstds(:,k) = NaN; popcweights(:,k) = NaN; continue; end;
    % find model fit:
    if numel(data),
        tmp_PEM = PlasmidExpressionModel(data(:,k),CFP_af,CFP_noise,getERFPerPlasmid(analysisParams),drop_threshold,n_components);
        dist = get_fp_dist(tmp_PEM);
        if ~isempty(dist)
            sortable = zeros(n_components,3);
            sortable(:,1) = dist.mu; sortable(:,2) = dist.Sigma(:); sortable(:,3) = dist.weight;
            sortrows(sortable);
            popcmeans(:,k) = sortable(:,1); popcstds(:,k) = sortable(:,2);  popcweights(:,k) = sortable(:,3);
        end
        if k==c_index, PEM = tmp_PEM; end;
    end
end
nonexpressing = sum(nonexpressing_set);

% Horrible kludge: needs to be repaired for release 9.0
analysis_mode = 'geometric';
bins = getBins(analysisParams);
if(strcmp(getUnits(getChannel(colorModel,c_index)),'Boolean')), 
    analysis_mode = 'arithmetic'; 
    bins = BinSequence(log10(min(get_bin_edges(bins))),log10(get_bin_widths(bins)),log10(max(get_bin_edges(bins))),'arithmetic');
end;

% Create (possibly filtered) subpopulation statistics [note: ignore excluded, as already calculated above]
[counts,means,stds] = subpopulation_statistics(bins, data, c_index, analysis_mode, drop_thresholds);
if numel(data)
    plasmid_counts = estimate_plasmids(PEM,means(:,c_index));
    fraction_active = estimate_fraction_active(PEM, means(:,c_index));
else
    plasmid_counts = [];
    fraction_active = [];
    PEM = [];
end
n_out_of_range = size(data,1)-sum(counts);

sampleresults = SampleResults(experiment, samplename, analysisParams, counts, means, stds, ...
    fraction_active, plasmid_counts, PEM, ...
    histograms, popmeans, popstds, poppeaks, excluded, nonexpressing, ...
    popcmeans, popcstds, popcweights);
sampleresults.OutOfRange = n_out_of_range;

% PER_COLOR_CONSTITUTIVE analyzes the inputted batch_description for batch
% analysis and outputs the results and sampleresults for further plotting. 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [results, sampleresults] = per_color_constitutive_analysis(colorModel,batch_description,colors,AP, cloudNames)
% The 'results' here is not a standard ExperimentResults, but a similar scratch structure

TASBESession.warn('TASBE:Analysis','UpdateNeeded','Need to update per_color_constitutive_analysis to use new samplestatistics');
n_conditions = size(batch_description,1);

% check to make sure batch_file has the correct dimensions
if size(batch_description, 2) ~= 2
    TASBESession.error('TASBE:Analysis', 'DimensionMismatch', 'Batch analysis invoked with incorrect number of columns. Make sure batch_file is a n X 2 matrix.');
end
% ensure all files referenced are DataFiles
batch_description(:,2) = cellfun(@(entry)({cellfun(@(f)({ensureDataFile(f)}),entry)}),batch_description(:,2));

data = cell(n_conditions,1);
n_removed = data; n_events = data;
csv_filenames = {};
for i = 1:n_conditions
    condition_name = batch_description{i,1};
    fileset = batch_description{i,2};
    experiment = Experiment(condition_name,'', {0,fileset});
    [data{i},n_removed_sub] = read_data(colorModel, experiment, AP);
    n_removed{i} = [n_removed_sub{1}{:}];
    for j=1:numel(fileset), n_events{i}(j) = size(data{i}{1}{j},1); end;
    if exist('cloudNames', 'var')
        filenames = {cloudNames{i}};
    else
        datafiles = getInducerLevelsToFiles(experiment); % returns array of datafiles
        % convert datafiles to filenames
        filenames = {};
        for k=1:numel(datafiles)
            perInducerDataFiles = datafiles{k};
            perInducerFiles = {};
            for p=1:numel(perInducerDataFiles)
                perInducerFiles{end+1} = getFile(perInducerDataFiles{p});
            end
            filenames{end+1} = perInducerFiles;
        end
    end
    csv_filename = writeFcsPointCloudCSV(colorModel, filenames, data{i});
    if ~isempty(csv_filename)
        csv_filenames = [csv_filenames csv_filename];
    end
end

% Create point cloud header
if TASBEConfig.get('flow.outputPointCloud')
    writePointCloudHeader(colorModel, csv_filenames);
end

% first do all the processing
rawresults = cell(size(colors));
for i=1:numel(colors)
    fprintf(['Processing for color ' colors{i} '...\n']);
    AP = setChannelLabels(AP,{'constitutive',channel_named(colorModel,colors{i})});
    rawresults{i} = process_constitutive_batch( colorModel, batch_description, AP, data);
end

bincenters = get_bin_centers(getBins(AP));

results = cell(n_conditions,1); sampleresults = results;
threshold = TASBEConfig.get('flow.onThreshold');
for i=1:n_conditions
    replicatecounts = numel(rawresults{1}{i,2});
    samplebincounts = cell(replicatecounts,1);
    sample_gmm_means = samplebincounts; sample_gmm_stds = samplebincounts; sample_gmm_weights = samplebincounts;
    samplemeans = zeros(replicatecounts,numel(colors)); samplestds = samplemeans;
    results{i}.condition = batch_description{i,1};
    results{i}.bincenters = bincenters;
    for j=1:numel(colors)
        % kludge needs generalization
        derived_channel = strcmp(getUnits(getChannel(colorModel,j)),'Boolean');
        ER = rawresults{j}{i,1};
        SR = rawresults{j}{i,2};
        rawbincounts = getBinCounts(ER);
        results{i}.bincounts(:,j) = rawbincounts;
        if(~derived_channel)
            results{i}.means(j) = geomean(bincenters',rawbincounts);
            results{i}.stds(j) =  geostd(bincenters',rawbincounts);
        else
            results{i}.means(j) = log10(get_bin_widths(getBins(AP)))/2+wmean(log10(bincenters'),rawbincounts); % KLUDGE FOR BOOLEAN
            results{i}.stds(j) =  wstd(log10(bincenters'),rawbincounts);
        end
        [results{i}.gmm_means(:,j), results{i}.gmm_stds(:,j), results{i}.gmm_weights(:,j)] = get_channel_gmm_results(rawresults{j}{i,1},'constitutive');
        % per-sample histograms
        for k=1:replicatecounts
            samplebincounts{k}(:,j) = SR{k}.BinCounts;
            if(~derived_channel)
                samplemeans(k,j) = geomean(bincenters',SR{k}.BinCounts);
                samplestds(k,j) = geostd(bincenters',SR{k}.BinCounts);
            else
                samplemeans(k,j) = log10(get_bin_widths(getBins(AP)))/2+wmean(log10(bincenters'),SR{k}.BinCounts);
                samplestds(k,j) = wstd(log10(bincenters'),SR{k}.BinCounts);
            end
            color_column = find(colorModel,channel_named(colorModel,colors{j}));
            sample_gmm_means{k}(:,j) = SR{k}.PopComponentMeans(:,color_column);
            sample_gmm_stds{k}(:,j) = SR{k}.PopComponentStandardDevs(:,color_column);
            sample_gmm_weights{k}(:,j) = SR{k}.PopComponentWeights(:,color_column);
            results{i}.n_events_used(k,j) = sum(SR{k}.BinCounts);
            results{i}.n_events_outofrange(k,j) = SR{k}.OutOfRange;
        end
        if(~derived_channel)
            results{i}.stdofmeans(j) = geostd(samplemeans(:,j));
        else
            results{i}.stdofmeans(j) = std(samplemeans(:,j));
        end
        results{i}.stdofstds(j) = mean(samplestds(:,j));
    end
    on_fracs = {};
    off_fracs = {};
    for k=1:replicatecounts
        on_counts = sum(data{i}{1}{k} >= threshold);
        off_counts = sum(data{i}{1}{k} < threshold);
        on_total_count = sum(on_counts);
        off_total_count = sum(off_counts);
        on_frac = on_total_count/(on_total_count+off_total_count);
        off_frac = off_total_count/(on_total_count+off_total_count);
        on_fracs{end+1} = on_frac;
        off_fracs{end+1} = off_frac;
        sampleresults{i}{k} = SampleResults([], batch_description{i,2}{k}, setChannelLabels(AP,colors'), samplebincounts{k}, samplemeans(k,:), samplestds(k,:), ...
            [], [], [], [], [], [], [], [], [], sample_gmm_means{k}, sample_gmm_stds{k}, sample_gmm_weights{k}, on_frac, off_frac);
    end
    results{i}.on_fracMean = mean(cell2mat(on_fracs));
    results{i}.on_fracStd = std(cell2mat(on_fracs));
    results{i}.off_fracMean = mean(cell2mat(off_fracs));
    results{i}.off_fracStd = std(cell2mat(off_fracs));
    results{i}.n_events = n_events{i};
    results{i}.n_events_removed = n_removed{i};
end

%%%%%%%%%%%%%%%%%%
% walk through all results and see if there's reason for statistical concern
max_events = 0;
for i=1:n_conditions
    if ~isempty(results{i}.n_events)
        max_events = max(max_events,max(results{i}.n_events));
    end
end
for i=1:n_conditions
    cur_events = min(results{i}.n_events);
    if ~isempty(cur_events) && (max_events/cur_events > TASBEConfig.get('flow.conditionEventRatioWarning'))
        TASBESession.warn('TASBE:Analysis','HighConditionSizeVariation','High variation in events per condition:\n  max=%i, Condition "%s" = %i',...
            max_events, results{i}.condition, cur_events);
    end
end

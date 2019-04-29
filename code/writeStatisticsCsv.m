% WRITESTATISTICSCSV outputs the stats file with summary data from batch
% analysis.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function statisticsFile = writeStatisticsCsv(colorModel, channels, sampleIds, sampleresults)
    if TASBEConfig.get('flow.outputStatisticsFile')
        baseName = sanitize_filename(TASBEConfig.get('OutputSettings.StemName'));

        path = TASBEConfig.get('flow.dataCSVPath');
        path = end_with_slash(path);
        if ~isdir(path)
            TASBESession.notify('TASBE:Utilities','MakeDirectory','Directory does not exist, attempting to create it: %s',path);
            mkdir(path);
        end

        % First create the default output filename.
        statisticsFile = [path baseName '_statisticsFile.csv'];

        numConditions = numel(sampleIds);
        numComponents = size(sampleresults{1}{1}.PopComponentMeans,1);

        totalCounts = cell(numConditions, 1);
        geoMeans = totalCounts; geoStdDev = totalCounts;
        gmmMeans = totalCounts; gmmStds = totalCounts; gmmWeights = totalCounts;
        onFracs = totalCounts; offFracs = totalCounts;

        replicates = zeros(numConditions, 1);

        for i=1:numConditions
            replicates(i) = numel(sampleresults{i});
            totalCounts{i} = cell(1,replicates(i));
            geoMeans{i} = totalCounts{i}; geoStdDev{i} = totalCounts{i};
            gmmMeans{i} = totalCounts{i}; gmmStds{i} = totalCounts{i}; gmmWeights{i} = totalCounts{i}; 
            onFracs{i} = totalCounts{i}; offFracs{i} = totalCounts{i};
            for j=1:replicates(i)
                totalCounts{i}{j} = sum(sampleresults{i}{j}.BinCounts);
                geoMeans{i}{j} = limitPrecision(sampleresults{i}{j}.Means,4);
                geoStdDev{i}{j} = limitPrecision(sampleresults{i}{j}.StandardDevs,4);
                geoStdDev{i}{j} = limitPrecision(sampleresults{i}{j}.StandardDevs,4);
                gmmMeans{i}{j} = limitPrecision(sampleresults{i}{j}.PopComponentMeans,4);
                gmmStds{i}{j} = limitPrecision(sampleresults{i}{j}.PopComponentStandardDevs,4);
                gmmWeights{i}{j} = limitPrecision(sampleresults{i}{j}.PopComponentWeights,4);
                onFracs{i}{j} = limitPrecision(sampleresults{i}{j}.on_frac,4);
                offFracs{i}{j} = limitPrecision(sampleresults{i}{j}.off_frac,4);
            end
        end

        effective_channels = cell(size(channels));
        for i=1:numel(channels)
            colors = getChannelNames(sampleresults{1}{1}.AnalysisParameters); % assume batch analysis AP
            color_column = find(colorModel,channel_named(colorModel,colors{i}));
            effective_channels{i} = channels{color_column};
        end

        columnNames = buildDefaultStatsFileHeader(effective_channels, numComponents);
        numColumns = numel(columnNames);
        totalReplicates = sum(replicates);

        statsTable = cell(totalReplicates+1, numColumns);
        statsTable(1, 1:numColumns) = columnNames;
        endingRow = 1;  % Because the column labels are in the first row.

        % Put everything in a cell array for Octave
        for i=1:numConditions
            startingRow = endingRow + 1;
            endingRow = startingRow + replicates(i) - 1;
            statsTable(startingRow:endingRow,1:numColumns) = formatDataPerSampleIndivdualColumns(numel(channels), sampleIds{i}, totalCounts{i}, geoMeans{i}, geoStdDev{i}, gmmMeans{i}, gmmStds{i}, gmmWeights{i}, onFracs{i}, offFracs{i});
        end

        % Needed to add column names when I created the tables due to conflicts
        % with the default names.  For a table, the column names must be valid
        % matlab variable names so I filtered out spaces and hypens and
        % replaced them with underscores.
        if (is_octave)
            cell2csv(statisticsFile, statsTable);
        else
            t = table(statsTable);
            writetable(t, statisticsFile, 'WriteVariableNames', false);
        end
    else
        statisticsFile = 'none';
    end
end

function perSampleTable = formatDataPerSampleIndivdualColumns(numChannels, sampleId, totalCounts, means, stddevs, gmm_means, gmm_stds, gmm_weights, on_fracs, off_fracs)
    % SampleId should just be a string. Means and stddevs should be a 1 by
    % number of channels matrix.  TotalCounts should be a 1 by number of
    % channels matrix.
    % Place replicates on separate lines. Padding will be necessary in
    % order to build a table.  Separate into individual columns..
    numReplicates = numel(totalCounts);
    numComponents = size(gmm_means{1},1);
    
    % Number of rows to pad
    rowsOfPadding = numReplicates-1;
    
    % Need to pad with a column vector
    sampleIdPadding = cell(rowsOfPadding, 1);
    
    % Split by the channels so the table will have the correct column labels.
    geoMeans = cell(numReplicates, numChannels);
    geoStdDevs = geoMeans; counts = geoMeans;
    gmmMeans = cell(numReplicates, numChannels * numComponents);
    gmmStds = gmmMeans; gmmWeights = gmmMeans;
    
    for i=1:numChannels        
        for j=1:numReplicates
            counts{j,i} = totalCounts{j}(i);
            geoMeans{j,i} = means{j}(i);
            geoStdDevs{j,i} = stddevs{j}(i);
            for k=1:numComponents
                gmmMeans{j,(i-1)*numComponents + k} = 10.^gmm_means{j}(k,i);
                gmmStds{j,(i-1)*numComponents + k} = 10.^gmm_stds{j}(k,i);
                gmmWeights{j,(i-1)*numComponents + k} = gmm_weights{j}(k,i);
            end
        end
    end
    
    % Pad the sampleId
    ID = [{sampleId}; sampleIdPadding];
    
    perSampleTable = [ID, counts, geoMeans, geoStdDevs, gmmMeans, gmmStds, gmmWeights, transpose(on_fracs), transpose(off_fracs)];
    
end

function fileHeader = buildDefaultStatsFileHeader(channels, numComponents)
    % Default file header to match the default file format.
    numChannels = numel(channels);
    
    binNames = cell(1,numChannels);
    meanNames = binNames; stdDevNames = binNames;
    gmmMeanNames = cell(1,numChannels*numComponents);
    gmmStdNames = gmmMeanNames; gmmWeightNames = gmmMeanNames;
    
    % Not elegant, but it gets the job done.
    for i=1:numChannels
        channelName = sanitizeColumnName([getPrintName(channels{i}) '_' getUnits(channels{i})]);
        binNames{i} = ['BinCount_' channelName];
        meanNames{i} = ['GeoMean_' channelName];
        stdDevNames{i} = ['GeoStdDev_' channelName];
        for j=1:numComponents
            gmmMeanNames{(i-1)*numComponents + j} = sprintf('GMM_Component%i_Mean_%s',j,channelName);
            gmmStdNames{(i-1)*numComponents + j} = sprintf('GMM_Component%i_Std_%s',j,channelName);
            gmmWeightNames{(i-1)*numComponents + j} = sprintf('GMM_Component%i_Weight_%s',j,channelName);
        end
    end
    
    % Don't separate with commas. We want all the column names in a cell
    % array so we can pass them to a table.
    fileHeader = {'Id', binNames{:}, meanNames{:}, stdDevNames{:}, gmmMeanNames{:}, gmmStdNames{:}, gmmWeightNames{:}, 'On Fraction', 'Off Fraction'};
end


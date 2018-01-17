% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function statisticsFile = writeStatisticsCsv(channels, sampleIds, sampleresults, baseName)

    % First create the default output filename.
    statisticsFile = [baseName '_statisticsFile.csv'];
    
    numConditions = numel(sampleIds);
    
    totalCounts = cell(numConditions, 1);
    geoMeans = cell(numConditions, 1);
    geoStdDev = cell(numConditions, 1);
    
    replicates = zeros(numConditions, 1);
    
    for i=1:numConditions
        replicates(i) = numel(sampleresults{i});
        totalCounts{i} = cell(1,replicates(i));
        geoMeans{i} = cell(1,replicates(i));
        geoStdDev{i} = cell(1,replicates(i));
        for j=1:replicates(i)
            totalCounts{i}{j} = sum(sampleresults{i}{j}.BinCounts);
            geoMeans{i}{j} = sampleresults{i}{j}.Means;
            geoStdDev{i}{j} = sampleresults{i}{j}.StandardDevs;
        end
    end
    
    columnNames = buildDefaultStatsFileHeader(channels);
    numColumns = numel(columnNames);
    totalReplicates = sum(replicates);
    
    statsTable = cell(totalReplicates+1, numColumns);
    statsTable(1, 1:numColumns) = columnNames;
    endingRow = 1;  % Because the column labels are in the first row.
    
    % Put everything in a cell array for Octave
    for i=1:numConditions
        startingRow = endingRow + 1;
        endingRow = startingRow + replicates(i) - 1;
        statsTable(startingRow:endingRow,1:numColumns) = formatDataPerSampleIndivdualColumns(channels, sampleIds{i}, totalCounts{i}, geoMeans{i}, geoStdDev{i});
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
end

function perSampleTable = formatDataPerSampleIndivdualColumns(channels, sampleId, totalCounts, means, stddevs)
    % SampleId should just be a string. Means and stddevs should be a 1 by
    % number of channels matrix.  TotalCounts should be a 1 by number of
    % channels matrix.
    % Place replicates on separate lines. Padding will be necessary in
    % order to build a table.  Separate into individual columns..
    numChannels = numel(channels);
    numReplicates = numel(totalCounts);
    
    % Number of rows to pad
    rowsOfPadding = numReplicates-1;
    
    % Need to pad with a column vector
    sampleIdPadding = cell(rowsOfPadding, 1);
    
    % Split by the channels so the table will have the correct column labels.
    geoMeans = cell(numReplicates, numChannels);
    geoStdDevs = cell(numReplicates, numChannels);
    counts = cell(numReplicates, numChannels);
    
    for i=1:numChannels        
        for j=1:numReplicates
            counts{j,i} = totalCounts{j}(i);
            geoMeans{j,i} = means{j}(i);
            geoStdDevs{j,i} = stddevs{j}(i);
        end
    end
    
    % Pad the sampleId
    ID = [{sampleId}; sampleIdPadding];
    
    perSampleTable = [ID, counts, geoMeans, geoStdDevs];
    
end

function fileHeader = buildDefaultStatsFileHeader(channels)
    % Default file header to match the default file format.
    numChannels = numel(channels);
    
    binNames = cell(1,numChannels);
    meanNames = cell(1,numChannels);
    stdDevNames = cell(1,numChannels);
    
    % Not elegant, but it gets the job done.
    for i=1:numChannels
        channelName = getName(channels{i});
        invalidChars = '-|\s';  % Matlab does not like hypens or whitespace in variable names.
        matlabValidVariableNameChannelName = regexprep(channelName,invalidChars,'_');
        binNames{i} = ['BinCount_' matlabValidVariableNameChannelName];
        meanNames{i} = ['GeoMean_' matlabValidVariableNameChannelName];
        stdDevNames{i} = ['GeoStdDev_' matlabValidVariableNameChannelName];
    end
    
    % Don't separate with commas. We want all the column names in a cell
    % array so we can pass them to a table.
    fileHeader = {'ID', binNames{:}, meanNames{:}, stdDevNames{:}};
end


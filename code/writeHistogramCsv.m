% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function histogramFile = writeHistogramCsv(channels, sampleIds, sampleresults, binCenters, baseName)

    % First create the default output filename.
    histogramFile = [baseName '_histogramFile.csv'];
    
    numConditions = numel(sampleIds);
    replicates = zeros(numConditions, 1);
    numBinsPerChannel = numel(binCenters);
    
    % Pull histogram data out of sample results.
    binCounts = cell(numConditions, 1);
    for i=1:numConditions
        replicates(i) = numel(sampleresults{i});
        binCounts{i} = cell(1,replicates(i));
        for j=1:replicates(i)
            binCounts{i}{j} = sampleresults{i}{j}.BinCounts;
        end
    end
    
    columnNames = buildDefaultHistFileHeader(channels);
    numColumns = numel(columnNames);
    totalReplicates = sum(replicates);
    
    histTable = cell(totalReplicates*numBinsPerChannel+1, numColumns);
    histTable(1, 1:numColumns) = columnNames;
    endingRow = 1;  % Because the column labels are in the first row.
    
    for i=1:numConditions
        startingRow = endingRow + 1;
        endingRow = startingRow + replicates(i)*numBinsPerChannel - 1;
        histTable(startingRow:endingRow,1:numColumns) = formatDataPerSample(channels, sampleIds{i}, binCenters, binCounts{i});
    end
    
    % Needed to add column names when I created the tables due to conflicts
    % with the default names.  For a table, the column names must be valid
    % matlab variable names so I filtered out spaces and hypens and
    % replaced them with underscores.
    if (is_octave)
        cell2csv(histogramFile, histTable);
    else
        t = table(histTable);
        writetable(t, histogramFile, 'WriteVariableNames', false);
    end
        
end

function perSampleTable = formatDataPerSample(channels, sampleId, binCenters, counts)
    % The channels are actually the column labels for the data and the
    % binCenters are actually the row labels.  To make writing the data to
    % a CSV file easier, I'm going to include the binCenters in the table
    % like they are data, not just row names.
    
    % File formatted as follows:
    % SampleId, BinCenters, Channel_1_counts, Channel_2_counts, ...
    
    % Only the sampleId needs to be padded to create the table.
    numChannels = numel(channels);
    numReplicates = numel(counts);
    numBinsPerChannel = numel(binCenters);
    
    % Number of rows to pad sample id
    rowsOfPadding = numBinsPerChannel * numReplicates - 1;
    
    % Need to pad the sampleId with a column vector
    sampleIdPadding = cell(rowsOfPadding, 1);
    
    % Pad the sampleId
    sampleIdPadded = [{sampleId}; sampleIdPadding];
    
    % Make sure all rows match up with a bin center
    binCentersForAllReplicates = repmat(binCenters', numReplicates, 1);
    
    % Split by the channels so the table will have the correct column labels.
    binCounts = zeros(numReplicates*numBinsPerChannel, numChannels);
    for i=1:numChannels        
        for j=1:numReplicates
            startRow = (j-1) * numBinsPerChannel +1;
            stopRow = j * numBinsPerChannel;
            binCounts(startRow:stopRow, i) = counts{j}(:,i);
        end
    end

    perSampleTable = [sampleIdPadded, num2cell(binCentersForAllReplicates), num2cell(binCounts)];
end

function fileHeader = buildDefaultHistFileHeader(channels)
    % Default file header to match the default file format.
    numChannels = numel(channels);
    binHeaders = cell(1,numChannels);
    
    % Not elegant, but it gets the job done.
    for i=1:numChannels
        channelName = getName(channels{i});
        invalidChars = '-|\s';  % Matlab does not like hypens or whitespace in variable names.
        matlabValidVariableNameChannelName = regexprep(channelName,invalidChars,'_');
        binHeaders{i} = ['BinCount_' matlabValidVariableNameChannelName];
    end
    
    % Don't separate with commas. We want all the column names in a cell
    % array so we can pass them to a table.
    fileHeader = {'ID', 'BinCenters', binHeaders{:}};
end

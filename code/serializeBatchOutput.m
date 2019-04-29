% SERIALIZEBATCHOUTPUT grabs all the data in separate data structures. Then formats it for output
% files used in batch analysis.    
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [statisticsFile, histogramFile] = serializeBatchOutput(file_pairs, CM, AP, sampleresults)
    channel_names = getChannelNames(sampleresults{1}{1}.AnalysisParameters); % channel names are same across conditions and replicates
    channels = cell(numel(channel_names),1);
    for i=1:numel(channel_names)
        channels{i} = channel_named(CM, channel_names{i});
    end
    sampleIds = file_pairs(:,1);
    binCenters = get_bin_centers(getBins(AP));
    
    % Formats and writes the output to the Statistics file.
    statisticsFile = writeStatisticsCsv(CM, channels, sampleIds, sampleresults);
    
    % Formats and writes the output to the Histogram file.
    histogramFile = writeHistogramCsv(channels, sampleIds, sampleresults, binCenters);
end


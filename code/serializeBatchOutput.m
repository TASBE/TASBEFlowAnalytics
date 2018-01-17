% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [statisticsFile, histogramFile] = serializeBatchOutput(file_pairs, CM, AP, sampleresults, baseName)

    % Grab all the data in separate data structures. Then format for output
    % files.    
    channels = getChannels(CM);
    sampleIds = file_pairs(:,1);
    binCenters = get_bin_centers(getBins(AP));
    
    % Formats and writes the output to the Statistics file.
    statisticsFile = writeStatisticsCsv(channels, sampleIds, sampleresults, baseName);
    
    % Formats and writes the output to the Histogram file.
    histogramFile = writeHistogramCsv(channels, sampleIds, sampleresults, binCenters, baseName);

end


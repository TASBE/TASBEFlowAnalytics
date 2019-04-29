% PROCESS_DATA processes the data from the fcs files into useful arrays. Here
% the structure of processed data mirrors that of the filenames in the
% experiment. 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function sampleresults = process_data( colorModel, experiment, analysisParams, data)
filenames = getInducerLevelsToFiles(experiment); % array of file names
n_conditions = numel(filenames);

% Process each file for each condition in turn, computing results
% incrementally
sampleresults = cell(size(filenames));
for i=1:n_conditions
    perInducerFiles = filenames{i};
    numberOfPerInducerFiles = numel(perInducerFiles);
    if (numberOfPerInducerFiles == 0), TASBESession.warn('TASBE:Analysis','MissingDataFile','An inducer level is missing a data file'); end;
    for j = 1:numberOfPerInducerFiles
        fileName = perInducerFiles{j};
        % Extract statistics        
        sampleresults{i}{j} = compute_sample_statistics(colorModel,experiment,fileName,analysisParams,data{i}{j});
    end
end

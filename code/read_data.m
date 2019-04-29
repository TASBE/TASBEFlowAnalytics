% In READ_DATA, the structure of processed data mirrors that of the filenames in the
% experiment. 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [data,n_removed] = read_data( colorModel, experiment, analysisParams)
datafiles = getInducerLevelsToFiles(experiment); % array of file names
n_conditions = numel(datafiles);

% Process each file for each condition in turn, computing results
% incrementally
data = cell(size(datafiles));
n_removed = data;
for i=1:n_conditions
    perInducerFiles = datafiles{i};
    numberOfPerInducerFiles = numel(perInducerFiles);
    if (numberOfPerInducerFiles == 0), TASBESession.warn('TASBE:Analysis','MissingDataFile','An inducer level is missing a data file'); end;
    for j = 1:numberOfPerInducerFiles
        datafile = perInducerFiles{j};
        datafile = ensureDataFile(datafile);
        % Read data and extract statistics
        [data{i}{j},n_removed{i}{j}] = readfcs_compensated_ERF(colorModel, datafile, getUseAutoFluorescence(analysisParams), true);
    end
end

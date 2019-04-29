% PROCESS_PLUSMINUS_BATCH analyzes the batch_description for plusminus
% analysis, conducts comparisons between batches and plots some of the
% results
%
% batch_description is a cell-array of: {condition_name, inducer_name, batch_names, plus_level_file_pairs, minus_level_file_pairs}
% pm_results is a cell-array of PlusMinusResults
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [pm_results, pm_sampleresults] = process_plusminus_batch(colorModel, batch_description, analysisParams)
batch_size = numel(batch_description);

% check to make sure batch_description has the correct dimensions
category_size = size(batch_description, 1);

for i=1:category_size
    section_size = numel(batch_description{i});
    batch_names_size = numel(batch_description{i}{3});
    if section_size ~= batch_names_size+3
        TASBESession.error('TASBE:Analysis', 'SetDimensionMismatch', 'Plus Minus analysis invoked with incorrect number of sets. Make sure batch_description is a n X 2 matrix with the correct number of sets (size of batch_names including condition_name and inducer_name).');
    end
    for j=4:section_size
        if size(batch_description{i}{j}, 2) ~= 2
            TASBESession.error('TASBE:Analysis', 'ColumnDimensionMismatch', 'Plus Minus analysis invoked with incorrect number of columns. Make sure batch_description is a n X 2 matrix.');
        end
    end
end

% verify that all files exist
% Begin by scanning to make sure all expected files are present
fprintf('Confirming files are present...\n');
for i=1:batch_size
    batch_names_size = numel(batch_description{i}{3});
    level_file_pairs = cell(batch_names_size);
    for j=1:batch_names_size
        level_file_pairs{j} = batch_description{i}{j+3};
    end
    for j=1:size(level_file_pairs,1)
        fileset = level_file_pairs{j,2};
        for k=1:size(fileset,1)
            if ~exist(char(fileset{k,2}),'file')
                TASBESession.error('TASBE:Analysis','MissingFile','Could not find file: %s',char(fileset{k,2}));
            end
        end
    end 
end

pm_results = cell(batch_size,1);
pm_sampleresults = cell(batch_size,1);
for i = 1:batch_size
    batch_names_size = numel(batch_description{i}{3});
    condition_name = batch_description{i}{1};
    inducer_name = batch_description{i}{2}; 
    level_file_pairs = cell(batch_names_size,1);
    experiment_names = cell(batch_names_size,1);
    experiments = cell(batch_names_size,1);
    sampleresults = cell(batch_names_size,1);
    results = cell(batch_names_size,1);
    % Go through batch_names
    for j = 1:batch_names_size
        level_file_pair = batch_description{i}{j+3};
        experiment_name = [condition_name ': ' inducer_name ' ' batch_description{i}{3}{j}];
        experiment = Experiment(experiment_name,{inducer_name}, level_file_pair);
        data{i} = read_data(colorModel, experiment, analysisParams);
        fprintf(['Starting analysis of ' experiment_name '...\n']);
        sampleresult = process_data(colorModel,experiment,analysisParams, data{i});
        result = summarize_data(colorModel,experiment,analysisParams,sampleresult);
        experiment_names{j} = experiment_name;
        experiments{j} = experiment;
        sampleresults{j} = sampleresult;
        results{j} = result;
        level_file_pairs{j} = level_file_pair;
    end
    
    % Comparisons between batches
    fprintf(['Computing comparison for ' condition_name '...\n']);
    comp_results = cell(numel(results)-1,1);
    for j = 1:numel(results)-1
        comp_results{j} = compare_plusminus(results{j},results{end}); % compare against last batch
    end
    pm_sampleresults{i} = sampleresults; % array of arrays of sampleresults
    pm_results{i} = comp_results; % array of arrays of comparisons
    % dump bincounts files
    % get/set should be replaced by push/pop of output settings
    % stemName = TASBEConfig.getexact('OutputSettings.StemName',[]);
    ERROR = []; % This is an ugly kludge, and we should figure out a way to avoid it
    try
        for j = 1:batch_names_size
            TASBEConfig.set('OutputSettings.StemName', [condition_name '-' batch_description{i}{3}{j}]);
            plot_bin_statistics(sampleresults{j}, getInducerLevelsToFiles(experiments{j},1));
        end
        
        TASBEConfig.set('OutputSettings.StemName', condition_name);
        % Plot bin counts on constitutive channel
        constitutive_index = find(colorModel,getChannel(analysisParams, 'constitutive'));
        plot_bin_statistics_all(sampleresults, batch_description{i}{3}, getInducerLevelsToFiles(experiments{1},1), constitutive_index, 'Constitutive');
        % Plot bin counts on output channel
        output_index = find(colorModel,getChannel(analysisParams, 'output'));
        plot_bin_statistics_all(sampleresults, batch_description{i}{3}, getInducerLevelsToFiles(experiments{1},1), output_index, 'Output');
    catch ERROR
    end
    % TASBEConfig.set('OutputSettings.StemName', stemName);
    if ~isempty(ERROR)
        rethrow(ERROR);
    end
end

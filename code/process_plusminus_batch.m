% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [pm_results pm_sampleresults] = process_plusminus_batch( colorModel, batch_description, analysisParams, OSbin)
% batch_description is a cell-array of: {condition_name, inducer_name, plus_level_file_pairs, minus_level_file_pairs}
% pm_results is a cell-array of PlusMinusResults

batch_size = numel(batch_description);

% verify that all files exist
% Begin by scanning to make sure all expected files are present
fprintf('Confirming files are present...\n');
for i=1:batch_size
    level_file_pairs = [batch_description{i}{3}; batch_description{i}{4}];
    for j=1:size(level_file_pairs,1),
        fileset = level_file_pairs{j,2};
        for k=1:numel(fileset),
            if ~exist(fileset{k},'file'),
                error('Could not find file: %s',fileset{k});
            end
        end
    end
end

pm_results = cell(batch_size,1);
pm_sampleresults = cell(batch_size,2);
for i = 1:batch_size
    condition_name = batch_description{i}{1};
    inducer_name = batch_description{i}{2};
    p_level_file_pairs = batch_description{i}{3};
    m_level_file_pairs = batch_description{i}{4};
    
    experiment_name = [condition_name ': ' inducer_name ' +'];
    p_experiment = Experiment(experiment_name,{inducer_name}, p_level_file_pairs);
    fprintf(['Starting analysis of ' experiment_name '...\n']);
    p_sampleresults = process_data(colorModel,p_experiment,analysisParams);
    p_results = summarize_data(colorModel,p_experiment,analysisParams,p_sampleresults);
    
    experiment_name = [condition_name ': ' inducer_name ' -'];
    m_experiment = Experiment(experiment_name,{inducer_name}, m_level_file_pairs);
    fprintf(['Starting analysis of ' experiment_name '...\n']);
    m_sampleresults = process_data(colorModel,m_experiment,analysisParams);
    m_results = summarize_data(colorModel,m_experiment,analysisParams,m_sampleresults);
    
    fprintf(['Computing comparison for ' condition_name '...\n']);
    pm_results{i} = compare_plusminus(p_results,m_results);
    pm_sampleresults{i,1} = p_sampleresults;
    pm_sampleresults{i,2} = m_sampleresults;
    
    if nargin>3, % if supplied an output setting, dump bincounts files
        OSp = OSbin; OSp.StemName = [OSp.StemName condition_name '-plus'];
        plot_bin_statistics(p_sampleresults,OSp);
        OSm = OSbin; OSm.StemName = [OSm.StemName condition_name '-minus'];
        plot_bin_statistics(m_sampleresults,OSm);
    end
end

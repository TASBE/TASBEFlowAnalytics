% Function that runs plusminus analysis given a template spreadsheet. A
% TemplateExtraction object and optional Color Model are inputs. Outputs
% results and batch descriptions used for test functions. 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
function [all_results, all_batch_descrips] = plusminus_analysis_excel(extractor, CM)
    % Reset and update TASBEConfig and get exp name
    extractor.TASBEConfig_updates();
    path = extractor.path;
    TASBEConfig.set('template.displayErrors', 1);
    experimentName = extractor.getExcelValue('experimentName', 'char'); 
    TASBEConfig.set('template.displayErrors', 0);
    % Find preference_row
    preference_row = 0;
    col_num = extractor.getColNum('first_compGroup_PM');
    sh_num = extractor.getSheetNum('first_compGroup_PM');
    for i=1:size(extractor.sheets{sh_num},1)
        try
            value = extractor.getExcelValuePos(sh_num, i, col_num, 'char');
            if strcmp(value, 'Required: Instructions to use button below')
                preference_row = i + 2;
                break
            end
        catch
            continue
        end
    end
    if preference_row == 0
        TASBESession.error('plusminus_analysis_excel', 'MissingPreference', 'No end row number found for plusminus analysis. Make sure run button and instructions are in column A.')
    end

    % Load the color model
    if nargin < 2
        % Obtain the CM_name
        try
            coords = extractor.getExcelCoordinates('inputName_CM', 2);
            CM_name = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'char'); 
            [~,name,~] = fileparts(CM_name);
            CM_name = [name '.mat'];
        catch
            TASBESession.warn('plusminus_analysis_excel', 'MissingPreference', 'Missing CM Filename in "Comparative Analysis" sheet. Looking in "Calibration" sheet.');
            try
                CM_name = extractor.getExcelValue('outputName_CM', 'char');
                [~,name,~] = fileparts(CM_name);
                CM_name = [name '.mat'];
            catch
                TASBESession.warn('plusminus_analysis_excel', 'MissingPreference', 'Missing Output Filename in "Calibration" sheet. Defaulting to exp name.');
                CM_name = [experimentName '-ColorModel.mat'];
            end
        end
        
        % Obtain the CM_path
        try
            coords = extractor.getExcelCoordinates('inputPath_CM', 2);
            CM_path = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'char');
            CM_path = make_filename_absolute(CM_path, path);
        catch
            TASBESession.warn('plusminus_analysis_excel', 'MissingPreference', 'Missing CM Filepath in "Comparative Analysis" sheet. Looking in "Calibration" sheet.'); 
            try
                CM_path = extractor.getExcelValue('outputPath_CM', 'char');
                CM_path = make_filename_absolute(CM_path, path);
            catch
                TASBESession.warn('plusminus_analysis_excel', 'MissingPreference', 'Missing Output Filepath in "Calibration" sheet. Defaulting to template path.'); 
                CM_path = path;
            end
        end
        CM_file = [CM_path CM_name];

        try 
            load(CM_file);
        catch
            TASBESession.warn('plusminus_analysis_excel', 'MissingPreference', 'Could not load CM file, creating a new one.');
            CM = make_color_model_excel(extractor);
        end
    end

    % Set TASBEConfigs and create variables needed to run plusminus analysis
    % Analyze on a histogram of 10^[first] to 10^[third] ERF, with bins every 10^[second]
    try
        coords = extractor.getExcelCoordinates('binseq_min', 2);
        binseq_min = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
        coords = extractor.getExcelCoordinates('binseq_max', 2);
        binseq_max = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
        coords = extractor.getExcelCoordinates('binseq_pdecade', 2);
        binseq_pdecade = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
        bins = BinSequence(binseq_min, (1/binseq_pdecade), binseq_max, 'log_bins');
    catch
        bins = BinSequence();
    end

    % Designate which channels have which roles
    [channel_roles, ~] = getChannelRoles(CM, extractor);
    
    if isempty(channel_roles)
        TASBESession.warn('plusminus_analysis_excel', 'MissingPreference', 'Missing constitutive, input, output in "Calibration" sheet');
        APs = {AnalysisParameters(bins,{})};
    else
        APs = {};
        for j=1:numel(channel_roles)
            APs{end+1} = AnalysisParameters(bins,channel_roles{j});
        end
    end
    
    % Obtain the necessary sample filenames and print names
    sample_names = {};
    file_names = {};
    sh_num3 = extractor.getSheetNum('first_compGroup_PM');
    sh_num2 = extractor.getSheetNum('first_sample_num');
    first_sample_row = extractor.getRowNum('first_sample_num');
    sample_num_col = extractor.getColNum('first_sample_num');
    col_names = {};
    row_nums = {};
    col_nums = {};
    last_sampleColName_row = preference_row - 5;
    % Determine the number of plusminus analysis to run
    for i=extractor.getRowNum('primary_sampleColName_PM'): last_sampleColName_row
        col_name = {};
        try
            value = extractor.getExcelValuePos(sh_num3, i, extractor.getColNum('primary_sampleColName_PM'), 'char');
            col_name{end+1} = value;
            row_nums{end+1} = i;
        catch
            try
                value = num2str(extractor.getExcelValuePos(sh_num3, i, extractor.getColNum('primary_sampleColName_PM'), 'numeric'));
                if ~isempty(value)
                    col_name{end+1} = value;
                    row_nums{end+1} = i;
                end
            catch 
                continue
            end
        end
        try
            col_name{end+1} = extractor.getExcelValuePos(sh_num3, i, extractor.getColNum('secondary_sampleColName_PM'), 'char');
            if ~isempty(col_name)
                col_names{end+1} = col_name;
            end
        catch
            try
                value = num2str(extractor.getExcelValuePos(sh_num3, i, extractor.getColNum('secondary_sampleColName_PM'), 'numeric'));
                if ~isempty(value)
                    col_name{end+1} = value;
                end
                if ~isempty(col_name)
                    col_names{end+1} = col_name;
                end
            catch
                % no secondary column add col_name to col_names and
                % continue
                if ~isempty(col_name)
                    col_names{end+1} = col_name;
                end
                continue
            end
        end
    end

    % Go through columns to find column numbers of primary and secondary col names
    for i=1:numel(col_names)
        col_name = col_names{i};
        col_num = {};
        for j=1:numel(col_name)
            pos = find(ismember(extractor.col_names, col_name{j}), 1);
            if ~isempty(pos)
                col_num{j} = pos;
            end
        end
        if numel(col_num) ~= numel(col_name)
            TASBESession.error('plusminus_analysis_excel', 'InvalidColumnName', 'Primary/secondary column names in "Comparative Analysis" does not match with any column name in "Samples".');
        end
        col_nums{end+1} = col_num;
    end
    
    % Go though comparison groups and store values and column numbers in
    % cell array
    comp_groups = {};
    comp_group_names = {};
    first_group_col = extractor.getColNum('first_compGroup_PM');
    outputNames = {};
    outputPaths = {};
    plotPaths = {};
    inducer_names = {};
    for i=1:numel(row_nums)
        if i == numel(row_nums)
            end_row = last_sampleColName_row;
        else
            end_row = row_nums{i+1}-1; 
        end
        comp_group = {};
        for j=row_nums{i}:end_row
            try
                group = extractor.getExcelValuePos(sh_num3, j, first_group_col, 'char');
            catch
                try
                    group = num2str(extractor.getExcelValuePos(sh_num3, j, first_group_col, 'numeric'));
                    if isempty(group)
                        continue
                    end
                catch 
                    continue
                end
            end
            [group_names, values] = getCompGroups(group);
            % Go through group_names and find column numbers
            pos = {};
            for k=1:numel(group_names)
                temp_pos = find(ismember(extractor.col_names, group_names{k}), 1);
                if isempty(temp_pos)
                    TASBESession.error('plusminus_analysis_excel', 'InvalidColumnName', 'Sample column name, %s, under Comparison Groups in "Comparative Analysis" does not match with any column name in "Samples".', col_name);
                else
                    pos{end+1} = temp_pos;
                end
            end
            comp_group{end+1} = {pos, values};
            comp_group_names{end+1} = group;
        end
        if isempty(comp_group)
            comp_group{end+1} = {};
        end
        comp_groups{end+1} = comp_group;
        % Get unique preferences
        % Obtain output name
        try 
            outputName = extractor.getExcelValuePos(sh_num3, row_nums{i}, extractor.getColNum('outputName_PM'), 'char');
            [~, name, ~] = fileparts(outputName);
            outputNames{end+1} = [name '.mat'];
        catch
            TASBESession.warn('plusminus_analysis_excel', 'MissingPreference', 'Missing Output File Name for Plusminus Analysis %s in "Comparative Analysis" sheet', num2str(i));
            if i > 1
                outputNames{end+1} = [experimentName '-CompAnalysis' num2str(i) '.mat'];
            else
                outputNames{end+1} = [experimentName '-CompAnalysis.mat'];
            end
        end
        
        % Obtain output path
        try
            outputPath = extractor.getExcelValuePos(sh_num3, row_nums{i}, extractor.getColNum('outputPath_PM'), 'char');
            outputPath = make_filename_absolute(outputPath, path);
            outputPaths{end+1} = outputPath;
        catch
            TASBESession.warn('plusminus_analysis_excel', 'MissingPreference', 'Missing Output File Path for Plusminus Analysis %s in "Comparative Analysis" sheet', num2str(i));
            outputPaths{end+1} = path;
        end
        
        try
            plotPath_coord = extractor.getExcelCoordinates('plots.plotPath');
            plot_path = extractor.getExcelValuePos(sh_num3, row_nums{i}, plotPath_coord{3}{3}, 'char');
            plot_path = make_filename_absolute(plot_path, path);
            plotPaths{end+1} = plot_path;
        catch
            TASBESession.warn('plusminus_analysis_excel', 'MissingPreference', 'Missing plot path for Plusminus Analysis %s in "Comparative Analysis" sheet', num2str(i));
            plot_path = end_with_slash(fullfile(path, 'plots/'));
            plotPaths{end+1} = plot_path;
        end
        
        inducer_names{end+1} = col_names{i}{1};
    end

    % Get list of distinct values from primary sample column names using
    % condition key col_nums{i}. Also get list of distinct values (in
    % order) from secondary sample column names
    condition_col = extractor.getColNum('first_condition_key');
    condition_sh = extractor.getSheetNum('first_condition_key');
    first_condition_row = extractor.getRowNum('first_condition_key');
    all_keys = {};
    for i=1:numel(col_names)
        col_name = col_names{i};
        keys = {};
        for j=first_condition_row:size(extractor.sheets{condition_sh}, 1)
            try
                value = extractor.getExcelValuePos(condition_sh, j, condition_col, 'char');
                if ~isempty(strfind(value, 'Sample Column Name'))
                    try
                        column_name = extractor.getExcelValuePos(condition_sh, j, condition_col+1, 'char');
                        ind = find(ismember(col_name, column_name), 1);
                        if ~isempty(ind)
                            % get keys
                            key = {};
                            for k=j+2:size(extractor.sheets{condition_sh}, 1)
                                try
                                    key{end+1} = extractor.getExcelValuePos(condition_sh, k, condition_col, 'char');
                                catch
                                    try
                                        part = extractor.getExcelValuePos(condition_sh, k, condition_col, 'numeric');
                                        if isempty(part)
                                            break
                                        end
                                        key{end+1} = part;
                                    catch
                                        break
                                    end
                                end
                            end
                            keys{end+1} = key;
                        end
                    catch
                        try
                            column_name = num2str(extractor.getExcelValuePos(condition_sh, j, condition_col+1, 'numeric'));
                            if isempty(column_name)
                                continue
                            end
                            ind = find(ismember(col_name, column_name), 1);
                            if ~isempty(ind)
                                % get keys
                                key = {};
                                for k=j+2:size(extractor.sheets{condition_sh}, 1)
                                    try
                                        key{end+1} = extractor.getExcelValuePos(condition_sh, k, condition_col, 'char');
                                    catch
                                        try
                                            part = extractor.getExcelValuePos(condition_sh, k, condition_col, 'numeric');
                                            if isempty(part)
                                                break
                                            end
                                            key{end+1} = part;
                                        catch
                                            break
                                        end
                                    end
                                end
                                keys{end+1} = key;
                            end
                        catch
                            continue
                        end
                    end
                end
            catch
                continue
            end
        end
        all_keys{end+1} = keys;
    end
    
    all_results = {};
    all_batch_descrips = {};
    for i=1:numel(col_names)
        TASBEConfig.set('plots.plotPath', plotPaths{i});
        outputName = outputNames{i};
        outputPath = outputPaths{i};
        keys = all_keys{i};
        comp_group = comp_groups{i};
        col_num = col_nums{i};
        col_name = col_names{i};
        batch_description = {}; % contains all the sample row numbers for groups
        % Find the different sets for each group 
        for k=1:numel(comp_group)
            group = {};
            for j=first_sample_row:extractor.getRowNum('last_sample_num')
                try
                    if isempty(comp_group{k})
                        group{end+1} = j;
                    else
                        %Go through all of the sample column names and make sure all values equal
                        equal = true;
                        for x=1:numel(comp_group{k}{1})
                            value = extractor.getExcelValuePos(sh_num2, j, comp_group{k}{1}{x});
                            if isa(value, 'numeric')
                                value = num2str(value);
                            end
                            if ~strcmp(value, comp_group{k}{2}{x})
                                equal = false;
                            end
                        end
                        if equal
                            group{end+1} = j;
                        end
                    end
                catch
                    continue
                end
            end
            if ~isempty(comp_group{k})
                if ~isa(comp_group{k}{2}, 'char')
                    batch_description{end+1} = {[comp_group_names{k}]; col_name{1}; keys{1}; group};
                else
                    batch_description{end+1} = {[comp_group_names{k}]; col_name{1}; keys{1}; group};
                end
            else
                batch_description{end+1} = {experimentName; col_name{1}; keys{1}; group};
            end
        end

        % Go through row numbers in groups and extract and organize by keys{1}
        for k=1:numel(batch_description)
            sets = {};
            for j=1:numel(batch_description{k}{4})
                try
                    value = extractor.getExcelValuePos(sh_num2, batch_description{k}{4}{j}, col_num{1});
                    ind = find(ismember(keys{1}, value), 1);
                    if ~isempty(ind)
                        try 
                            isempty(sets{ind})
                            sets{ind}{end+1} = batch_description{k}{4}{j};
                        catch
                            sets{ind} = {batch_description{k}{4}{j}};
                        end
                    end
                catch
                    continue
                end
            end
            for j=1:numel(sets)
                batch_description{k}{3+j} = sets{j};
            end
        end

        % Reorder the sample within the sets if applicable to match with the
        % order of the secondary sample column name (else just keep old order)
        if numel(col_name) > 1 && numel(keys) > 1
            for z=1:numel(batch_description)
                for j=4:numel(batch_description{z})
                    set = batch_description{z}{j};
                    ordered_set = {};
                    for k=1:numel(set)
                        % Get value at col_num{2} and compare with keys{2}
                        try
                            value = extractor.getExcelValuePos(sh_num2, set{k}, col_num{2});
                            ind = zeros(0,0);
                            if isa(value, 'numeric')
                                for c=1:numel(keys{2})
                                    if keys{2}{c} == value
                                        ind = c;
                                        break
                                    end
                                end
                            else
                                ind = find(ismember(keys{2}, value), 1);
                            end
                            if ~isempty(ind)
                                ordered_set{ind,1} = value;
                                ordered_set{ind,2} = getExcelFilename(extractor, set{k});
                            end
                        catch
                            continue
                        end
                    end
                    batch_description{z}{j} = ordered_set;
                end
            end
        else
            % do the same but no reordering
            for z=1:numel(batch_description)
                for j=4:numel(batch_description{z})
                    set = batch_description{z}{j};
                    ordered_set = {};
                    for k=1:numel(set)
                        % Get value at col_num{2}}
                        if numel(col_name) > 1
                            try
                                value = extractor.getExcelValuePos(sh_num2, set{k}, col_num{2});
                                ordered_set{end+1,1} = value;
                                ordered_set{end,2} = getExcelFilename(extractor, set{k});
                            catch
                                continue
                            end
                        else
                            try
                                value = extractor.getExcelValuePos(sh_num2, set{k}, col_num{1});
                                ordered_set{end+1,1} = k; % default to just index
                                ordered_set{end,2} = getExcelFilename(extractor, set{k});
                            catch
                                continue
                            end
                        end
                    end
                    batch_description{z}{j} = ordered_set;
                end
            end
        end
        
        assignin('base','batch_description',batch_description);
        
        % Go through all possible APs
        for j=1:numel(APs)
            AP = APs{j};
            % Ignore any bins with less than valid count as noise
            try
                coords = extractor.getExcelCoordinates('minValidCount', 2);
                minValidCount = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
                if isempty(minValidCount)
                    error('empty preference');
                end
                AP=setMinValidCount(AP,minValidCount);
            catch
                TASBESession.warn('plusminus_analysis_excel', 'ImportantMissingPreference', 'Missing Min Valid Count in "Comparative Analysis" sheet');
            end
            
            % Corresponds to pem_drop_threshold for histogram computation
            try
                coords = extractor.getExcelCoordinates('minValidau', 2);
                minValidau = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
                if isempty(minValidau)
                    error('empty preference');
                end
                AP=setPemDropThreshold(AP,minValidau);
            catch
                TASBESession.warn('plusminus_analysis_excel', 'ImportantMissingPreference', 'Missing Min Valid au in "Samples" sheet');
            end

            % Add autofluorescence back in after removing for compensation?
            try
                coords = extractor.getExcelCoordinates('autofluorescence', 2);
                autofluorescence = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
                if isempty(autofluorescence)
                    error('empty preference');
                end
                AP=setUseAutoFluorescence(AP,autofluorescence);
            catch
                TASBESession.warn('plusminus_analysis_excel', 'ImportantMissingPreference', 'Missing Use Auto Fluorescence in "Comparative Analysis" sheet');
            end

            try
                coords = extractor.getExcelCoordinates('minFracActive', 2);
                minFracActive = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
                if isempty(minFracActive)
                    error('empty preference');
                end
                AP=setMinFractionActive(AP,minFracActive);
            catch
                TASBESession.warn('plusminus_analysis_excel', 'ImportantMissingPreference', 'Missing Min Fraction Active in "Comparative Analysis" sheet');
            end
            
            if j > 1
                outputName_parts = strtrim(strsplit(outputName, '.'));
                outputName = [outputName_parts{1} '-' num2str(j) '.' outputName_parts{2}];
                for z=1:numel(batch_description)
                    batch_description{z}{1} = [batch_description{z}{1} '-' num2str(j)];
                end
            end

            % Execute the actual analysis
            results = process_plusminus_batch(CM, batch_description, AP);
            all_results{end+1} = results;
            all_batch_descrips{end+1} = batch_description;
            % Make additional output plots
            for k=1:numel(results)
                TASBEConfig.set('OutputSettings.StemName', batch_description{k}{1});
                plot_plusminus_comparison(results{k}, batch_description{k}{3});
            end
            
            if ~isdir(outputPath)
                sanitized_path = strrep(outputPath, '/', '&#47;');
                sanitized_path = strrep(sanitized_path, '\', '&#92;');
                sanitized_path = strrep(sanitized_path, ':', '&#58;');
                TASBESession.notify('OutputFig','MakeDirectory','Directory does not exist, attempting to create it: %s',sanitized_path);
                mkdir(outputPath);
            end
            save('-V7',[outputPath outputName],'batch_description','AP','results');
        end
    end
end

% Helper function that returns cell array of sample column names and their corresponding values
% for row in Comparison Groups
function [group_names, values] = getCompGroups(group)
    group_names = {};
    values = {};
    pairs = strtrim(strsplit(group, ','));
    for i=1:numel(pairs)
        sections = strtrim(strsplit(pairs{i}, '='));
        if numel(sections) ~= 2
            TASBESession.error('plusminus_analysis_excel', 'InvalidCompGroup', 'Comparison groups must come in "sample column = value" pairs. %s not valid.', group);
        end
        group_names{end+1} = sections{1};
        values{end+1} = sections{2};
    end
    if numel(group_names) ~= numel(values)
        TASBESession.error('plusminus_analysis_excel', 'InvalidCompGroup', 'Comparison groups must come in "sample column = value" pairs. %s not valid.', group);
    end
end

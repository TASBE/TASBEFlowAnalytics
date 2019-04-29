% BATCH_ANALYSIS_EXCEL runs batch analysis given a template spreadsheet. A
% TemplateExtraction object and optional Color Model are inputs. Outputs
% results, stats file, and histogram file used in test functions.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.
function [results, statisticsFile, histogramFile] = batch_analysis_excel(extractor, CM)
    % Reset and update TASBEConfig and obtain experiment name
    extractor.TASBEConfig_updates();
    TASBEConfig.set('template.displayErrors', 1);
    experimentName = extractor.getExcelValue('experimentName', 'char');
    preference_row = extractor.getRowNum('last_sample_num') + 5;
    TASBEConfig.set('template.displayErrors', 0);
    
    path = extractor.path;
    
    % Double checking the preference_row
    preference_row2 = 0;
    col_num = extractor.getColNum('last_sample_num');
    sh_num = extractor.getSheetNum('last_sample_num');
    for i=1:size(extractor.sheets{sh_num},1)
        try
            value = extractor.getExcelValuePos(sh_num, i, col_num, 'char');
            if strcmp(value, 'Required: Instructions to use button below')
                preference_row2 = i + 2;
                break
            end
        catch
            continue
        end
    end
    if preference_row ~= preference_row2
        TASBESession.error('batch_analysis_excel', 'MissingPreference', 'Make sure there is an empty row between last sample and Preferences for Batch Analysis')
    end
   
    % Load the color model
    if nargin < 2
        % Obtain the CM_name
        try
            coords = extractor.getExcelCoordinates('inputName_CM', 1);
            CM_name = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'char'); 
            [~,name,~] = fileparts(CM_name);
            CM_name = [name '.mat'];
        catch
            TASBESession.warn('batch_analysis_excel', 'MissingPreference', 'Missing CM Filename in "Samples" sheet. Looking in "Calibration" sheet.');
            try
                CM_name = extractor.getExcelValue('outputName_CM', 'char');
                [~,name,~] = fileparts(CM_name);
                CM_name = [name '.mat'];
            catch
                TASBESession.warn('batch_analysis_excel', 'MissingPreference', 'Missing Output Filename in "Calibration" sheet. Defaulting to exp name.');
                CM_name = [experimentName '-ColorModel.mat'];
            end
        end
        
        % Obtain the CM_path
        try
            coords = extractor.getExcelCoordinates('inputPath_CM', 1);
            CM_path = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'char');
            CM_path = make_filename_absolute(CM_path, path);
        catch
            TASBESession.warn('batch_analysis_excel', 'MissingPreference', 'Missing CM Filepath in "Samples" sheet. Looking in "Calibration" sheet.'); 
            try
                CM_path = extractor.getExcelValue('outputPath_CM', 'char');
                CM_path = make_filename_absolute(CM_path, path);
            catch
                TASBESession.warn('batch_analysis_excel', 'MissingPreference', 'Missing Output Filepath in "Calibration" sheet. Defaulting to template path.'); 
                CM_path = path;
            end
        end
        CM_file = [CM_path CM_name];

        try 
            load(CM_file);
        catch
            TASBESession.warn('batch_analysis_excel', 'MissingPreference', 'Could not load CM file, creating a new one.');
            CM = make_color_model_excel(extractor);
        end
    end

    % Set TASBEConfigs and create variables needed to run batch analysis
    try
        coords = extractor.getExcelCoordinates('plots.plotPath', 2);
        plot_path = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'char');
        plot_path = make_filename_absolute(plot_path, path);
        TASBEConfig.set('plots.plotPath', plot_path);
    catch
        TASBESession.warn('batch_analysis_excel', 'MissingPreference', 'Missing plot path in "Samples" sheet');
        plot_path = end_with_slash(fullfile(path, 'plots/'));
        TASBEConfig.set('plots.plotPath', plot_path);
    end

    % Analyze on a histogram of 10^[first] to 10^[third] ERF, with bins every 10^[second]
    try
        coords = extractor.getExcelCoordinates('binseq_min', 1);
        binseq_min = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
        coords = extractor.getExcelCoordinates('binseq_max', 1);
        binseq_max = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
        coords = extractor.getExcelCoordinates('binseq_pdecade', 1);
        binseq_pdecade = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
        bins = BinSequence(binseq_min, (1/binseq_pdecade), binseq_max, 'log_bins');
    catch
        bins = BinSequence();
    end

    % Designate which channels have which roles
    [channel_roles, print_names] = getChannelRoles(CM, extractor);
    
    if isempty(channel_roles)
        TASBESession.warn('batch_analysis_excel', 'MissingPreference', 'Missing constitutive, input, output in "Calibration" sheet');
        APs = {AnalysisParameters(bins,{})};
    else
        APs = {};
        for j=1:numel(channel_roles)
            APs{end+1} = AnalysisParameters(bins,channel_roles{j});
        end
    end
    
    % Obtain template number for cloud names
    cloudNames = {};
    try
        coords = extractor.getExcelCoordinates('cloudName_BA'); 
        template_num = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
    catch
        TASBESession.warn('batch_analysis_excel', 'MissingPreference', 'Missing filename template number Point Cloud Filename in "Samples" sheet');
        template_num = 0;
    end
    
    % Obtain the necessary sample filenames and print names
    sample_names = {};
    file_names = {};
    sh_num2 = extractor.getSheetNum('first_sample_num');
    first_sample_row = extractor.getRowNum('first_sample_num');
    sample_num_col = extractor.getColNum('first_sample_num');
    sample_name_col = find(ismember(extractor.col_names, 'SAMPLE NAME'), 1);
    if isempty(sample_name_col)
        TASBESession.error('batch_analysis_excel', 'InvalidHeaderName', 'The header, SAMPLE NAME, does not match with any column titles in "Samples" sheet.');
    end
    sample_exclude_col = find(ismember(extractor.col_names, 'Exclude from Batch Analysis'), 1);
    if isempty(sample_exclude_col)
        TASBESession.error('batch_analysis_excel', 'InvalidHeaderName', 'The header, Exclude from Batch Analysis, does not match with any column titles in "Samples" sheet.');
    end
    for i=first_sample_row:size(extractor.sheets{sh_num2},1)
        try
            num = extractor.getExcelValuePos(sh_num2, i, sample_num_col, 'numeric');
            if isempty(num)
                break
            end
        catch
            break
        end
        % check if sample should be included in batch analysis
        try
            extractor.getExcelValuePos(sh_num2, i, sample_exclude_col, 'char');
        catch
            sample_names{end+1} = extractor.getExcelValuePos(sh_num2, i, sample_name_col, 'char');
            datafiles = getExcelFilename(extractor, i);
            file_names{end+1} = datafiles;
            % Obtain point cloud name
            if template_num ~= 0
                name = getCloudName(extractor, i, template_num);
                cloudNames{end+1} = name;
            end
        end
    end

    % Make a map of condition names to file sets
    file_pairs = {};
    file_pairs(:,1) = sample_names;
    file_pairs(:,2) = file_names;
    
    for i=1:numel(APs)
        AP = APs{i};
        try 
            coords = extractor.getExcelCoordinates('OutputSettings.StemName', 1);   
            stemName = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'char');
        catch
            TASBESession.warn('batch_analysis_excel', 'MissingPreference', 'Missing Stem Name in "Samples" sheet');
            stemName = experimentName;
        end
        
        % Obtain output name
        try
            coords = extractor.getExcelCoordinates('outputName_BA'); 
            outputName = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'char');
            [~, name, ~] = fileparts(outputName);
            outputName = [name '.mat'];
        catch
            TASBESession.warn('batch_analysis_excel', 'MissingPreference', 'Missing Output File Name for Batch Analysis in "Samples" sheet');
            outputName = [experimentName '-BatchAnalysis.mat'];
        end
        
        % Obtain output path
        try
            coords = extractor.getExcelCoordinates('outputPath_BA'); 
            outputPath = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'char');
            outputPath = make_filename_absolute(outputPath, path);
        catch
            TASBESession.warn('batch_analysis_excel', 'MissingPreference', 'Missing Output File Path in "Samples" sheet');
            outputPath = path;
        end
        
        % Obtain stat name
        try
            coords = extractor.getExcelCoordinates('statName_BA'); 
            statName = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'char');
        catch
            TASBESession.warn('batch_analysis_excel', 'MissingPreference', 'Missing Statistics Filename for Batch Analysis in "Samples" sheet');
            statName = stemName;
        end
        
        % Obtain stat path
        try
            coords = extractor.getExcelCoordinates('statPath_BA'); 
            statPath = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'char');
            statPath = make_filename_absolute(statPath, path);
        catch
            TASBESession.warn('batch_analysis_excel', 'MissingPreference', 'Missing Statistics Filepath in "Samples" sheet');
            statPath = path;
        end
        
        % Obtain cloud path
        try
            coords = extractor.getExcelCoordinates('cloudPath_BA'); 
            cloudPath = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'char');
            cloudPath = make_filename_absolute(cloudPath, path);
        catch
            TASBESession.warn('batch_analysis_excel', 'MissingPreference', 'Missing Point Cloud Filepath in "Samples" sheet');
            cloudPath = path;
        end
        
        if i > 1
            outputName_parts = strtrim(strsplit(outputName, '.'));
            outputName = [outputName_parts{1} num2str(i) '.' outputName_parts{2}];
            stemName = [stemName num2str(i)];
            statName = [statName num2str(i)];
        end
        
        TASBEConfig.set('OutputSettings.StemName', stemName);
        TASBEConfig.set('flow.pointCloudPath',cloudPath);
        TASBEConfig.set('flow.dataCSVPath',statPath);
        
        % Ignore any bins with less than valid count as noise
        try
            coords = extractor.getExcelCoordinates('minValidCount', 1);
            minValidCount = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
            if isempty(minValidCount)
                error('empty preference');
            end
            AP=setMinValidCount(AP,minValidCount);
        catch
            TASBESession.warn('batch_analysis_excel', 'ImportantMissingPreference', 'Missing Min Valid Count in "Samples" sheet');
        end
        
        % Corresponds to pem_drop_threshold for histogram computation
        try
            coords = extractor.getExcelCoordinates('minValidau', 1);
            minValidau = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
            if isempty(minValidau)
                error('empty preference');
            end
            AP=setPemDropThreshold(AP,minValidau);
        catch
            TASBESession.warn('batch_analysis_excel', 'ImportantMissingPreference', 'Missing Min Valid au in "Samples" sheet');
        end

        % Add autofluorescence back in after removing for compensation?
        try
            coords = extractor.getExcelCoordinates('autofluorescence', 1);
            autofluorescence = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
            if isempty(autofluorescence)
                error('empty preference');
            end
            AP=setUseAutoFluorescence(AP,autofluorescence);
        catch
            TASBESession.warn('batch_analysis_excel', 'ImportantMissingPreference', 'Missing Use Auto Fluorescence in "Samples" sheet');
        end

        try
            coords = extractor.getExcelCoordinates('minFracActive', 1);
            minFracActive = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
            if isempty(minFracActive)
                error('empty preference');
            end
            AP=setMinFractionActive(AP,minFracActive);
        catch
            TASBESession.warn('batch_analysis_excel', 'ImportantMissingPreference', 'Missing Min Fraction Active in "Samples" sheet');
        end
        
        if ~isempty(cloudNames)
            [results, sampleresults] = per_color_constitutive_analysis(CM,file_pairs,print_names,AP, cloudNames);
        else
            [results, sampleresults] = per_color_constitutive_analysis(CM,file_pairs,print_names,AP);
        end
       
        % Make output plots
        plot_batch_histograms(results,sampleresults,CM);
        
        TASBEConfig.set('OutputSettings.StemName', statName);
        [statisticsFile, histogramFile] = serializeBatchOutput(file_pairs, CM, AP, sampleresults);
        TASBEConfig.set('OutputSettings.StemName', stemName);
        
        if ~isdir(outputPath)
            sanitized_path = strrep(outputPath, '/', '&#47;');
            sanitized_path = strrep(sanitized_path, '\', '&#92;');
            sanitized_path = strrep(sanitized_path, ':', '&#58;');
            TASBESession.notify('OutputFig','MakeDirectory','Directory does not exist, attempting to create it: %s',sanitized_path);
            mkdir(outputPath);
        end

        save([outputPath outputName],'AP','bins','file_pairs','results','sampleresults');
    end
end

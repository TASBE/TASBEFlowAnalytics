% Function that runs transfer curve analysis given a template spreadsheet. An Excel
% object and optional Color Model are inputs
function transfercurve_analysis_excel(path, extractor, CM)
    % Reset and update TASBEConfig and get exp name
    extractor.TASBEConfig_updates();
    TASBEConfig.set('template.displayErrors', 1);
    experimentName = extractor.getExcelValue('experimentName', 'char');
    TASBEConfig.set('template.displayErrors', 0);
    
    % Determine the number of transfer curve analysis to run
    sh_num3 = extractor.getSheetNum('first_sampleColName_TC');
    sh_num2 = extractor.getSheetNum('first_sample_num');
    first_sample_row = extractor.getRowNum('first_sample_num');
    sample_num_col = extractor.getColNum('first_sample_num');
    col_names = {};
    col_nums = {};
    row_nums = {};
    num_runs = 0;
    for i=extractor.getRowNum('sampleColName_TC'):size(extractor.sheets{sh_num3},1)
        try
            value = extractor.getExcelValuePos(sh_num3, i, extractor.getColNum('sampleColName_TC'), 'char');
            col_names{end+1} = value;
            row_nums{end+1} = i;
            num_runs = num_runs + 1;
        catch
            try
                value = num2str(extractor.getExcelValuePos(sh_num3, i, extractor.getColNum('sampleColName_TC'), 'numeric'));
                if isempty(value)
                    break
                end
                col_names{end+1} = value;
                row_nums{end+1} = i;
                num_runs = num_runs + 1;
            catch 
                break
            end
        end
    end
    
    if num_runs > 4
        preference_row = num_runs + 5;
    else
        preference_row = 13;
    end

    % Load the color model
    if nargin < 3
        % Obtain the CM_name
        try
            coords = extractor.getExcelCoordinates('inputName_CM', 3);
            CM_name = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'char'); 
            [~,name,~] = fileparts(CM_name);
            CM_name = [name '.mat'];
        catch
            TASBESession.warn('transfercurve_analysis_excel', 'MissingPreference', 'Missing CM Filename in "Transfer Curve Analysis" sheet. Looking in "Calibration" sheet.');
            try
                CM_name = extractor.getExcelValue('outputName_CM', 'char');
                [~,name,~] = fileparts(CM_name);
                CM_name = [name '.mat'];
            catch
                TASBESession.warn('transfercurve_analysis_excel', 'MissingPreference', 'Missing Output Filename in "Calibration" sheet. Defaulting to exp name.');
                CM_name = [experimentName '-ColorModel.mat'];
            end
        end
        
        % Obtain the CM_path
        try
            coords = extractor.getExcelCoordinates('inputPath_CM', 3);
            CM_path = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'char');
            javaFileObj = javaObject("java.io.File", end_with_slash(CM_path));
            if javaFileObj.isAbsolute()
                CM_path = end_with_slash(CM_path);
            else
                CM_path = end_with_slash(fullfile(path, CM_path));
            end
        catch
            TASBESession.warn('transfercurve_analysis_excel', 'MissingPreference', 'Missing CM Filepath in "Transfer Curve Analysis" sheet. Looking in "Calibration" sheet.'); 
            try
                CM_path = extractor.getExcelValue('outputPath_CM', 'char');
                javaFileObj = javaObject("java.io.File", end_with_slash(CM_path));
                if javaFileObj.isAbsolute()
                    CM_path = end_with_slash(CM_path);
                else
                    CM_path = end_with_slash(fullfile(path, CM_path));
                end
            catch
                TASBESession.warn('transfercurve_analysis_excel', 'MissingPreference', 'Missing Output Filepath in "Calibration" sheet. Defaulting to template path.'); 
                CM_path = path;
            end
        end
        CM_file = [CM_path CM_name];

        try 
            load(CM_file);
        catch
            TASBESession.warn('transfercurve_analysis_excel', 'MissingPreference', 'Could not load CM file, creating a new one.');
            CM = make_color_model_excel(extractor);
        end
    end

    % Set TASBEConfigs and create variables needed to run transfer curve analysis
    % Analyze on a histogram of 10^[first] to 10^[third] ERF, with bins every 10^[second]
    try
        coords = extractor.getExcelCoordinates('binseq_min', 3);
        binseq_min = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
        coords = extractor.getExcelCoordinates('binseq_max', 3);
        binseq_max = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
        coords = extractor.getExcelCoordinates('binseq_pdecade', 3);
        binseq_pdecade = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
        bins = BinSequence(binseq_min, (1/binseq_pdecade), binseq_max, 'log_bins');
    catch
        bins = BinSequence();
    end

    % Designate which channels have which roles
    [channel_roles, ~] = getChannelRoles(CM, extractor);
    
    if isempty(channel_roles)
        TASBESession.warn('transfercurve_analysis_excel', 'MissingPreference', 'Missing constitutive, input, output in "Calibration" sheet');
        APs = {AnalysisParameters(bins,{})};
    else
        APs = {};
        for j=1:numel(channel_roles)
            APs{end+1} = AnalysisParameters(bins,channel_roles{j});
        end
    end
   
    % Obtain the necessary sample filenames and print names
    % Go through columns to find column number
    for i=1:numel(col_names)
        col_name = col_names{i};
        for j=sample_num_col:size(extractor.sheets{sh_num2},2)
            try 
                ref_header = extractor.getExcelValuePos(sh_num2, first_sample_row-1, j, 'char');
            catch
                try
                    ref_header = num2str(extractor.getExcelValuePos(sh_num2, first_sample_row-1, j, 'numeric'));
                    if isempty(ref_header)
                        continue
                    end
                catch 
                    continue
                end
            end
            if strcmp(col_name, ref_header)
                col_nums{end+1} = j;
                break
            end 
        end
    end
    if numel(col_nums) ~= numel(col_names)
        TASBESession.error('transfercurve_analysis_excel', 'InvalidColumnName', 'Not all Sample Column Names in "Transfer Curve Analysis" are matched with column names in "Samples".');
    end
    % Go though comparison groups and store values and column numbers in
    % cell array
    comp_groups = {};
    comp_group_names = {};
    first_group_col = extractor.getColNum('first_sampleColName_TC');
    outputNames = {};
    outputPaths = {};
    stemNames = {};
    plotPaths = {};
    inducer_names = {};
    
    for i=1:numel(row_nums)
        try
            col_name = extractor.getExcelValuePos(sh_num3, row_nums{i}, first_group_col, 'char');
            checkError = true;
            % Get column number of col_name
            for j=sample_num_col:size(extractor.sheets{sh_num2},2)
                try 
                    ref_header = extractor.getExcelValuePos(sh_num2, first_sample_row-1, j, 'char');
                catch
                    try
                        ref_header = num2str(extractor.getExcelValuePos(sh_num2, first_sample_row-1, j, 'numeric'));
                        if isempty(ref_header)
                            continue
                        end
                    catch 
                        continue
                    end
                end
                if strcmp(col_name, ref_header)
                    checkError = false;
                    comp_groups{end+1} = {j, extractor.getExcelValuePos(sh_num3, row_nums{i}, extractor.getColNum('first_sampleVal_TC'))};
                    comp_group_names{end+1} = col_name;
                    break
                end 
            end
            if checkError
                TASBESession.error('transfercurve_analysis_excel', 'InvalidColumnName', 'Sample column name, %s, under Comparison Groups in "Transfer Curve Analysis" does not match with any column name in "Samples".', col_name);
            end
        catch
            try
                col_name = num2str(extractor.getExcelValuePos(sh_num3, row_nums{i}, first_group_col, 'numeric'));
                if isempty(col_name)
                    % Add empty col_name
                    comp_groups{end+1} = {};
                else
                    checkError = true;
                    % Get column number of col_name
                    for j=sample_num_col:size(extractor.sheets{sh_num2},2)
                        try 
                            ref_header = extractor.getExcelValuePos(sh_num2, first_sample_row-1, j, 'char');
                        catch
                            try
                                ref_header = num2str(extractor.getExcelValuePos(sh_num2, first_sample_row-1, j, 'numeric'));
                                if isempty(ref_header)
                                    continue
                                end
                            catch 
                                continue
                            end
                        end
                        if strcmp(col_name, ref_header)
                            checkError = false;
                            comp_groups{end+1} = {j, extractor.getExcelValuePos(sh_num3, row_nums{i}, extractor.getColNum('first_sampleVal_TC'))};
                            comp_group_names{end+1} = col_name;
                            break
                        end 
                    end
                    if checkError
                        TASBESession.error('transfercurve_analysis_excel', 'InvalidColumnName', 'Sample column name, %s, under Comparison Groups in "Transfer Curve Analysis" does not match with any column name in "Samples".', col_name);
                    end
                end
            catch 
                % Add empty col_name
                comp_groups{end+1} = {};
            end
        end
        % Get unique preferences 
        % Obtain output name
        try 
            outputName = extractor.getExcelValuePos(sh_num3, row_nums{i}, extractor.getColNum('outputName_TC'), 'char');
            [~, name, ~] = fileparts(outputName);
            outputNames{end+1} = [name '.mat'];
        catch
            TASBESession.warn('transfercurve_analysis_excel', 'MissingPreference', 'Missing Output File Name for Transfer Curve Analysis %s in "Transfer Curve Analysis" sheet', num2str(i));
            if i > 1
                outputNames{end+1} = [experimentName '-TransAnalysis' num2str(i) '.mat'];
            else
                outputNames{end+1} = [experimentName '-TransAnalysis.mat'];
            end
        end
        
        % Obtain output path
        try
            outputPath = extractor.getExcelValuePos(sh_num3, row_nums{i}, extractor.getColNum('outputPath_TC'), 'char');
            javaFileObj = javaObject("java.io.File", end_with_slash(outputPath));
            if javaFileObj.isAbsolute()
                outputPath = end_with_slash(outputPath);
            else
                outputPath = end_with_slash(fullfile(path, outputPath));
            end
            outputPaths{end+1} = outputPath;
        catch
            TASBESession.warn('transfercurve_analysis_excel', 'MissingPreference', 'Missing Output File Path for Transfer Curve Analysis %s in "Transfer Curve" sheet', num2str(i));
            outputPaths{end+1} = path;
        end

        try 
            stemName_coord = extractor.getExcelCoordinates('OutputSettings.StemName');
            stemNames{end+1} = extractor.getExcelValuePos(sh_num3, row_nums{i}, stemName_coord{3}{3}, 'char');
        catch
            TASBESession.warn('transfercurve_analysis_excel', 'MissingPreference', 'Missing Stem Name for Transfer Curve Analysis %s in "Transfer Curve Analysis" sheet. Defaulting to Comparison Groups', num2str(i));
            try
                if isa(comp_groups{i}{2}, 'numeric')
                    stemNames{end+1} = [comp_group_names{i} '=' num2str(comp_groups{i}{2})];
                else
                    stemNames{end+1} = [comp_group_names{i} '=' comp_groups{i}{2}];
                end
            catch
                TASBESession.warn('transfercurve_analysis_excel', 'MissingPreference', 'Stem Name for Transfer Curve Analysis %s in "Transfer Curve Analysis" sheet defaulting to exp name.', num2str(i));
                if i > 1
                    stemNames{end+1} = [experimentName num2str(i)];
                else
                    stemNames{end+1} = experimentName;
                end
            end
        end
        
        try
            plotPath_coord = extractor.getExcelCoordinates('plots.plotPath');
            plot_path = extractor.getExcelValuePos(sh_num3, row_nums{i}, plotPath_coord{4}{3}, 'char');
            javaFileObj = javaObject("java.io.File", end_with_slash(plot_path));
            if javaFileObj.isAbsolute()
                plot_path = end_with_slash(plot_path);
            else
                plot_path = end_with_slash(fullfile(path, plot_path));
            end
            plotPaths{end+1} = plot_path;
        catch
            TASBESession.warn('transfercurve_analysis_excel', 'MissingPreference', 'Missing plot path for Transfer Curve Analysis %s in "Transfer Curve Analysis" sheet', num2str(i));            
            plot_path = end_with_slash(fullfile(path, 'plots/'));
            plotPaths{end+1} = plot_path;
        end
        
        inducer_names{end+1} = col_names{i};
    end
    
    for i=1:numel(col_names)
        sample_names = {};
        file_names = {};
        inducer_name = inducer_names{i};
        stemName = stemNames{i};
        TASBEConfig.set('plots.plotPath', plotPaths{i});
        outputName = outputNames{i};
        outputPath = outputPaths{i};
        % Go though sample rows of selected column and add to cell arrays
        for j=first_sample_row:extractor.getRowNum('last_sample_num')
            try
                value = extractor.getExcelValuePos(sh_num2, j, col_nums{i}, 'numeric');
                if isempty(comp_groups{i})
                    sample_names{end+1} = value;
                    file = getExcelFilename(extractor, j, path);
                    file_names{end+1} = file;
                elseif isa(extractor.getExcelValuePos(sh_num2, j, comp_groups{i}{1}), 'numeric') && (extractor.getExcelValuePos(sh_num2, j, comp_groups{i}{1}) == comp_groups{i}{2})
                    sample_names{end+1} = value;
                    file = getExcelFilename(extractor, j, path);
                    file_names{end+1} = file;
                elseif isa(extractor.getExcelValuePos(sh_num2, j, comp_groups{i}{1}), 'char') && strcmp(extractor.getExcelValuePos(sh_num2, j, comp_groups{i}{1}), comp_groups{i}{2})
                    sample_names{end+1} = value;
                    file = getExcelFilename(extractor, j, path);
                    file_names{end+1} = file;
                end
            catch
                continue
            end
        end

        % Make a map of condition names to file sets
        level_file_pairs = {};
        level_file_pairs(:,1) = sample_names;
        level_file_pairs(:,2) = file_names;
        
        % Go through all possible APs
        for j=1:numel(APs)
            AP = APs{j};
            % Ignore any bins with less than valid count as noise
            try
                coords = extractor.getExcelCoordinates('minValidCount', 3);
                minValidCount = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
                AP=setMinValidCount(AP,minValidCount);
            catch
                TASBESession.warn('transfercurve_analysis_excel', 'ImportantMissingPreference', 'Missing Min Valid Count in "Transfer Curve Analysis" sheet');
            end

            % Add autofluorescence back in after removing for compensation?
            try
                coords = extractor.getExcelCoordinates('autofluorescence', 3);
                autofluorescence = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
                AP=setUseAutoFluorescence(AP,autofluorescence);
            catch
                TASBESession.warn('transfercurve_analysis_excel', 'ImportantMissingPreference', 'Missing Use Auto Fluorescence in "Transfer Curve Analysis" sheet');
            end

            try
                coords = extractor.getExcelCoordinates('minFracActive', 3);
                minFracActive = extractor.getExcelValuePos(coords{1}, preference_row, coords{3}, 'numeric');
                AP=setMinFractionActive(AP,minFracActive);
            catch
                TASBESession.warn('transfercurve_analysis_excel', 'ImportantMissingPreference', 'Missing Min Fraction Active in "Transfer Curve Analysis" sheet');
            end
            
            if j > 1
                stemName = [stemName '-' num2str(j)];
                outputName_parts = strsplit(outputName, '.');
                outputName = [outputName_parts{1} '-' num2str(j) '.' outputName_parts{2}];
            end
            
            TASBEConfig.set('OutputSettings.StemName', stemName);
            experiment = Experiment(experimentName,{inducer_name}, level_file_pairs);
    
            % Execute the actual analysis
            fprintf('Starting analysis...\n');
            [results, sampleresults] = process_transfer_curve( CM, experiment, AP);

            % Plot how the constitutive fluorescence was distributed
            TASBEConfig.set('histogram.displayLegend',false);
            plot_bin_statistics(sampleresults, getInducerLevelsToFiles(experiment,1));

            % Plot the relation between inducer and input fluorescence
            TASBEConfig.set('OutputSettings.DeviceName',inducer_name);
            plot_inducer_characterization(results);

            % Plot the relation between input and output fluorescence
            plot_IO_characterization(results);

            % Save the results of computation
            save('-V7',[outputPath outputName],'experiment','AP','sampleresults','results');
        end
    end
end

% Function that runs batch analysis given a template spreadsheet. An Excel
% object and optional Color Model are inputs
function batch_analysis_excel(extractor, CM)
    % Reset and update TASBEConfig and obtain experiment name
    extractor.TASBEConfig_updates();
    experimentName = extractor.getExcelValue('experimentName', 'char');

    % Load the color model
    if nargin < 2
        try
            CM_file = extractor.getExcelValue('inputName_CM', 'char', 1); 
        catch
            try
                CM_file = extractor.getExcelValue('outputName_CM', 'char');
            catch
                CM_file = [experimentName '-ColorModel.mat'];
            end
        end

        try 
            load(CM_file);
        catch
            CM = make_color_model_excel();
        end
    end

    % Set TASBEConfigs and create variables needed to run batch analysis
    try
        outputName = extractor.getExcelValue('outputName_BA', 'char');
    catch
        TASBESession.warn('make_color_model_excel', 'MissingPreference', 'Missing Output File Name for Batch Analysis in "Samples" sheet');
        outputName = [experimentName '-BatchAnalysis.mat'];
    end

    try 
        stemName = extractor.getExcelValue('OutputSettings.StemName', 'char', 1);
        TASBEConfig.set('OutputSettings.StemName', stemName);
    catch
        TASBESession.warn('make_color_model', 'MissingPreference', 'Missing Stem Name in "Samples" sheet');
        TASBEConfig.set('OutputSettings.StemName', experimentName);
    end

    extractor.setTASBEConfig('plots.plotPath', 'char', 2);
    % Analyze on a histogram of 10^[first] to 10^[third] ERF, with bins every 10^[second]
    try
        binseq_min = extractor.getExcelValue('binseq_min', 'numeric', 1);
        binseq_max = extractor.getExcelValue('binseq_max', 'numeric', 1);
        binseq_pdecade = extractor.getExcelValue('binseq_pdecade', 'numeric', 1);
        bins = BinSequence(binseq_min, (1/binseq_pdecade), binseq_max, 'log_bins');
    catch
        bins = BinSequence();
    end

    % Designate which channels have which roles
    ref_channels = {'constitutive', 'input', 'output'};
    outputs = {};
    print_names = {};
    sh_num1 = extractor.getSheetNum('first_flchrome_name');
    first_flchrome_row = extractor.getRowNum('first_flchrome_name');
    flchrome_name_col = extractor.getColNum('first_flchrome_name');
    flchrome_type_col = extractor.getColNum('first_flchrome_type');
    for i=first_flchrome_row:size(extractor.sheets{sh_num1},1)
        try
            print_name = extractor.getExcelValuePos(sh_num1, i, flchrome_name_col, 'char');
        catch
            break
        end
        print_names{end+1} = print_name;
        try
            channel_type = extractor.getExcelValuePos(sh_num1, i, flchrome_type_col, 'char');
        catch
            continue
        end
        for j=1:numel(ref_channels)
            if strcmpi(ref_channels{j}, channel_type)
                outputs{j} = channel_named(CM, print_name);
            end
        end
    end

    if numel(outputs) == 3
        AP = AnalysisParameters(bins,{'input',outputs{2}; 'output',outputs{3}; 'constitutive',outputs{1}});
    else
        TASBESession.warn('make_color_model', 'ImportantMissingPreference', 'Missing constitutive, input, output in "Calibration" sheet');
        AP = AnalysisParameters(bins,{});
    end

    % Ignore any bins with less than valid count as noise
    try
        minValidCount = extractor.getExcelValue('minValidCount', 'numeric', 1);
        AP=setMinValidCount(AP,minValidCount);
    catch
        TASBESession.warn('make_color_model', 'ImportantMissingPreference', 'Missing Min Valid Count in "Samples" sheet');
    end

    % Add autofluorescence back in after removing for compensation?
    try
        autofluorescence = extractor.getExcelValue('autofluorescence', 'numeric', 1);
        AP=setUseAutoFluorescence(AP,autofluorescence);
    catch
        TASBESession.warn('make_color_model', 'ImportantMissingPreference', 'Missing Use Auto Fluorescence in "Samples" sheet');
    end

    try
        minFracActive = extractor.getExcelValue('minFracActive', 'numeric', 1);
        AP=setMinFractionActive(AP,minFracActive);
    catch
        TASBESession.warn('make_color_model', 'ImportantMissingPreference', 'Missing Min Fraction Active in "Samples" sheet');
    end
    
    % Obtain the necessary sample filenames and print names
    sample_names = {};
    file_names = {};
    sh_num2 = extractor.getSheetNum('first_sample_num');
    first_sample_row = extractor.getRowNum('first_sample_num');
    sample_num_col = extractor.getColNum('first_sample_num');
    sample_name_col = extractor.getColNum('first_sample_name');
    sample_exclude_col = extractor.getColNum('first_sample_exclude');
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
            file = getFilename(extractor, i);
            file_names{end+1} = file;
        end
    end

    % Make a map of condition names to file sets
    file_pairs = {};
    file_pairs(:,1) = sample_names;
    file_pairs(:,2) = file_names;

    [results, sampleresults] = per_color_constitutive_analysis(CM,file_pairs,print_names,AP);

    % Make output plots
    plot_batch_histograms(results,sampleresults,CM);

    [statisticsFile, histogramFile] = serializeBatchOutput(file_pairs, CM, AP, sampleresults);

    save(outputName,'AP','bins','file_pairs','results','sampleresults');
end

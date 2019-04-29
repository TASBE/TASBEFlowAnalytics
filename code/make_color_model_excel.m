% Function that uses a template spreadsheet to create and save a Color Model. 
% Input includes a TemplateExtraction object.  
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics

function [CM] = make_color_model_excel(extractor)
    % Reset and update TASBEConfig 
    extractor.TASBEConfig_updates();
    path = extractor.path;
    % Set TASBEConfigs and create variables needed to generate the CM
    TASBEConfig.set('template.displayErrors', 1);
    experimentName = extractor.getExcelValue('experimentName', 'char');
    TASBEConfig.set('template.displayErrors', 0);
    try
        outputName = extractor.getExcelValue('outputName_CM', 'char');
        [~, name, ~] = fileparts(outputName);
        outputName = [name '.mat'];
    catch
        TASBESession.warn('make_color_model_excel', 'MissingPreference', 'Missing Output File Name in "Calibration" sheet');
        outputName = [experimentName '-ColorModel.mat'];
    end
    
    try
        outputPath = extractor.getExcelValue('outputPath_CM', 'char');
        outputPath = make_filename_absolute(outputPath, path);
    catch
        TASBESession.warn('make_color_model_excel', 'MissingPreference', 'Missing Output File Path in "Calibration" sheet');
        outputPath = path;
    end
    
    try
        transChannelMin = extractor.getExcelValue('transChannelMin', 'cell');
    catch
        TASBESession.warn('make_color_model_excel', 'MissingPreference', 'Missing Translation Channel Min in "Calibration" sheet');
        transChannelMin = {};
    end
    
    try
        plot_path = extractor.getExcelValue('plots.plotPath', 'char', 1);
        display(plot_path);
        plot_path = make_filename_absolute(plot_path, path);
        TASBEConfig.set('plots.plotPath', plot_path);
    catch
        TASBESession.warn('make_color_model_excel', 'MissingPreference', 'Missing plot path in "Calibration" sheet');
        plot_path = end_with_slash(fullfile(path, 'plots/'));
        TASBEConfig.set('plots.plotPath', plot_path);
    end
    
    extractor.setTASBEConfig('beads.beadModel', 'char');
    extractor.setTASBEConfig('beads.beadBatch', 'char');
    extractor.setTASBEConfig('beads.rangeMin', 'numeric');
    extractor.setTASBEConfig('beads.rangeMax', 'numeric');
    extractor.setTASBEConfig('beads.peakThreshold', 'numeric');
    extractor.setTASBEConfig('beads.beadChannel', 'char');
    extractor.setTASBEConfig('beads.secondaryBeadChannel', 'char');
    
    % Extract bead, blank, and all files (and size bead file if applicable)
    % ref_filenames = {'blank','beads','all','sizebead'};
    sizebeadfile = [];
    TASBEConfig.set('template.displayErrors', 1);
    ref_filenames = {extractor.getExcelValue('blank_name', 'char'), extractor.getExcelValue('all_name', 'char')};
    TASBEConfig.set('template.displayErrors', 0);
    size_bead = false;
    try
        size_bead_name = extractor.getExcelValue('size_bead_name', 'char');
        ref_filenames{end+1} = size_bead_name;
        size_bead = true;
    catch
        TASBESession.notify('make_color_model_excel', 'NoSizeBeads', 'Size bead feature not being used.');
    end
    output_filenames = {};
    sh_num1 = extractor.getSheetNum('first_sample_num');
    first_sample_row = extractor.getRowNum('first_sample_num');
    sample_num_col = extractor.getColNum('first_sample_num');
    sample_name_col = find(ismember(extractor.col_names, 'SAMPLE NAME'), 1);
    if isempty(sample_name_col)
        TASBESession.error('make_color_model_excel', 'InvalidHeaderName', 'The header, SAMPLE NAME, does not match with any column titles in "Samples" sheet.');
    end
    % Go through samples in "Samples" sheet and look for matches in name to
    % elements in ref_filenames
    for i=first_sample_row:size(extractor.sheets{sh_num1},1)
        try
            num = extractor.getExcelValuePos(sh_num1, i, sample_num_col, 'numeric');
            name = extractor.getExcelValuePos(sh_num1, i, sample_name_col, 'char');
            if isempty(num)
                break
            end
        catch
            break
        end
        for j=1:numel(ref_filenames)
            if strcmpi(name, ref_filenames{j})
                datafile = getExcelFilename(extractor, i);
                output_filenames{j} = datafile{1};
            end
        end
    end
    bead_files = getBeadFile(extractor);
    beads_file = bead_files{1};
    blank_file = output_filenames{1};
    all_file = output_filenames{2};
    if size_bead
        sizebeadfile = output_filenames{3};
    end
    
    % Autodetect gating with an N-dimensional gaussian-mixture-model
    autogate = GMMGating(blank_file);
    
    % Dealing with channels 
    print_names = {};
    channel_names = {};
    sample_ids = {};
    all_channels = false;
    % Determine relevant channels
    try 
        rel_channels = extractor.getExcelValue('relevant_channels', 'cell');
    catch
        % this means that all channels are relevant
        all_channels = true;
    end
    % Create one channel for each color
    channels = {}; 
    sh_num2 = extractor.getSheetNum('first_flchrome_name');
    first_flchrome_row = extractor.getRowNum('first_flchrome_name');
    flchrome_name_col = extractor.getColNum('first_flchrome_name');
    flchrome_channel_col = extractor.getColNum('first_flchrome_channel');
    flchrome_wavlen_col = extractor.getColNum('first_flchrome_wavlen');
    flchrome_filter_col = extractor.getColNum('first_flchrome_filter');
    flchrome_color_col = extractor.getColNum('first_flchrome_color');
    flchrome_id_col = extractor.getColNum('first_flchrome_id');
    for i=first_flchrome_row:size(extractor.sheets{sh_num2},1)
        try
            print_name = extractor.getExcelValuePos(sh_num2, i, flchrome_name_col, 'char');
        catch
            break
        end
        ind = i - (first_flchrome_row-1);
        if all_channels || ~isempty(find([rel_channels{:}] == ind, 1))
            % Extract the rest of the information for the channel
            channel_name = extractor.getExcelValuePos(sh_num2, i, flchrome_channel_col, 'char');
            excit_wavelen = extractor.getExcelValuePos(sh_num2, i, flchrome_wavlen_col, 'numeric');
            filter = strtrim(strsplit(extractor.getExcelValuePos(sh_num2, i, flchrome_filter_col, 'char'), '/'));
            color = extractor.getExcelValuePos(sh_num2, i, flchrome_color_col, 'char');
            print_names{ind} = print_name;
            channel_names{ind} = channel_name;
            channels{ind} = Channel(channel_name, excit_wavelen, str2double(filter{1}), str2double(filter{2}));
            channels{ind} = setPrintName(channels{ind}, print_name); % Name to print on charts
            channels{ind} = setLineSpec(channels{ind}, color); % Color for lines, when needed
            try
                id = extractor.getExcelValuePos(sh_num2, i, flchrome_id_col, 'char');
            catch
                id = print_name;
            end
            sample_ids{ind} = id;
        end
    end
    
    % Create all non-fluorescence channels
    num_nonflr = 0;
    first_nonflr_row = extractor.getRowNum('first_nonflr_name');
    nonflr_name_col = extractor.getColNum('first_nonflr_name');
    nonflr_channel_col = extractor.getColNum('first_nonflr_channel');
    nonflr_wavlen_col = extractor.getColNum('first_nonflr_wavlen');
    nonflr_filter_col = extractor.getColNum('first_nonflr_filter');
    nonflr_color_col = extractor.getColNum('first_nonflr_color');
    for i=first_nonflr_row:size(extractor.sheets{sh_num2},1)
        try
            print_name = extractor.getExcelValuePos(sh_num2, i, nonflr_name_col, 'char');
        catch
            break
        end
        % Extract the rest of the information for the channel
        channel_name = extractor.getExcelValuePos(sh_num2, i, nonflr_channel_col, 'char');
        excit_wavelen = extractor.getExcelValuePos(sh_num2, i, nonflr_wavlen_col, 'numeric');
        filter = strtrim(strsplit(extractor.getExcelValuePos(sh_num2, i, nonflr_filter_col, 'char'), '/'));
        color = extractor.getExcelValuePos(sh_num2, i, nonflr_color_col, 'char');
        % FSC and SSC channels can be added to be read unprocessed to MEFL
        channels{end+1} = Channel(channel_name, excit_wavelen, str2double(filter{1}), str2double(filter{2}));
        channels{end} = setPrintName(channels{end}, print_name); % Name to print on charts
        channels{end} = setLineSpec(channels{end}, color); % Color for lines, when needed
        % If the name is FSC or SSC (or one of those with '-A', '-H', or '-W') it will automatically be unprocessed; otherwise, set it 
        channels{end} = setIsUnprocessed(channels{end}, true);
        num_nonflr = num_nonflr + 1;
    end
    
    % Make sure non-fluorescence channel made if size_bead is true
    if size_bead && num_nonflr == 0
        TASBESession.warn('make_color_model_excel', 'NoSizeBeadChannel', 'Size bead feature requires at least one non-fluorescence channel.');
    end
    
    % Obtain channel filenames using sample_ids
    colorfiles = {};
    for i=first_sample_row:size(extractor.sheets{sh_num1},1)
        try
            num = extractor.getExcelValuePos(sh_num1, i, sample_num_col, 'numeric');
            name = extractor.getExcelValuePos(sh_num1, i, sample_name_col, 'char');
            if isempty(num)
                break
            end
        catch
            break
        end
        for j=1:numel(sample_ids)
            if strcmpi(name, sample_ids{j})
                datafile = getExcelFilename(extractor, i);
                colorfiles{j} = datafile{1};
            end
        end
    end
    
    % Included a check to make sure that the number of channels matches with
    % the number in the template
    if numel(channels) ~= extractor.getExcelValue('num_channels', 'numeric')
        TASBESession.warn('make_color_model_excel', 'DimensionMismatch', 'Channel dimensions do not match with number of channels from template');
    end

    % Multi-color controls are used for converting other colors into FITC units
    colorpairfiles = {};
    % Entries are: channel1, channel2, constitutive channel, filename
    % This allows channel1 and channel2 to be converted into one another.
    % If you only have two colors, you can set consitutive-channel to equal channel1 or channel2
    if size_bead
        n_channels = numel(channels) - num_nonflr;
    else
        n_channels = numel(channels);
    end
    
    if n_channels == 2
        colorpairfiles{1} = {channels{1}, channels{2}, channels{2}, all_file};
    else
        for i=1:n_channels-1
            if i >= 2
                colorpairfiles{end+1} = {channels{1}, channels{i+1}, channels{i-1}, all_file};
            else
                colorpairfiles{end+1} = {channels{1}, channels{i+1}, channels{i+2}, all_file};
            end
        end
    end

    % Making the color model
    CM = ColorModel(beads_file, blank_file, channels, colorfiles, colorpairfiles, sizebeadfile);
    CM = set_ERF_channel_name(CM, channel_names{1});
    CM = add_prefilter(CM,autogate);
    
    if ~isempty(transChannelMin) 
        TASBEConfig.set('colortranslation.channelMinimum',cell2mat(transChannelMin));
    end
    
    if size_bead
        % Setting size bead configs
        CM=set_um_channel_name(CM, 'FSC-A');
        extractor.setTASBEConfig('sizebeads.beadModel', 'char');
        extractor.setTASBEConfig('sizebeads.beadBatch', 'char');
        extractor.setTASBEConfig('sizebeads.rangeMin', 'numeric');
        extractor.setTASBEConfig('sizebeads.rangeMax', 'numeric');
        extractor.setTASBEConfig('sizebeads.peakThreshold', 'numeric');
        extractor.setTASBEConfig('sizebeads.beadChannel', 'char');
    end
    
    % Execute and save the model
    CM = resolve(CM);
    if ~isdir(outputPath)
        sanitized_path = strrep(outputPath, '/', '&#47;');
        sanitized_path = strrep(sanitized_path, '\', '&#92;');
        sanitized_path = strrep(sanitized_path, ':', '&#58;');
        TASBESession.notify('OutputFig','MakeDirectory','Directory does not exist, attempting to create it: %s',sanitized_path);
        mkdir(outputPath);
    end
    save('-V7',[outputPath, outputName],'CM');
    
    % Conduct bead comparisons if applicable
    % If size of beadfiles is greater than 1, then run bead comparisons
    if numel(bead_files) > 1
        try
            tolerance = extractor.getExcelValue('bead_tolerance', 'numeric');
        catch
            tolerance = 0.5;
            TASBESession.notify('make_color_model_excel', 'MissingPreference', 'Bead comparison tolerance defaulting to 0.5');
        end
        
        for i=2:numel(bead_files)
            [ok, ratios] = check_beads_identical(CM, bead_files{i}, tolerance);
            TASBESession.notify('make_color_model_excel', 'BeadCompRatio', 'Ratios between peaks: %.2f \n', ratios);
            % A warning is thrown if two beadfiles are not considered identical
            if ok
                TASBESession.succeed('make_color_model_excel', 'Identical', 'Bead files are sufficiently identical to each other');
            else
                TASBESession.warn('make_color_model_excel', 'NotIdentical', 'Bead files are not identical to each other!');
            end
        end
    end
end

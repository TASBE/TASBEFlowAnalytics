% Excel class with objects representing the information for a given template
% spreadsheet
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
classdef TemplateExtraction
    properties
        % Coordinates is a cell array of different template variable names
        % and their coordinates in the form of {sheet number, row, col}. (A
        % few variables have multiple coordinates.) 
        coordinates;
        % Sheets is an array of raw data from four main sheets in template
        sheets;
        % Save the inputted template filepath and filename
        file;
        % The path of the template used to for relative filepath inputs
        path;
        % Cell array of sample column names 
        col_names;
    end
    methods
        % Constuctor with filepath of template and optional coordinates
        % property as inputs
        function obj = TemplateExtraction(file, path, coords)
            obj.file = file;
            if exist('path', 'var')
                obj.path = path;
            else
                [obj.path, ~, ~] = fileparts(file);
            end
            % Read in Excel for information, Experiment sheet
            [~,~,s1] = xlsread(file, 'Experiment');
            % Read in Excel for information, Samples sheet
            [~,~,s2] = xlsread(file, 'Samples');
            % Read in Excel for information, Calibration sheet
            [~,~,s3] = xlsread(file, 'Calibration');
            % Read in Excel for information, Comparative Analysis sheet
            try
                [~,~,s4] = xlsread(file, 'Comparative Analysis');
            catch
                TASBESession.warn('TemplateExtraction', 'MissingAnalysisSheet', 'Comparative Analysis sheet is missing. Add if want to run comparative analysis.');
                s4 = {};
            end
            % Read in Excel for information, Transfer Curve Analysis sheet
            try
                [~,~,s5] = xlsread(file, 'Transfer Curve Analysis');
            catch
                TASBESession.warn('TemplateExtraction', 'MissingAnalysisSheet', 'Transfer Curve Analysis sheet is missing. Add if want to run transfer curve analysis.');
                s5 = {};
            end
            % Read in Excel for information, Additional Settings sheet
            [~,~,s6] = xlsread(file, 'Optional Settings');
            obj.sheets = {s1, s2, s3, s4, s5, s6};
            % Find Calibration rows
            rowCoordinates = obj.findCalibrationRows(3);
            if nargin < 3
                obj.coordinates = {...
                    % Coords for variables in "Experiment"
                    {'experimentName'; {1, 5, 1}}; % uses all three vals 
                    {'first_filename_template'; {1, 16, 5}}; % uses sh & col vals
                    {'first_condition_key'; {1, 16, 1}}; % uses sh & col vals
                    % Coords for variables in "Samples"
                    {'plots.plotPath'; {{3, rowCoordinates{7}, 2}, {2, 30, 2}, {4, 5, 15}, {5, 5, 16}}}; % uses sh & col vals
                    {'first_sample_num'; {2, 5, 1}}; % uses all three vals 
                    {'inputName_CM'; {{2, 30, 3}, {4, 13, 2}, {5, 13, 2}}}; % uses sh & col vals
                    {'inputPath_CM'; {{2, 30, 4}, {4, 13, 3}, {5, 13, 3}}}; % uses sh & col vals
                    {'OutputSettings.StemName'; {{2, 30, 5}, {4, 5, 13}, {5, 5, 13}}}; % uses sh & col vals
                    {'binseq_min'; {{2, 30, 16}, {4, 13, 8}, {5, 13, 8}}}; % uses sh & col vals
                    {'binseq_pdecade'; {{2, 30, 17}, {4, 13, 9}, {5, 13, 9}}}; % uses sh & col vals
                    {'binseq_max'; {{2, 30, 18}, {4, 13, 10}, {5, 13, 10}}}; % uses sh & col vals
                    {'minValidCount'; {{2, 30, 12}, {4, 13, 4}, {5, 13, 4}}}; % uses sh & col vals
                    {'minValidau'; {{2, 30, 13}, {4, 13, 5}, {5, 13, 5}}}; % uses sh & col vals
                    {'autofluorescence'; {{2, 30, 14}, {4, 13, 6}, {5, 13, 6}}}; % uses sh & col vals
                    {'minFracActive'; {{2, 30, 15}, {4, 13, 7}, {5, 13, 7}}}; % uses sh & col vals
                    {'outputName_BA'; {2, 30, 6}}; % uses sh & col vals
                    {'outputPath_BA'; {2, 30, 7}}; % uses sh & col vals
                    {'statName_BA'; {2, 30, 8}}; % uses sh & col vals
                    {'statPath_BA'; {2, 30, 9}}; % uses sh & col vals
                    {'cloudName_BA'; {2, 30, 10}}; % uses sh & col vals
                    {'cloudPath_BA'; {2, 30, 11}}; % uses sh & col vals
                    % Coords for variables in "Calibration"
                    {'beads.beadModel'; {3, rowCoordinates{1}, 2}}; % uses all three vals 
                    {'beads.beadBatch'; {3, rowCoordinates{1}, 1}}; % uses all three vals 
                    {'beads.rangeMin'; {3, rowCoordinates{1}, 3}}; % uses all three vals 
                    {'beads.rangeMax'; {3, rowCoordinates{1}, 4}}; % uses all three vals 
                    {'beads.peakThreshold'; {3, rowCoordinates{1}, 5}}; % uses all three vals 
                    {'beads.beadChannel'; {3, rowCoordinates{1}, 6}}; % uses all three vals 
                    {'bead_name'; {3, rowCoordinates{1}, 7}}; % uses all three vals 
                    {'bead_tolerance'; {3, rowCoordinates{1}, 8}}; % uses all three vals 
                    {'blank_name'; {3, rowCoordinates{2}, 2}}; % uses all three vals
                    {'first_flchrome_name'; {3, rowCoordinates{3}, 2}}; % uses all three vals 
                    {'first_flchrome_channel'; {3, rowCoordinates{3}, 3}}; % uses all three vals 
                    {'first_flchrome_type'; {3, rowCoordinates{3}, 4}}; % uses all three vals, whether constitutive or input or output
                    {'first_flchrome_wavlen'; {3, rowCoordinates{3}, 5}}; % uses all three vals 
                    {'first_flchrome_filter'; {3, rowCoordinates{3}, 6}}; % uses all three vals 
                    {'first_flchrome_color'; {3, rowCoordinates{3}, 7}}; % uses all three vals 
                    {'first_flchrome_id'; {3, rowCoordinates{3}, 8}}; % uses all three vals
                    {'first_nonflr_name'; {3, rowCoordinates{4}, 2}}; % uses all three vals 
                    {'first_nonflr_channel'; {3, rowCoordinates{4}, 3}}; % uses all three vals 
                    {'first_nonflr_wavlen'; {3, rowCoordinates{4}, 4}}; % uses all three vals 
                    {'first_nonflr_filter'; {3, rowCoordinates{4}, 5}}; % uses all three vals 
                    {'first_nonflr_color'; {3, rowCoordinates{4}, 6}}; % uses all three vals
                    {'relevant_channels'; {3, rowCoordinates{5}, 2}}; % uses all three vals 
                    {'transChannelMin'; {3, rowCoordinates{5}, 3}}; % uses all three vals
                    {'num_channels'; {3, rowCoordinates{5}, 1}}; % uses all three vals  
                    {'all_name'; {3, rowCoordinates{5}, 4}}; % uses all three vals 
                    {'size_bead_name'; {3, rowCoordinates{6}, 7}}; % uses all three vals
                    {'sizebeads.beadModel'; {3, rowCoordinates{6}, 2}}; % uses all three vals
                    {'sizebeads.rangeMin'; {3, rowCoordinates{6}, 3}}; % uses all three vals
                    {'sizebeads.rangeMax'; {3, rowCoordinates{6}, 4}}; % uses all three vals
                    {'sizebeads.peakThreshold'; {3, rowCoordinates{6}, 5}}; % uses all three vals
                    {'sizebeads.beadChannel'; {3, rowCoordinates{6}, 6}}; % uses all three vals
                    {'sizebeads.beadBatch'; {3, rowCoordinates{6}, 1}}; % uses all three vals
                    {'beads.secondaryBeadChannel'; {3, rowCoordinates{7}, 3}}; % uses all three vals 
                    {'outputName_CM'; {3, rowCoordinates{7}, 4}}; % uses all three vals 
                    {'outputPath_CM'; {3, rowCoordinates{7}, 5}}; % uses all three vals    
                    % Coords for variables in "Comparative Analysis"
                    {'outputName_PM'; {4, 5, 13}}; % uses sh & col vals
                    {'outputPath_PM'; {4, 5, 14}}; % uses sh & col vals
                    {'primary_sampleColName_PM'; {4, 5, 7}}; % uses sh & col vals
                    {'secondary_sampleColName_PM'; {4, 5, 10}}; % uses sh & col vals
                    {'first_compGroup_PM'; {4, 5, 1}}; % uses sh & col vals
                    % Coords for variables in "Transfer Curve Analysis"
                    {'outputName_TC'; {5, 5, 14}}; % uses sh & col vals
                    {'outputPath_TC'; {5, 5, 15}}; % uses sh & col vals
                    {'sampleColName_TC'; {5, 5, 7}}; % uses sh & col vals
                    {'first_compGroup_TC'; {5, 5, 1}}; % uses sh & col vals
                    % Coords for variables in "Optional Settings"
                    {'first_preference_name'; {6, 3, 1}}; % uses all three vals 
                    {'first_preference_value'; {6, 3, 3}}; % uses sh & col vals
                    };
            else
                obj.coordinates = coords;
            end
            % Find the number of templates and update coordinates with info
            obj.coordinates = obj.findTemplates();
            % Find sample column names and save into col_names property
            obj.col_names = obj.findSampleCols();
            % Find last sample row
            obj.coordinates = obj.findLastSampleRow();
            % Check on condition keys
            obj.checkConditions();
        end

        function TASBEConfig_updates(obj)
        % Update any relevant TASBEConfig preferences from the Optional Settings sheet in the
        % template spreadsheet
            TASBEConfig.checkpoint(TASBEConfig.checkpoints());
            raw = obj.sheets{obj.getSheetNum('first_preference_name')};
            name_col = obj.getColNum('first_preference_name');
            val_col = obj.getColNum('first_preference_value');
            for i=obj.getRowNum('first_preference_name'):size(raw,1)
                % Set the TASBEConfig if value column not empty for given
                % row
                if ~isnan(cell2mat(raw(i,val_col)))
                    if ~isempty(strfind(char(cell2mat(raw(i,val_col))), ','))
                        bounds = strtrim(strsplit(char(cell2mat(raw(i,val_col))), ','));
                        if ~isempty(strfind(char(cell2mat(raw(i,name_col))), 'Size'))
                            TASBEConfig.set(char(cell2mat(raw(i,name_col))), [str2double(bounds{1}), str2double(bounds{2})]);
                        else
                            TASBEConfig.set(char(cell2mat(raw(i,name_col))), bounds);
                        end
                    elseif ~isempty(strfind(char(cell2mat(raw(i,val_col))), '^'))
                        val = str2num(char(cell2mat(raw(i,val_col))));
                        TASBEConfig.set(char(cell2mat(raw(i,name_col))), val);
                    else
                        TASBEConfig.set(char(cell2mat(raw(i,name_col))), cell2mat(raw(i,val_col)));
                    end
                end
            end
        end
        
        function new_value = sanitizeFromExcel(obj, value)
        % Cleans up the Excel cell values by calling trim function is value
        % is a string
            try 
                new_value = strtrim(value);
            catch
                new_value = value;
            end    
        end
        
        function position = getExcelCoordinates(obj, name, index)
        % Returns the ExcelCoordinates stored within obj.coordinates with
        % name of variable as input and optional index if multiple
        % coordinates are stored for a single variable
            for i=1:numel(obj.coordinates)
                if strcmp(name, obj.coordinates{i}{1})
                    position = obj.coordinates{i}{2};
                    if exist('index', 'var')
                        position = position{index};
                    end
                    return
                end
            end
            TASBESession.error('TemplateExtraction','CoordNotFound','Inputted name, %s, not valid. No match found in coordinates.', name);
        end
        
        function sheet_num = getSheetNum(obj, name, index)
        % Returns the first coordinate value (sheet number) of a given variable
            if exist('index', 'var')
                pos = obj.getExcelCoordinates(name, index);
            else
                pos = obj.getExcelCoordinates(name);
            end
            sheet_num = pos{1};
        end
        
        function row = getRowNum(obj, name, index)
        % Returns the second coordinate value (row number) of a given variable
            if exist('index', 'var')
                pos = obj.getExcelCoordinates(name, index);
            else
                pos = obj.getExcelCoordinates(name);
            end
            row = pos{2};
        end
        
        function col = getColNum(obj, name, index)
        % Returns the third coordinate value (col number) of a given variable
            if exist('index', 'var')
                pos = obj.getExcelCoordinates(name, index);
            else
                pos = obj.getExcelCoordinates(name);
            end
            col = pos{3};
        end
        
        function new_obj = setExcelCoordinates(obj, name, coords)
        % Returns a new Excel object with updated coordinates. Takes the
        % name of the variable to change and new coords as inputs.
            new_coords = obj.coordinates;
            for i=1:numel(obj.coordinates)
                if strcmp(name, obj.coordinates{i}{1})
                    new_coords{i}{2} = coords;
                    new_obj = TemplateExtraction(obj.file, obj.path, new_coords);
                    return
                end
            end
            TASBESession.error('TemplateExtraction','CoordNotFound','Inputted name, %s, not valid. No match found in coordinates.', name);
        end
        
        function new_coords = addExcelCoordinates(obj, names, coords)
        % Returns a new obj.coordinates with updated coordinates. Takes the
        % name of the variable to add and its coords as inputs.
            new_coords = obj.coordinates;
            for i=1:numel(names)
                new_coords{end+1} = {names{i}; coords{i}};
            end
        end
        
        function new_coords = findTemplates(obj)
        % Finds the coordinates for all of the filename templates and adds
        % to coordinates property
            template_col = obj.getColNum('first_filename_template');
            template_sh = obj.getSheetNum('first_filename_template');
            first_template_row = 1;
            coords = {};
            for i=first_template_row:size(obj.sheets{template_sh}, 1)
                try
                    value = obj.getExcelValuePos(template_sh, i, template_col, 'char');
                    if ~isempty(strfind(value, 'Filename Template'))
                        coords{end+1} = {template_sh, i+1, template_col};
                    end
                catch
                    continue
                end
            end
            new_coords = obj.addExcelCoordinates({'filename_templates'}, {coords});
        end
        
        function row_coords = findCalibrationRows(obj, sheet_num)
        % Finds row numbers of Calibration parameters 
            row_headers = {'Rainbow Beads', 'Autofluorescence', 'Fluorochromes', 'Non-fluorescence Channels (unprocessed)', 'Color Translation', 'Size Beads', 'Additional Preferences for Color Model'};
            row_coords = {0};
            for i=1:size(obj.sheets{sheet_num}, 1)
                try
                    value = obj.getExcelValuePos(sheet_num, i, 1, 'char');
                    if ~isempty(find(strcmp(row_headers, value)))
                        row_coords{find(strcmp(row_headers, value))} = i+3;
                    end
                catch
                    continue
                end
            end
            tf = cellfun('isempty', row_coords);
            row_coords(tf) = {0};
        end
        
        function checkConditions(obj)
        % Looks through the condition keys in "Experiment" and raises a
        % warning if a value in "Samples" has a value that is not in the key
            condition_col = obj.getColNum('first_condition_key');
            condition_sh = obj.getSheetNum('first_condition_key');
            first_condition_row = 1;
            for i=first_condition_row:size(obj.sheets{condition_sh}, 1)
                try
                    value = obj.getExcelValuePos(condition_sh, i, condition_col, 'char');
                    if ~isempty(strfind(value, 'Sample Column Name'))
                        % Found a condition key, save its keys in an array
                        % and compare with column in "Samples"
                        try
                            column_name = obj.getExcelValuePos(condition_sh, i, condition_col+1, 'char');
                            obj.checkConditions_helper(i, column_name);
                            
                        catch
                            try
                                column_name = num2str(obj.getExcelValuePos(condition_sh, i, condition_col+1, 'numeric'));
                                if isempty(column_name)
                                    continue
                                end
                                obj.checkConditions_helper(i, column_name);
                            catch
                                continue
                            end
                        end
                    end
                catch
                    continue
                end
            end
        end
        
        function checkConditions_helper(obj, i, column_name)
        % Helper function for checkConditions
            condition_col = obj.getColNum('first_condition_key');
            condition_sh = obj.getSheetNum('first_condition_key');
            sh_num1 = obj.getSheetNum('first_sample_num');
            sample_start_col = obj.getColNum('first_sample_num');
            sample_start_row = obj.getRowNum('first_sample_num') - 1;
            keys = {};
            for j=i+2:size(obj.sheets{condition_sh}, 1)
                try
                    keys{end+1} = obj.getExcelValuePos(condition_sh, j, condition_col, 'char');
                catch
                    try
                        key = num2str(obj.getExcelValuePos(condition_sh, j, condition_col, 'numeric'));
                        if isempty(key)
                            break
                        end
                        keys{end+1} = key;
                    catch
                        break
                    end
                end
            end
            % look at column in "Samples"
            % Find the matching section name in filename template
            % and column header in "Samples"
            pos = find(ismember(obj.col_names, column_name), 1);
            if isempty(pos)
                TASBESession.error('TemplateExtraction', 'InvalidHeaderName', 'The header, %s, does not match with any column titles in "Samples" sheet.', column_name);
            else
                % go through all the rows and make sure
                % matches with one of the keys, raise
                % warning if not
                for k=sample_start_row+1:obj.getRowNum('last_sample_num')
                    try
                        value = obj.getExcelValuePos(sh_num1, k, pos, 'char');
                        ind = find(ismember(keys, value), 1);
                        if isempty(ind)
                            TASBESession.warn('TemplateExtraction', 'InvalidValue', 'The value of %s at row %s col %s does not match with listed keys.', value, num2str(k), column_name);
                        end
                    catch
                        try
                            value = num2str(obj.getExcelValuePos(sh_num1, k, pos, 'numeric'));
                            if isempty(value)
                                continue
                            end
                            ind = find(ismember(keys, value), 1);
                            if isempty(ind)
                                TASBESession.warn('TemplateExtraction', 'InvalidValue', 'The value of %s at row %s col %s does not match with listed keys.', value, num2str(k), column_name);
                            end
                        catch
                            continue
                        end
                    end
                end
            end
        end
        
        function col_names = findSampleCols(obj)
        % Find column names in "Samples" sheet and add to cell array
            sh_num1 = obj.getSheetNum('first_sample_num');
            sample_start_col = obj.getColNum('first_sample_num');
            sample_start_row = obj.getRowNum('first_sample_num') - 1;
            col_names = {};
            % look through columns in "Samples"
            for i=sample_start_col:size(obj.sheets{sh_num1},2)
                try 
                    ref_header = obj.getExcelValuePos(sh_num1, sample_start_row, i, 'char');
                catch 
                    try
                        ref_header = num2str(obj.getExcelValuePos(sh_num1, sample_start_row, i, 'numeric'));
                        if isempty(ref_header)
                            continue
                        end
                    catch 
                        continue
                    end
                end
                col_names{i} = ref_header;
            end
        end
        
        function new_coords = findLastSampleRow(obj)
        % Finds the last sample row in Samples sheet and adds info to
        % coordinates property 
            sh_num1 = obj.getSheetNum('first_sample_num');
            sample_start_col = obj.getColNum('first_sample_num');
            sample_start_row = obj.getRowNum('first_sample_num');
            coords = {};
            row_num = 0;
            for i=sample_start_row:size(obj.sheets{sh_num1},1)
                try 
                    num = num2str(obj.getExcelValuePos(sh_num1, i, sample_start_col, 'numeric'));
                    if isempty(num)
                        % found last
                        row_num = i-1;
                        break
                    end
                catch
                    % found last
                    row_num = i-1;
                    break
                end 
            end
            coords{end+1} = {sh_num1, row_num, sample_start_col};
            new_coords = obj.addExcelCoordinates({'last_sample_num'}, coords);
        end
        
        function value = getExcelValuePos(obj, sheet_num, row, col, type)
        % Returns the value at an inputted position. Error checks make sure
        % that the value is of the correct type
            sheet = obj.sheets{sheet_num};
            value = cell2mat(sheet(row,col));
            if and(isnan(value), TASBEConfig.get('template.displayErrors')) 
                TASBESession.error('TemplateExtraction','ValueNotFound','No value at position (%s, %s, %s).', num2str(sheet_num), num2str(row), num2str(col));
            elseif isnan(value)
                error('ValueNotFound');
            end
            if exist('type', 'var')
                if strcmp(type, 'cell') 
                    bounds = strtrim(strsplit(char(value), ','));
                    if and(isnan(str2double(bounds)), TASBEConfig.get('template.displayErrors')) 
                        TASBESession.error('TemplateExtraction','IncorrectType','Value at (%s, %s, %s) does not make a numeric array. Make sure value is in the form of #,#,#.', num2str(sheet_num), num2str(row), num2str(col));
                    elseif isnan(str2double(bounds))
                        error('IncorrectType');
                    else
                        value = num2cell(str2double(bounds));
                    end
                end
                if and(~isa(value, type), TASBEConfig.get('template.displayErrors') )
                    TASBESession.error('TemplateExtraction','IncorrectType','Value at (%s, %s, %s) of type %s, does not match with required type, %s.', num2str(sheet_num), num2str(row), num2str(col), class(value), type);
                elseif ~isa(value, type)
                    error('IncorrectType');
                end
            end
            value = obj.sanitizeFromExcel(value);
        end
        
        function value = getExcelValue(obj, name, type, index)
        % Returns the value of an inputted variable name using
        % getExcelValuePos. Also takes optional index for variables with
        % multiple coordinates.
            pos = obj.getExcelCoordinates(name);
            % index is considered for variables with more than one
            % coordinates
            if exist('index', 'var')
                pos = pos{index};
            end
            if exist('type', 'var')
                value = obj.getExcelValuePos(pos{1}, pos{2}, pos{3}, type);
            else
                value = obj.getExcelValuePos(pos{1}, pos{2}, pos{3});
            end
        end
        
        function setTASBEConfig(obj, name, type, index)
        % Sets a TASBEConfig given a variable name. getExcelValue is used
        % extensively in this function. Warnings/ errors are placed to make
        % sure there is a value found and that the variable is a
        % TASBEConfig. 
            try
                if exist('index', 'var')
                    value = obj.getExcelValue(name, type, index);
                else
                    value = obj.getExcelValue(name, type);
                end
                
            catch
                TASBESession.warn('TemplateExtraction','ValueNotFound','Name, %s, has no value at recorded position.', name);
                return
            end
            
            try
                TASBEConfig.isSet(name)
            catch
                TASBESession.error('TemplateExtraction', 'NotTASBEConfig', 'Could not get any preference for: %s', name);
            end
            TASBEConfig.set(name, value);
        end
    end
end

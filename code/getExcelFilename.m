% Helper function that obtains the correct filename for a given row in the
% "Samples" sheet with a TemplateExtraction object and row number as inputs. 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
function [filename] = getExcelFilename(extractor, row)
    % Take filenames from Excel and throw error if no filename is found
    path = extractor.path;
    sh_num1 = extractor.getSheetNum('first_sample_num');
    exclude_col = find(ismember(extractor.col_names, 'Exclude from Batch Analysis'), 1);
    if isempty(exclude_col)
        TASBESession.error('getExcelFilename', 'InvalidHeaderName', 'The header, Exclude from Batch Analysis, does not match with any column titles in "Samples" sheet.');
    end
    template_col = find(ismember(extractor.col_names, 'Template #'), 1);
    if isempty(template_col)
        TASBESession.error('getExcelFilename', 'InvalidHeaderName', 'The header, Template #, does not match with any column titles in "Samples" sheet.');
    end
    TASBEConfig.set('template.displayErrors', 1);
    template_num = extractor.getExcelValuePos(sh_num1, row, template_col, 'numeric');
    filename_templates = extractor.getExcelCoordinates('filename_templates');
    template_pos = filename_templates{template_num};
    TASBEConfig.set('template.displayErrors', 0);
    % Extracting the data stem path 
    try
        stem = extractor.getExcelValuePos(template_pos{1}, template_pos{2}+1, template_pos{3}+4, 'char');
        stem = make_filename_absolute(stem, path);
    catch
        TASBESession.warn('getExcelFilename','ValueNotFound','Template %s has no data stem.', num2str(template_num));
        stem = '';
    end
    
    % Get csv header (if applicable)
    try
        header = extractor.getExcelValuePos(template_pos{1}, template_pos{2}+1, template_pos{3}+5, 'char');
        [filepath,name,ext] = fileparts(header);
        filepath = make_filename_absolute(stem, filepath);
        header = [filepath name ext];
    catch
        header = '';
    end
    
    names = {};
    
    % Extracting the filenames
    for i=exclude_col+1:size(extractor.sheets{sh_num1},2)
        try
            name = extractor.getExcelValuePos(sh_num1, row, i);
            if isempty(name) && isempty(names)
                TASBESession.error('getExcelFilename', 'FilenameNotFound', 'No filename found at row %s. Make sure "Click to Update Filenames" button was pressed.', num2str(row));
            elseif ~isempty(name)
                names{end+1} = name;
            end
        catch
            if isempty(names)
                TASBESession.error('getExcelFilename', 'FilenameNotFound', 'No filename found at row %s. Make sure "Click to Update Filenames" button was pressed.', num2str(row));
            end
            break
        end
    end
    
    % Add stem name
    for i=1:numel(names)
        names{i} = [stem names{i}];
    end
    
    datafiles = {};
    for p=1:numel(names)
        if ~isempty(header)
            datafile = DataFile('csv', names{p}, header);
        else
            datafile = DataFile('fcs', names{p});
        end
        datafiles{end+1} = datafile;
    end
    
    filename = datafiles;
end
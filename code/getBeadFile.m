% GETBEADFILE returns FCS filenames for all bead files listed in the Sample Name cell
% for the Rainbow Bead section in the Calibration sheet
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.
function beadfiles = getBeadFile(extractor)
    TASBEConfig.set('template.displayErrors', 1);    
    bead_name = {extractor.getExcelValue('bead_name', 'char')};
    TASBEConfig.set('template.displayErrors', 0);
    if ~isempty(strfind(bead_name{1}, ',')) 
        bead_name = strtrim(strsplit(bead_name{1}, ','));
    end

    % Name contains a cell array of bead sample names
    % find their corresponding row numbers
    beadfiles = {};
    sh_num1 = extractor.getSheetNum('first_sample_num');
    first_sample_row = extractor.getRowNum('first_sample_num');
    sample_num_col = extractor.getColNum('first_sample_num');
    sample_name_col = find(ismember(extractor.col_names, 'SAMPLE NAME'), 1);
    if isempty(sample_name_col)
        TASBESession.error('getBeadFile', 'InvalidHeaderName', 'The header, SAMPLE NAME, does not match with any column titles in "Samples" sheet.');
    end
    % Go through samples in "Samples" sheet and look for matches in name to
    % elements in ref_filenames and add filenames to beadfiles array
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
        ind = find(ismember(bead_name, name), 1);
        if ~isempty(ind)
            datafile = getExcelFilename(extractor, i);
            beadfiles{ind} = datafile{1};
        end
    end
end
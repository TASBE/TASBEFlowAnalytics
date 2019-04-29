% Writes TASBEConfig preferences into a given sheet number and filepath of a
% batch_template spreadsheet
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
function TASBEConfig_preferences(filename, sheet_num)
    TASBEConfig.checkpoint('init');
    % Obtain all of the TASBEConfig preferences 
    [settings, ~, documentation] = TASBEConfig.list();
    % Extract just the names of the structure
    names = fieldnames(documentation);

    preferences = [];
    % Start from 2 to ignore the overall doc
    for i=2:numel(names)
        name = char(names{i});
        % Extract a section of TASBEConfig (i.e. beads)
        doc_section = getfield(documentation, name);
        set_section = getfield(settings, name);
        doc_names = fieldnames(doc_section);
        about_note = getfield(doc_section, doc_names{1});
        row = {name, char(about_note), ''};
        preferences = [preferences; row];
        % Go through each section and add preferences, their default values,
        % and documentation
        for j=2:numel(doc_names)
            sub_name = doc_names{j};
            note = getfield(doc_section, doc_names{j});
            value = getfield(set_section, doc_names{j});
            if iscell(value)
                value = strjoin(value);
                row = {[name '.' char(sub_name)], char(note), value};
            elseif ischar(value)
                row = {[name '.' char(sub_name)], char(note), value};
            else
                row = {[name '.' char(sub_name)], char(note), num2str(value)};
            end
            
            preferences = [preferences; row];
        end
    end
    % Write the preferences into the template spreadsheet. Formatting is
    % done by hand afterwards
    xlswrite(filename, preferences, sheet_num);
end


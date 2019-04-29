% ANALYZEFROMEXCEL: wrapper function that creates a TemplateExtraction object and calls the correct analysis. 
% The only needed inputs are the template file path and type of analysis that needs to be run.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.
function analyzeFromExcel(file, type)
    try
        % Setting up TASBESession log key
        TASBESession.warn('analyzeFromExcel', 'ExampleWarning', 'This is what a warning looks like.');
        TASBESession.succeed('analyzeFromExcel', 'ExampleSuccess', 'This is what a success/ notification looks like.');
        TASBESession.skip('analyzeFromExcel', 'ExampleSkip', 'A skip is also green.');
        try
            TASBESession.error('analyzeFromExcel', 'ExampleError', 'This is what an error looks like.');
        catch
            % Just makes sure that the program does not terminate after
            % establishing key for TASBESession log
        end
        
        % Editing file to get path, used for relative paths
        [filepath, name, ext] = fileparts(file);
        extractor = TemplateExtraction([end_with_slash(filepath) name ext]);
        
        % Running the actual analysis
        switch type
            case {'colormodel', 'CM', 'Colormodel'}
                % Make color model
                make_color_model_excel(extractor);
            case {'batch', 'BA', 'Batch'}
                % Run batch analysis
                batch_analysis_excel(extractor);
            case {'plusminus', 'PM', 'Plusminus', 'comparativeanalysis', 'companalysis', 'comparative'}
                % Run plus minus analysis
                plusminus_analysis_excel(extractor);
            case {'transfercurve', 'TC', 'Transfercurve'}
                % Run transfer curve analysis 
                transfercurve_analysis_excel(extractor);
            otherwise
                TASBESession.error('analyzeFromExcel', 'InvalidType', 'Input type of %s is invalid. The choices are colormodel, batch, plusminus, and transfercurve.', type);
        end
    catch exception
        % Turn MATLAB error into a TASBESession error
        if isempty(exception.identifier)
            TASBESession.error('analyzeFromExcel', 'NoIdentifier', exception.message);
        else
            if is_octave()
                msg = exception.message;
            else
                msg = strrep(sprintf(getReport(exception, 'extended', 'hyperlinks', 'off')), '\n', '');
            end
            id = exception.identifier;
            id_parts = strtrim(strsplit(id, ':'));
            if numel(id_parts) > 1
                name = '';
                for i=2:numel(id_parts)
                    if i > 2
                        name =[name ':' id_parts{i}];
                    else
                        name = id_parts{i};
                    end
                end
                TASBESession.error(id_parts{1}, name, msg);
            else
                TASBESession.error('analyzeFromExcel', exception.identifier, msg);
            end
        end
    end
end
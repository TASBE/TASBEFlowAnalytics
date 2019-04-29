% Constuctor for DataFile class with properties:
% type ('fcs' and 'csv')
% header
% file
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function DF = DataFile(type, file, header)
    % map of file types to CSVs
    key = {{'csv',{'.csv'}}, {'fcs',{'.fcs','.lmd'}}};
    
    if nargin == 0 % blank initialization
        DF.type = '';
        DF.file = '';
        DF.header = '';
    elseif nargin==1 % file name only
            %Create obj with one arguement if file has an extension of 'fcs'
            DF.type = 'fcs';
            DF.file = type;
            DF.header = '';
            [~, ~, fext] = fileparts(DF.file);
            for i=1:numel(key)
                sub_key = key{i};
                if strcmp(sub_key{1}, DF.type)
                    if isempty(find(strcmpi(sub_key{2}, fext), 1))
                        TASBESession.error('TASBE:DataFile','Underspecified','Only files with a fcs extension can create a DataFile obj with no supplied type');
                    end
                end
            end
    else % full constructor
        DF.type = type;
        DF.file = file;
        DF.header = '';

        % Make sure csv has header
        if strcmp(type, 'csv') && (nargin < 3)
            TASBESession.error('TASBE:DataFile','Underspecified','DataFile object with csv type needs a header');
        elseif strcmp(type, 'csv')
            DF.header = header;
        end

        % Make sure fcs doesn't have header
        if strcmp(type, 'fcs') && (nargin > 2)
            TASBESession.error('TASBE:DataFile','Overspecified','DataFile object with fcs type should not have a header');
        end

        % Check to make sure extension of file agrees with type
        [~, ~, fext] = fileparts(DF.file);
        error = 'InvalidType';
        for i=1:numel(key)
            sub_key = key{i};
            if strcmp(sub_key{1}, type)
                error = '';
                if isempty(find(strcmpi(sub_key{2}, fext), 1))
                    error = 'InvalidExtension';
                end
            end
        end

        if ~isempty(error)
            if strcmp(error, 'InvalidType')
                TASBESession.error('TASBE:DataFile',error,'File is not of type csv or fcs');
            elseif strcmp(error, 'InvalidExtension')
                TASBESession.error('TASBE:DataFile',error,'File extension does not match with possible file extensions');
            end
        end
    end
    
    DF=class(DF,'DataFile');
    
end
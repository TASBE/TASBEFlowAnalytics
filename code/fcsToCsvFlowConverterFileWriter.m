% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

% Call readfsc_compensated_MEFL to convert the data. Write the data to a
% CSV file and return the data.  This will overwrite any existing data.

function data = fcsToCsvFlowConverterFileWriter(CM, filename, with_AF, floor)
    % process the file to obtain point cloud
    data = readfcs_compensated_ERF(CM, filename, with_AF, floor);
    
    % create output filename for cloud
    [filepath,name,ext] = fileparts(filename);
    path = TASBEConfig.get('flow.pointCloudPath');
    if ~isdir(path),
        warning('TASBE:Utilities','Directory does not exist, attempting to create it: %s',path);
        mkdir(path);
    end
    
    csvName = [path '/' sanitize_name(name) '_PointCloud.csv'];
    
    % sanitize the channel names
    channels = getChannels(CM);
    sanitizedChannelName = cell(1, numel(channels));
    
    for i=1:numel(channels)
        channelName = getName(channels{i});
        invalidChars = '-|\s';  % Matlab does not like hypens or whitespace in variable names.
        sanitizedChannelName{i} = regexprep(channelName,invalidChars,'_');
    end
    
    % Use the channel names as the column labels
    columnLabels = strjoin(sanitizedChannelName, ',');
    
    % Write column labels to file
    fprintf('Writing Point Cloud CSV file: %s\n', csvName);
    fid = fopen(csvName,'w');
    fprintf(fid, '%s\n', columnLabels);
    fclose(fid);
    
    % Write the data to the file
    dlmwrite(csvName, data, '-append','precision','%.2f');
end


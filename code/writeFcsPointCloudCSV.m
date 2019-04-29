% WRITEFCSPOINTCLOUDCSV calls readfsc_compensated_MEFL to convert the data. Write the data to a
% CSV file and return the data.  This will overwrite any existing data.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function csv_filenames = writeFcsPointCloudCSV(CM, filenames, data)
    csv_filenames = {};
    if TASBEConfig.get('flow.outputPointCloud')
        n_conditions = numel(filenames);
        % Write each file for each condition
        for i=1:n_conditions
            perInducerFiles = filenames{i};
            numberOfPerInducerFiles = numel(perInducerFiles);
            for j = 1:numberOfPerInducerFiles
                fileName = perInducerFiles{j};
                % Write data
                filename = writeIndividualPointCloud(CM, fileName, data{i}{j});
                csv_filenames{end+1} = filename;
            end
        end
    end
end

function csv_filename = writeIndividualPointCloud(CM, filename, data)
    % create output filename for cloud
    [~,name,~] = fileparts(filename);
    path = TASBEConfig.get('flow.pointCloudPath');
    path = end_with_slash(path);
    if ~isdir(path)
        TASBESession.notify('TASBE:Utilities','MakeDirectory','Directory does not exist, attempting to create it: %s',path);
        mkdir(path);
    end
    name = sanitize_filename(name);
    % If filename already has point cloud extension don't need to add it
    if numel(name)>10 && strcmp(name(end-10:end), '_PointCloud')
        csvName = [path name '.csv'];
    else
        csvName = [path name '_PointCloud.csv'];
    end
    if TASBEConfig.get('flow.pointCloudFileType')
        % Store relative filename in JSON header
        if numel(name)>10 && strcmp(name(end-10:end), '_PointCloud')
            csv_filename = [name '.csv'];
        else
            csv_filename = [name '_PointCloud.csv'];
        end
    else
        csv_filename = csvName;
    end
    
    % sanitize the channel names
    channels = getChannels(CM);
    sanitizedChannelName = cell(1, numel(channels));

    for i=1:numel(channels)
        channelName = [getPrintName(channels{i}) '_' getUnits(channels{i})];
        sanitizedChannelName{i} = sanitizeColumnName(channelName);
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


function writePointCloudHeader(CM, filenames)
    % create JSON file of point cloud header (one for each batch)
    % header consists of the channel names
    path = TASBEConfig.get('flow.pointCloudPath');
    path = end_with_slash(path);
    if ~isdir(path)
        TASBESession.notify('TASBE:Utilities','MakeDirectory','Directory does not exist, attempting to create it: %s',path);
        mkdir(path);
    end
    channels = getChannels(CM);
    channel_structs = {};
    for i=1:numel(channels)
        channel_name = getName(channels{i});
        print_name = getPrintName(channels{i});
        unit = getUnits(channels{i});
        temp_struct = struct('name', channel_name, 'print_name', print_name, 'unit', unit);
        channel_structs{end+1} = temp_struct;
    end
    s = struct;
    s.channels = channel_structs;
    s.filenames = filenames;
    string = savejson('',s);
    stem_name = TASBEConfig.get('OutputSettings.StemName');
    if isempty(stem_name)
        stem_name = 'PointCloudHeader';
    end
    filename = strcat(path, stem_name, '.json');
    fid = fopen(filename, 'w');
    if fid == -1, error('Cannot create JSON file'); end
    fwrite(fid, string, 'char');
    fclose(fid);
end
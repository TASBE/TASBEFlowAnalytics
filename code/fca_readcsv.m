function [fcsdat, fcshdr, fcsdatscaled] = fca_readcsv(filename, headername, clip_events)
% [fcsdat, fcshdr, fcsdatscaled] = fca_readcsv(filename);
%
% Read CSV of flow cytometry data file and put the list mode  
% parameters to the fcsdat array with size of [NumOfPar TotalEvents]. 
% Some important header data are stored in the fcshdr structure:
% TotalEvents, NumOfPar, name

if isempty(filename)
    TASBESession.warn('CSV:Read','NoFile','No file provided! Returning empty dataset.'); 
    fcsdat = []; fcshdr = []; fcsdatscaled = [];
    return;
end
filecheck = dir(filename);
if size(filecheck,1) == 0
    TASBESession.warn('CSV:Read','NoFile',[filename,': The file does not exist! Returning empty dataset.']); 
    fcsdat = []; fcshdr = []; fcsdatscaled = [];
    return;
end
filecheck = dir(headername);
if size(filecheck,1) == 0
    TASBESession.warn('CSV:Read','NoHeaderFile',[headername,': The header file does not exist! Returning empty dataset.']); 
    fcsdat = []; fcshdr = []; fcsdatscaled = [];
    return;
end
if nargin<3, clip_events = 1e6; end

% If filename arg. only contain PATH, set the default dir to this
% before issuing the uigetfile command. This is an option for the "fca"
% tool
[FilePath, FileNameMain, fext] = fileparts(filename);
FilePath = [FilePath filesep];
FileName = [FileNameMain, fext];
if  isempty(FileNameMain)
    currend_dir = cd;
    cd(FilePath);
    [FileName, FilePath] = uigetfile('*.*','Select CSV file');
     filename = [FilePath,FileName];
     if FileName == 0
          fcsdat = []; fcshdr = [];
          return;
     end
     cd(currend_dir);
end

% Read in the data
if is_octave()
    T = csv2cell(filename);
    VarNames = T(1,:);
    % TotalEvents is number of rows
    fcshdr.TotalEvents = size(T(2:end,:),1);
else
    T = readtable(filename);
    VarNames = T.Properties.VariableNames;
    % TotalEvents is number of rows
    fcshdr.TotalEvents = height(T);
end

% Read in header info
fcshdr.fcstype = 'CSV1.0';
fcshdr.filename = FileName;
fcshdr.filepath = FilePath;
% NumOfPar is number of cols
fcshdr.NumOfPar = size(VarNames,2);
 
% Read in JSON header info to get fcshdr names
if nargin > 1
    fid = fopen(headername); 
    raw = fread(fid,inf); 
    string = char(raw'); 
    fclose(fid); 
    header = loadjson(string);
    channels = header.channels;
    filenames = header.filenames; 
    units = cell(fcshdr.NumOfPar,1);
    if numel(channels) ~= fcshdr.NumOfPar
        TASBESession.error('fca_readcsv', 'NumParameterMismatch', 'Number of columns in CSV %s not equal to number of channels specified in JSON header file %s',filename,headername);
    else
        for i=1:fcshdr.NumOfPar
            channel = channels{i};
            channel_fields = fieldnames(channel);
            % ensure all expected channel fields are populated, at least by defaults
            fcshdr.par(i) = fcs_channel();
            % make sure required fields are present in channel
            required_fields = {'name','unit','print_name'};
            for f=1:numel(required_fields)
                if isempty(find(cellfun(@(cf)(strcmp(required_fields{f},cf)),channel_fields), 1))
                    TASBESession.error('fca_readcsv', 'MissingRequiredHeaderField', 'Channel %i in %s missing required field %s',i,headername,required_fields{f});
                end
            end
            % copy channel fields, making sure they match target
            target_fields = fieldnames(fcshdr.par(i));
            for f=1:numel(channel_fields)
                % ensure the field is known
                if ~isempty(find(cellfun(@(tf)(strcmp(channel_fields{f},tf)),target_fields), 1))
                    fcshdr.par(i).(channel_fields{f}) = channel.(channel_fields{f});
                else
                    TASBESession.warn('fca_readcsv', 'UnknownHeaderField', 'Channel %i in %s contains unrecognized field %s',i,headername,channel_fields{f});
                end
            end
            % TODO: consider just using the values in the hdr.par fields?
            units{i} = channel.unit;
        end
    end
    
    % Double check units
    allowed_pattern = {'a\.u\.','ERF','Eum','M\w*','Boolean'};
    calibrated_unit_pattern = {'ERF','Eum','M\w*'};
    non_calibrated = 0;
    fcshdr.non_au = 1;
    for i = 1:numel(units)
        unit = units{i};
        valid_unit = 0;
        for j = 1:numel(allowed_pattern)
            matches = regexp(unit,allowed_pattern{j},'match');
            if(numel(matches)==1 && strcmp(unit,matches{1}))
                valid_unit = 1;
                if isempty(find(cellfun(@(p)(~isempty(regexp(unit,p,'match'))),calibrated_unit_pattern), 1));
                    non_calibrated = non_calibrated + 1;
                end
                break
            end
        end
        if valid_unit == 0
            TASBESession.error('fca_readcsv','UnknownUnits','Unit named %s is not a known permitted type',unit);
        end
    end
    % The file is a non-a.u. file if _any_ element is not a.u.
    fcshdr.non_au = (non_calibrated < numel(units));
    
    % consider both absolute and relative in comparing with filenames
    [HdrPath, ~,~] = fileparts(headername);
    % Read in filenames
    file_match = 0;
    filename_to_compare = strrep(filename, '\', '/');
    for i=1:numel(filenames)
        temp_filename = filenames{i};
        temp_filename = strrep(temp_filename, '\', '/');
        if strcmp(temp_filename, filename_to_compare) || strcmp(strcat(HdrPath, '/', temp_filename),filename_to_compare)
            file_match = 1;
            break
        end
    end
    
    if file_match ~= 1
        TASBESession.warn('fca_readcsv','FilenameMismatch','CSV file %s is not listed in JSON header %s',filename,headername);
    end
end
 
% Optionally truncate events to avoid memory problems with extremely large FCS files -JSB
if fcshdr.TotalEvents>clip_events
    TASBESession.warn('FCS:Read','TooManyEvents','FCS file has more than %i events; truncating to avoid memory problems',clip_events);
    fcshdr.TotalEvents = clip_events;
end

% Reading the events by setting fcsdat
if is_octave()
    T2 = cell2mat(T(2:end,:));
    fcsdat = double(T2);
else
    fcsdat = double(table2array(T));
end

% I don't believe we need fcsdatscaled because we don't have any log scales
fcsdatscaled = fcsdat;

end


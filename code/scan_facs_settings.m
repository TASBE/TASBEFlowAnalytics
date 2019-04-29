% SCAN_FACS_SETTINGS checks settings in FACS files, to allow validation of lasers, etc.
% - fileset is a cell array: {{filename [, nickname]}, ... } of files to compare
% - channels is a cell array: {{channel [, nickname]}, ... } of channels to investigate
% - fields is a cell array of fields to report values for
%   default is all non-null fields except: bit, range, name, decade, log, logzero
% - csvout is an optional name for a .csv file to dump a report to
% Returns a cell array of the channel descriptors from all files
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function descs = scan_facs_settings(fileset,channels,fields,csvout)
% if no fields given, take fields from first file's first channel
if nargin<3 || numel(fields) == 0
    [dat hdr] = fca_read(fileset{1}{1});
    [ignored desc] = get_fcs_color(dat,hdr,channels{1}{1},true);
    allfields = fieldnames(desc);
    fields = cell(0);
    ignorable = {'name','range','bit','decade','log','logzero'};
    for i=1:numel(allfields)
        if ~sum(cellfun(@(v)(strcmp(v,allfields{i})),ignorable)),
            fields{end+1} = allfields{i};
        end
    end
    nonnull_only = true;
else
    nonnull_only = false;
end

fprintf('Scanning FCS files');
nonnull = zeros(size(fields));
descs = cell(numel(fileset),numel(channels));
for f=1:numel(fileset)
    fprintf('.');
    [dat hdr] = fca_read(fileset{f}{1});
    for c=1:numel(channels)
        [ignored desc] = get_fcs_color(dat,hdr,channels{c}{1},true);
        descs{f,c} = desc;
        % check for contents
        if(numel(desc))
            for i=1:numel(fields)
                value = desc.(fields{i});
                if numel(value), nonnull(i) = 1; end;
            end
        end
    end
end
fprintf('\n');

% if we're only emitting non-null, get rid of the null fields
if nonnull_only,
    fprintf('Removing null fields\n');
    newfields = {};
    for i=1:numel(fields)
        if nonnull(i), newfields{end+1} = fields{i}; end;
    end
    fields = newfields;
end

% dump optional CSV file
if nargin >= 4
    fprintf('Outputting CSV file: %s\n',csvout);
    fid = fopen(csvout,'w+');
    % print header
    fprintf(fid,'File');
    for c=1:numel(channels)
        fprintf(fid,',%s',channels{c}{end});
        for i=2:numel(fields), fprintf(fid,','); end
    end
    fprintf(fid,'\n');
    for c=1:numel(channels)
        for i=1:numel(fields), fprintf(fid,',%s',fields{i}); end
    end
    
    for f=1:numel(fileset)
        fprintf(fid,'\n%s',fileset{f}{end});
        for c=1:numel(channels)
            desc = descs{f,c};
            for i=1:numel(fields),
                if(numel(desc))
                    value = desc.(fields{i});
                    if ischar(value), fprintf(fid,',%s',value);
                    else fprintf(fid,',%d',value);
                    end
                else
                    fprintf(fid,',---');
                end
            end
        end
    end
    fclose(fid);
end

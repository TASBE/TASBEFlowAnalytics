% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.
%
% Returns a list of possible beads that the user can reference when finding
% peaks

function [models] = get_bead_models()
    catalog = getBeadCatalog();
    models = cell(numel(catalog), 1);
    % go through all of the models and add to indices list
    for i=1:numel(catalog)
        modelEntry = catalog{i};
        models{i} =[modelEntry{1}];
    end
    return;
end

% The catalog has a 4-layer cell structure:
% { {'model', {'batch', {Channel, [peak, peak, ...]}, ... }, ...}, ...}
function returned = getBeadCatalog(forceReload)
    persistent catalog;
    if nargin<1, forceReload=false; end;
    if isempty(catalog) || forceReload,
        % ISSUE-81: The cell range of the spreadsheet must be updated whenever BeadCatalog.xlsx is updated.
        [nums txts combo] = xlsread(TASBEConfig.get('beads.catalogFileName'), 1, 'A1:M114');
        catalog = parseCatalog(combo);
    end
    returned = catalog;
end

function x = emptyOrNaN(v)
    if isempty(v),
        x = true; return;
    elseif isnumeric(v) && isnan(v(1)),
        x = true; return;
    else
        x = false; return;
    end
end

function catalog = parseCatalog(entries)
    currentLine = 1;
    catalog = {};
    while currentLine<=size(entries,1)
        if emptyOrNaN(entries{currentLine,1}) % skip lines that start with blanks
            currentLine = currentLine+1;
        else % otherwise, try to parse as a bead model entry
            [catalog{end+1} currentLine] = parseModel(entries,currentLine);
            % ensure there are no duplicates
            assert(isempty(find(cellfun(@(x)(strcmp(catalog{end}{1},x{1})),catalog(1:end-1)),1)),'Line %i: Multiple bead catalog entries for model %s',currentLine,catalog{end}{1});
        end
    end
end

function [model currentLine] = parseModel(entries,currentLine)
    % from first line, take first cell to be name; ignore rest
    name = entries{currentLine,1};
    assert(ischar(name),'Line %i: Expected bead model name, but first cell in row is not a string',currentLine);
    model = {name};
    currentLine = currentLine+1; % advance to start reading batches
    % now parse batches until we find one not starting with a string
    while currentLine<=size(entries,1) && ischar(entries{currentLine,1})
        [model{end+1} currentLine] = parseBatch(entries,currentLine);
        % ensure there are no duplicates
        assert(isempty(find(cellfun(@(x)(strcmp(model{end}{1},x{1})),model(2:end-1)),1)),'Line %i: Multiple bead catalog entries for batch %s in %s',currentLine,model{end}{1},name);
    end
end

function [batch currentLine] = parseBatch(entries,currentLine)
    % from first line, take first cell to be batch name; rest will be parsed amongst
    name = entries{currentLine,1};
    assert(ischar(name),'Line %i: Expected bead batch name, but first cell in row is not a string',currentLine);
    batch = {name};
    % read read of first line as channel, then advance
    batch{2} = parseChannel(entries,currentLine);
    currentLine = currentLine+1;
    % now parse batches until we hit the end or find a non-nan line or a blank link
    while currentLine<=size(entries,1) && ~(ischar(entries{currentLine,1}) || ~ischar(entries{currentLine,2}))
        batch{end+1} = parseChannel(entries,currentLine);
        currentLine = currentLine+1;
    end
end

function channelEntry = parseChannel(entries,currentLine)
    row = entries(currentLine,2:end);
    name = row{1};
    laser = row{2};
    if ischar(row{3}), filter = row{3};
    elseif emptyOrNaN(row{3}), filter = 'Unspecified';
    else error('Line %i: filter must be either a string or blank',currentLine);
    end
    units = row{4};
    try
        lastPeak = find(~cellfun(@emptyOrNaN,row(5:end)),1,'last');
    catch e
        error('Line %i: couldn''t interpret peak specifications',currentLine);
    end
    peaks = [row{4+(1:lastPeak)}];
    if isempty(peaks)
        error('Line %i: bead peak specification must contain at least one peak.',currentLine);
    end

    % Code added to make sure Octave doesn't throw away the empty peak
    % cells.
    cellPeaks = {row{4+(1:lastPeak)}};
    empties = cellfun('isempty', cellPeaks);
    cellPeaks(empties) = {NaN};
    peaks = cell2mat(cellPeaks);

    channelEntry = {name, units, peaks};
end

% READ_FILTERED_AU reads in an inputted FCS filename and outputs the raw
% data and header directory. 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [data,fcshdr,n_removed] = read_filtered_au(CM,datafile)
    % Get read FCS file and select channels of interest
    [~,fcshdr,rawfcs] = fca_read(datafile);
    if (isempty(fcshdr))
        TASBESession.error('TASBE:ReadFCS','CannotReadFile','Could not process FCS file %s', getFile(datafile));
    end
    
    data = rawfcs;
    n_raw_events = size(data,1);
    if n_raw_events<TASBEConfig.get('flow.smallFileWarning')
        TASBESession.warn('TASBE:ReadFCS','UnusuallySmallFile','FCS file "%s" is unusually small: only %i events', getFile(datafile), n_raw_events);
    end
    
    try
        non_au = fcshdr.non_au;
    catch
        non_au = 0;
    end
    
    if non_au ~= 1
        % optional discarding of filtered data (e.g., debris, time contamination)
        for i=1:numel(CM.prefilters)
            data = applyFilter(CM.prefilters{i},fcshdr,data);
        end
        % make sure we didn't throw away huge amounts...
        if numel(data)<numel(rawfcs)*TASBEConfig.get('flow.preGateDiscardsWarning')
            TASBESession.warn('TASBE:ReadFCS','TooMuchDataDiscarded','a.u. (pre)filters may be discarding too much data: only %d%% retained in %s',numel(data)/numel(rawfcs)*100,getFile(datafile));
        end
    end
    
    % if requested to dequantize, add a random value in [-0.5, 0.5]
    if(CM.dequantize), data = data + rand(size(data)) - 0.5; end

    % count how many have been removed, all told
    n_removed = n_raw_events - size(data,1);
    
% READFCS_COMPENSATED_AU reads in an inputted FCS filename by calling the
% read_filtered_au function with additional corrections done regarding
% autofluorescence and values below 1. 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [data,n_removed, non_au] = readfcs_compensated_au(CM,datafile,with_AF,floor)
    if(CM.initialized<0.5), TASBESession.error('TASBE:ReadFCS','Unresolved','Cannot read a.u.: ColorModel not yet resolved'); end % ensure initted

    % Acquire initial data, discarding likely contaminated portions
    [rawfcs, fcshdr, n_preremoved] = read_filtered_au(CM,datafile);
    try
        non_au = fcshdr.non_au;
    catch
        non_au = 0;
    end

    [rawdata, channel_desc] = select_channels(CM.Channels,rawfcs,fcshdr);
    
    if non_au ~= 1
        % check if voltages are identical to those of the color model, warn if otherwise
        ok = true;
        for i=1:numel(CM.Channels)
            ok = ok & confirm_channel(CM.Channels{i},channel_desc{i});
        end
        if(~ok), TASBESession.warn('TASBE:ReadFCS','ColorModelMismatch','File %s does not match color model',getFile(datafile)); end

        % Check to make sure not too many negative values:
        for i=1:numel(CM.Channels)
            frac_neg = sum(rawdata(:,i)<0)/size(rawdata,1);
            if(frac_neg>0.60) % more than 60% subzero
                TASBESession.warn('TASBE:ReadFCS','MajorityNegativeData','More than 60%% of channel %s negative (%.1f%%) in ''%s''',getName(CM.Channels{i}),100*frac_neg,getFile(datafile));
            end
        end

        % Remove autofluorescence from processed channels
        no_AF_data = zeros(size(rawdata));
        for i=1:numel(CM.Channels)
            if(~isUnprocessed(CM.Channels{i}))
                no_AF_data(:,i) = rawdata(:,i)-getMean(CM.autofluorescence_model{i});
            else
                no_AF_data(:,i) = rawdata(:,i);
            end
        end
        % make sure nothing's below 1, for compensation and geometric statistics
        % (compensation can be badly thrown off by negative values)
        if(floor && numel(CM.Channels)>1), no_AF_data(no_AF_data<1) = 1; end
        % Compensate for spectral bleed
        data = color_compensate(CM.compensation_model,no_AF_data);
        % Return autofluorescence, if desired
        if(with_AF)
            for i=1:numel(CM.Channels)
                if(~isUnprocessed(CM.Channels{i})) % only adjust processed channels
                    data(:,i) = data(:,i)+getMean(CM.autofluorescence_model{i});
                end
            end
        end
        % make sure nothing's below 1, for geometric statistics
        if(floor), data(data<1) = 1; end
        n_removed = n_preremoved + (size(rawfcs,1)-size(data,1));
    else
        data = rawdata;
        n_removed = n_preremoved + (size(rawfcs,1)-size(data,1));
    end
    

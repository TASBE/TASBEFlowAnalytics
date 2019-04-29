% READFCS_COMPENSATED_ERF reads in inputted FCS filename by calling the readfcs_compensated_au
% function and converts the outputted data to ERF units
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [data,n_removed] = readfcs_compensated_ERF(CM,datafile,with_AF,floor)
    if(CM.initialized<1), TASBESession.error('TASBE:ReadFCS','Unresolved','Cannot read ERF: ColorModel not yet resolved'); end % ensure initted
    
    % Read to arbitrary units
    [audata,n_preremoved, non_au] = readfcs_compensated_au(CM,datafile,with_AF,floor);
    
    if non_au ~= 1 
        % Translate each (processed) channel to ERF channel 
        ERF_channel_data = zeros(size(audata));
        for i=1:numel(CM.Channels)
            if(~isUnprocessed(CM.Channels{i}))
                ERF_channel_data(:,i) = translate(CM.color_translation_model,audata(:,i),CM.Channels{i},CM.ERF_channel);
            else
                ERF_channel_data(:,i) = audata(:,i);
            end
        end
        % Translate ERF AU to ERFs
        k_ERF= getK_ERF(CM.unit_translation);
        for i=1:numel(CM.Channels)
            if(~isUnprocessed(CM.Channels{i})) % only for processed channels
                data(:,i) = ERF_channel_data(:,i)*k_ERF;
            else
                data(:,i) = ERF_channel_data(:,i);
            end
        end
        % if possible, translate um channel AU to um
        if ~isempty(CM.size_unit_translation)
            i_um = find(CM, CM.um_channel);
            data(:,i_um) = um_channel_AU_to_um(CM.size_unit_translation,data(:,i_um));
        end
    else
        data = audata;
    end
    
    % apply post-filters
    % optional discarding of filtered data (e.g., poorly transfected cells)
    for i=1:numel(CM.postfilters)
        data = applyFilter(CM.postfilters{i},CM.Channels,data);
    end
    % make sure we didn't throw away huge amounts...
    if numel(data)<numel(audata)*TASBEConfig.get('flow.postGateDiscardsWarning')
        TASBESession.warn('TASBE:ReadFCS','TooMuchDataDiscarded','ERF (post)filters may be discarding too much data: only %d%% retained in %s',numel(data)/numel(audata)*100,getFile(datafile));
    end
    n_removed = n_preremoved + (size(audata,1)-size(data,1));

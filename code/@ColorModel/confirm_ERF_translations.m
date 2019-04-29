% CONFIRM_ERF_TRANSLATIONS checks to see that all channels can be translated to ERF
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [ok CM] = confirm_ERF_translations(CM)
    scales = getScales(CM.color_translation_model);
    ok = true;
    fi = indexof(CM.Channels,CM.ERF_channel);
    for i=1:numel(CM.Channels)
        if(i==fi || isUnprocessed(CM.Channels{i})) continue; end;
        if(isnan(scales(i,fi))), 
            ok = false;
            TASBESession.warn('TASBE:ColorModel','MissingTranslation','No pairwise translation for %s to %s; using pseudoERF',getPrintName(CM.Channels{i}),getPrintName(CM.ERF_channel));
            scales(i,fi) = 1;
            CM.Channels{i} = setIsPseudo(CM.Channels{i},1);
            
            % pseudo-ERFize channel
            AFMi = CM.autofluorescence_model{i};
            k_ERF=getK_ERF(CM.unit_translation);
            CM.autofluorescence_model{i}=ERFize(AFMi,scales(i,fi),k_ERF);
        else
            CM.Channels{i} = setUnits(CM.Channels{i},getUnits(CM.Channels{fi}));
        end
    end
    
    % if ERF channel is pseudo (e.g., becaus of missing beads), then make all channels pseudo
    if(isPseudo(CM.ERF_channel)),
        TASBESession.warn('TASBE:ColorModel','ERFChannelUncalibrated','ERF channel is pseudo, so all other channels are pseudo as well');
        for i=1:numel(CM.Channels),
            CM.Channels{i} = setIsPseudo(CM.Channels{i},1);
        end
    end
    
    if ~ok,
        TASBESession.warn('TASBE:ColorModel','UncalibratedChannels','Not all channels can be translated to standard ERF units.');
        CM.color_translation_model = setScales(CM.color_translation_model,scales);
    end

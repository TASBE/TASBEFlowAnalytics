% TRANSLATE scales inputted data from an original channel to a final
% channel for a ColorTranslationModel object.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function translated = translate(CTM,data,original,final)
    if(original==final), translated=data; return; end;
    % Color_j = matrix(
    i = indexof(CTM.Channels,original);
    if(i==-1), TASBESession.error('TASBE:ColorTranslation','MissingTranslationChannel','Missing translation channel: %s',getPrintName(original)); end;
    j = indexof(CTM.Channels,final);
    if(j==-1), TASBESession.error('TASBE:ColorTranslation','MissingTranslationChannel','Missing translation channel: %s',getPrintName(final)); end;
    if(isnan(CTM.scales(i,j))), TASBESession.error('TASBE:ColorTranslation','MissingTranslation','No pairwise translation for %s to %s',getPrintName(original),getPrintName(final)); end;

    translated = CTM.scales(i,j) * data;

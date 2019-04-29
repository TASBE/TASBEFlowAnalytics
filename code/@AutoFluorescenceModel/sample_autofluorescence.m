% SAMPLE_AUTOFLUORESCENCE computes autofluorescence from an
% AutoFluorescenceModel object. Also truncates values less than 1 to 1 if
% specified.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function AF = sample_autofluorescence(AFM,n_samples,truncate)

if nargin<3, truncate = 0; end;

AF = getMeanERF(AFM) + getStdERF(AFM)*randn(n_samples,1);

if truncate
    AF(AF<=1) = 1;
end

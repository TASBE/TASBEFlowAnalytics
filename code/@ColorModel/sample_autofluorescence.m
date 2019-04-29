% SAMPLE_AUTOFLUORESCENCE returns the sample autofluorescence for each
% channel in the inputted ColorModel object.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function AF = sample_autofluorescence(CM,n_samples,truncate,channels)

n_channels = numel(getChannels(CM));

if nargin<3, truncate = 0; end;
if nargin<4, channels = 1:n_channels; end;

AF = zeros(n_samples,numel(channels));
for i=1:numel(channels),
    AFM = get_autofluorescence_model(CM,channels(i));
    AF(:,i) = sample_autofluorescence(AFM,n_samples,truncate);
end

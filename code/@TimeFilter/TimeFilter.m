% Constructor of the TimeFilter class
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function TF = TimeFilter(exclusion)
% Optional discarding of early (possibly contaminated) data

if nargin == 0
    TF.early_data_exclusion = 25; % Discard first 25/100 seconds of data as likely contaminated
else
    TF.early_data_exclusion = exclusion;
end

TF = class(TF,'TimeFilter',Filter());

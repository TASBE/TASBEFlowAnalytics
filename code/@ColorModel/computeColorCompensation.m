% COMPUTECOLORCOMPENSATION produces an NxN matrix of linear compensation models
% Row j, column i is the fraction of channel i that bleeds into channel j
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function compensation_model = computeColorCompensation(CM)
n = numel(CM.Channels);
matrix = zeros(n,n); error = zeros(n,n);

for i=1:n
    if(isUnprocessed(CM.Channels{i}))
        TASBESession.notify('TASBE:Compensation','UnprocessedChannel','Skipping compensation computation for unprocessed channel %s',getName(CM.Channels{i}));
        matrix(i,i) = 1; % let channel through uncompensated
        continue;
    end
    for j=1:n
        % Don't need to compensate for self or for unprocessed
        if (i==j), matrix(j,i) = 1; error(j,i) = 0; continue; end
        if isUnprocessed(CM.Channels{j}), matrix(j,i) = 0; error(j,i) = 0; continue; end
        % Compute model
        [b_ij, b_ij_err] = make_linear_compensation_model(CM, CM.ColorFiles{i}, i, j);
        matrix(j,i) = b_ij; error(j,i) = b_ij_err;
    end
end

compensation_model = LinearCompensationModel(matrix,error);

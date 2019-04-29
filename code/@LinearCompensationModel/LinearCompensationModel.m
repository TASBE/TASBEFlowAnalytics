% Constructor of LinearCompensationModel class.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function LCM = LinearCompensationModel(matrix, error)
    if nargin == 0
        LCM.matrix = zeros(3,3);
        LCM.error = zeros(3,3);
    elseif nargin == 2
        LCM.matrix = matrix;
        LCM.error = error;
    end
    LCM= class(LCM,'LinearCompensationModel');
        

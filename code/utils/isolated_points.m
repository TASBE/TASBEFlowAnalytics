% ISOLATED_POINTS returns an indicator set that is true for all locations in the sequence 
% that will not plot properly without an explicit marker, being isolated 
% by a NaN, infinity, or complex number
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function which = isolated_points(sequence,log_scale)

if nargin>1 && log_scale, sequence(sequence<=0) = NaN; end;

if ~isvector(sequence), TASBESession.error('TASBE:Utilities','IsolatedPointVector','Sequence expected to be vector, but is not'); end;
if size(sequence,1)>1, sequence = sequence'; end; % transpose if needed

nonplotting = isnan(sequence) | isinf(sequence) | ~isreal(sequence);

which = ~nonplotting & [1 nonplotting(1:(end-1))] & [nonplotting(2:end) 1];

% tested by test_utilities

% CLEAR_PREFILTERS sets the prefilters property of an inputted ColorModel
% object to an empty array.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function CM=clear_prefilters(CM)
   CM.prefilters = {}; 

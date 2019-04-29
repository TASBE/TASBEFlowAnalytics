% CLEAR_FILTERS clears the prefilters and postfilters properties for an
% inputted ColorModel object
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function CM=clear_filters(CM)
   CM = clear_prefilters(CM); 
   CM = clear_postfilters(CM); 

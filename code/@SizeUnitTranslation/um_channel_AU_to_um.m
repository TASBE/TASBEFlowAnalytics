% ERF_CHANNEL_AU_TO_ERF converts ERF channel arbitrary units into standard ERF units
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.
    
function um = um_channel_AU_to_um(UT, data)
  um = 10.^polyval(UT.um_poly,log10(data));

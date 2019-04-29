% APPLYFILTER overrides applyFilter method in Filter class. Corresponds to
% a TimeFilter object, excludes data that occurred before a certain time threshold.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function data = applyFilter(TF,fcshdr,rawfcs)
% Optional discarding of early (possibly contaminated) data

if(numel(rawfcs)==0), data = rawfcs; return; end

found = false;
for j=1:numel(fcshdr.par)
    if(strcmp('Time',fcshdr.par(j).name))
        timechannel=rawfcs(:,j); 
        found=true; continue;
    end
end
if(~found), 
    TASBESession.warn('TASBE:TimeFilter','NoTimeChannel','Could not find Time channel'); 
    data = rawfcs;
else
    not_early = timechannel>TF.early_data_exclusion;
    data = rawfcs(not_early,:);
end

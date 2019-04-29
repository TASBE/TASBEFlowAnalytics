% GET_FCS_COLOR: search through the fcshdr file to determine which column
%   contains the color of interest, then return that column's data
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [data hdr] = get_fcs_color(fcsdat,fcshdr,color,suppress_errors)
if nargin<4, suppress_errors = 0; end;

data = []; hdr = [];
for i=1:size(fcshdr.par,2)
  if(strcmp(fcshdr.par(i).name,color))
      data = fcsdat(:,i); 
      hdr = fcshdr.par(i);
      return;
  elseif (strcmp(fcshdr.par(i).rawname,color))
      data = fcsdat(:,i); 
      hdr = fcshdr.par(i);
      return;
  end;
end;

if ~suppress_errors
    TASBESession.error('FCS:Select','MissingColor','Could not find color %s in FCS data',color);
end

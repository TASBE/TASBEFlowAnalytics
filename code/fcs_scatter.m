function [data figh] = fcs_scatter(filename,xcolor,ycolor,density,range,visible,largeoutliers,linear,filters)
% FCS_SCATTER(filename,xcolor,ycolor,density,range,visible): 
%   Plot a log-log scatter graph of positive points in FCS file
%   Defaults to non-visible, density 10, range = []

% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

if nargin < 6, visible = false; end;
if nargin < 7, largeoutliers = false; end;
if nargin < 8 || isempty(linear), linear = [0 0]; end;
if nargin < 9, filters = {}; end;

[fcsraw,fcshdr,fcsdat] = fca_readfcs(strtrim(filename));
% optional discarding of filtered data (e.g., debris, time contamination)
for i=1:numel(filters)
    fcsdat = applyFilter(filters{i},fcshdr,fcsdat);
end

xc = get_fcs_color(fcsdat,fcshdr,xcolor);
yc = get_fcs_color(fcsdat,fcshdr,ycolor);

pos = xc>0 & yc>0;
if linear(1), xv = xc; else xv = log10(xc(pos)); end;
if linear(2), yv = yc; else yv = log10(yc(pos)); end;

if nargin >= 4 && density
    smoothing = 10;
    h = figure('PaperPosition',[1 1 5 5]);
    if ~visible, set(h,'visible','off'); end;
    if nargin < 5, range = []; end;
    if density >= 1, type = 'image'; else type = 'contour'; end
    smoothhist2D([xv yv],smoothing,[200, 200],[],type,range,largeoutliers);
else
    h = figure('PaperPosition',[1 1 5 5]);
    if ~visible, set(h,'visible','off'); end;
    plot(xv,yv,'.','MarkerSize',1);
    if (nargin >= 5 && ~isempty(range)), xlim(range(:,1)); ylim(range(:,2)); end;
end
xlabel(clean_for_latex(['log_{10} ' xcolor ' a.u.'])); ylabel(clean_for_latex(['log_{10} ' ycolor ' a.u.']));
title(clean_for_latex(sprintf('Scatter of %s vs. %s for %s',xcolor,ycolor,filename)));

data = [xc yc];
figh = h;

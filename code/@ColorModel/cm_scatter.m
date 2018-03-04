function [data figh] = cm_scatter(CM,filename,xcolor,ycolor,ERF,density,range,visible,largeoutliers)
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

if nargin < 8, visible = false; end;
if nargin < 9, largeoutliers = false; end;

if ERF,
    data = readfcs_compensated_ERF(CM,strtrim(filename),false,false);
else
    data = readfcs_compensated_au(CM,strtrim(filename),false,false);
end

xchan = channel_named(CM,xcolor);
ychan = channel_named(CM,ycolor);
xc = data(:,find(CM,xchan));
yc = data(:,find(CM,ychan));

if nargin >= 6 && density
    h = figure('PaperPosition',[1 1 5 5]);
    if ~visible, set(h,'visible','off'); end;
    pos = xc>=1e7 & yc>=1;
    if nargin < 7, range = []; end;
    if density > 1, type = 'contour'; else type = 'image'; end
    smoothhist2D(log10([xc(pos) yc(pos)]),10,[200, 200],[],type,range,largeoutliers);
else
    h = figure('PaperPosition',[1 1 5 5]);
    if ~visible, set(h,'visible','off'); end;
    pos = xc>=1 & yc>=1;
    loglog(xc(pos),yc(pos),'.','MarkerSize',1);
    if (nargin >= 7 && ~isempty(range)), xlim(10.^range(:,1)); ylim(10.^(range(:,2))); end;
end
if ERF
    xlabel(['log_{10} ' xcolor ' ' clean_for_latex(getUnits(xchan))]); ylabel(['log_{10} ' ycolor ' ' clean_for_latex(getUnits(ychan))]);
else
    xlabel(clean_for_latex(['log_{10} ' xcolor ' a.u.'])); ylabel(clean_for_latex(['log_{10} ' ycolor ' a.u.']));
end
title(clean_for_latex(sprintf('Scatter of %s vs. %s for %s',xcolor,ycolor,filename)));

data = [xc yc];
figh = h;

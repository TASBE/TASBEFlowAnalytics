function h=ZeroOnLog(zero,start)

% Created by Jacob Beal in 2014
% based on BreakXAxis by Julie Haas (BSD license, copyright 2004), after Michael Robbins
% Assumes an already existing plot and axes, just waiting for relabeling and marking


yrange = ylim;
t1=text((zero+start)/2,yrange(1),'//','fontsize',15);
% t2=text((zero+start)/2,yrange(2),'//','fontsize',15);
% For y-axis breaks, use set(t1,'rotation',270);

% remap tick marks, and 'erase' them in the gap
xrange = xlim;

% Can't control minor range, so turn off and replace
set(gca,'XMinorTick','off')

mintick = floor(log10(start));
maxtick = floor(log10(xrange(2)));
zerotick = floor(log10(xrange(1)));

xtick = zeros((maxtick-mintick+1)*9,1);
for i=mintick:maxtick
    xtick((i-zerotick)*9+(1:9)) = (1:9)*10^i;
end
xtick = [zero; xtick(xtick>=start&xtick<=xrange(2))];

labels = cell(numel(xtick),1);
labels{1} = '0';

offset = 0.02;
for i=2:numel(xtick)
    if log10(xtick(i))==floor(log10(xtick(i)))
        label = sprintf('10^%i',log10(xtick(i)));
        text(xtick(i),yrange(1),...%(yrange(1)-offset*(yrange(2)-yrange(1))),...
            label,'HorizontalAlignment','center',...
            'VerticalAlignment','top','rotation',0);%,'interpreter','LaTeX');
    end
    labels{i} = '';
end

set(gca,'XTick',xtick);
set(gca,'XTickLabel',labels);

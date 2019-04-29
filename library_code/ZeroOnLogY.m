function h=ZeroOnLogY(zero,start)

% Create by Jacob Beal in 2014
% based on BreakXAxis by Julie Haas (BSD license, copyright 2004), after Michael Robbins
% Assumes an already existing plot and axes, just waiting for relabeling and marking


xrange = xlim;
t1=text(xrange(1),(zero+start)/2,'//','fontsize',15);
t2=text(xrange(2),(zero+start)/2,'//','fontsize',15);
set(t1,'rotation',270);
set(t2,'rotation',270);

% remap tick marks, and 'erase' them in the gap
yrange = ylim;

% Can't control minor range, so turn off and replace
set(gca,'YMinorTick','off')

mintick = floor(log10(start));
maxtick = floor(log10(yrange(2)));
zerotick = floor(log10(yrange(1)));

ytick = zeros((maxtick-mintick+1)*9,1);
for i=mintick:maxtick
    ytick((i-zerotick)*9+(1:9)) = (1:9)*10^i;
end
ytick = [zero; ytick(ytick>=start&ytick<=yrange(2))];

labels = cell(numel(ytick),1);
labels{1} = '0';

offset = 0.02;
for i=2:numel(ytick)
    if log10(ytick(i))==floor(log10(ytick(i)))
        label = sprintf('10^{%i}',log10(ytick(i)));
        text(xrange(1),ytick(i),...%(yrange(1)-offset*(yrange(2)-yrange(1))),...
            label,'HorizontalAlignment','right',...
            'VerticalAlignment','middle','rotation',0);%,'interpreter','LaTeX');
    end
    labels{i} = '';
end

set(gca,'YTick',ytick);
set(gca,'YTickLabel',labels);

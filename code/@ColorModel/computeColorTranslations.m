% Copyright (C) 2010-2017, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [colorTranslationModel CM] = computeColorTranslations(CM,settings)
% COMPUTECOLORTRANSLATIONS generates a model to use for color
% translation. The model relates pairs of colors to a fit. The
% transformation is of the form Color_j = scales(i,j)*Color_i 

n = numel(CM.Channels);
scales = zeros(n,n)*NaN;

for i=1:numel(CM.ColorPairFiles)
    cp = CM.ColorPairFiles{i};
    cX = indexof(CM.Channels,cp{1});
    if(cX==-1), error('Missing channel %s',getPrintName(cp{1})); end
    cY = indexof(CM.Channels,cp{2});
    if(cY==-1), error('Missing channel %s',getPrintName(cp{2})); end
    cCtrl = indexof(CM.Channels,cp{3});
    if(cCtrl==-1), error('Missing channel %s',getPrintName(cp{3})); end
    
    data = readfcs_compensated_au(CM,cp{4},false,true); % Leave out AF, use floor
    if(cX==cCtrl || cY==cCtrl),
        [scales(cX,cY) CM] = compute_two_color_translation_scale(CM,data,cX,cY,settings);
        [scales(cY,cX) CM] = compute_two_color_translation_scale(CM,data,cY,cX,settings);
    else
        [scales(cX,cY) CM] = compute_translation_scale(CM,data,cX,cY,cCtrl,settings);
        [scales(cY,cX) CM] = compute_translation_scale(CM,data,cY,cX,cCtrl,settings);
    end
    transerror = 10^(abs(log10(scales(cX,cY)*scales(cY,cX))));
    if(transerror > 1.05)
        warning('Model:Color','Translation from %s to %s not invertible (round trip error = %.2f)',getPrintName(cp{1}),getPrintName(cp{2}),transerror);
    end
end

colorTranslationModel = ColorTranslationModel(CM.Channels,scales);

end

function plot_translation_graph(CM,data,i,j,settings,scale,means,stds,which)
    h = figure('PaperPosition',[1 1 6 4]);
    set(h,'visible','off');
    %loglog(data(:,i),data(:,j),'b.','MarkerSize',1); hold on;
    %plot(means(which,i),means(which,j),'g*-');
    %plot([1e0 1e6],scale*[1e0 1e6],'r-');
    pos = data(:,i)>1 & data(:,j)>1;
    smoothhist2D(log10([data(pos,i) data(pos,j)]),10,[200, 200],[],[],[0 0; 6 6]); hold on;
    plot(log10(means(which,i)),log10(means(which,j)),'k*-');
    plot(log10(means(which,i)./stds(which,i)),log10(means(which,j).*stds(which,j)),'k:');
    plot(log10(means(which,i).*stds(which,i)),log10(means(which,j)./stds(which,j)),'k:');
    plot([0 6],log10(scale)+[0 6],'r-');
    xlim([0 6]); ylim([0 6]);
    xlabel(sprintf('%s a.u.',getName(CM.Channels{i})));
    ylabel(sprintf('%s a.u.',getName(CM.Channels{j})));
    title('Color Translation Model');
    path = getSetting(settings, 'path', './');
    outputfig(h,sprintf('color-translation-%s-to-%s', getPrintName(CM.Channels{i}),getPrintName(CM.Channels{j})), path);
end

function [scale CM] = compute_translation_scale(CM,data,i,j,ctrl,settings)
    % Average subpopulations, then find the ratio between them.
    bins = BinSequence(1.0,0.1,5.5,'log_bins'); % previously defaulted to 2.5
    % If minimums have been set, filter data to exclude any point that
    % doesn't meet them.
    if(~isempty(CM.translation_channel_min))
        which = data(:,i)>=10^CM.translation_channel_min(i) & ...
                data(:,j)>=10^CM.translation_channel_min(j) & ...
                data(:,ctrl)>=10^CM.translation_channel_min(ctrl);
        data = data(which,:);
        minbin_i = 10^CM.translation_channel_min(i);
        minbin_j = 10^CM.translation_channel_min(j);
    else
        minbin_i = 1e3; minbin_j = 1e3;
    end
    
    [counts means stds] = subpopulation_statistics(bins,data,ctrl,'geometric');
    which = find(counts(:)>CM.translation_channel_min_samples & means(:,i)>minbin_i & means(:,j)>minbin_j & means(:,i)<1e5 & means(:,j)<1e5); % ignore points without significant support or that are near saturation
    scale = means(which,i) \ means(which,j); % find ratio of j/i
    
    % ERFize channel if possible
    AFMi = CM.autofluorescence_model{i};
    k_ERF=getK_ERF(CM.unit_translation);
    if(CM.Channels{j}==CM.ERF_channel), CM.autofluorescence_model{i}=ERFize(AFMi,scale,k_ERF); end

    if CM.translation_plot
        plot_translation_graph(CM,data,i,j,settings,scale,means,stds,which);
    
% Plots against control, if needed
%         figure('PaperPosition',[1 1 6 4]);
%         loglog(data(:,ctrl),data(:,i),'b.','MarkerSize',1); hold on;
%         plot(means(which,ctrl),means(which,i),'g*-');
%         xlabel(sprintf('%s a.u.',CM.Channels{ctrl}.Filter));
%         ylabel(sprintf('%s a.u.',CM.Channels{i}.Filter));
%         title('Color Translation Model');
% 
%         figure('PaperPosition',[1 1 6 4]);
%         loglog(data(:,ctrl),data(:,j),'b.','MarkerSize',1); hold on;
%         plot(means(which,ctrl),means(which,j),'g*-');
%         xlabel(sprintf('%s a.u.',CM.Channels{ctrl}.Filter));
%         ylabel(sprintf('%s a.u.',CM.Channels{j}.Filter));
%         title('Color Translation Model');
    end
end

function [scale CM] = compute_two_color_translation_scale(CM,data,i,j,settings)
    % Average subpopulations, then find the ratio between them.
    bins = BinSequence(1.0,0.1,5.5,'log_bins');% previously defaulted to 2.5
    % If minimums have been set, filter data to exclude any point that
    % doesn't meet them.
%     if(~isempty(CM.translation_channel_min))
%         which = data(:,i)>=10^CM.translation_channel_min(i) & ...
%                 data(:,j)>=10^CM.translation_channel_min(j);
%         data = data(which,:);
%     end
    if(~isempty(CM.translation_channel_min))
        minbin_i = 10^CM.translation_channel_min(i);
        minbin_j = 10^CM.translation_channel_min(j);
    else
        minbin_i = 1e3; minbin_j = 1e3;
    end
    
    spine = sqrt(data(:,i).*data(:,j));
    num_rows = size(data,2);
    data(:,num_rows+1) = spine;
    
    [counts means stds] = subpopulation_statistics(bins,data,num_rows+1,'geometric');
    which = find(counts(:)>CM.translation_channel_min_samples & means(:,i)>minbin_i & means(:,j)>minbin_j & means(:,i)<1e5 & means(:,j)<1e5); % ignore points without significant support or that are near saturation
    scale = means(which,i) \ means(which,j); % find ratio of j/i
    
    % ERFize channel if possible
    AFMi = CM.autofluorescence_model{i};
    k_ERF=getK_ERF(CM.unit_translation);
    if(CM.Channels{j}==CM.ERF_channel), CM.autofluorescence_model{i}=ERFize(AFMi,scale,k_ERF); end

    if CM.translation_plot
        plot_translation_graph(CM,data,i,j,settings,scale,means,stds,which);
    
% Plots against control, if needed
%         figure('PaperPosition',[1 1 6 4]);
%         loglog(data(:,ctrl),data(:,i),'b.','MarkerSize',1); hold on;
%         plot(means(which,ctrl),means(which,i),'g*-');
%         xlabel(sprintf('%s a.u.',CM.Channels{ctrl}.Filter));
%         ylabel(sprintf('%s a.u.',CM.Channels{i}.Filter));
%         title('Color Translation Model');
% 
%         figure('PaperPosition',[1 1 6 4]);
%         loglog(data(:,ctrl),data(:,j),'b.','MarkerSize',1); hold on;
%         plot(means(which,ctrl),means(which,j),'g*-');
%         xlabel(sprintf('%s a.u.',CM.Channels{ctrl}.Filter));
%         ylabel(sprintf('%s a.u.',CM.Channels{j}.Filter));
%         title('Color Translation Model');
    end
end

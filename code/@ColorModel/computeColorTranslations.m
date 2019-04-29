% COMPUTECOLORTRANSLATIONS generates a model to use for color
% translation. The model relates pairs of colors to a fit. The
% transformation is of the form Color_j = scales(i,j)*Color_i 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function [colorTranslationModel, CM] = computeColorTranslations(CM)
n = numel(CM.Channels);
scales = zeros(n,n)*NaN;

for i=1:numel(CM.ColorPairFiles)
    cp = CM.ColorPairFiles{i};
    cX = indexof(CM.Channels,cp{1});
    if(cX==-1), TASBESession.error('TASBE:ColorTranslation','MissingTranslationChannel','Missing channel %s',getPrintName(cp{1})); end
    cY = indexof(CM.Channels,cp{2});
    if(cY==-1), TASBESession.error('TASBE:ColorTranslation','MissingTranslationChannel','Missing channel %s',getPrintName(cp{2})); end
    cCtrl = indexof(CM.Channels,cp{3});
    if(cCtrl==-1), TASBESession.error('TASBE:ColorTranslation','MissingTranslationChannel','Missing channel %s',getPrintName(cp{3})); end
    
    if(isUnprocessed(CM.Channels{cX})), TASBESession.error('TASBE:ColorTranslation','UnprocessedChannel','Cannot translate unprocessed channel %s',getName(cX)); end;
    if(isUnprocessed(CM.Channels{cY})), TASBESession.error('TASBE:ColorTranslation','UnprocessedChannel','Cannot translate unprocessed channel %s',getName(cY)); end;
    
    data = readfcs_compensated_au(CM,cp{4},false,true); % Leave out AF, use floor
    if(cX==cCtrl || cY==cCtrl),
        [scales(cX,cY), CM] = compute_two_color_translation_scale(CM,data,cX,cY);
        [scales(cY,cX), CM] = compute_two_color_translation_scale(CM,data,cY,cX);
    else
        [scales(cX,cY), CM] = compute_translation_scale(CM,data,cX,cY,cCtrl);
        [scales(cY,cX), CM] = compute_translation_scale(CM,data,cY,cX,cCtrl);
    end
    transerror = 10^(abs(log10(scales(cX,cY)*scales(cY,cX))));
    if(transerror > 1.05)
        TASBESession.warn('TASBE:ColorTranslation','NotInvertible','Translation from %s to %s not invertible (round trip error = %.2f)',getPrintName(cp{1}),getPrintName(cp{2}),transerror);
    end
end

colorTranslationModel = ColorTranslationModel(CM.Channels,scales);

end

function plot_translation_graph(CM,data,i,j,scale,means,stds,which)
    visiblePlots = TASBEConfig.get('colortranslation.visiblePlots');
    plotPath = TASBEConfig.get('colortranslation.plotPath');
    plotSize = TASBEConfig.get('colortranslation.plotSize');
    
    h = figure('PaperPosition',[1 1 plotSize]);
    if(~visiblePlots), set(h,'visible','off'); end;
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
    xlabel(sprintf('%s a.u.',clean_for_latex(getName(CM.Channels{i}))));
    ylabel(sprintf('%s a.u.',clean_for_latex(getName(CM.Channels{j}))));
    title(sprintf('Color Translation Model: %s to %s',clean_for_latex(getName(CM.Channels{i})),clean_for_latex(getName(CM.Channels{j}))));
    outputfig(h,sprintf('color-translation-%s-to-%s', clean_for_latex(getPrintName(CM.Channels{i})),clean_for_latex(getPrintName(CM.Channels{j}))), plotPath);
end

function [scale, CM] = compute_translation_scale(CM,data,i,j,ctrl)
    rangeMin = TASBEConfig.get('colortranslation.rangeMin');
    rangeMax = TASBEConfig.get('colortranslation.rangeMax');
    binIncrement = TASBEConfig.get('colortranslation.binIncrement');
    minSamples = TASBEConfig.get('colortranslation.minSamples');
    channelMinimum = TASBEConfig.getexact('colortranslation.channelMinimum',{});
    
    % Average subpopulations, then find the ratio between them.
    bins = BinSequence(rangeMin,binIncrement,rangeMax,'log_bins');
    % If minimums have been set, filter data to exclude any point that
    % doesn't meet them.
    if(~isempty(channelMinimum))
        which = data(:,i)>=10^channelMinimum(i) & ...
                data(:,j)>=10^channelMinimum(j) & ...
                data(:,ctrl)>=10^channelMinimum(ctrl);
        data = data(which,:);
        minbin_i = 10^channelMinimum(i);
        minbin_j = 10^channelMinimum(j);
    else
        minbin_i = 1e3; minbin_j = 1e3;
    end
    
    [counts, means, stds] = subpopulation_statistics(bins,data,ctrl,'geometric');
    nearMax = 10^(rangeMax-0.5);
    which = find(counts(:)>minSamples & means(:,i)>minbin_i & means(:,j)>minbin_j & means(:,i)<nearMax & means(:,j)<nearMax); % ignore points without significant support or that are near saturation
    scale = means(which,i) \ means(which,j); % find ratio of j/i
    
    % ERFize channel if possible
    AFMi = CM.autofluorescence_model{i};
    k_ERF=getK_ERF(CM.unit_translation);
    if(CM.Channels{j}==CM.ERF_channel), CM.autofluorescence_model{i}=ERFize(AFMi,scale,k_ERF); end

    if TASBEConfig.get('colortranslation.plot')
        plot_translation_graph(CM,data,i,j,scale,means,stds,which);
    end
end

function [scale, CM] = compute_two_color_translation_scale(CM,data,i,j)
    rangeMin = TASBEConfig.get('colortranslation.rangeMin');
    rangeMax = TASBEConfig.get('colortranslation.rangeMax');
    binIncrement = TASBEConfig.get('colortranslation.binIncrement');
    minSamples = TASBEConfig.get('colortranslation.minSamples');
    channelMinimum = TASBEConfig.getexact('colortranslation.channelMinimum',{});
    
    % Average subpopulations, then find the ratio between them.
    bins = BinSequence(rangeMin,binIncrement,rangeMax,'log_bins');
    % If minimums have been set, filter data to exclude any point that
    % doesn't meet them.
    if(~isempty(channelMinimum))
        minbin_i = 10^channelMinimum(i);
        minbin_j = 10^channelMinimum(j);
    else
        minbin_i = 1e3; minbin_j = 1e3;
    end
    
    spine = sqrt(data(:,i).*data(:,j));
    num_rows = size(data,2);
    data(:,num_rows+1) = spine;
    
    [counts, means, stds] = subpopulation_statistics(bins,data,num_rows+1,'geometric');
    nearMax = 10^(rangeMax-0.5);
    which = find(counts(:)>minSamples & means(:,i)>minbin_i & means(:,j)>minbin_j & means(:,i)<nearMax & means(:,j)<nearMax); % ignore points without significant support or that are near saturation
    scale = means(which,i) \ means(which,j); % find ratio of j/i
    
    % ERFize channel if possible
    AFMi = CM.autofluorescence_model{i};
    k_ERF=getK_ERF(CM.unit_translation);
    if(CM.Channels{j}==CM.ERF_channel), CM.autofluorescence_model{i}=ERFize(AFMi,scale,k_ERF); end

    if TASBEConfig.get('colortranslation.plot')
        plot_translation_graph(CM,data,i,j,scale,means,stds,which);
    end
end

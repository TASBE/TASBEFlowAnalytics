% A specific script written to create curves for a 6 FP x 6 channel coarse experiment 
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

%% Create curves for 6 FP x 6 channel coarse experiment
sources = {'Cerulean','EBFP2','AmCyan','EYFP','EGFP','mKate'};
source_files = {'Cer 1000 Dox','Hef1a-EBFP','AmCyan 1000 Dox','EYFP 1000 Dox','EGFP 1000 Dox','mKate 1000 Dox'};
filestem = 'ND04132011-P3-events/04-13-11_%s_P3.fcs';
source_channels = {{'Pacific Blue-A','AmCyan-A'},{'Pacific Blue-A','AmCyan-A'}, ...
    {'Pacific Blue-A','AmCyan-A'},{'FITC-A'},{'FITC-A'},{'PE YG-A','PE-Cy5-5 YG-A','PE-TxRed YG-A'}};
channels = {'Pacific Blue-A','AmCyan-A','FITC-A','PE YG-A','PE-Cy5-5 YG-A','PE-TxRed YG-A'};
o=1; clear med_interference_pct med_interference_abs;
for i=1:size(sources,2)
    for k=1:size(source_channels{i},2)
        for j=1:size(channels,2)
            [level,fit,err] = make_compensation_model(sprintf(filestem,source_files{i}),source_channels{i}{k},channels{j},0);
            which = find(level(:)>3.8 & level(:)<=4.5 & isnan(fit(:))==0);
            the_ylabel{o} = sprintf('%s/%s',sources{i},source_channels{i}{k});
            med_interference_pct(o,j) = mean((10.^fit(which))./(10.^level(which)));
        end
        o=o+1;
    end
end

% Result: med_interference =
%         Channel   450/50    510/50    530/50    575/26    695/40    610/20
%           alias  (PacBlue) (AmCyan)   (FITC)    (PE YG)  (PE-Cy5)   (TxRed)
% FP/Filter
% Cer/PacBlue       1.0246    0.4275    0.0241    0.0047    0.0076    0.0015
% Cer/AmCyan        2.4315    1.0248    0.0575    0.0074    0.0119    0.0027
% EBFP/PacBlue      1.0234    0.0520    0.0013    0.0036    0.0068    0.0012
% EBFP/AmCyan      17.2314    1.0253    0.0024    0.0093    0.0151    0.0019
% AmCyan/PacBlue    1.0227    0.6952    0.3718    0.0061    0.0085    0.0020
% AmCyan/AmCyan     1.1068    1.0275    0.5650    0.0072    0.0101    0.0029
% EYFP/FITC         0.0069    0.0087    1.0228    0.0201    0.0081    0.0016
% EGFP/FITC         0.0156    0.0559    1.0226    0.0068    0.0110    0.0020
% mKate/PE-YG       0.0195    0.0103    0.0109    1.0266    9.1194    0.8357
% mKate/PE-Cy5      0.0079    0.0035    0.0020    0.1162    1.0231    0.0914
% mKate/PE-TxRed    0.0241    0.0134    0.0129    1.2499   10.9898    1.0203
% 
% Note that because we're measuring only up to 10^4.5, and autofluorescence
% is around 10^2, any values around 0.01 or under are essentially identical.
% So, what we see is the following:
% - Filters cannot separate any of the FPs measured by the same laser 
%   (e.g. EBFP and Cerulean)
% - EBFP does not interfere with the 530/50 read
% - Cerulean, AmCyan, and mKate all interfere with the 530/50 read.
% - EGFP interferes with the 510/50 read; EYFP with the 575/26
% - It's best to read blues with 450/50, reds with 610/20
%
% No four-color combinations work.
% Three color combinations: EBFP x {EGFP,EYFP} x mKate
% Two color combinations: 
%    - subsets of three-color
%    - {AmCyan,Cerulean} x mKate
%
% This also means that we can expect the most valuable color permutations
% to be: Constitutive = EBFP, Input = {EYFP/mKate}, Output = {mKate/EYFP}
% The reason is that we want the constitutive to have the lowest interference.

% Graphical representation:
figure('PaperPosition',[1 1 11 6]);
imagesc(1:6,1:11,log10(med_interference_pct));
% colormap(gray); %% uncomment for grayscale version
colorbar;
xlabel('Read Channel','FontSize',16);
ylabel('Source/Filter','FontSize',16);
title('Mid-range Bleed-over percent (log_{10})','FontSize',16);
set(gca,'XTickLabel',channels);
set(gca,'YTickLabel',the_ylabel);
set(gca,'FontSize',11);
% Print as PNG because eps is "smartly" anti-aliasing
print('-dpng','coarse-bleedover-matrix.png');

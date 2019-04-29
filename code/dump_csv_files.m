% DUMP_CSV_FILES: Turn some of the key experiment results into
%    comma-separated-value files, which are easy for biologists
%    to pull into Excel and play with.
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed 
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function dump_csv_files(filestem, results)
% Get Inductions
experiment = getExperiment(results);
inductions = getInducerLevelsToFiles(experiment,1);
AP = getAnalysisParameters(results);
bins = get_bin_centers(getBins(AP));
cfp_units = getChannelUnits(AP,'constitutive');

% Note: the output depends on dlmwrite defaulting the delimiter to ','

% Counts file
fid = fopen([filestem '-counts.csv'],'w'); fprintf(fid,'Induction Levels\n'); fclose(fid);
dlmwrite([filestem '-counts.csv'],inductions,'-append');
fid = fopen([filestem '-counts.csv'],'a'); fprintf(fid,['\nBin centers (log10 ' cfp_units ')\n']); fclose(fid);
dlmwrite([filestem '-counts.csv'],log10(bins),'-append');
fid = fopen([filestem '-counts.csv'],'a'); fprintf(fid,'\nBin Counts\n'); fclose(fid);
dlmwrite([filestem '-counts.csv'],getBinCounts(results),'-append');
fid = fopen([filestem '-counts.csv'],'a'); fprintf(fid,'\nEstimated Fraction Active\n'); fclose(fid);
dlmwrite([filestem '-counts.csv'],getFractionActive(results),'-append');

% Values file
fid = fopen([filestem '-values.csv'],'w'); fprintf(fid,'Per-Channel Results\n\nInduction Levels\n'); fclose(fid);
dlmwrite([filestem '-values.csv'],inductions,'-append');
fid = fopen([filestem '-values.csv'],'a'); fprintf(fid,['\nBin centers (log10 ' cfp_units ')\n']); fclose(fid);
dlmwrite([filestem '-values.csv'],log10(bins),'-append');
channel_names = getChannelNames(AP);
for i=1:numel(channel_names),
    fid = fopen([filestem '-values.csv'],'a'); fprintf(fid,['\n\nChannel: ' channel_names{i}]); fclose(fid);
    [means stds stdofmeans stdofstds] = get_channel_results(results,channel_names{i});
    fid = fopen([filestem '-values.csv'],'a'); fprintf(fid,'\ngeometric means\n'); fclose(fid);
    dlmwrite([filestem '-values.csv'],means,'-append');
    fid = fopen([filestem '-values.csv'],'a'); fprintf(fid,'\ngeometric standard deviations\n'); fclose(fid);
    dlmwrite([filestem '-values.csv'],stds,'-append');
    fid = fopen([filestem '-values.csv'],'a'); fprintf(fid,'\ngeometric standard deviation of means across replicates\n'); fclose(fid);
    dlmwrite([filestem '-values.csv'],stdofmeans,'-append');
    fid = fopen([filestem '-values.csv'],'a'); fprintf(fid,'\ngeometric standard deviation of standard deviations across replicates\n'); fclose(fid);
    dlmwrite([filestem '-values.csv'],stdofstds,'-append');
end

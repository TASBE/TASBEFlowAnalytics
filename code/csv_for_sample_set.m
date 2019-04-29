% CSV_FOR_SAMPLE_SET displays information of inputted results saved within
% a CSV file
%
% Copyright (C) 2010-2018, Raytheon BBN Technologies and contributors listed
% in the AUTHORS file in TASBE analytics package distribution's top directory.
%
% This file is part of the TASBE analytics package, and is distributed
% under the terms of the GNU General Public License, with a linking
% exception, as described in the file LICENSE in the TASBE analytics
% package distribution's top directory.

function csv_for_sample_set(batch_description,results,filename)

if nargin<3, filename = 'batch.csv'; end

batch_size = size(batch_description,1);

fprintf('Outputting summary data');

f = fopen(filename,'w');

AP = getAnalysisParameters(results{1,1});
channels = getChannelNames(AP);
fprintf(f,'Sample Name');
for i=1:numel(channels), fprintf(f,',%s,,',channels{i}); end
fprintf(f,'\n');
for i=1:numel(channels), fprintf(f,',Count,Mean (log10 ERF),S.D. of Events (fold)'); end
fprintf(f,'\n');

for i=1:batch_size,
    condition_name = batch_description{i,1};
    r = results{i,1};
    
    fprintf(f,'%s',condition_name);
    for j=1:numel(channels),
        [means stds histograms] = get_channel_population_results(r,channels{j});
        fprintf(f,',%i,%.2f,%.2f',sum(histograms),log10(means),stds);
    end
    fprintf(f,'\n');
    fprintf('.');
end
fprintf('\n');


for i=1:size(results,1),
end

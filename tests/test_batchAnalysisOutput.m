function test_suite = test_batchAnalysisOutput
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_batchAnalysisEndtoend

TASBEConfig.set('flow.outputPointCloud','true');
TASBEConfig.set('flow.pointCloudPath','/tmp/CSV/');
TASBEConfig.set('flow.onThreshold', 10^4);

CM = load_or_make_testing_colormodel();
stem1011 = '../TASBEFlowAnalytics-Tutorial/example_assay/LacI-CAGop_';

% set up metadata
experimentName = 'LacI Transfer Curve';

% create default filenames based on experiment name
baseName = ['/tmp/' regexprep(experimentName,' ','_')];

% Configure the analysis
% Analyze on a histogram of 10^[first] to 10^[third] ERF, with bins every 10^[second]
bins = BinSequence(4,0.1,10,'log_bins');

% Designate which channels have which roles
AP = AnalysisParameters(bins,{});
% Ignore any bins with less than valid count as noise
AP=setMinValidCount(AP,100');
% Ignore any raw fluorescence values less than this threshold as too contaminated by instrument noise
AP=setPemDropThreshold(AP,5');
% Add autofluorescence back in after removing for compensation?
AP=setUseAutoFluorescence(AP,false');

% Make a map of condition names to file sets
file_pairs = {...
  'Dox 0.1/0.2',    {DataFile('fcs', [stem1011 'B3_P3.fcs']), DataFile('fcs', [stem1011 'B4_P3.fcs'])}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
  'Dox 0.5/1.0',    {DataFile('fcs', [stem1011 'B5_P3.fcs']), DataFile('fcs', [stem1011 'B6_P3.fcs'])};
  'Dox 2.0/5.0',    {DataFile('fcs', [stem1011 'B7_P3.fcs']), DataFile('fcs', [stem1011 'B8_P3.fcs'])};
% Remove these to let it be faster:
%  'Dox 10.0/20.0',   {[stem1011 'B9_B09_P3.fcs'], [stem1011 'B10_B10_P3.fcs']};
%  'Dox 50.0/100.0',   {[stem1011 'B11_B11_P3.fcs'], [stem1011 'B12_B12_P3.fcs']};
%  'Dox 200.0/500.0',  {[stem1011 'C1_C01_P3.fcs'], [stem1011 'C2_C02_P3.fcs']};
%  'Dox 1000.0/2000.0', {[stem1011 'C3_C03_P3.fcs'], [stem1011 'C4_C04_P3.fcs']};
  };

n_conditions = size(file_pairs,1);

% Execute the actual analysis
TASBEConfig.set('OutputSettings.StemName','LacI-CAGop');
[results, sampleresults] = per_color_constitutive_analysis(CM,file_pairs,{'EBFP2','EYFP','mKate'},AP);

% Make output plots
TASBEConfig.set('plots.plotPath','/tmp/plots');
TASBEConfig.set('OutputSettings.FixedInputAxis',[1e4 1e10]);
plot_batch_histograms(results,sampleresults,CM,{'b','y','r'});

% Make output plots without linespecs
TASBEConfig.set('plots.plotPath','/tmp/plots2');
TASBEConfig.set('OutputSettings.FixedInputAxis',[1e4 1e10]);
plot_batch_histograms(results,sampleresults,CM);

save('/tmp/LacI-CAGop-batch.mat','AP','bins','file_pairs','results','sampleresults');

TASBEConfig.set('flow.outputPointCloud','false');

% Test serializing the output
[statisticsFile, histogramFile] = serializeBatchOutput(file_pairs, CM, AP, sampleresults);

% Read the files into matlab tables
if (is_octave)
    statsTable = csv2cell(statisticsFile);
    histTable = csv2cell(histogramFile);
    statsCell = statsTable(2:end,:);
    histCell = histTable(2:end,:);
else
    statsTable = readtable(statisticsFile);
    histTable = readtable(histogramFile);
    statsCell = table2cell(statsTable);
    histCell = table2cell(histTable);
end

% Split the stats table
geoMeans = statsCell(:,5:7);
geoStdDevs = statsCell(:,8:10);
gmmMeans = statsCell(:,11:16);
gmmStds = statsCell(:,17:22);
gmmWeights = statsCell(:,23:28);
onFracs = statsCell(:,29);
offFracs = statsCell(:,30);

% Check on/off frac mean and std values
assertElementsAlmostEqual(round(results{1}.on_fracMean.*10000)./10000, 0.5791, 'relative', 1e-2);
assertElementsAlmostEqual(round(results{1}.off_fracMean.*10000)./10000, 0.4209, 'relative', 1e-2);
assertElementsAlmostEqual(round(results{1}.on_fracStd.*10000)./10000, 0.0038, 'relative', 1e-2);
assertElementsAlmostEqual(round(results{1}.off_fracStd.*10000)./10000, 0.0038, 'relative', 1e-2);

% Split the hist table
binCounts = histCell(:,3:5);

% Strip out the padding put into the sampleIds, means, and stdDevs
sampleIDListWithPadding = statsCell(:,1);
sampleIDs = sampleIDListWithPadding(find(~cellfun(@isempty,sampleIDListWithPadding)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CSV files:

% The first five rows should be enough to verify writing the histogram file
% correctly.
expected_bincounts = [...
        6799         637          36; % clipped by the drop threshold
        8012        2732        2696;
        8780        3327        2638;
        8563        4637        2632;
        7622        4623        3741;
        ];
       
% Means and stddevs tests writing the statistics file correctly.
expected_means = [...
    22170	260800	429100
    22200	260700	426500
    22110	269300	444900
    22050	271400	454400
    22160	262600	450600
    22550	260800	448200
    22820	266800	445000
    25400	268800	460800
    37920	253000	482200
    48930	242800	492200
    69260	220000	486300
    109400	182500	592500
    159100	160100	697900
    194800	146800	772000
    ];

expected_stds = [...
    1.601	6.611	7.921
    1.600	6.702	7.942
    1.599	6.710	7.944
    1.603	6.746	8.014
    1.604	6.599	7.922
    1.643	6.641	8.080
    1.703	6.595	7.936
    1.992	6.613	8.106
    3.057	6.308	8.224
    3.688	6.028	8.242
    4.508	5.731	8.056
    5.283	5.163	8.538
    5.603	4.613	8.392
    5.578	4.301	8.263
    ];

expected_gmm_means = 10.^[...
    4.2543    4.2543    4.6200    5.9921    4.6817    6.0832; % row 1
    ];

expected_gmm_stds = 10.^[...
    0.0692    0.0692    0.0783    0.3112    0.0900    0.5154; % row 1
    ];

expected_gmm_weights = [...
    0.5000    0.5000    0.4195    0.5805    0.3215    0.6785; % row 1
    ];

expected_on_fracs = [...
    0.5818
    0.5764
    0.5793
    0.5770
    0.5815
    0.5838
    ];

expected_off_fracs = [...
    0.4182
    0.4236
    0.4207
    0.4230
    0.4185
    0.4162
    ];
    
assertEqual(numel(sampleIDs), 3);
%assertEqual(numel(sampleIDs), 7);

% spot-check names
assertEqual(sampleIDs{1}, 'Dox 0.1/0.2');
assertEqual(sampleIDs{3}, 'Dox 2.0/5.0');
%assertEqual(sampleIDs{7}, 'Dox 1000.0/2000.0');

% spot-check first five rows of binCounts
assertElementsAlmostEqual(cell2mat(binCounts(1:5,:)), expected_bincounts, 'relative', 2e-2);

% spot-check geo means and geo std devs.
for i=1:6, % was 7
    assertElementsAlmostEqual(cell2mat(geoMeans(i,:)), expected_means(i,:), 'relative', 1e-2);
    assertElementsAlmostEqual(cell2mat(geoStdDevs(i,:)), expected_stds(i,:),  'relative', 1e-2);
end

% spot check Gaussian Mixture Model materials:
assertElementsAlmostEqual(cell2mat(gmmMeans(1,:)), expected_gmm_means, 'relative', 1e-2);
assertElementsAlmostEqual(cell2mat(gmmStds(1,:)), expected_gmm_stds, 'relative', 1e-2);
assertElementsAlmostEqual(cell2mat(gmmWeights(1,:)), expected_gmm_weights, 'relative', 1e-2);
assertElementsAlmostEqual(cell2mat(onFracs), expected_on_fracs, 'relative', 1e-2);
assertElementsAlmostEqual(cell2mat(offFracs), expected_off_fracs, 'relative', 1e-2);

% Check the first five rows of the first point cloud file
expected_pointCloud = [...
    42801.34    40500.46    33567.67
    2456.10     42822.39    1039.11
    70903.34    68176.25    20623.25
    2830130.69  17561178.05   1039.11
    8742.07     2238.27     1039.11
    ];

% The first point cloud file: /tmp/LacI-CAGop_B3_B03_P3_PointCloud.csv
firstPointCloudFile = '/tmp/CSV/LacI-CAGop_B3_P3_PointCloud.csv';

% Read the point cloud into matlab tables
if (is_octave)
    cloudTable = csv2cell(firstPointCloudFile);
    cloudCell = cloudTable(2:end,:);
else
  cloudTable = readtable(firstPointCloudFile);
  cloudCell = table2cell(cloudTable);
end

% Split the cloud table
points = cloudCell(1:5,:);

% spot-check first five rows of binCounts
assertElementsAlmostEqual(cell2mat(points), expected_pointCloud, 'relative', 1e-2);
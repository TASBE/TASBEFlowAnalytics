function test_suite = test_batchAnalysisOutput
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_batchAnalysisEndtoend

TASBEConfig.set('flow.outputPointCloud','true');
TASBEConfig.set('flow.pointCloudPath','/tmp/CSV/');

load('../TASBEFlowAnalytics-Tutorial/template_colormodel/CM120312.mat');
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
  'Dox 0.1/0.2',    {[stem1011 'B3_B03_P3.fcs'], [stem1011 'B4_B04_P3.fcs']}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
  'Dox 0.5/1.0',    {[stem1011 'B5_B05_P3.fcs'], [stem1011 'B6_B06_P3.fcs']};
  'Dox 2.0/5.0',    {[stem1011 'B7_B07_P3.fcs'], [stem1011 'B8_B08_P3.fcs']};
  'Dox 10.0/20.0',   {[stem1011 'B9_B09_P3.fcs'], [stem1011 'B10_B10_P3.fcs']};
  'Dox 50.0/100.0',   {[stem1011 'B11_B11_P3.fcs'], [stem1011 'B12_B12_P3.fcs']};
  'Dox 200.0/500.0',  {[stem1011 'C1_C01_P3.fcs'], [stem1011 'C2_C02_P3.fcs']};
  'Dox 1000.0/2000.0', {[stem1011 'C3_C03_P3.fcs'], [stem1011 'C4_C04_P3.fcs']};
  };

n_conditions = size(file_pairs,1);

% Execute the actual analysis
[results, sampleresults] = per_color_constitutive_analysis(CM,file_pairs,{'EBFP2','EYFP','mKate'},AP);

% Make output plots
OS = OutputSettings('LacI-CAGop','','','/tmp/plots');
OS.FixedInputAxis = [1e4 1e10];
plot_batch_histograms(results,sampleresults,OS,{'b','y','r'},CM);

save('/tmp/LacI-CAGop-batch.mat','AP','bins','file_pairs','OS','results','sampleresults');

TASBEConfig.set('flow.outputPointCloud','false');

% Test serializing the output
[statisticsFile, histogramFile] = serializeBatchOutput(file_pairs, CM, AP, sampleresults, baseName);

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
        6799        2200        1383;
        8012        2732        2696;
        8780        3327        2638;
        8563        4637        2632;
        7622        4623        3741;
        ];
       
% Means and stddevs tests writing the statistics file correctly.
expected_means = 1e5 * [...
    0.2217    2.4948    4.1064
    0.2219    2.4891    4.0757
    0.2211    2.5766    4.2599
    0.2205    2.5874    4.3344
    0.2216    2.5099    4.3095
    0.2255    2.4862    4.2764
    0.2281    2.5457    4.2586
    0.2539    2.5739    4.4073
    0.3791    2.4218    4.6213
    0.4891    2.3266    4.7217
    0.6924    2.1068    4.6593
    1.0930    1.7513    5.6729
    1.5909    1.5451    6.7144
    1.9472    1.4175    7.4609
    ];

expected_stds = [...
    1.6006    6.7653    8.1000
    1.5990    6.8670    8.1306
    1.5981    6.8650    8.1230
    1.6036    6.9155    8.2135
    1.6035    6.7565    8.1069
    1.6427    6.8020    8.2742
    1.7030    6.7618    8.1220
    1.9914    6.7701    8.2937
    3.0568    6.4579    8.4052
    3.6868    6.1704    8.4187
    4.5068    5.8686    8.2393
    5.2819    5.2780    8.7369
    5.6018    4.7061    8.5892
    5.5773    4.3900    8.4391
    ];


assertEqual(numel(sampleIDs), 7);

% spot-check names
assertEqual(sampleIDs{1}, 'Dox 0.1/0.2');
assertEqual(sampleIDs{7}, 'Dox 1000.0/2000.0');

% spot-check first five rows of binCounts
assertElementsAlmostEqual(cell2mat(binCounts(1:5,:)), expected_bincounts, 'relative', 1e-2);

% spot-check geo means and geo std devs.
for i=1:7,
    assertElementsAlmostEqual(cell2mat(geoMeans(i,:)), expected_means(i,:), 'relative', 1e-2);
    assertElementsAlmostEqual(cell2mat(geoStdDevs(i,:)), expected_stds(i,:),  'relative', 1e-2);
end

% Check the first five rows of the first point cloud file
expected_pointCloud = [...
    42801.34    40500.46    33567.67
    2456.10     42822.39    1039.11
    70903.34    68176.25    20623.25
    2830130.69  17561178.05   1039.11
    8742.07     2238.27     1039.11
    ];

% The first point cloud file: /tmp/LacI-CAGop_B3_B03_P3_PointCloud.csv
firstPointCloudFile = '/tmp/CSV/LacI-CAGop_B3_B03_P3_PointCloud.csv';

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
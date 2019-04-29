function test_suite = test_batch_analysis
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_batch_analysis_endtoend

CM = load_or_make_testing_colormodel();
stem1011 = '../TASBEFlowAnalytics-Tutorial/example_assay/LacI-CAGop_';

% set up metadata
experimentName = 'LacI Transfer Curve';

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
  'Dox 0.1',    {[stem1011 'B3_P3.fcs']}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
  'Dox 0.2',    {[stem1011 'B4_P3.fcs']}; % first couple are just strings to test implicit conversion
  'Dox 0.5',    {DataFile('fcs', [stem1011 'B5_P3.fcs'])};
  'Dox 1.0',    {DataFile('fcs', [stem1011 'B6_P3.fcs'])};
  'Dox 2.0',    {DataFile('fcs', [stem1011 'B7_P3.fcs'])};
  'Dox 5.0',    {DataFile('fcs', [stem1011 'B8_P3.fcs'])};
  'Dox 10.0',   {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
  'Dox 20.0',   {DataFile('fcs', [stem1011 'B10_P3.fcs'])};
  'Dox 50.0',   {DataFile('fcs', [stem1011 'B11_P3.fcs'])};
  'Dox 100.0',  {DataFile('fcs', [stem1011 'B12_P3.fcs'])};
  'Dox 200.0',  {DataFile('fcs', [stem1011 'C1_P3.fcs'])};
  'Dox 500.0',  {DataFile('fcs', [stem1011 'C2_P3.fcs'])};
  'Dox 1000.0', {DataFile('fcs', [stem1011 'C3_P3.fcs'])};
  'Dox 2000.0', {DataFile('fcs', [stem1011 'C4_P3.fcs'])};
  };

n_conditions = size(file_pairs,1);

% Execute the actual analysis
TASBEConfig.set('OutputSettings.StemName','LacI-CAGop');
[results, sampleresults] = per_color_constitutive_analysis(CM,file_pairs,{'EBFP2','EYFP','mKate'},AP);

% Make output plots
TASBEConfig.set('plots.plotPath','/tmp/plots');
TASBEConfig.set('OutputSettings.FixedInputAxis',[1e4 1e10]);
plot_batch_histograms(results,sampleresults,CM,{'b','g','r'});

save('/tmp/LacI-CAGop-batch.mat','AP','bins','file_pairs','results','sampleresults');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in results:

result1_expected_bincounts = [...
        6806         637           0; % clipped by the drop threshold
        8017        2753        2706;
        8782        3323        2637;
        8558        4640        2623;
        7617        4624        3739;
        6343        5595        4714;
        3931        5937        5304;
        1817        6282        5434;
         511        5747        5801;
         124        4272        4683;
           0        3097        4012;
           0        2284        3469;
           0        2340        2917;
           0        2545        3200;
           0        2845        3612;
           0        3390        3985;
           0        3755        4034;
           0        4031        3985;
           0        4246        4135;
           0        4436        4179;
           0        4502        4199;
           0        4289        4095;
           0        4007        3890;
           0        3630        3817;
           0        3244        3685;
           0        2738        3509;
           0        2203        3248;
           0        1731        3032;
           0        1406        2598;
           0         989        2401;
           0         769        1920;
           0         493        1626;
           0         391        1353;
           0         214         995;
           0         150         808;
           0         101         634;
           0           0         428;
           0           0         272;
           0           0         176;
           0           0         122;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           ];
       
result_expected_means = 1e5 * [...
    0.2214    2.5920    4.2576
    0.2217    2.5773    4.2344
    0.2201    2.6618    4.4116
    0.2193    2.6768    4.4791
    0.2202    2.5979    4.4447
    0.2230    2.5794    4.4459
    0.2246    2.6410    4.3857
    0.2489    2.6600    4.5635
    0.3723    2.5019    4.7657
    0.4765    2.3962    4.8756
    0.6809    2.1736    4.8155
    1.0768    1.8055    5.8724
    1.5798    1.5830    6.9178
    1.9350    1.4546    7.6877
    ];

result_expected_stds = [...
    1.5964    6.5596    7.8289
    1.5955    6.6127    7.8521
    1.5856    6.6173    7.8504
    1.5881    6.6379    7.8693
    1.5876    6.5152    7.7858
    1.6027    6.5557    7.9894
    1.6348    6.5152    7.7879
    1.8831    6.5316    7.9963
    2.9448    6.2255    8.0919
    3.5238    5.9371    8.1359
    4.3923    5.6446    7.9458
    5.1748    5.0904    8.4391
    5.5515    4.5412    8.2867
    5.5330    4.2440    8.2005
    ];

% Blue, Yellow, Red
result_expected1_gmm_means = [...
    4.2543    4.6200    4.6817
    4.2543    5.9921    6.0832
    ];
result_expected1_gmm_stds = [...
    0.0692    0.0783    0.0900
    0.0692    0.3112    0.5154
    ];
result_expected1_gmm_weights = [...
    0.5000    0.4195    0.3215
    0.5000    0.5805    0.6785
    ];

result_expected14_gmm_means = [...
    4.3427    4.6252    4.7382
    5.6856    5.5610    6.2619
    ];
result_expected14_gmm_stds = [...
    0.0904    0.0804    0.1034
    0.2825    0.2650    0.5099
    ];
result_expected14_gmm_weights = [...
    0.3486    0.4211    0.2456
    0.6514    0.5789    0.7544
    ];
    
assertEqual(numel(results), 14);

% spot-check name, bincenter, bin-count
assertEqual(results{1}.condition, 'Dox 0.1');
assertElementsAlmostEqual(log10(results{1}.bincenters([1 10 40 end])), [4.0500    4.9500    7.9500    9.9500], 'relative', 1e-2);
assertElementsAlmostEqual(results{1}.bincounts, result1_expected_bincounts,     'relative', 1e-2);

assertEqual(results{14}.condition, 'Dox 2000.0');
assertElementsAlmostEqual(log10(results{14}.bincenters([1 10 40 end])), [4.0500    4.9500    7.9500    9.9500], 'relative', 1e-2);

for i=1:14,
    assertElementsAlmostEqual(results{i}.means, result_expected_means(i,:), 'relative', 1e-2);
    assertElementsAlmostEqual(results{i}.stds,  result_expected_stds(i,:),  'relative', 1e-2);
end

assertElementsAlmostEqual(results{1}.gmm_means,  result_expected1_gmm_means,  'relative', 1e-2);
assertElementsAlmostEqual(results{1}.gmm_stds,  result_expected1_gmm_stds,  'relative', 1e-2);
assertElementsAlmostEqual(results{1}.gmm_weights,  result_expected1_gmm_weights,  'relative', 1e-2);
assertElementsAlmostEqual(results{14}.gmm_means,  result_expected14_gmm_means,  'relative', 1e-2);
assertElementsAlmostEqual(results{14}.gmm_stds,  result_expected14_gmm_stds,  'relative', 1e-2);
assertElementsAlmostEqual(results{14}.gmm_weights,  result_expected14_gmm_weights,  'relative', 1e-2);

% raw, filtered
firstlast_event_counts = [220379 183753; 161222 145854];
firstlast_events_used = [52544 107765 118200; 102451 98339 118601];
firstlast_events_outofrange = [108678 53457 43022; 43403 47515 27253];
assertElementsAlmostEqual(results{1}.n_events,  firstlast_event_counts(2,1),  'relative', 1e-2);
assertElementsAlmostEqual(results{14}.n_events, firstlast_event_counts(2,2),  'relative', 1e-2);
assertElementsAlmostEqual(results{1}.n_events_removed,  firstlast_event_counts(1,1)-firstlast_event_counts(2,1),  'relative', 1e-2);
assertElementsAlmostEqual(results{14}.n_events_removed, firstlast_event_counts(1,2)-firstlast_event_counts(2,2),  'relative', 1e-2);
assertElementsAlmostEqual(results{1}.n_events_used,  firstlast_events_used(1,:),  'relative', 1e-2);
assertElementsAlmostEqual(results{14}.n_events_used, firstlast_events_used(2,:),  'relative', 1e-2);
assertElementsAlmostEqual(results{1}.n_events_outofrange,  firstlast_events_outofrange(1,:),  'relative', 1e-2);
assertElementsAlmostEqual(results{14}.n_events_outofrange, firstlast_events_outofrange(2,:),  'relative', 1e-2);

function test_batch_analysis_nodrops

CM = load_or_make_testing_colormodel();
stem1011 = '../TASBEFlowAnalytics-Tutorial/example_assay/LacI-CAGop_';

% Configure the analysis
% Analyze on a histogram of 10^[first] to 10^[third] ERF, with bins every 10^[second]
bins = BinSequence(4,0.1,10,'log_bins');

% Designate which channels have which roles
AP = AnalysisParameters(bins,{});
% Ignore any bins with less than valid count as noise
AP=setMinValidCount(AP,100');
% Ignore any raw fluorescence values less than this threshold as too contaminated by instrument noise
AP=setPemDropThreshold(AP,0);
% Add autofluorescence back in after removing for compensation?
AP=setUseAutoFluorescence(AP,false');

% Make a map of condition names to file sets
file_pairs = {...
  'Dox 0.1',    {DataFile([stem1011 'B3_P3.fcs'])};
  'Dox 2000.0', {DataFile([stem1011 'C4_P3.fcs'])};
  };

% Execute the actual analysis
TASBEConfig.set('OutputSettings.StemName','LacI-CAGop');
[results, sampleresults] = per_color_constitutive_analysis(CM,file_pairs,{'EBFP2','EYFP','mKate'},AP);

% Make output plots
TASBEConfig.set('plots.plotPath','/tmp/plots');
TASBEConfig.set('OutputSettings.FixedInputAxis',[1e4 1e10]);
plot_batch_histograms(results,sampleresults,CM,{'b','g','r'});

save('/tmp/LacI-CAGop-nodrop.mat','AP','bins','file_pairs','results','sampleresults');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in results:

result1_expected_bincounts = [...
        6806        2178        1373; % not clipped by the drop threshold
        8017        2753        2706;
        8782        3323        2637;
        8558        4640        2623;
        7617        4624        3739;
        6343        5595        4714;
        3931        5937        5304;
        1817        6282        5434;
         511        5747        5801;
         124        4272        4683;
           0        3097        4012;
           0        2284        3469;
           0        2340        2917;
           0        2545        3200;
           0        2845        3612;
           0        3390        3985;
           0        3755        4034;
           0        4031        3985;
           0        4246        4135;
           0        4436        4179;
           0        4502        4199;
           0        4289        4095;
           0        4007        3890;
           0        3630        3817;
           0        3244        3685;
           0        2738        3509;
           0        2203        3248;
           0        1731        3032;
           0        1406        2598;
           0         989        2401;
           0         769        1920;
           0         493        1626;
           0         391        1353;
           0         214         995;
           0         150         808;
           0         101         634;
           0           0         428;
           0           0         272;
           0           0         176;
           0           0         122;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           0           0           0;
           ];
       
% spot-check name, bincenter, bin-count
assertEqual(results{1}.condition, 'Dox 0.1');
assertElementsAlmostEqual(log10(results{1}.bincenters([1 10 40 end])), [4.0500    4.9500    7.9500    9.9500], 'relative', 1e-2);
assertElementsAlmostEqual(results{1}.bincounts, result1_expected_bincounts,     'relative', 1e-2);


function test_batch_analysis_plot_warnings
% Test for warnings in plot_batch_histograms
CM2 = load_or_make_testing_colormodel2();
stem1011 = '../TASBEFlowAnalytics-Tutorial/example_assay/LacI-CAGop_';

% set up metadata
experimentName = 'LacI Transfer Curve';

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
  'Dox 0.1',    {DataFile('fcs', [stem1011 'B3_P3.fcs'])}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
  'Dox 0.2',    {DataFile('fcs', [stem1011 'B4_P3.fcs'])};
%   'Dox 0.5',    {DataFile('fcs', [stem1011 'B5_P3.fcs'])};
%   'Dox 1.0',    {DataFile('fcs', [stem1011 'B6_P3.fcs'])};
%   'Dox 2.0',    {DataFile('fcs', [stem1011 'B7_P3.fcs'])};
%   'Dox 5.0',    {DataFile('fcs', [stem1011 'B8_P3.fcs'])};
%   'Dox 10.0',   {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
%   'Dox 20.0',   {DataFile('fcs', [stem1011 'B10_P3.fcs'])};
%   'Dox 50.0',   {DataFile('fcs', [stem1011 'B11_P3.fcs'])};
%   'Dox 100.0',  {DataFile('fcs', [stem1011 'B12_P3.fcs'])};
%   'Dox 200.0',  {DataFile('fcs', [stem1011 'C1_P3.fcs'])};
%   'Dox 500.0',  {DataFile('fcs', [stem1011 'C2_P3.fcs'])};
  'Dox 1000.0', {DataFile('fcs', [stem1011 'C3_P3.fcs'])};
  'Dox 2000.0', {DataFile('fcs', [stem1011 'C4_P3.fcs'])};
  };

n_conditions = size(file_pairs,1);

% Execute the actual analysis
TASBEConfig.set('OutputSettings.StemName','LacI-CAGop');
[results, sampleresults] = per_color_constitutive_analysis(CM2,file_pairs,{'EBFP2','EYFP','mKate'},AP);

% Make output plots
TASBEConfig.set('plots.plotPath','/tmp/plots');
TASBEConfig.set('OutputSettings.FixedInputAxis',[1e4 1e10]);
plot_batch_histograms(results,sampleresults,CM2);
[~, ~] = serializeBatchOutput(file_pairs, CM2, AP, sampleresults); % should not error
log = TASBESession.list();
assertEqual(log{end}.contents{end}.name, 'NoLineSpecs');
assertExceptionThrown(@()plot_batch_histograms(results,sampleresults,CM2,{'b'}), 'plot_batch_histograms:LineSpecDimensionMismatch', 'No error was raised.');

% Execute the actual analysis
TASBEConfig.set('OutputSettings.StemName','LacI-CAGop');
[results, sampleresults] = per_color_constitutive_analysis(CM2,file_pairs,{'EYFP','mKate'},AP); % shouldn't error, but should set one color to 'k'
% Make output plots
TASBEConfig.set('plots.plotPath','/tmp/plots');
TASBEConfig.set('OutputSettings.FixedInputAxis',[1e4 1e10]);
plot_batch_histograms(results,sampleresults,CM2);
log = TASBESession.list();
assertEqual(log{end}.contents{end}.name, 'NoLineSpecs');
[~, ~] = serializeBatchOutput(file_pairs, CM2, AP, sampleresults); % should not error


function test_singlecolor_batch_analysis

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';

beadfile = [stem0312 'Beads_P3.fcs'];
blankfile = [stem0312 'blank_P3.fcs'];

% Create one channel / colorfile pair for each color
channels = {}; colorfiles = {};
channels{1} = Channel('FITC-A', 488, 515, 20);
channels{1} = setPrintName(channels{1}, 'EYFP'); % Name to print on charts
channels{1} = setLineSpec(channels{1}, 'y'); % Color for lines, when needed

% Multi-color controls are used for converting other colors into ERF units
% Any channel without a control mapping it to ERF will be left in arbirary units.
colorpairfiles = {};

CM = ColorModel(beadfile, blankfile, channels, colorfiles, colorpairfiles);

TASBEConfig.set('beads.beadModel','SpheroTech RCP-30-5A'); % Entry from BeadCatalog.xls matching your beads
TASBEConfig.set('beads.beadBatch','Lot AA01, AA02, AA03, AA04, AB01, AB02, AC01, GAA01-R'); % Entry from BeadCatalog.xls containing your lot
TASBEConfig.set('beads.rangeMin', 2);
CM=set_ERF_channel_name(CM, 'FITC-A');
TASBEConfig.set('plots.plotPath', '/tmp/plots');
% Execute the model
CM=resolve(CM);

stem1011 = '../TASBEFlowAnalytics-Tutorial/example_assay/LacI-CAGop_';
% set up metadata
experimentName = 'LacI Transfer Curve';
% Configure the analysis
bins = BinSequence(4,0.1,10,'log_bins');
AP = AnalysisParameters(bins,{});
AP=setMinValidCount(AP,100');
AP=setPemDropThreshold(AP,5');
AP=setUseAutoFluorescence(AP,false');

% Make a map of condition names to file sets
file_pairs = {...
  'Dox 0.1',    {DataFile([stem1011 'B3_P3.fcs'])};
  'Dox 0.2',    {DataFile([stem1011 'B4_P3.fcs'])};
  };

% Execute the actual analysis
TASBEConfig.set('OutputSettings.StemName','LacI-CAGop');
[results, sampleresults] = per_color_constitutive_analysis(CM,file_pairs,{'EYFP'},AP);

% Make output plots
TASBEConfig.set('plots.plotPath','/tmp/plots');
TASBEConfig.set('OutputSettings.FixedInputAxis',[1e4 1e10]);
plot_batch_histograms(results,sampleresults,CM);

save('/tmp/LacI-CAGop-batch-single.mat','AP','bins','file_pairs','results','sampleresults');


function test_dataless_batch_analysis

CM = load_or_make_testing_colormodel();
stem1011 = '../TASBEFlowAnalytics-Tutorial/example_assay/LacI-CAGop_';

% configure to discard all data from bins
bins = BinSequence(14,0.1,20,'log_bins');

AP = AnalysisParameters(bins,{});
AP=setMinValidCount(AP,100');
AP=setPemDropThreshold(AP,5');
AP=setUseAutoFluorescence(AP,false');

% Make a map of condition names to file sets
file_pairs = {...
  'Dox 0.1',    {[stem1011 'B3_P3.fcs']};
  };

% Execute the actual analysis
TASBEConfig.set('OutputSettings.StemName','LacI-CAGop');
[results, sampleresults] = per_color_constitutive_analysis(CM,file_pairs,{'EBFP2','EYFP','mKate'},AP);

% Make output plots
TASBEConfig.set('plots.plotPath','/tmp/plots');
TASBEConfig.set('OutputSettings.FixedInputAxis',[1e4 1e10]);
plot_batch_histograms(results,sampleresults,CM,{'b','g','r'});

save('/tmp/LacI-dataless-batch.mat','AP','bins','file_pairs','results','sampleresults');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in results:

assertEqual(numel(results), 1);

% see that counts are all zeros
assertElementsAlmostEqual(results{1}.bincounts, zeros(size(results{1}.bincounts)),     'relative', 1e-2);

assertElementsAlmostEqual(results{1}.means, nan(1,3), 'relative', 1e-2);
assertElementsAlmostEqual(results{1}.stds,  nan(1,3),  'relative', 1e-2);

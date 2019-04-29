function test_suite = test_transfercurve_analysis
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_batch_analysis_endtoend

CM = load_or_make_testing_colormodel();



% set up metadata
experimentName = 'LacI Transfer Curve';
device_name = 'LacI-CAGop';
inducer_name = 'Dox';

% Configure the analysis
% Analyze on a histogram of 10^[first] to 10^[third] ERF, with bins every 10^[second]
bins = BinSequence(4,0.1,10,'log_bins');

% Designate which channels have which roles
input = channel_named(CM, 'EBFP2');
output = channel_named(CM, 'EYFP');
constitutive = channel_named(CM, 'mKate');
AP = AnalysisParameters(bins,{'input',input; 'output',output; 'constitutive' constitutive});
% Ignore any bins with less than valid count as noise
AP=setMinValidCount(AP,100');
% Ignore any raw fluorescence values less than this threshold as too contaminated by instrument noise
AP=setPemDropThreshold(AP,5');
% Add autofluorescence back in after removing for compensation?
AP=setUseAutoFluorescence(AP,false');

% Make a map of induction levels to file sets
stem1011 = '../TASBEFlowAnalytics-Tutorial/example_assay/LacI-CAGop_';
level_file_pairs = {...
  0.1,    {DataFile('fcs', [stem1011 'B3_P3.fcs'])}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
  0.2,    {DataFile('fcs', [stem1011 'B4_P3.fcs'])};
  0.5,    {DataFile('fcs', [stem1011 'B5_P3.fcs'])};
%   1.0,    {DataFile('fcs', [stem1011 'B6_P3.fcs'])};
%   2.0,    {DataFile([stem1011 'B7_P3.fcs'])};
%   5.0,    {DataFile([stem1011 'B8_P3.fcs'])};
%   10.0,   {DataFile([stem1011 'B9_P3.fcs'])};
%   20.0,   {DataFile([stem1011 'B10_P3.fcs'])};
%   50.0,   {[stem1011 'B11_P3.fcs']};
%   100.0,  {[stem1011 'B12_P3.fcs']};
%   200.0,  {[stem1011 'C1_P3.fcs']};
%   500.0,  {[stem1011 'C2_P3.fcs']};
  1000.0, {DataFile('fcs', [stem1011 'C3_P3.fcs'])};
  2000.0, {DataFile('fcs', [stem1011 'C4_P3.fcs'])};
  };
experiment = Experiment(experimentName,{inducer_name}, level_file_pairs);

% Execute the actual analysis
fprintf('Starting analysis...\n');
[results, sampleresults] = process_transfer_curve( CM, experiment, AP);

% Make output plots
TASBEConfig.set('OutputSettings.StemName',device_name);
TASBEConfig.set('plots.plotPath','/tmp/plots');
% Plot how the constitutive fluorescence was distributed
TASBEConfig.set('OutputSettings.DeviceName',device_name);
plot_bin_statistics(sampleresults, getInducerLevelsToFiles(experiment,1));
% Plot the relation between inducer and input fluorescence
TASBEConfig.set('OutputSettings.DeviceName',inducer_name);
plot_inducer_characterization(results);
% Plot the relation between input and output fluorescence
TASBEConfig.set('OutputSettings.DeviceName',device_name);
plot_IO_characterization(results);

% Save the results of computation
save('-V7','/tmp/LacI-CAGop-Fine.mat','experiment','AP','sampleresults','results');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in results:

ers = struct(results);

result123_expected_bincounts = [...
           0           0           0;
        2706        2113        2088;
        2637        2204        1992;
        2623        2132        1995;
        3739        3058        2782;
        4714        3846        3671;
        5304        4421        4029;
        5434        4451        4323;
        5801        4614        4321;
        4683        3701        3442;
        4012        3186        3177;
        3469        2640        2661;
        2917        2394        2364;
        3200        2539        2482;
        3612        2874        2838;
        3985        3115        3092;
        4034        3229        3132;
        3985        3331        3165;
        4135        3274        3240;
        4179        3287        3260;
        4199        3265        3313;
        4095        3280        3103;
        3890        3220        3199;
        3817        3126        3095;
        3685        3074        2994;
        3509        2822        2756;
        3248        2612        2523;
        3032        2351        2366;
        2598        2202        2176;
        2401        1811        1856;
        1920        1645        1542;
        1626        1265        1351;
        1353        1088        1034;
         995         775         806;
         808         633         635;
         634         490         487;
         428         327         361;
         272         219         225;
         176         151         158;
         122         114         102;
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
       
% order of expected channels taken from batch analysis: 'EBFP2','EYFP','mKate'
% order from color model: EYFP, mKate, EBFP2
result_expected_means = 1e5 * [...
    2.6092    4.2917    0.1797
    2.6073    4.2659    0.1801
    2.6933    4.4501    0.1791
%     2.7147    4.5464    0.1786
%     2.6258    4.5050    0.1791
%     2.6080    4.4814    0.1822
%     2.6687    4.4504    0.1837
%     2.6884    4.6085    0.2033
%     2.5303    4.8246    0.2982
%     2.4281    4.9244    0.3844
%     2.2010    4.8634    0.5469
%     1.8253    5.9254    0.8774
    1.6010    6.9794    1.3203
    1.4686    7.7212    1.6497
    ];

result_expected_stds = [...
    6.6025    7.9140    1.8335
    6.6956    7.9290    1.8315
    6.7029    7.9377    1.8321
%     6.7406    8.0077    1.8316
%     6.5894    7.9130    1.8376
%     6.6332    8.0719    1.8720
%     6.5871    7.9260    1.9225
%     6.6026    8.0970    2.1899
%     6.2995    8.2167    3.2610
%     6.0207    8.2315    3.9419
%     5.7252    8.0455    4.8530
%     5.1557    8.5278    5.7880
    4.6059    8.3820    6.1890
    4.2960    8.2537    6.1751
    ];

assertEqual(size(ers.BinCounts), [60 5]);
assertEqual(size(ers.Means), [3 2]);
assertEqual(size(ers.Means{1,2}), [60 5]);

% spot-check bin counts
for i=1:3
    assertElementsAlmostEqual(ers.BinCounts(:,i), result123_expected_bincounts(:,i),     'relative', 1e-2);
end

assertEqual(struct(ers.PopMeans{1,1}).PrintName,'EYFP');
assertEqual(struct(ers.PopMeans{2,1}).PrintName,'mKate');
assertEqual(struct(ers.PopMeans{3,1}).PrintName,'EBFP2');

% order of expected channels taken from batch analysis: 'EBFP2','EYFP','mKate'
% order from color model: EYFP, mKate, EBFP2
assertElementsAlmostEqual(ers.PopMeans{1,2}, result_expected_means(:,1), 'relative', 1e-2);
assertElementsAlmostEqual(ers.PopMeans{2,2}, result_expected_means(:,2), 'relative', 1e-2);
assertElementsAlmostEqual(ers.PopMeans{3,2}, result_expected_means(:,3), 'relative', 1e-2);
assertElementsAlmostEqual(ers.PopStandardDevs{1,2}, result_expected_stds(:,1), 'relative', 1e-2);
assertElementsAlmostEqual(ers.PopStandardDevs{2,2}, result_expected_stds(:,2), 'relative', 1e-2);
assertElementsAlmostEqual(ers.PopStandardDevs{3,2}, result_expected_stds(:,3), 'relative', 1e-2);


% TODO: add validation of finer stats as well

function test_suite = test_AnalysisError
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_AnalysisCellCountErrors
% testing the extra cells issue for batch analysis, transfer curve
% analysis, and plus minus analysis

% BATCH ANALYSIS
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

% Make a faulty map of condition names to file sets
bad_file_pairs_batch1 = {...
  'Dox 0.1',    {DataFile('fcs', [stem1011 'B3_P3.fcs'])}, DataFile('fcs', [stem1011 'B3_P3.fcs']); 
  'Dox 0.2',    {DataFile('fcs', [stem1011 'B4_P3.fcs'])}, DataFile('fcs', [stem1011 'B3_P3.fcs']);
  'Dox 0.5',    {DataFile('fcs', [stem1011 'B5_P3.fcs'])}, DataFile('fcs', [stem1011 'B3_P3.fcs']);
  'Dox 1.0',    {DataFile('fcs', [stem1011 'B6_P3.fcs'])}, DataFile('fcs', [stem1011 'B3_P3.fcs']);
  'Dox 2.0',    {DataFile('fcs', [stem1011 'B7_P3.fcs'])}, DataFile('fcs', [stem1011 'B3_P3.fcs']);
  'Dox 5.0',    {DataFile('fcs', [stem1011 'B8_P3.fcs'])}, DataFile('fcs', [stem1011 'B3_P3.fcs']);
  'Dox 10.0',   {DataFile('fcs', [stem1011 'B9_P3.fcs'])}, DataFile('fcs', [stem1011 'B3_P3.fcs']);
  };

bad_file_pairs_batch2 = {...
  {DataFile('fcs', [stem1011 'B3_P3.fcs'])}; 
  {DataFile('fcs', [stem1011 'B4_P3.fcs'])};
  {DataFile('fcs', [stem1011 'B5_P3.fcs'])};
  {DataFile('fcs', [stem1011 'B6_P3.fcs'])};
  {DataFile('fcs', [stem1011 'B7_P3.fcs'])};
  {DataFile('fcs', [stem1011 'B8_P3.fcs'])};
  {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
  };

good_file_pairs_batch = {...
  'Dox 0.1',    {DataFile('fcs', [stem1011 'B3_P3.fcs'])}; 
  'Dox 0.2',    {DataFile('fcs', [stem1011 'B4_P3.fcs'])};
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


% Execute the actual analysis to see if an error gets thrown
assertExceptionThrown(@()per_color_constitutive_analysis(CM,bad_file_pairs_batch1,{'EBFP2','EYFP','mKate'},AP), 'TASBE:Analysis:DimensionMismatch', 'No error was raised.');
assertExceptionThrown(@()per_color_constitutive_analysis(CM,bad_file_pairs_batch2,{'EBFP2','EYFP','mKate'},AP), 'TASBE:Analysis:DimensionMismatch', 'No error was raised.');

% TRANSFER CURVE ANALYSIS
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

% Make a faulty map of induction levels to file sets
stem1011 = '../TASBEFlowAnalytics-Tutorial/example_assay/LacI-CAGop_';
bad_file_pairs_tc1 = {...
  0.1,    {DataFile('fcs', [stem1011 'B3_P3.fcs'])}, DataFile('fcs', [stem1011 'C4_P3.fcs']); % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
  0.2,    {DataFile('fcs', [stem1011 'B4_P3.fcs'])}, DataFile('fcs', [stem1011 'C4_P3.fcs']);
  0.5,    {DataFile('fcs', [stem1011 'B5_P3.fcs'])}, DataFile('fcs', [stem1011 'C4_P3.fcs']);
  1.0,    {DataFile('fcs', [stem1011 'B6_P3.fcs'])}, DataFile('fcs', [stem1011 'C4_P3.fcs']);
  2.0,    {DataFile('fcs', [stem1011 'B7_P3.fcs'])}, DataFile('fcs', [stem1011 'C4_P3.fcs']);
  5.0,    {DataFile('fcs', [stem1011 'B8_P3.fcs'])}, DataFile('fcs', [stem1011 'C4_P3.fcs']);
  10.0,   {DataFile('fcs', [stem1011 'B9_P3.fcs'])}, DataFile('fcs', [stem1011 'C4_P3.fcs']);
  };

bad_file_pairs_tc2 = {...
  {DataFile('fcs', [stem1011 'B3_P3.fcs'])}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
  {DataFile('fcs', [stem1011 'B4_P3.fcs'])};
  {DataFile('fcs', [stem1011 'B5_P3.fcs'])};
  {DataFile('fcs', [stem1011 'B6_P3.fcs'])};
  {DataFile('fcs', [stem1011 'B7_P3.fcs'])};
  {DataFile('fcs', [stem1011 'B8_P3.fcs'])};
  {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
  };

good_file_pairs_tc = {...
  0.1,    {DataFile('fcs', [stem1011 'B3_P3.fcs'])}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
  0.2,    {DataFile('fcs', [stem1011 'B4_P3.fcs'])};
  0.5,    {DataFile('fcs', [stem1011 'B5_P3.fcs'])};
  1.0,    {DataFile('fcs', [stem1011 'B6_P3.fcs'])};
  2.0,    {DataFile('fcs', [stem1011 'B7_P3.fcs'])};
  5.0,    {DataFile('fcs', [stem1011 'B8_P3.fcs'])};
  10.0,   {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
  };

% Execute the actual analysis to see if an error gets thrown
assertExceptionThrown(@()Experiment(experimentName,{inducer_name}, bad_file_pairs_tc1), 'TASBE:Experiment:DimensionMismatch', 'No error was raised.');
assertExceptionThrown(@()Experiment(experimentName,{inducer_name}, bad_file_pairs_tc2), 'TASBE:Experiment:DimensionMismatch', 'No error was raised.');

% PLUS MINUS ANALYSIS
CM = load_or_make_testing_colormodel();

% set up metadata
experimentName = 'LacI Transfer Curve';
device_name = 'LacI-CAGop';
inducer_name = '100xDox';

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

% Make a faulty map of the batches of plus/minus comparisons to test
% This analysis supports two variables: a +/- variable and a "tuning" variable
stem1011 = '../TASBEFlowAnalytics-Tutorial/example_assay/LacI-CAGop_';
bad_file_pairs_pm1 = {...
 {'Lows';'BaseDox';{'+', '-'};
  % First set is the matching "plus" conditions
  {0.1,  {DataFile('fcs', [stem1011 'B9_P3.fcs'])}, DataFile('fcs', [stem1011 'B9_P3.fcs']); % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
   0.2,  {DataFile('fcs', [stem1011 'B10_P3.fcs'])}, DataFile('fcs', [stem1011 'B9_P3.fcs'])};
  % Second set is the matching "minus" conditions 
  {0.1,  {DataFile('fcs', [stem1011 'B3_P3.fcs'])}, DataFile('fcs', [stem1011 'B9_P3.fcs']);
   0.2,  {DataFile('fcs', [stem1011 'B4_P3.fcs'])}, DataFile('fcs', [stem1011 'B9_P3.fcs'])}};
 {'Highs';'BaseDox';{'+', '-'};
  {10,   {DataFile('fcs', [stem1011 'C3_P3.fcs'])}, DataFile('fcs', [stem1011 'B9_P3.fcs']);
   20,   {DataFile('fcs', [stem1011 'C4_P3.fcs'])}, DataFile('fcs', [stem1011 'B9_P3.fcs'])};
  {10,   {DataFile('fcs', [stem1011 'B9_P3.fcs'])}, DataFile('fcs', [stem1011 'B9_P3.fcs']);
   20,   {DataFile('fcs', [stem1011 'B10_P3.fcs'])}, DataFile('fcs', [stem1011 'B9_P3.fcs'])}};
 };

bad_file_pairs_pm2 = {...
 {'Lows';'BaseDox';{'+', '-'};
  % First set is the matching "plus" conditions
  {{DataFile('fcs', [stem1011 'B9_P3.fcs'])}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
   {DataFile('fcs', [stem1011 'B10_P3.fcs'])}};
  % Second set is the matching "minus" conditions 
  {{DataFile('fcs', [stem1011 'B3_P3.fcs'])};
   {DataFile('fcs', [stem1011 'B4_P3.fcs'])}}};
 {'Highs';'BaseDox';{'+', '-'};
  {{DataFile('fcs', [stem1011 'C3_P3.fcs'])};
   {DataFile('fcs', [stem1011 'C4_P3.fcs'])}};
  {{DataFile('fcs', [stem1011 'B9_P3.fcs'])};
   {DataFile('fcs', [stem1011 'B10_P3.fcs'])}}};
 };

bad_file_pairs_pm3 = {...
 {'Lows';'BaseDox';{'+', '-'};'Extra';
  % First set is the matching "plus" conditions
  {0.1,  {DataFile('fcs', [stem1011 'B9_P3.fcs'])}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
   0.2,  {DataFile('fcs', [stem1011 'B10_P3.fcs'])}};
  % Second set is the matching "minus" conditions 
  {0.1,  {DataFile('fcs', [stem1011 'B3_P3.fcs'])};
   0.2,  {DataFile('fcs', [stem1011 'B4_P3.fcs'])}}};
 {'Highs';'BaseDox';{'+', '-'};'Extra';
  {10,   {DataFile('fcs', [stem1011 'C3_P3.fcs'])};
   20,   {DataFile('fcs', [stem1011 'C4_P3.fcs'])}};
  {10,   {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
   20,   {DataFile('fcs', [stem1011 'B10_P3.fcs'])}}};
 };

% Can become good if batch_names match:
bad_file_pairs_pm4 = {...
 {'Lows';'BaseDox';{'+', '-'};
  % First set is the matching "plus" conditions
  {0.1,  {DataFile('fcs', [stem1011 'B9_P3.fcs'])}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
   0.2,  {DataFile('fcs', [stem1011 'B10_P3.fcs'])}};
  % extra set
  {0.1,  {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
   0.2,  {DataFile('fcs', [stem1011 'B10_P3.fcs'])}};
  % Second set is the matching "minus" conditions 
  {0.1,  {DataFile('fcs', [stem1011 'B3_P3.fcs'])};
   0.2,  {DataFile('fcs', [stem1011 'B4_P3.fcs'])}}};
 {'Highs';'BaseDox';{'+', '-'};
  {10,   {DataFile('fcs', [stem1011 'C3_P3.fcs'])};
   20,   {DataFile('fcs', [stem1011 'C4_P3.fcs'])}};
  {0.1,  {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
   0.2,  {DataFile('fcs', [stem1011 'B10_P3.fcs'])}};
  {10,   {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
   20,   {DataFile('fcs', [stem1011 'B10_P3.fcs'])}}};
 };

good_file_pairs_pm1 = {...
 {'Lows';'BaseDox';{'+', '-'};
  % First set is the matching "plus" conditions
  {0.1,  {DataFile('fcs', [stem1011 'B9_P3.fcs'])}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
   0.2,  {DataFile('fcs', [stem1011 'B10_P3.fcs'])}};
  % Second set is the matching "minus" conditions 
  {0.1,  {DataFile('fcs', [stem1011 'B3_P3.fcs'])};
   0.2,  {DataFile('fcs', [stem1011 'B4_P3.fcs'])}}};
 {'Highs';'BaseDox';{'+', '-'};
  {10,   {DataFile('fcs', [stem1011 'C3_P3.fcs'])};
   20,   {DataFile('fcs', [stem1011 'C4_P3.fcs'])}};
  {10,   {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
   20,   {DataFile('fcs', [stem1011 'B10_P3.fcs'])}}};
 };

good_file_pairs_pm2 = {...
 {'Lows';'BaseDox';{'+', '-'};
  % First set is the matching "plus" conditions
  {0.1,  {DataFile('fcs', [stem1011 'B9_P3.fcs'])}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
   0.2,  {DataFile('fcs', [stem1011 'B10_P3.fcs'])}};
  % Second set is the matching "minus" conditions 
  {0.1,  {DataFile('fcs', [stem1011 'B3_P3.fcs'])};
   0.2,  {DataFile('fcs', [stem1011 'B4_P3.fcs'])}}};
 {'Highs';'BaseDox';{'+', '-'};
  {10,   {DataFile('fcs', [stem1011 'C3_P3.fcs'])};
   20,   {DataFile('fcs', [stem1011 'C4_P3.fcs'])}};
  {10,   {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
   20,   {DataFile('fcs', [stem1011 'B10_P3.fcs'])}}};
 {'Extra';'BaseDox';{'+', '-'};
  {10,   {DataFile('fcs', [stem1011 'C3_P3.fcs'])};
   20,   {DataFile('fcs', [stem1011 'C4_P3.fcs'])}};
  {10,   {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
   20,   {DataFile('fcs', [stem1011 'B10_P3.fcs'])}}};
 };

% Good if batch_names contains three elements
good_file_pairs_pm3 = {...
 {'Lows';'BaseDox';{'+', '-', 'extra'};
  % First set is the matching "plus" conditions
  {0.1,  {DataFile('fcs', [stem1011 'B9_P3.fcs'])}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
   0.2,  {DataFile('fcs', [stem1011 'B10_P3.fcs'])}};
  % extra set
  {0.1,  {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
   0.2,  {DataFile('fcs', [stem1011 'B10_P3.fcs'])}};
  % Second set is the matching "minus" conditions 
  {0.1,  {DataFile('fcs', [stem1011 'B3_P3.fcs'])};
   0.2,  {DataFile('fcs', [stem1011 'B4_P3.fcs'])}}};
 {'Highs';'BaseDox';{'+', '-', 'extra'};
  {10,   {DataFile('fcs', [stem1011 'C3_P3.fcs'])};
   20,   {DataFile('fcs', [stem1011 'C4_P3.fcs'])}};
  {0.1,  {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
   0.2,  {DataFile('fcs', [stem1011 'B10_P3.fcs'])}};
  {10,   {DataFile('fcs', [stem1011 'B9_P3.fcs'])};
   20,   {DataFile('fcs', [stem1011 'B10_P3.fcs'])}}};
 };

% Execute the actual analysis to see if an error gets thrown
assertExceptionThrown(@()process_plusminus_batch( CM, bad_file_pairs_pm1, AP), 'TASBE:Analysis:ColumnDimensionMismatch', 'No error was raised.');
assertExceptionThrown(@()process_plusminus_batch( CM, bad_file_pairs_pm2, AP), 'TASBE:Analysis:ColumnDimensionMismatch', 'No error was raised.');
assertExceptionThrown(@()process_plusminus_batch( CM, bad_file_pairs_pm3, AP), 'TASBE:Analysis:SetDimensionMismatch', 'No error was raised.');
assertExceptionThrown(@()process_plusminus_batch( CM, bad_file_pairs_pm4, AP), 'TASBE:Analysis:SetDimensionMismatch', 'No error was raised.');

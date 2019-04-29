function test_suite = test_beadpeaks
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

%%%%%%%%%%%%
% Note: test_colormodel tests the one and many bead cast; we just need to test a few special cases

function [CM] = setupRedPeakCM()

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';

beadfile = DataFile('fcs', [stem0312 'Beads_P3.fcs']);
blankfile = DataFile('fcs', [stem0312 'blank_P3.fcs']);

% Create one channel / colorfile pair for each color
channels = {}; colorfiles = {};
channels{1} = Channel('PE-Tx-Red-YG-A', 561, 610, 20);
channels{1} = setPrintName(channels{1}, 'mKate');
channels{1} = setLineSpec(channels{1}, 'r');
colorfiles{1} = DataFile('fcs', [stem0312 'mkate_P3.fcs']);

% Multi-color controls are used for converting other colors into ERF units
% Any channel without a control mapping it to ERF will be left in arbirary units.
colorpairfiles = {};

CM = ColorModel(beadfile, blankfile, channels, colorfiles, colorpairfiles);

TASBEConfig.set('beads.beadModel','SpheroTech RCP-30-5A'); % Entry from BeadCatalog.xls matching your beads
TASBEConfig.set('beads.beadBatch','Lot AA01, AA02, AA03, AA04, AB01, AB02, AC01, GAA01-R'); % Entry from BeadCatalog.xls containing your lot
TASBEConfig.set('beads.beadChannel','PE-TR');

CM=set_ERF_channel_name(CM, 'PE-Tx-Red-YG-A');

TASBEConfig.set('plots.plotPath', '/tmp/plots');


function test_twopeaks
TASBEConfig.checkpoint('test');

[CM] = setupRedPeakCM();
% Ignore all bead data below 10^rangeMin as being too "smeared" with noise
TASBEConfig.set('beads.rangeMin', 2.7);
TASBEConfig.set('beads.peakThreshold', 600);
% Execute and save the model
CM=resolve(CM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CM:
CMS = struct(CM);
assertTrue(strcmp(CMS.standardUnits,'MEPTR'));
UT = struct(CMS.unit_translation);
assertElementsAlmostEqual(UT.k_ERF,        64.5559,  'relative', 1e-2);
assertElementsAlmostEqual(UT.first_peak,    7);
assertElementsAlmostEqual(UT.fit_error,     0.00,   'absolute', 0.002);
assertElementsAlmostEqual(UT.peak_sets{1},  [855.4849 2.4685e+03], 'relative', 1e-2);


function test_onepeak
TASBEConfig.checkpoint('test');

[CM] = setupRedPeakCM();
% Ignore all bead data below 10^rangeMin as being too "smeared" with noise
TASBEConfig.set('beads.rangeMin', 3.1);
TASBEConfig.set('beads.peakThreshold', 600);
% Execute and save the model
CM=resolve(CM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CM:
CMS = struct(CM);
assertTrue(strcmp(CMS.standardUnits,'MEPTR'));
UT = struct(CMS.unit_translation);
assertElementsAlmostEqual(UT.k_ERF,        68.971,  'relative', 1e-2);
assertElementsAlmostEqual(UT.first_peak,    8);
assertElementsAlmostEqual(UT.fit_error,     0.00,   'absolute', 0.002);
assertElementsAlmostEqual(UT.peak_sets{1},  [2.4685e+03], 'relative', 1e-2);


function test_toomanypeaks
TASBEConfig.checkpoint('test');

[CM] = setupRedPeakCM();
 % set threshold and min too low so that it should sweep up lots of noise, get too many peaks
TASBEConfig.set('beads.rangeMin', 1);
TASBEConfig.set('beads.peakThreshold', 300);
% Execute and save the model
CM=resolve(CM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CM:
CMS = struct(CM);
assertTrue(strcmp(CMS.standardUnits,'MEPTR'));
UT = struct(CMS.unit_translation);
assertElementsAlmostEqual(UT.k_ERF,        11.2510,  'relative', 1e-2);
assertElementsAlmostEqual(UT.first_peak,    2);
assertElementsAlmostEqual(UT.fit_error,     0.4658,   'absolute', 0.002);
expected_peaks = 1e3 .* [0.0104    0.0114    0.0123    0.0138    0.0175    0.0372    0.1095    0.2280    0.8523 1.3515    2.4685    3.8884];
assertElementsAlmostEqual(UT.peak_sets{1},  expected_peaks, 'absolute', 1);


function test_nopeaks
TASBEConfig.checkpoint('test');

[CM] = setupRedPeakCM();
TASBEConfig.set('beads.peakThreshold', 1e7); % set too high to see anything
% Execute and save the model
CM=resolve(CM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CM:
CMS = struct(CM);
assertTrue(strcmp(CMS.standardUnits,'arbitrary units'));
UT = struct(CMS.unit_translation);
assertElementsAlmostEqual(UT.k_ERF,        1);
assertElementsAlmostEqual(UT.first_peak,    NaN);
assertTrue(isinf(UT.fit_error));
assertTrue(isempty(UT.peak_sets{1}));


function [CM] = setupBV421CM()

beadfile = DataFile('fcs', '../TASBEFlowAnalytics-Tutorial/example_controls/171221_E1_p1_AJ02.fcs');
blankfile = [];

% Create one channel / colorfile pair for each color
channels = {}; colorfiles = {};
channels{1} = Channel('VL1-A', 405, 440, 50);
channels{1} = setPrintName(channels{1}, 'BV421');
channels{1} = setLineSpec(channels{1}, 'b');
colorfiles{1} = [];

% Multi-color controls are used for converting other colors into ERF units
% Any channel without a control mapping it to ERF will be left in arbirary units.
colorpairfiles = {};

CM = ColorModel(beadfile, blankfile, channels, colorfiles, colorpairfiles);

TASBEConfig.set('beads.beadModel','SpheroTech URCP-38-2K'); % Entry from BeadCatalog.xls matching your beads
TASBEConfig.set('beads.beadBatch','Lot AJ02'); % Entry from BeadCatalog.xls containing your lot
TASBEConfig.set('beads.beadChannel','BV421');

CM=set_ERF_channel_name(CM, 'VL1-A');


function test_rightpeaks
TASBEConfig.checkpoint('test');

[CM] = setupBV421CM();
% Execute and save the model
TASBEConfig.set('plots.plotPath', '/tmp/plots');
TASBEConfig.set('calibration.overrideAutofluorescence',true);
CM=resolve(CM);

% Reset TASBEConfig to not contaminate other tests.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CM:
CMS = struct(CM);
assertTrue(strcmp(CMS.standardUnits,'MEBV421'));
UT = struct(CMS.unit_translation);
assertElementsAlmostEqual(UT.k_ERF,       0.2966, 'relative', 1e-2);
assertElementsAlmostEqual(UT.first_peak,    2);
assertElementsAlmostEqual(UT.fit_error, 0.0363,   'absolute', 0.002);
expected_peaks = 1e5 .* [0.0100    0.0689    0.2023    0.5471    1.5223];
assertElementsAlmostEqual(UT.peak_sets{1},  expected_peaks, 'relative', 1e-2);


function test_forcepeaks
TASBEConfig.checkpoint('test');

[CM] = setupRedPeakCM();
% Ignore all bead data below 10^rangeMin as being too "smeared" with noise
TASBEConfig.set('beads.rangeMin', 2.7);
TASBEConfig.set('beads.peakThreshold', 600);
% Execute and save the model
TASBEConfig.set('beads.forceFirstPeak',3);
CM=resolve(CM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CM:
CMS = struct(CM);
assertTrue(strcmp(CMS.standardUnits,'MEPTR'));
UT = struct(CMS.unit_translation);
assertElementsAlmostEqual(UT.k_ERF,        0.8833,  'relative', 1e-2);
assertElementsAlmostEqual(UT.first_peak,    3);
assertElementsAlmostEqual(UT.fit_error,     0.00,   'absolute', 0.002);
assertElementsAlmostEqual(UT.peak_sets{1},  [855.4849 2.4685e+03], 'relative', 1e-2);



function [CM] = setupYellowPeakCM()

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';

beadfile = DataFile('fcs', [stem0312 'Beads_P3.fcs']);
blankfile = DataFile('fcs', [stem0312 'blank_P3.fcs']);

% Create one channel / colorfile pair for each color
channels = {}; colorfiles = {};
channels{1} = Channel('PE-Tx-Red-YG-A', 561, 610, 20);
channels{1} = setPrintName(channels{1}, 'mKate');
channels{1} = setLineSpec(channels{1}, 'r');
colorfiles{1} = DataFile('fcs', [stem0312 'mkate_P3.fcs']);

channels{2} = Channel('FITC-A', 488, 515, 20);
channels{2} = setPrintName(channels{2}, 'EYFP');
channels{2} = setLineSpec(channels{2}, 'r');
colorfiles{2} = DataFile('fcs', [stem0312 'EYFP_P3.fcs']);

% Multi-color controls are used for converting other colors into ERF units
% Any channel without a control mapping it to ERF will be left in arbirary units.
colorpairfiles = {};

CM = ColorModel(beadfile, blankfile, channels, colorfiles, colorpairfiles);

TASBEConfig.set('beads.beadModel','SpheroTech RCP-30-5A'); % Entry from BeadCatalog.xls matching your beads
TASBEConfig.set('beads.beadBatch','Lot AA01, AA02, AA03, AA04, AB01, AB02, AC01, GAA01-R'); % Entry from BeadCatalog.xls containing your lot
TASBEConfig.set('beads.beadChannel','FITC');
CM=set_ERF_channel_name(CM, 'FITC-A');

TASBEConfig.set('plots.plotPath', '/tmp/plots');

function test_secondarypeaks
TASBEConfig.checkpoint('test');

[CM] = setupYellowPeakCM();
% Ignore all bead data below 10^rangeMin as being too "smeared" with noise
TASBEConfig.set('beads.rangeMin', 2.7);
TASBEConfig.set('beads.peakThreshold', 600);
% Execute and save the model
TASBEConfig.set('beads.secondaryBeadChannel','PE-Tx-Red-YG-A');
CM=resolve(CM);
TASBEConfig.clear('beads.secondaryBeadChannel');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CM:
CMS = struct(CM);
assertTrue(strcmp(CMS.standardUnits,'MEFL'));
UT = struct(CMS.unit_translation);
assertElementsAlmostEqual(UT.k_ERF,        2.3166e+03,  'relative', 1e-2);
assertElementsAlmostEqual(UT.first_peak,    7);
assertElementsAlmostEqual(UT.fit_error,     0.00,   'absolute', 0.002);
assertElementsAlmostEqual(UT.peak_sets{1},  [0.8555e+03    2.4685e+03], 'relative', 1e-2);
assertElementsAlmostEqual(UT.peak_sets{2},  [53.6450  128.2913], 'relative', 1e-2);

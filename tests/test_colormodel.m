function test_suite = test_colormodel
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_colormodel_endtoend
TASBEConfig.checkpoint('test');

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';

beadfile = [stem0312 'Beads_P3.fcs'];
blankfile = [stem0312 'blank_P3.fcs'];

% Autodetect gating with an N-dimensional gaussian-mixture-model
autogate = GMMGating(blankfile);

% Create one channel / colorfile pair for each color
channels = {}; colorfiles = {};
% Channel takes FCS channel name, laser frequency (nm), filter center (nm), filter width (nm)
% Do not duplicate laser/filter information, as this may cause analysis collisions
channels{1} = Channel('FITC-A', 488, 515, 20);
channels{1} = setPrintName(channels{1}, 'EYFP'); % Name to print on charts
channels{1} = setLineSpec(channels{1}, 'y'); % Color for lines, when needed
colorfiles{1} = [stem0312 'EYFP_P3.fcs'];

channels{2} = Channel('PE-Tx-Red-YG-A', 561, 610, 20);
channels{2} = setPrintName(channels{2}, 'mKate');
channels{2} = setLineSpec(channels{2}, 'r');
colorfiles{2} = [stem0312 'mkate_P3.fcs'];

channels{3} = Channel('Pacific Blue-A', 405, 450, 50);
channels{3} = setPrintName(channels{3}, 'EBFP2');
channels{3} = setLineSpec(channels{3}, 'b');
colorfiles{3} = [stem0312 'ebfp2_P3.fcs'];

% Multi-color controls are used for converting other colors into ERF units
% Any channel without a control mapping it to ERF will be left in arbirary units.
colorpairfiles = {};
% Entries are: channel1, channel2, constitutive channel, filename
% This allows channel1 and channel2 to be converted into one another.
% If you only have two colors, you can set consitutive-channel to equal channel1 or channel2
colorpairfiles{1} = {channels{1}, channels{2}, channels{3}, [stem0312 'mkate_EBFP2_EYFP_P3.fcs']};
colorpairfiles{2} = {channels{1}, channels{3}, channels{2}, [stem0312 'mkate_EBFP2_EYFP_P3.fcs']};

CM = ColorModel(beadfile, blankfile, channels, colorfiles, colorpairfiles);

TASBEConfig.set('beads.beadModel','SpheroTech RCP-30-5A'); % Entry from BeadCatalog.xls matching your beads
TASBEConfig.set('beads.beadBatch','Lot AA01, AA02, AA03, AA04, AB01, AB02, AC01, GAA01-R'); % Entry from BeadCatalog.xls containing your lot

% Ignore all bead data below 10^rangeMin as being too "smeared" with noise
TASBEConfig.set('beads.rangeMin', 2);
% The peak threshold determines the minumum count per bin for something to
% be considered part of a peak.  Set if automated threshold finds too many or few peaks
%TASBEConfig.set('beads.peakThreshold', 200);
CM=set_ERF_channel_name(CM, 'FITC-A');
% Ignore channel data for ith channel if below 10^[value(i)]
TASBEConfig.set('colortranslation.channelMinimum',[2,2,2]);

TASBEConfig.set('plots.plotPath', '/tmp/plots');
CM = add_prefilter(CM,autogate);

% Execute and save the model
CM=resolve(CM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CM:
CMS = struct(CM);

channels = getChannels(CM);
assertEqual(getUnits(channels{1}),'MEFL');
assertEqual(getUnits(channels{2}),'MEFL');
assertEqual(getUnits(channels{3}),'MEFL');

UT = struct(CMS.unit_translation);
assertElementsAlmostEqual(UT.k_ERF,        2267.3,   'relative', 1e-2);
assertElementsAlmostEqual(UT.first_peak,    8);
assertElementsAlmostEqual(UT.fit_error,     0);
assertElementsAlmostEqual(UT.peak_sets{1},  [128.35], 'relative', 1e-2);

AFM_Y = struct(CMS.autofluorescence_model{1});
assertElementsAlmostEqual(AFM_Y.af_mean,    3.2600,  'absolute', 0.5);
assertElementsAlmostEqual(AFM_Y.af_std,     17.0788, 'absolute', 0.5);
AFM_R = struct(CMS.autofluorescence_model{2});
assertElementsAlmostEqual(AFM_R.af_mean,    3.6683,  'absolute', 0.5);
assertElementsAlmostEqual(AFM_R.af_std,     17.5621, 'absolute', 0.5);
AFM_B = struct(CMS.autofluorescence_model{3});
assertElementsAlmostEqual(AFM_B.af_mean,    5.8697,  'absolute', 0.5);
assertElementsAlmostEqual(AFM_B.af_std,     16.9709, 'absolute', 0.5);

COMP = struct(CMS.compensation_model);
expected_matrix = [...
    1.0000      0.0056      0.0004;
    0.0010      1.0000      0.0022;
         0      0.0006      1.0000];

assertElementsAlmostEqual(COMP.matrix,      expected_matrix, 'absolute', 1e-3);

CTM = struct(CMS.color_translation_model);
expected_scales = [...
       NaN    1.0097    2.1796;
    0.9872       NaN       NaN;
    0.45834      NaN       NaN];

assertElementsAlmostEqual(CTM.scales,       expected_scales, 'absolute', 0.02);



function test_colormodel_warnings
TASBEConfig.checkpoint('test');

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';

beadfile = DataFile([stem0312 'Beads_P3.fcs']);
blankfile = DataFile([stem0312 'blank_P3.fcs']);

% Autodetect gating with an N-dimensional gaussian-mixture-model
autogate = GMMGating(blankfile);

% Create one channel / colorfile pair for each color
channels = {}; colorfiles = {};
% Channel takes FCS channel name, laser frequency (nm), filter center (nm), filter width (nm)
% Do not duplicate laser/filter information, as this may cause analysis collisions
channels{1} = Channel('FITC-A', 488, 515, 20);
channels{1} = setPrintName(channels{1}, 'EYFP'); % Name to print on charts
channels{1} = setLineSpec(channels{1}, 'y'); % Color for lines, when needed
colorfiles{1} = DataFile([stem0312 'EYFP_P3.fcs']);

channels{2} = Channel('PE-Tx-Red-YG-A', 561, 610, 20);
channels{2} = setPrintName(channels{2}, 'mKate');
channels{2} = setLineSpec(channels{2}, 'r');
colorfiles{2} = DataFile([stem0312 'mkate_P3.fcs']);

channels{3} = Channel('Pacific Blue-A', 405, 450, 50);
channels{3} = setPrintName(channels{3}, 'EBFP2');
channels{3} = setLineSpec(channels{3}, 'b');
colorfiles{3} = DataFile([stem0312 'ebfp2_P3.fcs']);

% Multi-color controls are used for converting other colors into ERF units
% Any channel without a control mapping it to ERF will be left in arbirary units.
colorpairfiles = {};
% Entries are: channel1, channel2, constitutive channel, filename
% This allows channel1 and channel2 to be converted into one another.
% If you only have two colors, you can set consitutive-channel to equal channel1 or channel2
colorpairfiles{1} = {channels{1}, channels{2}, channels{3}, DataFile([stem0312 'mkate_EBFP2_EYFP_P3.fcs'])};
colorpairfiles{2} = {channels{1}, channels{3}, channels{2}, DataFile([stem0312 'mkate_EBFP2_EYFP_P3.fcs'])};

CM = ColorModel(beadfile, blankfile, channels, colorfiles, colorpairfiles);

TASBEConfig.set('beads.beadModel','SpheroTech RCP-30-5A'); % Entry from BeadCatalog.xls matching your beads
TASBEConfig.set('beads.beadBatch','Lot AA01, AA02, AA03, AA04, AB01, AB02, AC01, GAA01-R'); % Entry from BeadCatalog.xls containing your lot

% Ignore all bead data below 10^rangeMin as being too "smeared" with noise
TASBEConfig.set('beads.rangeMin', 2);
% The peak threshold determines the minumum count per bin for something to
% be considered part of a peak.  Set if automated threshold finds too many or few peaks
%TASBEConfig.set('beads.peakThreshold', 200);
CM=set_ERF_channel_name(CM, 'FITC-A');
% Ignore channel data for ith channel if below 10^[value(i)]
TASBEConfig.set('colortranslation.channelMinimum',[2,2,2]);

TASBEConfig.set('plots.plotPath', '/tmp/plots');
CM = add_prefilter(CM,autogate);

TASBESession.reset();
TASBEConfig.set('beads.peakThreshold', 100);
TASBEConfig.set('beads.rangeMin', 1);
CM = resolve(CM);
log = TASBESession.list();
assertEqual(log{end-1}.contents{3}.name, 'PeaksNotAscending');
assertEqual(log{end-1}.contents{4}.name, 'HighestPeakUndetected');
assertEqual(log{end-1}.contents{5}.name, 'QuestionableFirstPeak');
assertEqual(log{end-1}.contents{7}.name, 'PotentialBeadClump');



function test_colormodel_singlered
TASBEConfig.checkpoint('test');

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';

beadfile = [stem0312 'Beads_P3.fcs'];
blankfile = [stem0312 'blank_P3.fcs'];

% Create one channel / colorfile pair for each color
channels = {}; colorfiles = {};
channels{1} = Channel('PE-Tx-Red-YG-A', 561, 610, 20);
channels{1} = setPrintName(channels{1}, 'mKate');
channels{1} = setLineSpec(channels{1}, 'r');
colorfiles{1} = [stem0312 'mkate_P3.fcs'];

% Multi-color controls are used for converting other colors into ERF units
% Any channel without a control mapping it to ERF will be left in arbirary units.
colorpairfiles = {};

CM = ColorModel(beadfile, blankfile, channels, colorfiles, colorpairfiles);

TASBEConfig.set('beads.beadModel','SpheroTech RCP-30-5A'); % Entry from BeadCatalog.xls matching your beads
TASBEConfig.set('beads.beadBatch','Lot AA01, AA02, AA03, AA04, AB01, AB02, AC01, GAA01-R'); % Entry from BeadCatalog.xls containing your lot
TASBEConfig.set('beads.beadChannel','PE-TR');

% Ignore all bead data below 10^rangeMin as being too "smeared" with noise
TASBEConfig.set('beads.rangeMin', 1.8);
% The peak threshold determines the minumum count per bin for something to
% be considered part of a peak.  Set if automated threshold finds too many or few peaks
TASBEConfig.set('beads.peakThreshold', 600);
CM=set_ERF_channel_name(CM, 'PE-Tx-Red-YG-A');

TASBEConfig.set('plots.plotPath', '/tmp/plots');
% Execute and save the model
CM=resolve(CM);
save('-V7','/tmp/CM120312.mat','CM');

assertEqual(TASBESession.getLast('TASBE:Beads','NonMEFL').message,'MEFL units are recommended, rather than MEPTR');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CM:
CMS = struct(CM);

UT = struct(CMS.unit_translation);
assertElementsAlmostEqual(UT.k_ERF,        59.9971,  'relative', 1e-2);
assertElementsAlmostEqual(UT.first_peak,    5);
assertElementsAlmostEqual(UT.fit_error,     0.019232,   'absolute', 0.002);
assertElementsAlmostEqual(UT.peak_sets{1},  [110.3949 227.5807 855.4849 2.4685e+03], 'relative', 1e-2);

AFM_R = struct(CMS.autofluorescence_model{1});
assertElementsAlmostEqual(AFM_R.af_mean,    3.6683,  'absolute', 0.5);
assertElementsAlmostEqual(AFM_R.af_std,     17.5621, 'absolute', 0.5);

COMP = struct(CMS.compensation_model);
assertElementsAlmostEqual(COMP.matrix,      1.0000, 'absolute', 1e-3);

CTM = struct(CMS.color_translation_model);
assertElementsAlmostEqual(CTM.scales,   NaN);



function test_colormodel_singlered_nocolorfile
TASBEConfig.checkpoint('test');

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';

beadfile = DataFile('fcs', [stem0312 'Beads_P3.fcs']);
blankfile = DataFile('fcs', [stem0312 'blank_P3.fcs']);

% Create one channel / colorfile pair for each color
channels = {}; colorfiles = {};
channels{1} = Channel('PE-Tx-Red-YG-A', 561, 610, 20);
channels{1} = setPrintName(channels{1}, 'mKate');
channels{1} = setLineSpec(channels{1}, 'r');

% Multi-color controls are used for converting other colors into ERF units
% Any channel without a control mapping it to ERF will be left in arbirary units.
colorpairfiles = {};

CM = ColorModel(beadfile, blankfile, channels, colorfiles, colorpairfiles);

TASBEConfig.set('beads.beadModel','SpheroTech RCP-30-5A'); % Entry from BeadCatalog.xls matching your beads
TASBEConfig.set('beads.beadBatch','Lot AA01, AA02, AA03, AA04, AB01, AB02, AC01, GAA01-R'); % Entry from BeadCatalog.xls containing your lot
TASBEConfig.set('beads.beadChannel','PE-TR');

% Ignore all bead data below 10^rangeMin as being too "smeared" with noise
TASBEConfig.set('beads.rangeMin', 1.8);
% The peak threshold determines the minumum count per bin for something to
% be considered part of a peak.  Set if automated threshold finds too many or few peaks
TASBEConfig.set('beads.peakThreshold', 600);
CM=set_ERF_channel_name(CM, 'PE-Tx-Red-YG-A');

TASBEConfig.set('plots.plotPath', '/tmp/plots');
% Execute and save the model
CM=resolve(CM);
save('-V7','/tmp/CM120312.mat','CM');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CM:
CMS = struct(CM);

UT = struct(CMS.unit_translation);
assertElementsAlmostEqual(UT.k_ERF,        59.9971,  'relative', 1e-2);
assertElementsAlmostEqual(UT.first_peak,    5);
assertElementsAlmostEqual(UT.fit_error,     0.019232,   'absolute', 0.002);
assertElementsAlmostEqual(UT.peak_sets{1},  [110.3949 227.5807 855.4849 2.4685e+03], 'relative', 1e-2);

AFM_R = struct(CMS.autofluorescence_model{1});
assertElementsAlmostEqual(AFM_R.af_mean,    3.6683,  'absolute', 0.5);
assertElementsAlmostEqual(AFM_R.af_std,     17.5621, 'absolute', 0.5);

COMP = struct(CMS.compensation_model);
assertElementsAlmostEqual(COMP.matrix,      1.0000, 'absolute', 1e-3);

CTM = struct(CMS.color_translation_model);
assertElementsAlmostEqual(CTM.scales,   NaN);



function test_colormodel_fsc_ssc
TASBEConfig.checkpoint('test');

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';

beadfile = DataFile('fcs', [stem0312 'Beads_P3.fcs']);
blankfile = DataFile('fcs', [stem0312 'blank_P3.fcs']);

% Create one channel / colorfile pair for each color
channels = {}; colorfiles = {};
channels{1} = Channel('PE-Tx-Red-YG-A', 561, 610, 20);
channels{1} = setPrintName(channels{1}, 'mKate');
channels{1} = setLineSpec(channels{1}, 'r');
colorfiles{1} = DataFile('fcs', [stem0312 'mkate_P3.fcs']);

channels{2} = Channel('FSC-A', 488, 488, 10);
channels{2} = setPrintName(channels{2}, 'FSC');
channels{2} = setLineSpec(channels{2}, 'k');

channels{3} = Channel('SSC-A', 488, 488, 10);
channels{3} = setPrintName(channels{3}, 'SSC');
channels{3} = setLineSpec(channels{3}, 'b');

% Multi-color controls are used for converting other colors into ERF units
% Any channel without a control mapping it to ERF will be left in arbirary units.
colorpairfiles = {};

CM = ColorModel(beadfile, blankfile, channels, colorfiles, colorpairfiles);

TASBEConfig.set('beads.beadModel','SpheroTech RCP-30-5A'); % Entry from BeadCatalog.xls matching your beads
TASBEConfig.set('beads.beadBatch','Lot AA01, AA02, AA03, AA04, AB01, AB02, AC01, GAA01-R'); % Entry from BeadCatalog.xls containing your lot
TASBEConfig.set('beads.beadChannel','PE-TR');

% Ignore all bead data below 10^rangeMin as being too "smeared" with noise
TASBEConfig.set('beads.rangeMin', 1.8);
% The peak threshold determines the minumum count per bin for something to
% be considered part of a peak.  Set if automated threshold finds too many or few peaks
TASBEConfig.set('beads.peakThreshold', 600);
CM=set_ERF_channel_name(CM, 'PE-Tx-Red-YG-A');

TASBEConfig.set('plots.plotPath', '/tmp/plots');
% Execute and save the model
CM=resolve(CM);
save('-V7','/tmp/CM120312.mat','CM');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CM:
CMS = struct(CM);

UT = struct(CMS.unit_translation);
assertElementsAlmostEqual(UT.k_ERF,        59.9971,  'relative', 1e-2);
assertElementsAlmostEqual(UT.first_peak,    5);
assertElementsAlmostEqual(UT.fit_error,     0.019232,   'absolute', 0.002);
assertElementsAlmostEqual(UT.peak_sets{1},  [110.3949 227.5807 855.4849 2.4685e+03], 'relative', 1e-2);

assertEqual(numel(CMS.autofluorescence_model), 3);
assertEqual(sum(~isempty(CMS.autofluorescence_model)), 1);

AFM_R = struct(CMS.autofluorescence_model{1});
assertElementsAlmostEqual(AFM_R.af_mean,    3.6683,  'absolute', 0.5);
assertElementsAlmostEqual(AFM_R.af_std,     17.5621, 'absolute', 0.5);

COMP = struct(CMS.compensation_model);
assertElementsAlmostEqual(COMP.matrix,      [1 0 0; 0 1 0; 0 0 1], 'absolute', 1e-3);

CTM = struct(CMS.color_translation_model);
assertElementsAlmostEqual(CTM.scales,   [NaN NaN NaN; NaN NaN NaN; NaN NaN NaN]);



function test_colormodel_fsc_ssc_nocolorfile
TASBEConfig.checkpoint('test');

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';

beadfile = DataFile('fcs', [stem0312 'Beads_P3.fcs']);
blankfile = DataFile('fcs', [stem0312 'blank_P3.fcs']);

% Create one channel / colorfile pair for each color
channels = {}; colorfiles = {};
channels{1} = Channel('PE-Tx-Red-YG-A', 561, 610, 20);
channels{1} = setPrintName(channels{1}, 'mKate');
channels{1} = setLineSpec(channels{1}, 'r');

channels{2} = Channel('FSC-A', 488, 488, 10);
channels{2} = setPrintName(channels{2}, 'FSC');
channels{2} = setLineSpec(channels{2}, 'k');

channels{3} = Channel('SSC-A', 488, 488, 10);
channels{3} = setPrintName(channels{3}, 'SSC');
channels{3} = setLineSpec(channels{3}, 'b');

% Multi-color controls are used for converting other colors into ERF units
% Any channel without a control mapping it to ERF will be left in arbirary units.
colorpairfiles = {};

CM = ColorModel(beadfile, blankfile, channels, colorfiles, colorpairfiles);

TASBEConfig.set('beads.beadModel','SpheroTech RCP-30-5A'); % Entry from BeadCatalog.xls matching your beads
TASBEConfig.set('beads.beadBatch','Lot AA01, AA02, AA03, AA04, AB01, AB02, AC01, GAA01-R'); % Entry from BeadCatalog.xls containing your lot
TASBEConfig.set('beads.beadChannel','PE-TR');

% Ignore all bead data below 10^rangeMin as being too "smeared" with noise
TASBEConfig.set('beads.rangeMin', 1.8);
% The peak threshold determines the minumum count per bin for something to
% be considered part of a peak.  Set if automated threshold finds too many or few peaks
TASBEConfig.set('beads.peakThreshold', 600);
CM=set_ERF_channel_name(CM, 'PE-Tx-Red-YG-A');

TASBEConfig.set('plots.plotPath', '/tmp/plots');
% Execute and save the model
CM=resolve(CM);
save('-V7','/tmp/CM120312.mat','CM');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CM:
CMS = struct(CM);

UT = struct(CMS.unit_translation);
assertElementsAlmostEqual(UT.k_ERF,        59.9971,  'relative', 1e-2);
assertElementsAlmostEqual(UT.first_peak,    5);
assertElementsAlmostEqual(UT.fit_error,     0.019232,   'absolute', 0.002);
assertElementsAlmostEqual(UT.peak_sets{1},  [110.3949 227.5807 855.4849 2.4685e+03], 'relative', 1e-2);

assertEqual(numel(CMS.autofluorescence_model), 3);
assertEqual(sum(~isempty(CMS.autofluorescence_model)), 1);

AFM_R = struct(CMS.autofluorescence_model{1});
assertElementsAlmostEqual(AFM_R.af_mean,    3.6683,  'absolute', 0.5);
assertElementsAlmostEqual(AFM_R.af_std,     17.5621, 'absolute', 0.5);

COMP = struct(CMS.compensation_model);
assertElementsAlmostEqual(COMP.matrix,      [1 0 0; 0 1 0; 0 0 1], 'absolute', 1e-3);

CTM = struct(CMS.color_translation_model);
assertElementsAlmostEqual(CTM.scales,   [NaN NaN NaN; NaN NaN NaN; NaN NaN NaN]);



function test_colormodel_error_missing_colorfiles
TASBEConfig.checkpoint('test');

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';

beadfile = DataFile('fcs', [stem0312 'Beads_P3.fcs']);
blankfile = DataFile('fcs', [stem0312 'blank_P3.fcs']);

% Create one channel / colorfile pair for each color
channels = {}; colorfiles = {};
channels{1} = Channel('FITC-A', 488, 515, 20);
channels{1} = setPrintName(channels{1}, 'EYFP'); % Name to print on charts
channels{1} = setLineSpec(channels{1}, 'y'); % Color for lines, when needed

channels{2} = Channel('PE-Tx-Red-YG-A', 561, 610, 20);
channels{2} = setPrintName(channels{2}, 'mKate');
channels{2} = setLineSpec(channels{2}, 'r');
colorfiles{2} = DataFile('fcs', [stem0312 'mkate_P3.fcs']);

% Multi-color controls are used for converting other colors into ERF units
% Any channel without a control mapping it to ERF will be left in arbirary units.
colorpairfiles = {};

assertExceptionThrown(@()(ColorModel(beadfile, blankfile, channels, colorfiles, colorpairfiles)),...
    'TASBE:ColorModel:OneColorfilePerChannel', 'No error was raised');


% Make sure manual setting of template file converts to datafile correctly
function test_colormodel_channeltemplate
TASBEConfig.checkpoint('test');

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';

beadfile = [stem0312 'Beads_P3.fcs'];
blankfile = [stem0312 'blank_P3.fcs'];
TASBEConfig.set('flow.channel_template_file',beadfile);

% Create one channel / colorfile pair for each color
channels = {}; colorfiles = {};
channels{1} = Channel('PE-Tx-Red-YG-A', 561, 610, 20);
channels{1} = setPrintName(channels{1}, 'mKate');
channels{1} = setLineSpec(channels{1}, 'r');
colorfiles{1} = [stem0312 'mkate_P3.fcs'];

% Multi-color controls are used for converting other colors into ERF units
% Any channel without a control mapping it to ERF will be left in arbirary units.
colorpairfiles = {};

CM = ColorModel(beadfile, blankfile, channels, colorfiles, colorpairfiles);

TASBEConfig.set('beads.beadModel','SpheroTech RCP-30-5A'); % Entry from BeadCatalog.xls matching your beads
TASBEConfig.set('beads.beadBatch','Lot AA01, AA02, AA03, AA04, AB01, AB02, AC01, GAA01-R'); % Entry from BeadCatalog.xls containing your lot
TASBEConfig.set('beads.beadChannel','PE-TR');
CM=set_ERF_channel_name(CM, 'PE-Tx-Red-YG-A');

TASBEConfig.set('plots.plotPath', '/tmp/plots');
% Execute and save the model
CM=resolve(CM);
save('-V7','/tmp/CM120312.mat','CM');


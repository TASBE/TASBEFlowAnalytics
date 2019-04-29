function test_suite = test_sizebeads
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

%%%%%%%%%%%%
% Note: test_colormodel tests the one and many bead cast; we just need to test a few special cases

function [CM] = setupSizePeakCM()

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';

beadfile = DataFile('fcs',[stem0312 'Beads_P3.fcs']);
blankfile = DataFile('fcs',[stem0312 'blank_P3.fcs']);

% Create one channel / colorfile pair for each color
channels = {}; colorfiles = {};
channels{1} = Channel('FITC-A', 488, 515, 20);
channels{1} = setPrintName(channels{1}, 'EYFP'); % Name to print on charts
channels{1} = setLineSpec(channels{1}, 'y'); % Color for lines, when needed
colorfiles{1} = DataFile('fcs',[stem0312 'EYFP_P3.fcs']); % If there is only one channel, the color file is optional

channels{2} = Channel('Pacific Blue-A', 405, 450, 50);
channels{2} = setPrintName(channels{2}, 'EBFP2');
channels{2} = setLineSpec(channels{2}, 'b');
colorfiles{2} = DataFile('fcs',[stem0312 'ebfp2_P3.fcs']);

channels{3} = Channel('FSC-A', 488, 488, 10);
channels{3} = setPrintName(channels{3}, 'FSC');
channels{3} = setLineSpec(channels{3}, 'k');

channels{4} = Channel('SSC-A', 488, 488, 10);
channels{4} = setPrintName(channels{4}, 'SSC');
channels{4} = setLineSpec(channels{4}, 'r');

% Multi-color controls are used for converting other colors into ERF units
% Any channel without a control mapping it to ERF will be left in arbirary units.
colorpairfiles = {};
colorpairfiles{1} = {channels{1}, channels{2}, channels{2}, DataFile('fcs',[stem0312 'mkate_EBFP2_EYFP_P3.fcs'])};

sizebeadfile = '../TASBEFlowAnalytics-Tutorial/example_controls/180614_PPS6K_A02.fcs';
sizedatafile = DataFile('fcs',sizebeadfile);

CM = ColorModel(beadfile, blankfile, channels, colorfiles, colorpairfiles, sizedatafile);

TASBEConfig.set('beads.beadModel','SpheroTech RCP-30-5A'); % Entry from BeadCatalog.xls matching your beads
TASBEConfig.set('beads.beadBatch','Lot AA01, AA02, AA03, AA04, AB01, AB02, AC01, GAA01-R'); % Entry from BeadCatalog.xls containing your lot

CM=set_ERF_channel_name(CM, 'FITC-A');

% Configuration for size beads, if used
TASBEConfig.set('sizebeads.beadModel','SpheroTech PPS-6K'); % Entry from BeadCatalog.xls matching your beads
% Can also set bead channel or batch, if alternatives are available
% Ignore all size bead data below 10^[rangeMin] as being too "smeared" with noise
TASBEConfig.set('sizebeads.rangeMin', 2);
CM=set_um_channel_name(CM, 'FSC-A');

TASBEConfig.set('plots.plotPath', '/tmp/plots');


function test_sizepeaks
TASBEConfig.checkpoint('test');

CM = setupSizePeakCM();
% Execute and save the model
CM=resolve(CM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CM:
CMS = struct(CM);
assertTrue(strcmp(CMS.sizeUnits,'Eum'));
expected_peaks = 1e5*[0.1257    0.2574    0.6506    1.3908    1.9087    3.1596];
UT = struct(CMS.size_unit_translation);
assertElementsAlmostEqual(UT.um_poly,       [0.5865 -2.0798],  'relative', 1e-2);
assertElementsAlmostEqual(UT.first_peak,    1);
assertElementsAlmostEqual(UT.fit_error,     0.0784,   'absolute', 0.01);
assertElementsAlmostEqual(UT.peak_sets{1},  expected_peaks, 'relative', 1e-2);

channels = getChannels(CM);
% make sure units are right
assertEqual(getUnits(channels{1}),'MEFL');
assertEqual(getUnits(channels{2}),'MEFL');
assertEqual(getUnits(channels{3}),'Eum');
assertEqual(getUnits(channels{4}),'a.u.');
% make sure translations are right
assertElementsAlmostEqual(au_to_ERF(CM,getChannel(CM,1),1),2267.3, 'relative', 1e-2);
assertElementsAlmostEqual(au_to_ERF(CM,getChannel(CM,2),1),1163.3, 'relative', 1e-2);
assertElementsAlmostEqual(au_to_ERF(CM,getChannel(CM,3),1),0.0083, 'relative', 1e-2);
assertElementsAlmostEqual(au_to_ERF(CM,getChannel(CM,4),1),1, 'relative', 1e-2);


function test_size_bead_reading

CM = setupSizePeakCM();
% Execute and save the model
CM=resolve(CM);

% make sure size channel isn't messed up by PEM drop
stem1011 = '../TASBEFlowAnalytics-Tutorial/example_assay/LacI-CAGop_';
bins = BinSequence(-3,0.1,10,'log_bins');

AP = AnalysisParameters(bins,{});
AP=setMinValidCount(AP,100');
AP=setPemDropThreshold(AP,5');
AP=setUseAutoFluorescence(AP,false');

% Make a map of condition names to file sets
file_pairs = {...
  'Dox 0.1',    {DataFile('fcs',[stem1011 'B3_P3.fcs'])};
  'Dox 2000.0', {DataFile('fcs',[stem1011 'C4_P3.fcs'])};
    };

[results, sampleresults] = per_color_constitutive_analysis(CM,file_pairs,{'EBFP2','EYFP','FSC','SSC'},AP);

% Make output plots
TASBEConfig.set('OutputSettings.StemName','LacI-CAGop');
TASBEConfig.set('plots.plotPath','/tmp/plots');
TASBEConfig.set('OutputSettings.FixedInputAxis',[1e4 1e10]);
plot_batch_histograms(results,sampleresults,CM,{'b','g','r','k'});

save('/tmp/size-batch.mat','AP','bins','file_pairs','results','sampleresults');

%%%%%%%%%%%%%%%%%%%%%%
% Run all comparisons

expectedBinCounts = [...
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0         351           0
           0           0         988           0
           0           0        1369           0
           0           0        1558           0
           0           0        1287           0
           0           0        1808           0
           0           0        5069           0
           0           0       16126           0
           0           0       54323           0
           0           0      113083           0
           0           0       24407           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0         362
           0           0           0        2039
           0           0           0        6089
           0           0           0       10663
           0           0           0       19672
           0           0           0       33667
           0           0           0       43173
           0           0           0       43585
           0           0           0       32119
           0           0           0       17692
        6920           0           0        7481
        6645           0           0        2635
        9624           0           0         846
        8838        4595           0         242
       10886        4695           0           0
       13706        6742           0           0
       10660        6478           0           0
       11065        7951           0           0
        7683        8270           0           0
        4219        8006           0           0
        1843        7968           0           0
         568        6619           0           0
         135        4127           0           0
           0        3371           0           0
           0        3108           0           0
           0        3519           0           0
           0        3672           0           0
           0        4319           0           0
           0        4803           0           0
           0        5217           0           0
           0        5254           0           0
           0        5568           0           0
           0        5631           0           0
           0        5470           0           0
           0        5085           0           0
           0        4593           0           0
           0        4243           0           0
           0        3588           0           0
           0        2948           0           0
           0        2428           0           0
           0        1989           0           0
           0        1484           0           0
           0        1168           0           0
           0         782           0           0
           0         622           0           0
           0         387           0           0
           0         249           0           0
           0         191           0           0
           0         124           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
           0           0           0           0
    ];

expectedMeans = [...
    2.0661e+04 2.5131e+05 7.8537 3.0293e+03;
    1.7136e+05 1.4958e+05 8.2587 2.5810e+03;
    ];


for i=1:numel(results)
    assertElementsAlmostEqual(expectedMeans(i,:),results{i}.means,'relative',1e-2);
end

assertElementsAlmostEqual(expectedBinCounts,results{1}.bincounts,'relative',1e-2);


function test_size_peak_forcing

CM = setupSizePeakCM();
TASBEConfig.set('sizebeads.forceFirstPeak',1);
TASBEConfig.set('sizebeads.rangeMax', 5);
% Execute and save the model
CM=resolve(CM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in CM:
CMS = struct(CM);
assertTrue(strcmp(CMS.sizeUnits,'Eum'));
expected_peaks = 1e5*[0.1257    0.2574    0.6506];
UT = struct(CMS.size_unit_translation);
assertElementsAlmostEqual(UT.um_poly,       [0.5565 -1.9418],  'relative', 1e-2);
assertElementsAlmostEqual(UT.first_peak,    1);
assertElementsAlmostEqual(UT.fit_error,     0.0099,   'absolute', 0.01);
assertElementsAlmostEqual(UT.peak_sets{1},  expected_peaks, 'relative', 1e-2);

channels = getChannels(CM);
% make sure units are right
assertEqual(getUnits(channels{1}),'MEFL');
assertEqual(getUnits(channels{2}),'MEFL');
assertEqual(getUnits(channels{3}),'Eum');
assertEqual(getUnits(channels{4}),'a.u.');
% make sure translations are right
assertElementsAlmostEqual(au_to_ERF(CM,getChannel(CM,1),1),2267.3, 'relative', 1e-2);
assertElementsAlmostEqual(au_to_ERF(CM,getChannel(CM,2),1),1163.3, 'relative', 1e-2);
assertElementsAlmostEqual(au_to_ERF(CM,getChannel(CM,3),1),0.0114, 'relative', 1e-2);
assertElementsAlmostEqual(au_to_ERF(CM,getChannel(CM,4),1),1, 'relative', 1e-2);



function test_size_peak_multiple_required

CM = setupSizePeakCM();
TASBEConfig.set('sizebeads.rangeMin', 4.5);
TASBEConfig.set('sizebeads.rangeMax', 5);
% Execute and save the model
CM=resolve(CM);

log = TASBESession.list();
assertEqual(log{end-3}.contents{end-1}.name, 'SinglePeak');

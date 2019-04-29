function test_suite = test_autogating
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_autogate

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';
blankfile = DataFile('fcs', [stem0312 'blank_P3.fcs']);

% Autodetect gating with an N-dimensional gaussian-mixture-model
TASBEConfig.set('gating.channelNames',{'FSC-A','SSC-A'});
gate = GMMGating(blankfile);

gate = struct(gate);

assertEqual(gate.selected_components, 1);

expected_mu = [5.1117    3.4087;    4.5555    3.3427];
GDS = struct(gate.distribution);
assertElementsAlmostEqual(GDS.mu,expected_mu,'absolute',0.01);

expected_sigma(:,:,1) = [0.0130    0.0149;    0.0149    0.0328];
expected_sigma(:,:,2) = [0.2053    0.0548;    0.0548    0.0515];
assertElementsAlmostEqual(GDS.Sigma,expected_sigma,'absolute',0.01);


function test_6D_autogate

stem0312 = '../TASBEFlowAnalytics-Tutorial/example_controls/2012-03-12_';
blankfile = DataFile('fcs', [stem0312 'blank_P3.fcs']);

% Autodetect gating with an N-dimensional gaussian-mixture-model
TASBEConfig.set('gating.channelNames',{'FSC-A','SSC-A','FSC-H','SSC-H','FSC-W','SSC-W'});
gate = GMMGating(blankfile);

[~,fcshdr,data] = fca_readfcs([stem0312 'EYFP_P3.fcs']);
% make sure it can apply properly, without errors
gated_data = applyFilter(gate,fcshdr,data);

assertElementsAlmostEqual(size(data,1),52093);
assertElementsAlmostEqual(size(gated_data,1),24161,'relative',0.01);

gate = struct(gate);

assertEqual(gate.selected_components, 1);
assertEqual(gate.channel_names,{'FSC-A'  'SSC-A'  'FSC-H'  'SSC-H'  'FSC-W'  'SSC-W'});

expected_mu = [...
    5.1022    3.3902    5.1184    3.2685    4.8002    4.9382;
    4.7724    3.4450    4.7893    3.2905    4.7996    4.9710;
    ];
GDS = struct(gate.distribution);
assertElementsAlmostEqual(GDS.mu,expected_mu,'absolute',0.01);

expected_sigma(:,:,1) = [...
    0.0148    0.0161    0.0101    0.0126    0.0047    0.0034
    0.0161    0.0319    0.0099    0.0271    0.0062    0.0047
    0.0101    0.0099    0.0076    0.0080    0.0025    0.0019
    0.0126    0.0271    0.0080    0.0240    0.0046    0.0031
    0.0047    0.0062    0.0025    0.0046    0.0021    0.0016
    0.0034    0.0047    0.0019    0.0031    0.0016    0.0016
    ];
expected_sigma(:,:,2) = [...
    0.2258    0.0714    0.1925    0.0457    0.0333    0.0257
    0.0714    0.0488    0.0595    0.0402    0.0119    0.0086
    0.1925    0.0595    0.1673    0.0398    0.0252    0.0197
    0.0457    0.0402    0.0398    0.0380    0.0059    0.0022
    0.0333    0.0119    0.0252    0.0059    0.0081    0.0060
    0.0257    0.0086    0.0197    0.0022    0.0060    0.0064
    ];
assertElementsAlmostEqual(GDS.Sigma,expected_sigma,'absolute',0.01);


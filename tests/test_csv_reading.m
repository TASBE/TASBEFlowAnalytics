function test_suite = test_csv_reading
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_fca_readcsv
f1 = '../TASBEFlowAnalytics-Tutorial/template_analysis/csv/LacI-CAGop_Dox01_PointCloud.csv';
header = 'tests/LacI-CAGop.json';

datafile = DataFile('csv', f1, header);

[data, hdr] = fca_read(datafile);

PACIFIC_BLUE_CHANNEL = 3;
NUM_CHANNELS = 3;
assert(strcmp(hdr.par(PACIFIC_BLUE_CHANNEL).name,'Pacific Blue-A'));
assert(all(size(data) == [161608 NUM_CHANNELS]));
assertElementsAlmostEqual(data(1,:),[4.2471e4 4.0352e4 3.7367e4],'absolute',1e0);

function test_cm_read_csv

CM = load_or_make_testing_colormodel();
f1 = '../TASBEFlowAnalytics-Tutorial/template_analysis/csv/LacI-CAGop_Dox01_PointCloud.csv';
% kludge to read as a.u. rather than MEFL
header = 'tests/LacI-CAGop-au.json';

datafile = DataFile('csv', f1, header);

CM = clear_filters(CM); % got to drop the filters because it's a.u. and there's no FSC and SSC
data = readfcs_compensated_ERF(CM,datafile,false,true);
NUM_CHANNELS = 3;
assert(all(size(data) == [161608 NUM_CHANNELS]));
% numbers should change because all a.u. is treated as uncalibrated
assertElementsAlmostEqual(data(1,:),[9.5741e7 9.0444e7 3.8822e7],'absolute',1e3);


function test_cm_read_mixed_csv

CM = load_or_make_testing_colormodel();
f1 = '../TASBEFlowAnalytics-Tutorial/template_analysis/csv/LacI-CAGop_Dox01_PointCloud.csv';
% kludge to read as a mix of a.u. and MEFL rather than all MEFL
header = 'tests/LacI-CAGop-mixed.json';

datafile = DataFile('csv', f1, header);

data = readfcs_compensated_ERF(CM,datafile,false,true);
NUM_CHANNELS = 3;
assert(all(size(data) == [161608 NUM_CHANNELS]));
% there is a non-a.u. channel, so it should all be treated as calibrated
assertElementsAlmostEqual(data(1,:),[4.2471e4 4.0352e4 3.7367e4],'absolute',1e0);


function test_cm_read_broken_csv

f1 = '../TASBEFlowAnalytics-Tutorial/template_analysis/csv/LacI-CAGop_Dox01_PointCloud.csv';
header = 'tests/LacI-CAGop-broken.json';

datafile = DataFile('csv', f1, header);

assertExceptionThrown(@()fca_read(datafile), 'fca_readcsv:NumParameterMismatch', 'No error was raised.');


function test_cm_read_broken_csv2

f1 = '../TASBEFlowAnalytics-Tutorial/template_analysis/csv/LacI-CAGop_Dox01_PointCloud.csv';
header = 'tests/LacI-CAGop-broken2.json';

datafile = DataFile('csv', f1, header);

assertExceptionThrown(@()fca_read(datafile), 'fca_readcsv:MissingRequiredHeaderField', 'No error was raised.');
log = TASBESession.list();
assertEqual(log{end}.contents{end-1}.name, 'UnknownHeaderField');


function test_cm_read_broken_csv3

f1 = '../TASBEFlowAnalytics-Tutorial/template_analysis/csv/LacI-CAGop_Dox01_PointCloud.csv';
header = 'tests/LacI-CAGop-broken3.json';

datafile = DataFile('csv', f1, header);

fca_read(datafile);
log = TASBESession.list();
assertEqual(log{end}.contents{end-1}.name, 'UnknownHeaderField');
assertEqual(log{end}.contents{end}.name, 'FilenameMismatch');

function test_fca_extracolumn
    f1 = 'tests/LacI-extra-column.csv';
    header = 'tests/LacI-extra.json';
    
    datafile = DataFile('csv', f1, header);

    [data, hdr] = fca_read(datafile);

    DERIVED_CHANNEL = 4;
    NUM_CHANNELS = 4;
    assert(strcmp(hdr.par(DERIVED_CHANNEL).name,'Derived_Value'));
    assert(all(size(data) == [50 NUM_CHANNELS]));
    assert(hdr.non_au);
    assertElementsAlmostEqual(data(1,:),[4.2471e4 4.0352e4 3.7367e4 1],'absolute',1e0);

function test_fca_extracolumn_au
    f1 = 'tests/LacI-extra-column.csv';
    header = 'tests/LacI-extra-au.json';
    
    datafile = DataFile('csv', f1, header);

    [data, hdr] = fca_read(datafile);

    DERIVED_CHANNEL = 4;
    NUM_CHANNELS = 4;
    assert(strcmp(hdr.par(DERIVED_CHANNEL).name,'Derived_Value'));
    assert(all(size(data) == [50 NUM_CHANNELS]));
    assert(~hdr.non_au);
    assertElementsAlmostEqual(data(1,:),[4.2471e4 4.0352e4 3.7367e4 1],'absolute',1e0);


function test_cm_extracolumn_ignored
    CM = load_or_make_testing_colormodel();
    
    f1 = 'tests/LacI-extra-column.csv';
    header = 'tests/LacI-extra.json';
    
    datafile = DataFile('csv', f1, header);

    data = readfcs_compensated_ERF(CM,datafile,false,true);
    NUM_CHANNELS = 3;
    assert(all(size(data) == [50 NUM_CHANNELS]));
    % there is a non-a.u. channel, so it should all be treated as calibrated
    assertElementsAlmostEqual(data(1,:),[4.2471e4 4.0352e4 3.7367e4],'absolute',1e0);

function test_cm_extracolumn_used
    CM = load_or_make_testing_colormodel();
    
    f1 = 'tests/LacI-extra-column.csv';
    header = 'tests/LacI-extra.json';
    
    datafile = DataFile('csv', f1, header);

    CM = add_derived_channel(CM,'Derived_Value','Derived','Boolean');
    data = readfcs_compensated_ERF(CM,datafile,false,true);
    NUM_CHANNELS = 4;
    assert(all(size(data) == [50 NUM_CHANNELS]));
    % there is a non-a.u. channel, so it should all be treated as calibrated
    assertElementsAlmostEqual(data(1,:),[4.2471e4 4.0352e4 3.7367e4 1],'absolute',1e0);

function test_cm_extrafiltered
    CM = load_or_make_testing_colormodel();
    CM = add_derived_channel(CM,'Derived_Value','Derived','Boolean');
    filter = RangeFilter('Derived_Value',[1 1]);
    CM = add_postfilter(CM,filter);

    f1 = 'tests/LacI-extra-column.csv';
    header = 'tests/LacI-extra.json';
    datafile = DataFile('csv', f1, header);
    
    data = readfcs_compensated_ERF(CM,datafile,false,true);
    NUM_CHANNELS = 4;
    assert(all(size(data) == [15 NUM_CHANNELS]));
    % there is a non-a.u. channel, so it should all be treated as calibrated
    assertElementsAlmostEqual(data(1,:),[4.2471e4 4.0352e4 3.7367e4 1],'absolute',1e0);
    assertElementsAlmostEqual(data(4,:),[0.2249e4 7.2319e4 0.1156e4 1],'absolute',1e0);
    

function test_suite = test_csv_writing
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_fca_writecsv

TASBEConfig.set('flow.outputPointCloud','true');
TASBEConfig.set('flow.pointCloudPath','/tmp/CSV/');

CM = load_or_make_testing_colormodel();

bins = BinSequence(4,0.1,10,'log_bins');
AP = AnalysisParameters(bins,{});
stem1011 = '../TASBEFlowAnalytics-Tutorial/example_assay/LacI-CAGop_';
experiment = Experiment('test','', {0,{DataFile('fcs', [stem1011 'B3_P3.fcs'])}});
data = read_data(CM, experiment, AP);
writeFcsPointCloudCSV(CM, {{'filename'}}, data);

% confirm that write is correct
written_file = '/tmp/CSV/filename_PointCloud.csv';
written_lines = {'EYFP_MEFL,mKate_MEFL,EBFP2_MEFL';
    '42808.53,40694.79,33597.78'};

fid = fopen(written_file);
for i=1:numel(written_lines),
    assert(strcmp(fgetl(fid),written_lines{i}));
end
fclose(fid);

function test_suite = test_batch_csv_excel
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_batch_csv_excel_endtoend
    % Create TemplateExtraction object
    [filepath, ~, ~] = fileparts(mfilename('fullpath'));
    extractor = TemplateExtraction('test_templates/test_batch_template_csv.xlsx', [end_with_slash(filepath) '../']);
    CM = load_or_make_testing_colormodel();
    [~, ~, ~] = batch_analysis_excel(extractor, CM);
    
    % Check the first five rows of the first point cloud file
    expected_pointCloud = [...
    42471.29    40352.2    37366.76
    2452.02     42665.24    1156.15
    70346.12    67919.83    22964.29
    2807178.51  17492541.02   1156.15
    8687.51     2229.52     1156.15
    ];

    % The first point cloud file: /tmp/LacI-CAGop_B3_B03_P3_PointCloud.csv
    firstPointCloudFile = '/tmp/CSV/LacI-CAGop_Dox01_PointCloud.csv';

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
    
    % Compare JSON headers (esp channel struct)
    real_header = 'LacI-CAGop.json';
    test_header = '/tmp/CSV/LacI Transfer Curve.json';

    fid1 = fopen(real_header); 
    raw1 = fread(fid1,inf); 
    string1 = char(raw1'); 
    fclose(fid1); 
    real = loadjson(string1);
    real_channels = real.channels;
    real_channel = real_channels{1};

    fid2 = fopen(test_header); 
    raw2 = fread(fid2,inf); 
    string2 = char(raw2'); 
    fclose(fid2); 
    test = loadjson(string2);
    test_channels = test.channels;
    test_channel = test_channels{1};

    assertEqual(numel(real_channels), numel(test_channels));
    assertEqual(real_channel.name, test_channel.name);
    assertEqual(real_channel.print_name, test_channel.print_name);
    assertEqual(real_channel.unit, test_channel.unit);
function test_suite = test_batch_output_excel
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_batch_output_excel_endtoend
    % Create TemplateExtraction object
    [filepath, ~, ~] = fileparts(mfilename('fullpath'));
    extractor = TemplateExtraction('test_templates/test_batch_template4.xlsx', [end_with_slash(filepath) '../']);
    CM = load_or_make_testing_colormodel();
    [results, statisticsFile, histogramFile] = batch_analysis_excel(extractor, CM);
    
    % Read the files into matlab tables
    if (is_octave)
        statsTable = csv2cell(statisticsFile);
        histTable = csv2cell(histogramFile);
        statsCell = statsTable(2:end,:);
        histCell = histTable(2:end,:);
    else
        statsTable = readtable(statisticsFile);
        histTable = readtable(histogramFile);
        statsCell = table2cell(statsTable);
        histCell = table2cell(histTable);
    end

    % Split the stats table
    geoMeans = statsCell(:,5:7);
    geoStdDevs = statsCell(:,8:10);
    gmmMeans = statsCell(:,11:16);
    gmmStds = statsCell(:,17:22);
    gmmWeights = statsCell(:,23:28);
    onFracs = statsCell(:,29);
    offFracs = statsCell(:,30);
    
    % Check on/off frac mean and std values
    assertElementsAlmostEqual(round(results{1}.on_fracMean.*10000)./10000, 0.5791, 'relative', 1e-2);
    assertElementsAlmostEqual(round(results{1}.off_fracMean.*10000)./10000, 0.4209, 'relative', 1e-2);
    assertElementsAlmostEqual(round(results{1}.on_fracStd.*10000)./10000, 0.0038, 'relative', 1e-2);
    assertElementsAlmostEqual(round(results{1}.off_fracStd.*10000)./10000, 0.0038, 'relative', 1e-2);

    % Split the hist table
    binCounts = histCell(:,3:5);

    % Strip out the padding put into the sampleIds, means, and stdDevs
    sampleIDListWithPadding = statsCell(:,1);
    sampleIDs = sampleIDListWithPadding(find(~cellfun(@isempty,sampleIDListWithPadding)));


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Check results in CSV files:

    % The first five rows should be enough to verify writing the histogram file
    % correctly.
    expected_bincounts = [...
             637          36       6799; % clipped by the drop threshold    
            2732        2696       8012;        
            3327        2638       8780;        
            4637        2632       8563;        
            4623        3741       7622;        
            ];

    % Means and stddevs tests writing the statistics file correctly.
    expected_means = [...
        260800	429100	22170
        260700	426500	22200
        269300	444900	22110
        271400	454400	22050
        262600	450600	22160
        260800	448200	22550
        266800	445000	22820
        268800	460800	25400
        253000	482200	37920
        242800	492200	48930
        220000	486300	69260
        182500	592500	109400
        160100	697900	159100
        146800	772000	194800
        ];

    expected_stds = [...
        6.611	7.921    1.601
        6.702	7.942    1.600
        6.710	7.944    1.599
        6.746	8.014    1.603
        6.599	7.922    1.604
        6.641	8.080    1.643
        6.595	7.936    1.703
        6.613	8.106    1.992
        6.308	8.224    3.057
        6.028	8.242    3.688
        5.731	8.056    4.508
        5.163	8.538    5.283
        4.613	8.392    5.603
        4.301	8.263    5.578
        ];

    expected_gmm_means = 10.^[...
        4.6200    5.9921    4.6817    6.0832   4.2543    4.2543;    % row 1
        ];

    expected_gmm_stds = 10.^[...
        0.0783    0.3112    0.0900    0.5154    0.0692    0.0692;     % row 1
        ];

    expected_gmm_weights = [...
        0.4195    0.5805    0.3215    0.6785    0.5000    0.5000;     % row 1
        ];
    
    expected_on_fracs = [...
        0.5818
        0.5764
        0.5793
        0.5770
        0.5815
        0.5838
        ];

    expected_off_fracs = [...
        0.4182
        0.4236
        0.4207
        0.4230
        0.4185
        0.4162
        ];

    assertEqual(numel(sampleIDs), 3);
    %assertEqual(numel(sampleIDs), 7);

    % spot-check names
    assertEqual(sampleIDs{1}, 'Dox 0.1/0.2');
    assertEqual(sampleIDs{3}, 'Dox 2.0/5.0');
    %assertEqual(sampleIDs{7}, 'Dox 1000.0/2000.0');

    % spot-check first five rows of binCounts
    assertElementsAlmostEqual(cell2mat(binCounts(1:5,:)), expected_bincounts, 'relative', 2e-2);

    % spot-check geo means and geo std devs.
    for i=1:6, % was 7
        assertElementsAlmostEqual(cell2mat(geoMeans(i,:)), expected_means(i,:), 'relative', 1e-2);
        assertElementsAlmostEqual(cell2mat(geoStdDevs(i,:)), expected_stds(i,:),  'relative', 1e-2);
    end

    % spot check Gaussian Mixture Model materials:
    assertElementsAlmostEqual(cell2mat(gmmMeans(1,:)), expected_gmm_means, 'relative', 1e-2);
    assertElementsAlmostEqual(cell2mat(gmmStds(1,:)), expected_gmm_stds, 'relative', 1e-2);
    assertElementsAlmostEqual(cell2mat(gmmWeights(1,:)), expected_gmm_weights, 'relative', 1e-2);
    assertElementsAlmostEqual(cell2mat(onFracs), expected_on_fracs, 'relative', 1e-2);
    assertElementsAlmostEqual(cell2mat(offFracs), expected_off_fracs, 'relative', 1e-2);
    
    % Check the first five rows of the first point cloud file
    expected_pointCloud = [...
        42801.34    40500.46    33567.67
        2456.10     42822.39    1039.11
        70903.34    68176.25    20623.25
        2830130.69  17561178.05   1039.11
        8742.07     2238.27     1039.11
        ];

    % The first point cloud file: /tmp/LacI-CAGop_B3_B03_P3_PointCloud.csv
    firstPointCloudFile = '/tmp/CSV/LacI-CAGop_B3_P3_PointCloud.csv';

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
    
    
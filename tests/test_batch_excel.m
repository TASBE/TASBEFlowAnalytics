function test_suite = test_batch_excel
    TASBEConfig.checkpoint('test');
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_batch_excel_endtoend
    % Create TemplateExtraction object
    [filepath, ~, ~] = fileparts(mfilename('fullpath'));
    extractor = TemplateExtraction('test_templates/test_batch_template1.xlsx', [end_with_slash(filepath) '../']);
    CM = load_or_make_testing_colormodel();
    [results, statisticsFile, histogramFile] = batch_analysis_excel(extractor, CM);
    
    % Make sure statistics and histogram files are 'none'
    assertEqual(statisticsFile, 'none');
    assertEqual(histogramFile, 'none');
    
    % Yellow, Red, Blue
    result1_expected_bincounts = [...
             637          36       6806; % clipped by the drop threshold        
            2753        2706       8017;        
            3323        2637       8782;        
            4640        2623       8558;        
            4624        3739       7617;        
            5595        4714       6343;        
            5937        5304       3931;        
            6282        5434       1817;        
            5747        5801       511;        
            4272        4683       124;        
            3097        4012       0;
            2284        3469       0;
            2340        2917       0;
            2545        3200       0;
            2845        3612       0;
            3390        3985       0;
            3755        4034       0;
            4031        3985       0;
            4246        4135       0;
            4436        4179       0;
            4502        4199       0;
            4289        4095       0;
            4007        3890       0;
            3630        3817       0;
            3244        3685       0;
            2738        3509       0;
            2203        3248       0;
            1731        3032       0;
            1406        2598       0;
            989        2401        0;
            769        1920        0;
            493        1626        0;
            391        1353        0;
            214         995        0;
            150         808        0;
            101         634        0;
               0           428        0;
               0           272        0;
               0           176        0; 
               0           122        0; 
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

    result_expected_means = 1e5 * [...
        2.5920    4.2576    0.2214
        2.5773    4.2344    0.2217
        2.6618    4.4116    0.2201
%         2.6768    4.4791    0.2193
%         2.5979    4.4447    0.2202
%         2.5794    4.4459    0.2230
%         2.6410    4.3857    0.2246
%         2.6600    4.5635    0.2489
%         2.5019    4.7657    0.3723
%         2.3962    4.8756    0.4765
%         2.1736    4.8155    0.6809
%         1.8055    5.8724    1.0768
        1.5830    6.9178    1.5798
        1.4546    7.6877    1.9350
        ];

    result_expected_stds = [...
        6.5596    7.8289    1.5964
        6.6127    7.8521    1.5955
        6.6173    7.8504    1.5856
%         6.6379    7.8693    1.5881
%         6.5152    7.7858    1.5876
%         6.5557    7.9894    1.6027
%         6.5152    7.7879    1.6348
%         6.5316    7.9963    1.8831
%         6.2255    8.0919    2.9448
%         5.9371    8.1359    3.5238
%         5.6446    7.9458    4.3923
%         5.0904    8.4391    5.1748
        4.5412    8.2867    5.5515
        4.2440    8.2005    5.5330
        ];

    result_expected1_gmm_means = [...
        4.6200    4.6817    4.2543    
        5.9921    6.0832    4.2543    
        ];
    result_expected1_gmm_stds = [...
        0.0783    0.0900    0.0692    
        0.3112    0.5154    0.0692    
        ];
    result_expected1_gmm_weights = [...
        0.4195    0.3215    0.5000    
        0.5805    0.6785    0.5000    
        ];

    result_expected14_gmm_means = [...
        4.6252    4.7382    4.3427    
        5.5610    6.2619    5.6856    
        ];
    result_expected14_gmm_stds = [...
        0.0804    0.1034    0.0904    
        0.2650    0.5099    0.2825    
        ];
    result_expected14_gmm_weights = [...
        0.4211    0.2456    0.3486    
        0.5789    0.7544    0.6514    
        ];

    assertEqual(numel(results), 5);

    % spot-check name, bincenter, bin-count
    assertEqual(results{1}.condition, 'Dox 0.1');
    assertElementsAlmostEqual(log10(results{1}.bincenters([1 10 40 end])), [4.0500    4.9500    7.9500    9.9500], 'relative', 1e-2);
    assertElementsAlmostEqual(results{1}.bincounts, result1_expected_bincounts,     'relative', 1e-2, 50);

    assertEqual(results{5}.condition, 'Dox 2000');
    assertElementsAlmostEqual(log10(results{5}.bincenters([1 10 40 end])), [4.0500    4.9500    7.9500    9.9500], 'relative', 1e-2);

    for i=1:5
        assertElementsAlmostEqual(results{i}.means, result_expected_means(i,:), 'relative', 1e-2);
        assertElementsAlmostEqual(results{i}.stds,  result_expected_stds(i,:),  'relative', 1e-2);
    end

    assertElementsAlmostEqual(results{1}.gmm_means,  result_expected1_gmm_means,  'relative', 1e-2);
    assertElementsAlmostEqual(results{1}.gmm_stds,  result_expected1_gmm_stds,  'relative', 1e-2);
    assertElementsAlmostEqual(results{1}.gmm_weights,  result_expected1_gmm_weights,  'relative', 1e-2);
    assertElementsAlmostEqual(results{5}.gmm_means,  result_expected14_gmm_means,  'relative', 1e-2);
    assertElementsAlmostEqual(results{5}.gmm_stds,  result_expected14_gmm_stds,  'relative', 1e-2);
    assertElementsAlmostEqual(results{5}.gmm_weights,  result_expected14_gmm_weights,  'relative', 1e-2);
    
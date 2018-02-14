function test_suite = test_plusminus_analysis
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;

function test_plusminus_analysis_endtoend

load('../TASBEFlowAnalytics-Tutorial/template_colormodel/CM120312.mat');


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

% Make a map of the batches of plus/minus comparisons to test
% This analysis supports two variables: a +/- variable and a "tuning" variable
stem1011 = '../TASBEFlowAnalytics-Tutorial/example_assay/LacI-CAGop_';
batch_description = {...
 {'Lows';'BaseDox';
  % First set is the matching "plus" conditions
  {0.1,  {[stem1011 'B9_B09_P3.fcs']}; % Replicates go here, e.g., {[rep1], [rep2], [rep3]}
   0.2,  {[stem1011 'B10_B10_P3.fcs']}};
  % Second set is the matching "minus" conditions 
  {0.1,  {[stem1011 'B3_B03_P3.fcs']};
   0.2,  {[stem1011 'B4_B04_P3.fcs']}}};
 {'Highs';'BaseDox';
  {10,   {[stem1011 'C3_C03_P3.fcs']};
   20,   {[stem1011 'C4_C04_P3.fcs']}};
  {10,   {[stem1011 'B9_B09_P3.fcs']};
   20,   {[stem1011 'B10_B10_P3.fcs']}}};
 };

% Execute the actual analysis
TASBEConfig.set('OutputSettings.DeviceName',device_name);
TASBEConfig.set('plots.plotPath','/tmp/plots');
results = process_plusminus_batch( CM, batch_description, AP);

% Make additional output plots
for i=1:numel(results)
    TASBEConfig.set('OutputSettings.StemName',batch_description{i}{1});
    TASBEConfig.set('OutputSettings.DeviceName',device_name);
    TASBEConfig.set('OutputSettings.PlotTickMarks',1);
    plot_plusminus_comparison(results{i})
end

save('-V7','/tmp/LacI-CAGop-plus-minus.mat','batch_description','AP','results');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check results in results:

expected_ratios1 = [...
    0.9703    1.0099    0.9604    0.9722    0.9359    0.9679    1.0177    ...
    1.0083    0.9675    1.0050    0.9423    1.0312    0.9417    0.9768    ...
    0.9417    0.9044    1.0601    0.9296    0.9193    0.9737    1.0260    ...
    0.9831    0.9670    1.0218    1.0087    1.0295    1.0014    0.9998    ...
    1.0122    0.9909    0.9744    1.0635    0.9223    0.9971    1.0246    ...
    0.9027    1.0512    0.8515    0.8538       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN;
    0.8953    0.9746    0.8948    0.9471    1.0188    0.9822    1.0057    ...
    0.9749    0.9391    0.9675    0.9057    0.9005    0.9873    1.0269    ...
    0.9367    1.0088    0.9284    0.8360    0.9268    0.9916    0.9157    ...
    0.9587    0.9869    0.9430    0.9397    0.9336    0.9287    0.9507    ...
    0.9031    0.9122    0.8376    0.9243    0.8898    0.7727    0.8641    ...
    0.8246    0.8080    0.8401    0.6543    0.6727       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN]';

expected_InSNR1 = [...
  -38.4867  -46.7731  -35.5302  -39.1779  -55.0857  -40.6221  -48.9905    ...
  -34.9919  -39.5790  -35.4469  -61.7799  -44.4814  -43.4279  -37.7208    ...
  -35.5624  -39.8118  -42.7565  -56.6078  -39.7806  -34.6803  -38.0924    ...
  -35.3169  -38.0582  -60.3043  -47.4826  -34.6696  -31.7850  -28.7535    ...
  -23.7945  -21.7045  -21.3608  -18.1324  -17.6059  -16.3149  -15.7278    ...
  -13.5991  -13.4115  -10.8743   -9.6963       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN;
  -32.0887  -33.3616  -49.8738  -45.4783  -45.4899  -42.9094  -34.9467    ...
  -41.8559  -51.3976  -41.1789  -38.7564  -29.0290  -44.9261  -35.9305    ...
  -39.7765  -54.8425  -35.1435  -40.6337  -61.2315  -32.1843  -28.9194    ...
  -40.3561  -27.9151  -23.7359  -24.7296  -19.7681  -18.9716  -16.3981    ...
  -13.8112  -13.3603  -12.6749  -11.4734   -9.4093   -9.0514   -9.0620    ...
   -8.4075   -8.1013   -8.8451   -4.7220   -5.3573       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN]';

expected_OutSNR1 = [...
  -41.3183  -50.2429  -38.2845  -41.7744  -34.1039  -40.3652  -45.7040    ...
  -52.2658  -40.3087  -57.0695  -35.6331  -41.9971  -36.8299  -45.2937    ...
  -37.2859  -32.8181  -37.6122  -35.5721  -34.2859  -43.8896  -44.0104    ...
  -47.1029  -40.7640  -43.6501  -51.1173  -39.6137  -65.0618  -67.2762    ...
  -44.7766  -47.2209  -37.4857  -29.5938  -27.4533  -56.8501  -37.8766    ...
  -25.6679  -32.1429  -22.4133  -23.3224       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN;
  -30.0115  -41.7135  -29.6292  -35.8754  -45.2217  -45.4710  -55.3789    ...
  -42.5061  -34.7584  -40.5219  -31.2667  -31.4094  -50.3513  -44.2578    ...
  -36.5406  -54.0465  -35.5705  -27.8890  -35.1778  -54.0854  -33.5441    ...
  -39.2779  -48.9068  -35.0830  -33.6274  -32.2079  -31.0901  -33.3312    ...
  -26.8562  -27.1486  -20.5093  -27.7771  -24.3766  -17.7613  -22.5291    ...
  -20.8888  -19.8080  -21.6260  -14.7413  -15.6473       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN]';

expected_ratios2 = [...
    0.9799    0.8603    0.9429    0.9002    0.8812    0.8419    0.8800    ...
    0.8202    0.7711    0.7017    0.7147    0.5842    0.6314    0.5934    ...
    0.6218    0.6014    0.5043    0.4882    0.4597    0.3598    0.3504    ...
    0.3199    0.2944    0.2396    0.2121    0.1934    0.1670    0.1706    ...
    0.1462    0.1249    0.1059    0.0964    0.0962    0.0751    0.0747    ...
    0.0835    0.0619    0.0547    0.0467       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN;
    0.9233    0.9879    0.8883    0.9438    0.8910    0.8466    0.8693    ...
    0.7810    0.7649    0.7215    0.6672    0.6337    0.5535    0.5742    ...
    0.6300    0.4653    0.4909    0.4391    0.3832    0.3361    0.3213    ...
    0.2669    0.2164    0.2147    0.1594    0.1556    0.1482    0.1182    ...
    0.1075    0.0924    0.0854    0.0756    0.0669    0.0627    0.0559    ...
    0.0515    0.0499    0.0385    0.0344    0.0295       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN]';

expected_InSNR2 = [...
  -21.0814  -26.6197  -21.9729  -17.9886  -20.6502  -19.1464  -18.4565    ...
  -17.5775  -17.2778  -15.1048  -12.8709  -10.2661   -7.6997   -8.1373    ...
   -7.0984   -6.2757   -5.2072   -4.4266   -3.4554   -2.2106   -1.0780    ...
    0.1714    0.6332    1.7853    2.5452    2.9177    3.4539    4.0354    ...
    4.0777    4.1819    4.5797    4.6572    4.1730    4.2073    3.9740    ...
    3.4139    3.4947    3.3711    2.7504       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN;
  -15.8142  -17.0562  -17.9996  -15.6329  -16.7130  -16.9782  -16.5726    ...
  -16.0637  -14.9390  -12.8546  -11.9615   -7.6430   -7.0551   -6.1516    ...
   -5.3003   -5.1714   -3.8763   -2.6068   -1.6076   -0.8768   -0.0232    ...
    1.2525    1.8627    2.2800    2.9933    3.3034    3.5991    3.5277    ...
    3.5822    3.5864    3.5078    3.3310    3.0951    2.9206    2.4236    ...
    2.2913    1.8108    2.4632    0.7311    0.5650       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN]';

expected_OutSNR2 = [...
  -43.9820  -26.3210  -34.9311  -30.3720  -28.2903  -25.6038  -28.3594    ...
  -24.6769  -22.3283  -20.0162  -20.4524  -16.9533  -18.9065  -18.1297    ...
  -19.0362  -18.4659  -15.9746  -15.5365  -14.8917  -12.3331  -11.8811    ...
  -10.8467   -9.8820   -7.9581   -6.9954   -6.0076   -4.9042   -4.5633    ...
   -3.6290   -2.9281   -2.0903   -1.5454   -1.7276   -1.3129   -1.0680    ...
   -2.2127   -1.3963   -0.8617   -1.0800       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN;
  -32.6122  -51.7024  -28.6628  -35.2323  -29.5257  -26.0661  -27.6671    ...
  -22.6600  -22.0884  -20.5135  -18.8750  -18.3714  -16.6315  -17.4932    ...
  -19.2210  -14.8864  -15.6098  -14.4152  -13.0392  -11.7019  -11.2292    ...
   -9.5699   -7.9860   -7.5146   -5.5817   -5.0231   -4.5080   -3.0905    ...
   -2.6013   -1.9559   -1.1906   -1.0204   -0.6503   -0.8013   -0.5238    ...
   -0.6172   -0.5548   -0.4493   -0.8065   -0.9529       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN       NaN       NaN       NaN    ...
       NaN       NaN       NaN       NaN]';


assertEqual(numel(results),2);
assertElementsAlmostEqual(results{1}.MeanRatio, [0.9764; 0.8880],   'relative', 0.01);
assertElementsAlmostEqual(results{1}.Ratios, expected_ratios1,      'relative', 0.01);
assertElementsAlmostEqual(results{1}.InputSNR, expected_InSNR1,     'relative', 0.1);
assertElementsAlmostEqual(results{1}.OutputSNR, expected_OutSNR1,   'relative', 0.1);

assertElementsAlmostEqual(results{2}.MeanRatio, [0.2379; 0.1786],   'relative', 0.01);
assertElementsAlmostEqual(results{2}.Ratios, expected_ratios2,      'relative', 0.01);
assertElementsAlmostEqual(results{2}.InputSNR, expected_InSNR2,     'relative', 0.1);
assertElementsAlmostEqual(results{2}.OutputSNR, expected_OutSNR2,   'relative', 0.1);


% <CopyRight 2019> Haoyin Zhou, Jagadeesan Jayender
% Surgical Planning Laboratory
% Brigham and Women's Hospital
% Harvard Medical School
% Zhou, Haoyin, and Jagadeesan Jayender. "Smooth Deformation Field-based Mismatch Removal in Real-time."
% arXiv preprint arXiv:2007.08553 (2020).

clear 
clc
close all
IniToolbox;

%% the input pair of images, we prepared some pairs of images for this demo
% ImgPath1 = [ '../data/DUT_Set001_Img025_05.bmp'];
% ImgPath2 = [ '../data/DUT_Set001_Img065_05.bmp'];
% ImgPath1 = [ '../data/liver1.jpg'];
% ImgPath2 = [ '../data/liver2.jpg'];
% ImgPath1 = [ '../data/HamlynData1.jpg'];
% ImgPath2 = [ '../data/HamlynData2.jpg'];
ImgPath1 = [ '../data/church1.jpg'];
ImgPath2 = [ '../data/church2.jpg'];

%% standard matlab way to obtain SURF matches, it could be replaced with other ways
I1 = imread(ImgPath1);
I2 = imread(ImgPath2);
I1_gray = rgb2gray(I1);
I2_gray = rgb2gray(I2);

points1 = detectSURFFeatures(I1_gray, 'MetricThreshold', 0.01, 'NumScaleLevels', 4);
[f1,vpts1] = extractFeatures(I1_gray,points1);
points2 = detectSURFFeatures(I2_gray, 'MetricThreshold', 0.01, 'NumScaleLevels', 4);
[f2,vpts2] = extractFeatures(I2_gray,points2);  

indexPairs = matchFeatures(f1,f2);
matchedPoints1 = vpts1(indexPairs(:,1));
matchedPoints2 = vpts2(indexPairs(:,2));            

X1 = matchedPoints1.Location';
X2 = matchedPoints2.Location';

%% EMDQ
% the initialization of EMDQ
[EMDQscale] = GetDataScaleForEMDQ_2D( size(I2,2),size(I2,1));
EMDQparams = EMDQ_Initialization_2D(EMDQscale);

% the main algorithm
tic
[inliersMask_RANSAC, inliersMask_final, dq_points, mu_points] = EMDQ(X1, X2, EMDQparams);
toc

% visulization 1
ShowMatchOutliersRemovalResults( I1, X1, X2, inliersMask_final, 'EMDQ' );
plot(X1(1,inliersMask_final),X1(2,inliersMask_final),'bo', 'MarkerSize',3,'LineWidth', 1.5);

% visulization 2
ShowMatchOutliersRemovalResults( I2, X1, X2, inliersMask_final, 'EMDQ' );
plot(X2(1,inliersMask_final),X2(2,inliersMask_final),'bo', 'MarkerSize',3,'LineWidth', 1.5);

% generation and visulization of the deformation field 
% note that the slow computaiton here is caused by Matlab drawing arrows.
[gridcoord, gridcoord_dq, dq_grid, mu_grid, gridmask] = GenerateGridDeformationField( X2(:,inliersMask_final), dq_points(inliersMask_final,:), mu_points(inliersMask_final,:), 50, EMDQparams, size(I1,2), size(I1,1), 0.0);
ShowGridDeformationField(I2, gridcoord, gridcoord_dq, gridmask);

%% VFC
if ~exist('conf', 'var'), conf = []; end
conf = VFC_init(conf);

[nX, nY, normal]=norm2(X1',X2');
VecFld=FastVFC(nX, nY-nX, conf); % could also try FastVFC or SparseVFC
VFCmask = zeros(1,size(X1,2));
VFCmask(VecFld.VFCIndex) = 1;
VFCmask = logical(VFCmask);
ShowMatchOutliersRemovalResults( I1, X1, X2, VFCmask, 'VFC' );
  
%% AHC
AHCidx = AHC([X1;ones(1,size(X1,2))],[X2;ones(1,size(X1,2))],45);
AHCmask = zeros(1,size(X1,2));
AHCmask(AHCidx) = 1;
AHCmask = logical(AHCmask);
ShowMatchOutliersRemovalResults( I1, X1, X2, AHCmask, 'AHC' );
 
%% LPM
LPMidx = LPM_func(X1',X2');
LPMmask = zeros(1,size(X1,2));
LPMmask(LPMidx) = 1;
LPMmask = logical(LPMmask);
ShowMatchOutliersRemovalResults( I1, X1, X2, LPMmask, 'LPM' );

%% LMR
% tested on Matlab 2014a, may have some problems on other versions, such as 2017b
LMRnet = load('./LMR/Trained_Model/Net.mat'); % used by the LMR algorithm

LMRidx = LMR_func(X1',X2',LMRnet.net);
LMRmask = zeros(1,size(X1,2));
LMRmask(LMRidx) = 1;
LMRmask = logical(LMRmask);
ShowMatchOutliersRemovalResults( I1, X1, X2, LMRmask, 'LMR' );
 
%% SIM
nX=[zscore(X1); ones(1,size(X1,2))];
nY=[zscore(X2); ones(1,size(X1,2))];
SIMidx=SIM(nX,nY);
SIMmask = zeros(1,size(X1,2));
SIMmask(SIMidx) = 1;
SIMmask = logical(SIMmask);
ShowMatchOutliersRemovalResults( I1, X1, X2, SIMmask, 'SIM');

%% save the results to files
% saveas(1, '../results/EMDQ_1.png');
% saveas(2, '../results/EMDQ_2.png');
% saveas(3, '../results/EMDQ_deformaiton_field.png');
% saveas(4, '../results/VFC.png');
% saveas(5, '../results/AHC.png');
% saveas(6, '../results/LPM.png');
% saveas(7, '../results/LMR.png');
% saveas(8, '../results/SIM.png');


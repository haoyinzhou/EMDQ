%   This is a demo for removing outliers.

% Authors: Jiayi Ma (jyma2010@gmail.com)
% Date:    04/17/2012

clear; 
% close all; 
% clc

% Read images
% ImgName1 = 'D:\DTUDataset\SET002\Img005_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET002\Img025_05.bmp';
% ImgName1 = 'D:\DTUDataset\SET012\Img049_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET012\Img025_05.bmp';
% ImgName1 = 'D:\DTUDataset\SET015\Img049_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET015\Img025_05.bmp';
% ImgName1 = 'D:\DTUDataset\SET018\Img049_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET018\Img025_05.bmp';
% ImgName1 = 'D:\DTUDataset\SET003\Img011_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET003\Img025_05.bmp';
% ImgName1 = 'D:\DTUDataset\SET010\Img065_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET010\Img025_05.bmp';
% ImgName1 = 'D:\DTUDataset\SET012\Img045_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET012\Img025_05.bmp';
% ImgName1 = 'D:\DTUDataset\SET008\Img071_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET008\Img025_05.bmp';
% ImgName1 = 'D:\DTUDataset\SET009\Img071_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET009\Img025_05.bmp';
% ImgName1 = 'D:\DTUDataset\SET025\Img071_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET025\Img025_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET008\Img059_05.bmp';
% ImgName1 = 'D:\DTUDataset\SET008\Img005_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET008\Img055_05.bmp';
% ImgName1 = 'D:\DTUDataset\SET008\Img005_05.bmp';
% ImgName1 = 'D:\DTUDataset\SET007\Img045_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET007\Img025_05.bmp';
% ImgName1 = 'C:\matlabworkplace\FeatureDQField\sampledata\Set002_Img085_05.bmp';
% ImgName2 = 'C:\matlabworkplace\FeatureDQField\sampledata\Set002_Img025_05.bmp';
% ImgName1 = 'C:\matlabworkplace\FeatureDQField\sampledata\Set010_Img075_05.bmp';
% ImgName2 = 'C:\matlabworkplace\FeatureDQField\sampledata\Set010_Img025_05.bmp';
% ImgName1 = 'D:\DTUDataset\SET016\Img001_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET016\Img025_05.bmp';
% ImgName1 = 'D:\DTUDataset\SET002\Img001_05.bmp';
% ImgName2 = 'D:\DTUDataset\SET002\Img025_05.bmp';
ImgName1 = 'D:\DTUDataset\SET017\Img105_05.bmp';
ImgName2 = 'D:\DTUDataset\SET017\Img025_05.bmp';

I1 = imread(ImgName1);
I2 = imread(ImgName2);

TH = 0.1;
I1_gray = rgb2gray(I1);
I2_gray = rgb2gray(I2);
points1 = detectSURFFeatures(I1_gray, 'MetricThreshold', TH);
[f1,vpts1] = extractFeatures(I1_gray,points1);
points2 = detectSURFFeatures(I2_gray, 'MetricThreshold', TH);
[f2,vpts2] = extractFeatures(I2_gray,points2);  
indexPairs = matchFeatures(f1,f2, 'MatchThreshold', 100.0);

% matchedPoints1 = vpts1(indexPairs(:, 1));
% matchedPoints2 = vpts2(indexPairs(:, 2));
% figure; 
% showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2,'montage');

vpts1_matched.Scale  = vpts1.Scale(indexPairs(:,1));
vpts1_matched.SignOfLaplacian  = vpts1.SignOfLaplacian(indexPairs(:,1));
vpts1_matched.Orientation  = vpts1.Orientation(indexPairs(:,1));
vpts1_matched.Location  = vpts1.Location(indexPairs(:,1),:);
vpts1_matched.Metric  = vpts1.Metric(indexPairs(:,1));
vpts1_matched.Count = size(indexPairs(:,1));
vpts2_matched.Scale  = vpts2.Scale(indexPairs(:,2));
vpts2_matched.SignOfLaplacian  = vpts2.SignOfLaplacian(indexPairs(:,2));
vpts2_matched.Orientation  = vpts2.Orientation(indexPairs(:,2));
vpts2_matched.Location  = vpts2.Location(indexPairs(:,2),:);
vpts2_matched.Metric  = vpts2.Metric(indexPairs(:,2));
vpts2_matched.Count = size(indexPairs(:,2));

X = points1.Location(indexPairs(:,1),:);
Y = points2.Location(indexPairs(:,2),:);

% tempoffset = 200.0 * randn(size(X,1),2);
% X = X + tempoffset;
% Y = Y + tempoffset;


% Load data: initial correspondences and ground truth
%load('church.mat');

% Data normalization
[nX, nY, normal]=norm2(X,Y);

tic;
% Initialization
conf.method = 'VFC';
if ~exist('conf', 'var'), conf = []; end
conf = VFC_init(conf);

% Ourlier removal
switch conf.method
    case 'VFC'
        VecFld=VFC(nX, nY-nX, conf);
    case 'FastVFC'
        VecFld=FastVFC(nX, nY-nX, conf);
    case 'SparseVFC'
        VecFld=SparseVFC(nX, nY-nX, conf);
end
toc;

% Denormalization
VecFld.V=(VecFld.V+nX)*normal.yscale+repmat(normal.ym,size(Y,1),1)-X;

% Evaluation
if ~exist('CorrectIndex', 'var'), CorrectIndex = VecFld.VFCIndex; end
[precise, recall, corrRate] = evaluate(CorrectIndex, VecFld.VFCIndex, size(X,1));

% Plot results
% plot_matches(I1, I2, X, Y, VecFld.VFCIndex, CorrectIndex);

A = X(VecFld.VFCIndex,:);
B = Y(VecFld.VFCIndex,:);


figure
imshow(I2);
hold on
for i = 1:size(VecFld.VFCIndex)
    temp = [A(i,:);B(i,:)];
    plot(temp(:,1),temp(:,2),'y-','LineWidth', 2); 
    hold on
end
%axis([0 800 0 600]);

% vpts1_matched.Scale  = vpts1.Scale(indexPairs(:,1));
% vpts1_matched.SignOfLaplacian  = vpts1.SignOfLaplacian(indexPairs(:,1));
% vpts1_matched.Orientation  = vpts1.Orientation(indexPairs(:,1));
% vpts1_matched.Location  = vpts1.Location(indexPairs(:,1),:);
% vpts1_matched.Metric  = vpts1.Metric(indexPairs(:,1));
% vpts1_matched.Count = size(indexPairs(:,1));
% figure; 
% showMatchedFeatures(I1,I2,matchedPoints1,matchedPoints2,'montage');


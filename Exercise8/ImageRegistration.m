function [Tx_RGB Ty_RGB]= ImageRegistration
% *************************************************************************
% Wavelets and Applications Course - Dr. P.L. Dragotti
% MATLAB mini-project 'Sampling Signals with Finite Rate of Innovation'
% Exercice 6
% *************************************************************************
% 
% FOR STUDENTS
%
% This function registers the set of 40 low-resolution images
% 'LR_Tiger_xx.tif' and returns the shifts for each image and each layer
% Red, Green and Blue. The shifts are calculated relatively to the first
% image 'LR_Tiger_01.tif'. Each low-resolution image is 64 x64 pixels.
%
%
% OUTPUT:   Tx_RGB: horizontal shifts, a 40x3 matrix
%           Ty_RGB: vertical shifts, a 40x3 matrix
%
% NOTE: _Tx_RGB(1,:) = Ty_RGB(1,:) = (0 0 0) by definition.
%       _Tx_RGB(20,2) is the horizontal shift of the Green layer of the
%       20th image relatively to the Green layer of the firs image.
%
%
% OUTLINE OF THE ALGORITHM:
%
% 1.The first step is to compute the continuous moments m_00, m_01 and m_10
% of each low-resolution image using the .mat file called:
% PolynomialReproduction_coef.mat. This file contains three matrices
% 'Coef_0_0', 'Coef_1_0' and 'Coef_0_1' used to calculate the continuous
% moments.
%
% 2.The second step consists in calculating the barycenters of the Red,
% Green and Blue layers of the low-resolution images.
%
% 3.By computing the difference between the barycenters of corresponding 
% layers between two images, the horizontal and vertical shifts can be 
% retrieved for each layer.
%
%
% Author:   Loic Baboulaz
% Date:     August 2006
%
% Imperial College London
% *************************************************************************

%clear all % delete this one
% Load the coefficients for polynomial reproduction

% -------- include your code here -----------

load('PolynomialReproduction_coef.mat','Coef_0_0','Coef_1_0','Coef_0_1');
%% Calculate xbar ybar for the first image
img1 = imread('LR_Tiger_01.tif');

%set threshold
thresh_denoi = 105; % threshold, the best is 105/255, its PSNR is 24.45dB
img1(img1<thresh_denoi) = 0; % Set to 0 if the value is below to the threshold
img1 = double(img1)/255;
%img1Rm00 = sum(sum(Coef_0_0 .* img1(:,:,1)));
%img1Gm00 = sum(sum(Coef_0_0 .* img1(:,:,1)));
img1RGBm = zeros(3, 3); %Row is RGB, Col is m00 -1 m01 -2 m10 -3
for i = 1: 3 % i refers to RGB respectively
    img1RGBm(i, 1) = sum(sum(Coef_0_0 .* img1(:,:,i)));
    img1RGBm(i, 2) = sum(sum(Coef_0_1 .* img1(:,:,i)));
    img1RGBm(i, 3) = sum(sum(Coef_1_0 .* img1(:,:,i)));
end
[xbar1, ybar1] = getXYbar(img1RGBm); % Calculate the xbar ybar respectively
%dx1 = zeros(1,3);
%dy1 = zeros(1,3);

Tx_RGB = zeros(40, 3);
Ty_RGB = zeros(40, 3);

%Tx_RGB(1, :) = dx1;
%Ty_RGB(1, :) = dy1;
%% FLow the similar routine above to calculate the dx dy of the other images

numFormat = '%02d';
for img_num = 2: 40
    file_name = ['LR_Tiger_'  num2str(img_num, numFormat)  '.tif']; % Set file name
    imgn = imread(file_name); % Load image
    imgn(imgn<thresh_denoi) = 0; % Use the threshold
    imgn = double(imgn)/255; % Normalise
    imgnRGBm = zeros(3, 3); %Row is RGB, Col is m00 -1 m01 -2 m10 -3
    for i = 1: 3 % i refers to RGB respectively
        imgnRGBm(i, 1) = sum(sum(Coef_0_0 .* imgn(:,:,i)));
        imgnRGBm(i, 2) = sum(sum(Coef_0_1 .* imgn(:,:,i)));
        imgnRGBm(i, 3) = sum(sum(Coef_1_0 .* imgn(:,:,i)));
    end
    [xbarn, ybarn] = getXYbar(imgnRGBm); % Calculate xbar ybar for image n
    % Calculate dx dy for image n
    Tx_RGB(img_num, :) = xbarn - xbar1; 
    Ty_RGB(img_num, :) = ybarn - ybar1;
end
end

%% The funciton to get Xbar and Ybar
function [xbar, ybar] = getXYbar(imgRGBm)
xbar = zeros(1, 3);
ybar = zeros(1, 3);
    for i = 1: 3
        xbar(i) = imgRGBm(i, 3)/imgRGBm(i, 1); % Get X bar
    end
    for i = 1: 3
        ybar(i) = imgRGBm(i, 2)/imgRGBm(i, 1); % Get Y bar
    end
        
end

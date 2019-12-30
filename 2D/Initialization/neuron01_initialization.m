function [InitialGuess] = neuron01_initialization(Img)

% parameter for detecting nucleus
kernelSize      = 80;               % size of LoG filter (unit: pixel) (depend on the scale of image)
kernelScale     = 20;               % scale of LoG filter (unit: pixel) (depend on the scale of image)
threRatio       = 0.5;              % threshod for detecting cell (not critical parameter)
levelsetR       = 40;               % parameter for level set (refer to "Distance Regularized Level Set Evolution and Its Application to Image Segmentation", in IEEE TRANSACTIONS ON IMAGE PROCESSING)

minSizeCell     = kernelSize/2;     % minimum size of cell
inputImg1 = double(Img);
inputImg = inputImg1/(max(inputImg1(:)));
[x, y]  = find_nucleus_center(inputImg, inputImg1, kernelSize, kernelScale, threRatio, minSizeCell, levelsetR);

[height, width] = size(inputImg);
% make kernel and convolution
LoG = LoG_kernel(kernelSize,kernelScale);
LoGImg = conv2(inputImg, LoG, 'same');

threshold = min(min(LoGImg));
threRatioImg = LoGImg/threshold;
threImg = (threRatioImg >= threRatio);
BW_soma = threImg;

%% fiber detection
BW = fibermetric(Img);
BW = imfill(BW,'holes');
BW_neurite = zeros(size(Img));
BW_neurite(BW>=0.07) = 1;

M = zeros(size(Img));
M(BW_neurite>0) = 1;
M(BW_soma>0) = 1;

%% clean up image
se = strel('disk',10);
M = imclose(M,se);
M = bwlabel(M,4);
M1 = zeros(size(M));
label =zeros(size(x,1));

for i=1:size(x,1)
    label(i) = M(y(i,1),x(i,1));
    M1(M==label(i)) = 1;
end


M = M1;

figure
imagesc(M)
colormap gray
set(gca,'position',[0 0 1 1],'units','normalized')

InitialGuess = M;
end



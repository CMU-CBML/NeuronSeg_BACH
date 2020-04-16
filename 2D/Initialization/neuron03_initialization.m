function [InitialGuess] = neuron03_initialization(Img)

%% Soma detection
% Source code for soma detection taken from: https://github.com/kilho/NIA
% Paper: Kim, K. M., Son, K., & Palmore, G. T. R. (2015). Neuron image analyzer: 
% Automated and accurate extraction of neuronal data from low quality images. Scientific reports, 5, 17062.
% Paper: Li, C., Xu, C., Gui, C., & Fox, M. D. (2010). Distance regularized level set evolution and its
% application to image segmentation. IEEE transactions on image processing, 19(12), 3243-3254.

%-- Input:

%Img: input image

%-- Output: 

%InitialGuess: binary image describing initial guess for level set function

% parameter for detecting nucleus
kernelSize      = 80;               % size of LoG filter (unit: pixel) (depend on the scale of image)
kernelScale     = 20;               % scale of LoG filter (unit: pixel) (depend on the scale of image)
threRatio       = 0.5;              % threshod for detecting cell (not critical parameter)
sigma = 3; %smooth image before intialization
neurite_threshold = 0.07; %threshold to detect neurite
disk_radius = 6;

minSizeCell     = kernelSize/2;     % minimum size of cell
Img = imgaussfilt(Img,sigma); %Apply Gaussian Filter on the image
inputImg1 = double(Img);
inputImg = inputImg1/(max(inputImg1(:)));

%Find the soma center using non-maximal suppression
[x, y]  = find_nucleus_center(inputImg, kernelSize, kernelScale, threRatio, minSizeCell);

% make kernel and convolution
LoG = LoG_kernel(kernelSize,kernelScale);
LoGImg = conv2(inputImg, LoG, 'same');

% Highlight soma using LoG filter
threshold = min(min(LoGImg));
threRatioImg = LoGImg/threshold;
threImg = (threRatioImg >= threRatio);

%Binary image with soma detected
BW_soma = threImg;

%% Neurite detection
BW = fibermetric(Img);
BW = imfill(BW,'holes');

%Binary image with neurites detected
BW_neurite = zeros(size(Img));
BW_neurite(BW>=neurite_threshold) = 1;

%Combine binary image of soma and neurites
M = zeros(size(Img));
M(BW_neurite>0) = 1;
M(BW_soma>0) = 1;

%% Clean initial guess using connected component analysis
se = strel('disk',disk_radius); %create a morphological disc component
M = imclose(M,se); %close the image before connected component analysis
M = bwlabel(M,4); %create a labelled image
M1 = zeros(size(M));
label =zeros(size(x,1)); %find label associated with soma center

%connect pixels with same label as cell center
for i=1:size(x,1)
    label(i) = M(y(i,1),x(i,1));
    M1(M==label(i)) = 1;
end

M = M1;
InitialGuess = M;
end



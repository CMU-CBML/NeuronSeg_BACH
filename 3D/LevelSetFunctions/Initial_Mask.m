function [ output_args ] = Initial_Mask( binary_image )
%This code is obtained from: https://ch.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/45082/versions/2/previews/Initial_Mask.m/index.html?access_key=
%This function compute the initial Signed Distance Function

% input: binary image
%output: image as signed distance function
N = 2 .* binary_image - 1; 
SDF = bwdist(N < 0) - bwdist(N >= 0);
output_args = SDF;
end
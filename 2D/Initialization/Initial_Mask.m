function [ output_args ] = Initial_Mask( binary_image )
%This function compute the initial Signed Distance Function
N = 2 .* binary_image - 1; 
SDF = bwdist(N < 0) - bwdist(N >= 0);
output_args = SDF;
end
% Main function
clc
clear all
close all

path = pwd;
addpath(strcat(path,'/readimg'));
addpath(strcat(path,'/images2D'));
addpath(strcat(path,'/images2D/neuron_data'));
addpath(strcat(path,'/bsplinelevelset'));
addpath(strcat(path,'/composition'));
addpath(strcat(path,'/setparameters'));
addpath(strcat(path,'/thbspline'));
addpath(strcat(path,'/iterationloop_funcs'));
addpath(strcat(path,'/levelsetseg_funcs'));

parameters = setparameters_neuron01_seg();

%% Read Image, Initialize
Img = imread('01.tif');
Img = double(Img(:,:,2));
Img = imresize(Img,[1040,1040]);
Img = (Img-min(Img(:)))./(max(Img(:))-min(Img(:)));
Img1 = Img;
Img = imgaussfilt(Img,3);
F = Img.*255;

[M] = neuron01_initialization(Img);
M = M.*255;

% %phi = -0.5.*ones(size(M));
% %phi(M==255) = 0.5;
% phi_temp = Initial_Mask(M);
% phi = -double(phi_temp);
% 
% F00 = F;
% M00 = M;
% phi00 = phi;
% 
% nx = size(F,1);
% ny = size(F,2);
% 
% %% Display intial config
% figure
% axis equal
% imagesc(M)
% colormap gray;
% set(gca,'position',[0 0 1 1],'units','normalized')
% 
% figure
% axis equal
% imagesc(Img1)
% colormap gray;
% set(gca,'position',[0 0 1 1],'units','normalized')
% 
% figure
% axis equal
% imagesc(F-M)
% colormap gray;
% set(gca,'position',[0 0 1 1],'units','normalized')
% 
% figure
% axis equal
% imagesc(F)
% colormap gray;
% set(gca,'position',[0 0 1 1],'units','normalized')
% phi_img = 255.*ones(nx,ny);
% phi_img(phi<=0) = -255;
% hold on
% contour(phi,[0,0],'g','Linewidth',3.24);
% hold off
% 
% tic
% [M,phi,VF] = MultipleResolution2D_neuron01_seg(M,F,phi,parameters);
% toc
% Main function for neuron segmentation

clc
clear all
close all

%add paths
path = pwd;
addpath(strcat(path,'/BsplineComposition'));
addpath(strcat(path,'/data'));
addpath(strcat(path,'/Initialization'));
addpath(strcat(path,'/IterationLoopFunctions'));
addpath(strcat(path,'/LevelSetFunctions'));
addpath(strcat(path,'/SetParameters'));
addpath(strcat(path,'/THBSpline'));

% set parameters
parameters = setparameters_neuron03_seg();

%% Read Neuron Cell Image

Img = imread('08.tif');
Img = double(Img(:,:,2));
Img = imresize(Img,[1040,1040]);
Img = (Img-min(Img(:)))./(max(Img(:))-min(Img(:)));
F = Img.*255;

%% Perform automatic initialization
[M] = neuron03_initialization(Img);
M = M.*255; %binary image after initialization

phi_temp = Initial_Mask(M);
phi = -double(phi_temp); %level set function

F00 = F;
M00 = M;
phi00 = phi;

nx = size(F,1);
ny = size(F,2);

%% Display intial image

figure
axis equal
imagesc(F)
colormap gray;
set(gca,'position',[0 0 1 1],'units','normalized')

figure
axis equal
imagesc(F)
colormap gray;
set(gca,'position',[0 0 1 1],'units','normalized')
phi_img = 255.*ones(nx,ny);
phi_img(phi<=0) = -255;
hold on
contour(phi,[0,0],'g','Linewidth',3.24);
hold off

[phi,VF] = MultipleResolution2D(F,phi,parameters);
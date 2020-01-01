% Main function
clc
clear all
close all

%% Add the paths of subfolders of the software
tic
path = pwd;
addpath(strcat(path,'/data'));
addpath(strcat(path,'/BsplineComposition'));
addpath(strcat(path,'/IterationLoopFunctions'));
addpath(strcat(path,'/LevelSetFunctions'));
addpath(strcat(path,'/PostProcessing'));
addpath(strcat(path,'/Initialization'));
addpath(strcat(path,'/SetParameters'));
addpath(strcat(path,'/THBSpline'));
disp('Added paths....');

%% Set the parameters for running registration process
disp('Setting parameters...');
param = setparameters_helix();

%% Read the image
disp('Reading image...');
load 'helix_img.mat';
Img = helix_img;
Img = permute(Img,[3,2,1]);
Img = (Img-min(Img(:)))./(max(Img(:))-min(Img(:)));

Img1 = Img;
Img = imnoise(Img,'salt & pepper',0.25);

Img = Img.*255; %Noise added to helix image
Img1 = Img1.*255; %Helix image without noise

Img = double(Img);
Img1 = double(Img1);

%Gaussain smoothing of image
Img = imgaussfilt(Img,1);

%Image size
[nx,ny,nz] = size(Img);

% disp('Initialization of contour...');
%Sphere center at (25,25,50) with radius 5, converted to signed distance
%function
phi = sdf3circle(size(Img,1), size(Img,2), size(Img,3),25,25,50,5);

%Original image
M = zeros(nx,ny,nz);
M(phi<=0) = 255;
F = Img;

figure
montage(Img,'DisplayRange',[0 64]);
colormap gray

F00 = F;
M00 = M;
phi00 = phi;

disp('Display images...');
figure
axis equal
imagesc(F(:,:,50))
colormap gray;
phi_img = 255.*ones(nx,ny,nz);
phi_img(phi<=0) = -255;
hold on
contour(phi_img(:,:,50),[0,0],'g','Linewidth',1.24);
set(gca,'position',[0 0 1 1],'units','normalized')
hold off

[X,Y,Z] = meshgrid(1:size(F,2),1:size(F,1),1:size(F,3));

figure
p = patch(isosurface(X,Y,Z,F,128));
isonormals(X,Y,Z,F,p);
p.FaceColor = 'red';
p.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud
hold on
p1 = patch(isosurface(X,Y,Z,phi,0));
isonormals(X,Y,Z,phi,p1);
p1.FaceColor = 'green';
p1.EdgeColor = 'none';
daspect([1 1 1])
view(3); 
axis tight
camlight 
lighting gouraud
hold on

[phi,VF] = MultipleResolution3D_heliximg(F,phi,param,Img1);


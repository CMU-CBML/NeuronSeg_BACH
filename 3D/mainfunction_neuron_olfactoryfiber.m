% Main function
clc
clear all
close all

%% Read Image, Initialize
tic
path = pwd;
addpath(strcat(path,'/data/OP_1'));
addpath(strcat(path,'/BsplineComposition'));
addpath(strcat(path,'/IterationLoopFunctions'));
addpath(strcat(path,'/LevelSetFunctions'));
addpath(strcat(path,'/PostProcessing'));
addpath(strcat(path,'/Initialization'));
addpath(strcat(path,'/SetParameters'));
addpath(strcat(path,'/THBSpline'));
addpath(strcat(path,'/initialization'));
p = genpath('OlfactoryFiber');
addpath(p);

disp('Added paths....');
%add the paths of subfolders of the software

%set the parameters for running registration process
disp('Setting parameters...');
param = setparameters_olfactoryfiber();

%read img
disp('Reading image...');
n_slices = 60;
Img = zeros(512,512,n_slices);
for k = 1:n_slices
    str1 = num2str(k);
    str2 = '.tif';
    str = strcat(str1,str2);
    Img1 = imread(str);
    Img(:,:,k) = Img1;
end

[tree, name, path] = load_tree ('OP_1.swc');
Img = double(Img);
Img = (Img-min(Img(:)))./(max(Img(:))-min(Img(:)));
Img = Img.*255;


disp('Initialization of contour...');
[M] = neuron_initialization(Img);
M = M.*255;
F = Img;

figure
montage(Img);
colormap gray

phi_temp = Initial_Mask(M);
phi = -double(phi_temp);

F00 = F;
phi00 = phi;

nx = size(F,1);
ny = size(F,2);
nz = size(F,3);

disp('Display images...');

figure

[X,Y,Z] = meshgrid(1:size(F,2),1:size(F,1),1:size(F,3));

figure
p = patch(isosurface(X,Y,Z,F,64));
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

[phi,V] = MultipleResolution3D_olfactoryfiber(F,phi,param,tree);


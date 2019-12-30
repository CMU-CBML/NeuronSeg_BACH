function boundaryNeclues = levelset_nucleus(centerPoint, r , img, height, width)

%% parameter setting
timestep=5;  % time step
mu=0.2/timestep;  % coefficient of the distance regularization term R(phi)
iter_inner=5;
iter_outer=50;
lambda=5; % coefficient of the weighted length term L(phi)
alfa=1.5;  % coefficient of the weighted area term A(phi)
epsilon=1.5; % papramater that specifies the width of the DiracDelta function

% initialize
maxX = centerPoint(1)+r;
minX = centerPoint(1)-r;
maxY = centerPoint(2)+r;
minY = centerPoint(2)-r;

if maxX > width
    maxX = width;
end
if minX < 1
    minX = 1;
end
if maxY > height
    maxY = height;
end
if minY < 1
    minY = 1;
end

% initialize
maxX2 = centerPoint(1)+2*r;
minX2 = centerPoint(1)-2*r;
maxY2 = centerPoint(2)+2*r;
minY2 = centerPoint(2)-2*r;

if maxX2 > width
    maxX2 = width;
end
if minX2 < 1
    minX2 = 1;
end
if maxY2 > height
    maxY2 = height;
end
if minY2 < 1
    minY2 = 1;
end

subImg = img(minY2:maxY2, minX2: maxX2);
% %%%%%%%%%%%%%%%%%%%%
% dImg = double(img)/max(img(:));
% drawSubImg = dImg(minY2:maxY2, minX2: maxX2);
% %%%%%%%%%%%%%%%%%%%%

%% run level set
tempSigma = 1.5;
G=fspecial('gaussian',15,tempSigma);
Img_smooth=conv2(subImg,G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.

% % initialize LSF as binary step function
c0=2;

initialMask = zeros(size(subImg));
centerPoint = ceil(size(initialMask)/2);
for i = 1:size(initialMask,1)
    for j = 1:size(initialMask,2)
        len = sqrt((i-centerPoint(1))^2+(j-centerPoint(2))^2);
        if len < r
            initialMask(i,j) = -c0;
        end
    end
end
% initialLSF=c0*ones(size(subImg));
% % generate the initial region R0 as a circle
% initialLSF(1+r:end-r, 1+r:end-r)=-c0;  
phi=initialMask;

% figure(2);
% %imagesc(imcomplement(drawSubImg),[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r', 'LineWidth',3);
% imshow(imcomplement(drawSubImg));hold on;  contour(phi, [0,0], 'r', 'LineWidth',3);hold off
% title('Initial zero level contour');
% pause;

potential=2;
if potential ==1
    potentialFunction = 'single-well';  % use single well potential p1(s)=0.5*(s-1)^2, which is good for region-based model 
elseif potential == 2
    potentialFunction = 'double-well';  % use double-well potential in Eq. (16), which is good for both edge and region based models
else
    potentialFunction = 'double-well';  % default choice of potential function
end


% start level set evolution
for n=1:iter_outer
    phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);
%     if mod(n,2)==0
%         figure(2);
% %         imagesc(imcomplement(drawSubImg),[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r', 'LineWidth', 3);
%         imshow(imcomplement(drawSubImg));hold on;  contour(phi, [0,0], 'r', 'LineWidth',3);hold off
%         pause
%     end
end
% refine the zero level contour by further level set evolution with alfa=0
alfa=0;
iter_refine = 10;
phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);


%% save the result
boundaryNeclues = c0*ones(size(img));
boundaryNeclues(minY2:maxY2, minX2: maxX2) = phi;

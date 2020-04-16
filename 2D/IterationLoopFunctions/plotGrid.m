function plotGrid(gridX, gridY)

%source code: https://github.com/stellaccl/cdmffd-image-registration
%paper: Chan, C. L., Anitescu, C., Zhang, Y., & Rabczuk, T. (2017). Two and three dimensional image registration based on B-spline composition and level sets. 
%Communications in Computational Physics, 21(2), 600-622.


%plots the deformation grid for visualization

%--Input Variable:
%gridX: the x-coordinate of the grid points
%gridY: the y-coordinate of the grid points

ngx = size(gridX,1);
ngy = size(gridX,2);

%plot the vertical lines
for i=1:ngx
    lineVectorX = zeros(1,ngy);
    lineVectorY = zeros(1,ngy);    
    for j=1:ngy
        lineVectorX(j)=gridX(i,j);
        lineVectorY(j)=gridY(i,j);
    end
    line(lineVectorX,lineVectorY);
    hold on
end

%plot the horizontal lines
for i=1:ngy
    lineVectorX = zeros(1,ngx);
    lineVectorY = zeros(1,ngx);    
    for j=1:ngx
        lineVectorX(j)=gridX(j,i);
        lineVectorY(j)=gridY(j,i);
    end
    line(lineVectorX,lineVectorY);
    hold on
end

        

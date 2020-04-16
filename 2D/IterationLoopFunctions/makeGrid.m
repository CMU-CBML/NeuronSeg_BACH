function [ gridX, gridY, pointsX, pointsY ] = makeGrid(ngx, ngy, nx, ny  )

%source code: https://github.com/stellaccl/cdmffd-image-registration
%paper: Chan, C. L., Anitescu, C., Zhang, Y., & Rabczuk, T. (2017). Two and three dimensional image registration based on B-spline composition and level sets. 
%Communications in Computational Physics, 21(2), 600-622.

%This function defines a grid of points

%--Input Variable:
%ngx, ngy: number of grid points in x and y direction
%nx,ny: image size nxXny

%--Output Variable:
%gridX, gridY: 2D mesh of grid points
%pointsX,pointsY: x and y coordinates of the points in the grid

pointsX = round(linspace(1,nx,ngx));
pointsY = round(linspace(1,ny,ngy));

[gridX, gridY] = meshgrid(pointsX, pointsY);

end

